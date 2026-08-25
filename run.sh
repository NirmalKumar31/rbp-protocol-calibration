#!/usr/bin/env bash
# The whole study, raw inputs to verified results, cloud only.
#
#   ./run.sh preflight          check GCP, spend nothing
#   ./run.sh preflight-modal    check Modal (run AFTER stage 1 creates the service account)
#   ./run.sh stage 3            run one stage
#   ./run.sh from 6             run stage 6 onward
#   ./run.sh all                everything (still stops at each paid gate)
#   ./run.sh status             where the artefacts are
#
# DESIGN RULES, each one paid for by a specific failure in the original build.
#
# 1. NOTHING RUNS BEFORE PREFLIGHT PASSES. Every failure of the first build was a missing
#    API, zero quota, a budget that reported $0, or an unauthenticated client -- all free to
#    detect, all discovered late, after spending.
#
# 2. PAID STAGES REQUIRE AN EXPLICIT CONFIRMATION. Stage 8 is ~$31 of real money and 95% of
#    the budget. It asks. Set RBP_YES=1 to run unattended once you have decided.
#
# 3. EVERY STAGE IS IDEMPOTENT AND RESUMABLE. Completion markers are written LAST, so a
#    stage killed midway redoes its work instead of being skipped. Rerunning a finished
#    stage costs seconds.
#
# 4. NO LOCAL COMPUTE. Every stage runs in a container on Batch or on Modal. The laptop
#    submits and reads; it never computes. That is the whole point of this rebuild.
#
# 5. STAGES 3 AND 10 NEED PUBLIC INTERNET (ENCODE, GENCODE, NCBI, UCSC phyloP). Workers have
#    Private Google Access only and no NAT by design, so those two run on a VM with an
#    external IP. Everything else stays sealed.

set -uo pipefail
cd "$(dirname "$0")" || exit 1

PY="${PY:-python3}"
export PYTHONPATH="$PWD/src:${PYTHONPATH:-}"

PROJECT_ID="${GOOGLE_CLOUD_PROJECT:-$($PY -c 'import sys;sys.path.insert(0,"src");from rbp.utils import cloud;print(cloud.project())' 2>/dev/null)}"
DERIVED="${DERIVED_BUCKET:-${PROJECT_ID}-derived}"
RAW="${RAW_BUCKET:-${PROJECT_ID}-raw}"
REGION="${REGION:-us-central1}"
EVERY="${EVERY:-2}"                 # panel: keep every Nth dataset by pair rank
export GOOGLE_CLOUD_PROJECT="$PROJECT_ID" DERIVED_BUCKET="$DERIVED" RAW_BUCKET="$RAW"

mkdir -p logs
LOG="logs/run-$(date +%Y%m%d).log"

say()  { printf '[%s] %s\n' "$(date '+%H:%M:%S')" "$*" | tee -a "$LOG"; }
die()  { say "STOP: $*"; exit 1; }

# A paid gate. Never silently spends.
confirm() {
  local what=$1 cost=$2
  say "PAID STAGE: $what  (estimated $cost)"
  if [ "${RBP_YES:-0}" = "1" ]; then say "RBP_YES=1, proceeding"; return 0; fi
  printf 'Type YES to spend %s on "%s": ' "$cost" "$what"
  read -r ans; [ "$ans" = "YES" ] || die "not confirmed"
}

# Has a stage already produced its marker? Used to make reruns cheap, never to skip work
# that only half happened -- markers are written after the payload, everywhere.
has() { gsutil -q stat "gs://${DERIVED}/$1" 2>/dev/null; }

gate_preflight() {
  [ -f .preflight-ok ] || die "run ./run.sh preflight first (it spends nothing)"
}

# THE MODAL GATE IS SEPARATE, AND IT HAS TO BE. The rbp-gcp secret holds a key for the
# rbp-modal service account, and that account does not exist until stage 1 has run. So a
# single all-or-nothing preflight is unsatisfiable on a fresh project: it would demand a
# secret that cannot yet be created, to permit the stage that creates its prerequisite.
# GCP stages need only the GCP gate; Modal stages need both.
gate_modal() {
  [ -f .preflight-modal-ok ] || die \
    "Modal is not verified yet. After stage 1:
       gcloud iam service-accounts keys create /tmp/k.json \\
         --iam-account=rbp-modal@${PROJECT_ID}.iam.gserviceaccount.com
       modal secret delete rbp-gcp 2>/dev/null || true
       modal secret create rbp-gcp SERVICE_ACCOUNT_JSON=\"\$(cat /tmp/k.json)\"
       rm /tmp/k.json
     then confirm your Modal balance and run:  ./run.sh preflight-modal"
}

# ---------------------------------------------------------------------------------------

s0_preflight() {
  say "stage 0: preflight (GCP only -- Modal is gated separately, see gate_modal)"
  $PY scripts/preflight.py --skip-modal ${PREFLIGHT_ARGS:-} \
    || die "preflight failed; fix and rerun"
  touch .preflight-ok
  say "GCP preflight recorded. GCP stages may run. Modal stages need ./run.sh preflight-modal"
}

preflight_modal() {
  say "verifying Modal: secret, auth, and a human-confirmed balance"
  $PY scripts/preflight.py --modal-credit-ok ${PREFLIGHT_ARGS:-} \
    || die "Modal preflight failed"
  touch .preflight-modal-ok
  say "Modal verified. Stages 9, 10 and 12 may run."
}

s1_terraform() {
  gate_preflight; say "stage 1: terraform -- buckets, service accounts, IAM, budget, killswitch"
  [ -f cloud/terraform/terraform.tfvars ] || die \
    "cloud/terraform/terraform.tfvars missing. Copy terraform.tfvars.example and set
     project_id and billing_account. It is gitignored on purpose."

  # The state bucket is per project and must exist before init. Creating it here rather than
  # in Terraform avoids the bootstrap paradox of storing state about the bucket that stores
  # the state.
  local TFSTATE="${PROJECT_ID}-tfstate"
  if ! gsutil ls -b "gs://${TFSTATE}" >/dev/null 2>&1; then
    say "creating state bucket gs://${TFSTATE}"
    gcloud storage buckets create "gs://${TFSTATE}" --project="$PROJECT_ID" \
      --location="$REGION" --uniform-bucket-level-access \
      || die "could not create the state bucket"
    gcloud storage buckets update "gs://${TFSTATE}" --versioning || true
  fi

  # -reconfigure, and the bucket passed explicitly. Without this, init reuses whatever
  # backend a previous checkout cached in .terraform/, which is how a plan ends up loading
  # another project's state and proposing to destroy it.
  ( cd cloud/terraform && terraform init -input=false -reconfigure \
      -backend-config="bucket=${TFSTATE}" ) || die "terraform init"

  ( cd cloud/terraform && terraform plan -input=false -no-color -out=tfplan.new ) \
    || die "terraform plan"

  # THE GUARD. A first apply on an empty project adds resources and destroys nothing. Any
  # destroy means the state does not describe this project, and applying would delete real
  # infrastructure somewhere else.
  local DESTROYS
  DESTROYS=$( cd cloud/terraform && terraform show -no-color tfplan.new \
              | grep -cE "^  # .* will be destroyed" || true )
  if [ "${DESTROYS:-0}" -gt 0 ]; then
    ( cd cloud/terraform && terraform show -no-color tfplan.new | grep -E "will be destroyed" | head -20 )
    die "plan contains ${DESTROYS} DESTROY actions. On a fresh project it must contain none.
     The state being read almost certainly belongs to a different project. Refusing to apply."
  fi
  say "plan is additive only (${DESTROYS} destroys). Applying."
  ( cd cloud/terraform && terraform apply -input=false -auto-approve tfplan.new ) \
    || die "terraform apply failed"
  ( cd cloud/terraform && rm -f tfplan.new )
}

s2_images() {
  gate_preflight; say "stage 2: build container images"
  confirm "Cloud Build, CPU + GPU images" "~\$0.50"
  gcloud builds submit --project="$PROJECT_ID" --config=docker/cloudbuild.cpu.yaml . || die "cpu image"
  gcloud builds submit --project="$PROJECT_ID" --config=docker/cloudbuild.gpu.yaml . || die "gpu image"
}

s3_ingest() {
  gate_preflight; say "stage 3: ingest raw inputs (NEEDS PUBLIC INTERNET -> external IP)"
  if has "raw-complete.json"; then say "already ingested, skipping"; return; fi
  confirm "ingest genome + GENCODE + ClinVar + ENCODE peaks" "~\$0.20"
  ./cloud/submit.sh ingest || die "ingest"
}

s4_panel() {
  gate_preflight; say "stage 4: build panel, both arms"
  ./cloud/submit.sh panel || die "panel"
}

s5_prep() {
  # PREP RUNS ON EVERY CANDIDATE, BEFORE THE PANEL EXISTS, AND THAT ORDER IS FORCED.
  # `pairs` counts the positives that could actually be matched to a negative, so it is a
  # RESULT of preprocessing, not an input. The study panel is a size-ranked sample, so it
  # cannot be chosen until prep has produced the counts. Prep is ~$2 for the full candidate
  # set, so nothing is saved by trying to invert this.
  gate_preflight; say "stage 5: preprocess ALL candidates, both arms, then finalize"
  confirm "preprocessing, both negative arms, full candidate set" "~\$2"
  $PY scripts/cloud_prep.py index    || die "prep index"
  $PY scripts/cloud_prep.py manifest || die "prep manifest"
  ./cloud/submit.sh prep             || die "prep"     # one job, both arms
  for arm in dinuc gc; do $PY scripts/cloud_prep.py finalize --arm "$arm" || die "finalize $arm"; done
}

s6_select() {
  gate_preflight; say "stage 6: define THE study panel (every=$EVERY)"
  $PY scripts/select_panel.py --every "$EVERY" || die "select_panel"
  $PY scripts/select_panel.py --show
  say "every later stage reads manifest/study_panel.tsv. Nothing else decides membership."
}

s7_rehearsal() {
  gate_preflight; say "stage 7: composition + k-mer, both arms  -> R1"
  confirm "rehearsal, both arms" "~\$0.60"
  for arm in dinuc gc; do
    ARM=$arm $PY scripts/cloud_rehearsal.py manifest --arm "$arm" || die "rehearsal manifest $arm"
    ARM=$arm ./cloud/submit.sh rehearsal || die "rehearsal $arm"
  done
  for arm in dinuc gc; do $PY scripts/cloud_rehearsal.py aggregate --arm "$arm"; done
}

s8_cnn() {
  gate_preflight; say "stage 8: CNN sweep  -> R2"
  confirm "CNN, 5 folds per dataset, Batch CPU" "~\$3"
  $PY scripts/cloud_train.py manifest --arm dinuc --models cnn --tag _cnn || die "cnn manifest"
  MODELS=cnn MANIFEST_TAG=_cnn ARM=dinuc ./cloud/submit.sh sweep || die "cnn"
  $PY scripts/cloud_train.py aggregate --arm dinuc --models cnn
}

s9_splicebert() {
  gate_preflight; gate_modal; say "stage 9: SpliceBERT sweep on Modal  -> R2"
  say "This is 95% of the money. Confirm your Modal balance is >= \$35 before continuing."
  confirm "SpliceBERT, 5 folds per dataset, Modal A10G" "~\$31 OUT OF POCKET"
  modal run cloud/modal/modal_sweep.py::sweep || die "splicebert"
  $PY scripts/cloud_train.py aggregate --arm dinuc --models splicebert
}

s10_locality() {
  gate_preflight; gate_modal; say "stage 10: ISM locality probe on Modal  -> R3"
  confirm "locality probe, Modal T4" "~\$0.30"
  modal run cloud/modal/modal_variants.py::locality_sweep || die "locality"
  $PY scripts/locality_probe.py --gather
}

s11_variants() {
  gate_preflight; say "stage 11: ClinVar assignment + phyloP (NEEDS PUBLIC INTERNET)"
  confirm "variant assignment and conservation fetch" "~\$0.30"
  ./cloud/submit.sh variants || die "variants"
}

s12_clinvar() {
  gate_preflight; gate_modal; say "stage 12: ClinVar scoring + mismatched-head control on Modal  -> R4"
  confirm "ClinVar scoring, matched and mismatched, Modal T4" "~\$0.60"
  $PY scripts/variant_splicebert.py --what tables || die "variant tables"
  modal run cloud/modal/modal_variants.py::sweep || die "clinvar matched"
  modal run cloud/modal/modal_variants.py::mismatch_sweep || die "clinvar mismatched"
  $PY scripts/variant_splicebert.py --what gather
  $PY scripts/variant_splicebert.py --what test
}

s13_analysis() {
  gate_preflight; say "stage 13: aggregate + figures"
  ./cloud/submit.sh analysis || die "analysis"
}

s14_verify() {
  say "stage 14: verify against golden numbers"
  $PY scripts/verify.py || die "THE SCIENCE DID NOT REPRODUCE -- see the failed claims above"
  say "reproduction verified"
}

STAGES=(s0_preflight s1_terraform s2_images s3_ingest s4_panel s5_prep s6_select \
        s7_rehearsal s8_cnn s9_splicebert s10_locality s11_variants s12_clinvar \
        s13_analysis s14_verify)

status() {
  say "project=$PROJECT_ID derived=$DERIVED raw=$RAW region=$REGION panel=every-$EVERY"
  for k in raw-complete.json manifest/study_panel.tsv \
           results/rehearsal_binding_dinuc.csv results/rehearsal_binding_gc.csv \
           results/tables/locality_ism.csv results/tables/variant_ladder.csv; do
    if has "$k"; then echo "  present  $k"; else echo "  MISSING  $k"; fi
  done
}

usage() { sed -n '2,30p' "$0"; }

case "${1:-}" in
  preflight)       s0_preflight ;;
  preflight-modal) preflight_modal ;;
  stage)     [ $# -ge 2 ] || die "which stage?"; "${STAGES[$2]}" ;;
  from)      [ $# -ge 2 ] || die "from which stage?"
             for i in $(seq "$2" $((${#STAGES[@]} - 1))); do "${STAGES[$i]}"; done ;;
  all)       for s in "${STAGES[@]}"; do "$s"; done ;;
  status)    status ;;
  *)         usage ;;
esac
