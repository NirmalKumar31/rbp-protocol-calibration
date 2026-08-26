#!/bin/bash
# Runs on the driver VM at boot. Sequences EVERY remaining stage, including Modal, so the
# operator's laptop is not part of the pipeline at any point.
#
# RESUMABILITY IS THE WHOLE DESIGN. Every stage is guarded by a completion marker in GCS, so
# this script can be killed at any instant and restarted with no loss and no double work:
#
#   variants   guarded by results/variants-complete.json
#   clinvar    guarded by per-dataset objects under variants/scores_sb/ and scores_mm/
#   analysis   guarded by results/analysis-complete.json
#
# The markers live under results/ rather than at the bucket root because rbp-analysis is
# restricted by an IAM condition to results/, variants/ and driver/. A root-level marker is
# a 403 at the end of an otherwise successful stage.
#
# If the VM is preempted, deleted, or the script crashes, relaunching it resumes from the
# first incomplete stage. Nothing here is idempotent by accident; each check is deliberate.
#
# It logs to GCS rather than local disk, because a log you can only read by sshing in is a
# log you will not read.

set -uo pipefail
META="http://metadata.google.internal/computeMetadata/v1/instance/attributes"
PROJECT=$(curl -s -H "Metadata-Flavor: Google" "$META/rbp-project")
DERIVED=$(curl -s -H "Metadata-Flavor: Google" "$META/rbp-derived")
REGION=us-central1
LOG=/tmp/driver.log
export GOOGLE_CLOUD_PROJECT="$PROJECT" DERIVED_BUCKET="$DERIVED" RAW_BUCKET="${PROJECT}-raw"
export RBP_YES=1 NO_WAIT=1

say() {
  echo "[$(date -u '+%H:%M:%S')] $*" | tee -a "$LOG"
  gsutil -q cp "$LOG" "gs://${DERIVED}/driver/driver.log" 2>/dev/null || true
}
die() { say "STOP: $*"; gsutil -q cp "$LOG" "gs://${DERIVED}/driver/driver.log"; exit 1; }

wait_job() {
  local job=$1 state
  while :; do
    state=$(gcloud batch jobs describe "$job" --project="$PROJECT" --location="$REGION" \
              --format="value(status.state)" 2>/dev/null)
    case "$state" in
      SUCCEEDED) say "  $job SUCCEEDED"; return 0 ;;
      FAILED|CANCELLED) say "  $job $state"; return 1 ;;
      # An unreadable state means the API call failed, NOT that the job is gone. Assuming
      # completion here is how a monitor reports success at 57% -- it happened once already.
      "") say "  cannot read $job state; assuming alive" ;;
      *) say "  $job $state" ;;
    esac
    sleep 60
  done
}

say "driver booting on $(hostname)"

# --- environment ------------------------------------------------------------------------
apt-get update -qq && apt-get install -y -qq python3-pip python3-venv >>"$LOG" 2>&1
mkdir -p /opt/rbp && cd /opt/rbp || die "cannot make workdir"
gsutil -q cp "gs://${DERIVED}/driver/repo.tgz" . && tar xzf repo.tgz || die "cannot fetch repo"
python3 -m venv /opt/venv >>"$LOG" 2>&1
/opt/venv/bin/pip install -q --upgrade pip >>"$LOG" 2>&1
/opt/venv/bin/pip install -q modal google-cloud-storage pandas numpy scipy pyyaml >>"$LOG" 2>&1 \
  || die "pip install failed"
export PY=/opt/venv/bin/python
export PATH="/opt/venv/bin:$PATH"
export PYTHONPATH=/opt/rbp/src
say "environment ready"

# --- wait for any in-flight image build --------------------------------------------------
# The variants stage needs the image containing pyBigWig built against libcurl. Submitting
# before the build lands would pin the OLD digest and reproduce the exact failure this rebuild
# exists to fix -- the same class as the stale-image bug that made a job report 189/189
# SUCCEEDED while doing nothing.
say "waiting for any in-flight Cloud Build"
for i in $(seq 1 60); do
  # NOT --ongoing: it returned zero while a build was plainly WORKING, so the driver used
  # the stale digest. Filter on the status field, which is what the console shows.
  W=$(gcloud builds list --project="$PROJECT" \
        --filter="status=WORKING OR status=QUEUED" --format="value(id)" 2>/dev/null | wc -l | tr -d ' ')
  [ "${W:-0}" -eq 0 ] && { say "  no builds in flight"; break; }
  say "  $W build(s) still running"
  sleep 60
done
say "image digest in use: $(gsutil cat gs://${PROJECT}-artifacts/images/cpu_digest.txt 2>/dev/null | head -c 22)"

# --- stage 11: variants -------------------------------------------------------------------
if gsutil -q stat "gs://${DERIVED}/results/variants-complete.json" 2>/dev/null; then
  say "stage 11 already complete, skipping"
else
  say "stage 11: variants (assign -> score -> phylop)"
  OUT=$(cd /opt/rbp && ./cloud/submit.sh variants 2>&1); echo "$OUT" >>"$LOG"
  JOB=$(echo "$OUT" | grep -oE "submitted [a-z0-9-]+" | awk '{print $2}')
  [ -n "$JOB" ] || die "variants did not submit: $(echo "$OUT" | tail -3)"
  say "  submitted $JOB"
  wait_job "$JOB" || die "variants failed; rerun the driver to resume"
fi

# --- stage 12a: the four-model table ------------------------------------------------------
# THIS RUNS BEFORE MODAL AND THAT ORDER IS THE POINT.
#
# The window cutter restricts to the datasets that have all four models, which it reads from
# results/tables/matched_four_models.csv. That table is written by cloud_analysis.four_models()
# -- a stage 13 function. So stage 12 depends on an artefact of stage 13, and the original
# order here (variants -> clinvar -> analysis) would have failed on a missing file after
# booting a VM and installing a toolchain.
#
# four_models() only needs results/rehearsal_binding_dinuc.csv and results/sweep_dinuc.csv,
# both of which stages 7 to 9 already produced, so it can run now. `--what tables` writes no
# completion marker, which is what makes running the analysis twice safe.
if gsutil -q stat "gs://${DERIVED}/results/tables/matched_four_models.csv" 2>/dev/null; then
  say "stage 12a: four-model table already present, skipping"
else
  say "stage 12a: four-model table (needed by the window cutter)"
  OUT=$(cd /opt/rbp && ANALYSIS_WHAT=tables ./cloud/submit.sh analysis 2>&1); echo "$OUT" >>"$LOG"
  JOB=$(echo "$OUT" | grep -oE "submitted [a-z0-9-]+" | awk '{print $2}')
  [ -n "$JOB" ] || die "analysis(tables) did not submit: $(echo "$OUT" | tail -3)"
  say "  submitted $JOB"
  wait_job "$JOB" || die "analysis(tables) failed"
fi

# --- stage 12b: cut the variant windows ---------------------------------------------------
# Runs as a Batch job rather than here, because it needs the 3.1 GB genome and pyfaidx, and
# cloud_variants.py already stages both. Doing it on this 2 GB VM would have meant a second
# genome download and a longer pip list.
if gsutil -q stat "gs://${DERIVED}/variants/variant_tasks.tsv" 2>/dev/null; then
  say "stage 12b: windows already cut, skipping"
else
  say "stage 12b: cutting ref/alt windows for every ClinVar variant"
  OUT=$(cd /opt/rbp && VARIANTS_WHAT=windows ./cloud/submit.sh variants 2>&1); echo "$OUT" >>"$LOG"
  JOB=$(echo "$OUT" | grep -oE "submitted [a-z0-9-]+" | awk '{print $2}')
  [ -n "$JOB" ] || die "windows did not submit: $(echo "$OUT" | tail -3)"
  say "  submitted $JOB"
  wait_job "$JOB" || die "window cutting failed"
fi

# --- stage 12c: ClinVar on Modal ----------------------------------------------------------
# The token comes from GCS, not from this script, and is removed from disk immediately after
# use. The bucket is private and the key is the one credential that leaves Google's network,
# so it should be rotated once the run finishes.
SB=$(gsutil ls "gs://${DERIVED}/variants/scores_sb/" 2>/dev/null | wc -l | tr -d ' ')
MM=$(gsutil ls "gs://${DERIVED}/variants/scores_mm/" 2>/dev/null | wc -l | tr -d ' ')
if [ "${SB:-0}" -ge 94 ] && [ "${MM:-0}" -ge 94 ]; then
  say "stage 12c already complete ($SB matched, $MM mismatched), skipping"
else
  say "stage 12c: ClinVar on Modal (matched + mismatched heads)"
  gsutil -q cp "gs://${DERIVED}/driver/modaltok" /tmp/tok || die "no modal token staged"
  modal token set --token-id "$(sed -n 1p /tmp/tok)" \
                  --token-secret "$(sed -n 2p /tmp/tok)" --profile=driver >>"$LOG" 2>&1
  modal profile activate driver >>"$LOG" 2>&1
  rm -f /tmp/tok
  cd /opt/rbp || die "workdir gone"
  # --detach: the app must outlive this shell, exactly as it must outlive a laptop.
  modal run --detach cloud/modal/modal_variants.py::sweep >>"$LOG" 2>&1 || say "  matched arm returned nonzero"
  modal run --detach cloud/modal/modal_variants.py::mismatch_sweep >>"$LOG" 2>&1 || say "  mismatched arm returned nonzero"
  say "  both Modal arms dispatched"
  for i in $(seq 1 60); do
    SB=$(gsutil ls "gs://${DERIVED}/variants/scores_sb/" 2>/dev/null | wc -l | tr -d ' ')
    MM=$(gsutil ls "gs://${DERIVED}/variants/scores_mm/" 2>/dev/null | wc -l | tr -d ' ')
    say "  matched $SB/94  mismatched $MM/94"
    [ "${SB:-0}" -ge 94 ] && [ "${MM:-0}" -ge 94 ] && break
    sleep 60
  done
fi

# --- stage 13: analysis --------------------------------------------------------------------
if gsutil -q stat "gs://${DERIVED}/results/analysis-complete.json" 2>/dev/null; then
  say "stage 13 already complete, skipping"
else
  say "stage 13: aggregate + figures"
  OUT=$(cd /opt/rbp && ./cloud/submit.sh analysis 2>&1); echo "$OUT" >>"$LOG"
  JOB=$(echo "$OUT" | grep -oE "submitted [a-z0-9-]+" | awk '{print $2}')
  [ -n "$JOB" ] && { say "  submitted $JOB"; wait_job "$JOB" || say "  analysis failed"; }
fi

say "DRIVER FINISHED. Remaining: stage 14 verify (local, reads GCS, costs nothing)."
gsutil -q cp /dev/null "gs://${DERIVED}/driver/DONE" 2>/dev/null || true
gsutil -q cp "$LOG" "gs://${DERIVED}/driver/driver.log" 2>/dev/null || true
