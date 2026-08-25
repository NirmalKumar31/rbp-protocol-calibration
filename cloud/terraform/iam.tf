# Service accounts, one per job type, each with the narrowest role set that works.
#
# WHY NOT ONE SERVICE ACCOUNT FOR EVERYTHING. It would work, and it is what most tutorials
# do. Three reasons not to:
#
#   * A bug in the preprocessing job cannot delete the metrics table if it has no BigQuery
#     write role. Least privilege turns a class of mistakes into an error message.
#   * When something writes to a bucket unexpectedly, the audit log names which identity did
#     it. With one account, every log line says the same thing.
#   * It is the difference between "I used GCP" and "I understand IAM", which is the point
#     of building this at all.
#
# The pattern below is: one account per stage, bound to specific buckets rather than
# project-wide storage roles.

locals {
  # Roles every job needs.
  #
  # `batch.agentReporter` is the non-obvious one and it cost a stuck job to discover. The
  # Cloud Batch agent runs ON the VM AS the job's service account, and it must call
  # batch.states.report to tell the Batch service that the task has started. Without that
  # permission the agent retries forever, the job sits in SCHEDULED, and the VM bills the
  # whole time while doing nothing. The failure mode is silent: the VM is healthy, the
  # container may even be running, and the only symptom is a state that never advances.
  #
  # A least-privilege design has to include the permissions the PLATFORM needs, not only
  # the ones the workload needs.
  common_roles = [
    "roles/logging.logWriter",
    "roles/monitoring.metricWriter",
    "roles/batch.agentReporter",
  ]
}

# ---------------------------------------------------------------------------------------
# Preprocessing (Cloud Batch, CPU)
# ---------------------------------------------------------------------------------------

resource "google_service_account" "prep" {
  account_id   = "rbp-prep"
  display_name = "RBP preprocessing jobs"
  depends_on   = [google_project_service.apis]
}

resource "google_project_iam_member" "prep_common" {
  for_each = toset(local.common_roles)
  project  = google_project.rbp.project_id
  role     = each.value
  member   = "serviceAccount:${google_service_account.prep.email}"
}

# Reads raw data, writes derived. Cannot touch artifacts or the metrics table.
resource "google_storage_bucket_iam_member" "prep_raw_read" {
  bucket = google_storage_bucket.raw.name
  role   = "roles/storage.objectViewer"
  member = "serviceAccount:${google_service_account.prep.email}"
}

resource "google_storage_bucket_iam_member" "prep_derived_write" {
  bucket = google_storage_bucket.derived.name
  role   = "roles/storage.objectAdmin"
  member = "serviceAccount:${google_service_account.prep.email}"
}

# ---------------------------------------------------------------------------------------
# Data acquisition -- the only identity allowed to WRITE raw data
# ---------------------------------------------------------------------------------------

resource "google_service_account" "ingest" {
  account_id   = "rbp-ingest"
  display_name = "RBP data acquisition"
  depends_on   = [google_project_service.apis]
}

resource "google_project_iam_member" "ingest_common" {
  for_each = toset(local.common_roles)
  project  = google_project.rbp.project_id
  role     = each.value
  member   = "serviceAccount:${google_service_account.ingest.email}"
}

# Separated from `prep` deliberately. Raw data is downloaded once and is expensive to
# replace; nothing that runs 374 times in parallel should be able to overwrite it.
resource "google_storage_bucket_iam_member" "ingest_raw_write" {
  bucket = google_storage_bucket.raw.name
  role   = "roles/storage.objectAdmin"
  member = "serviceAccount:${google_service_account.ingest.email}"
}

# ---------------------------------------------------------------------------------------
# Training (Vertex AI, GPU)
# ---------------------------------------------------------------------------------------

resource "google_service_account" "train" {
  account_id   = "rbp-train"
  display_name = "RBP training jobs"
  depends_on   = [google_project_service.apis]
}

resource "google_project_iam_member" "train_common" {
  for_each = toset(local.common_roles)
  project  = google_project.rbp.project_id
  role     = each.value
  member   = "serviceAccount:${google_service_account.train.email}"
}

# Reads datasets and panels from anywhere in the bucket.
resource "google_storage_bucket_iam_member" "train_derived_read" {
  bucket = google_storage_bucket.derived.name
  role   = "roles/storage.objectViewer"
  member = "serviceAccount:${google_service_account.train.email}"
}

# ...and writes ONLY under runs/ and ckpt/.
#
# The comment above the read grant used to claim it also wrote checkpoints and scores. It
# did not: objectViewer is read-only, so every upload in cloud_train.py would have failed
# with 403 on the first task of the sweep. Caught before the sweep ran, but the lesson is
# that a comment describing an intent is not the intent being implemented.
#
# The obvious fix -- objectAdmin on the whole bucket -- would let a training job overwrite
# the datasets it is being evaluated on, and those datasets are the thing whose byte-level
# reproducibility this project spent two days establishing. An IAM CONDITION scopes the
# write to the two prefixes the sweep actually produces. Conditions on a bucket require
# uniform bucket-level access, which is on for all four buckets.
#
# Read as: this principal may act as an object admin, but only on objects whose full
# resource name begins with one of these prefixes. Anything else, including every
# dataset.tsv, is still read-only to it.
resource "google_storage_bucket_iam_member" "train_derived_write_runs" {
  bucket = google_storage_bucket.derived.name
  role   = "roles/storage.objectAdmin"
  member = "serviceAccount:${google_service_account.train.email}"

  condition {
    title       = "runs-and-checkpoints-only"
    description = "Sweep outputs and in-flight checkpoints, nothing else in the bucket"
    expression  = <<-EOT
      resource.name.startsWith("projects/_/buckets/${google_storage_bucket.derived.name}/objects/runs/") ||
      resource.name.startsWith("projects/_/buckets/${google_storage_bucket.derived.name}/objects/ckpt/")
    EOT
  }
}

resource "google_storage_bucket_iam_member" "train_artifacts_write" {
  bucket = google_storage_bucket.artifacts.name
  role   = "roles/storage.objectAdmin"
  member = "serviceAccount:${google_service_account.train.email}"
}

# Appends run metrics. dataEditor, not dataOwner: it can insert rows, not drop the table.
resource "google_bigquery_dataset_iam_member" "train_bq" {
  dataset_id = google_bigquery_dataset.metrics.dataset_id
  role       = "roles/bigquery.dataEditor"
  member     = "serviceAccount:${google_service_account.train.email}"
}

resource "google_project_iam_member" "train_bq_jobs" {
  project = google_project.rbp.project_id
  role    = "roles/bigquery.jobUser"
  member  = "serviceAccount:${google_service_account.train.email}"
}

# ---------------------------------------------------------------------------------------
# Analysis
# ---------------------------------------------------------------------------------------

resource "google_service_account" "analysis" {
  account_id   = "rbp-analysis"
  display_name = "RBP analysis and figures"
  depends_on   = [google_project_service.apis]
}

resource "google_project_iam_member" "analysis_common" {
  for_each = toset(local.common_roles)
  project  = google_project.rbp.project_id
  role     = each.value
  member   = "serviceAccount:${google_service_account.analysis.email}"
}

resource "google_storage_bucket_iam_member" "analysis_derived_read" {
  bucket = google_storage_bucket.derived.name
  role   = "roles/storage.objectViewer"
  member = "serviceAccount:${google_service_account.analysis.email}"
}

resource "google_storage_bucket_iam_member" "analysis_artifacts" {
  bucket = google_storage_bucket.artifacts.name
  role   = "roles/storage.objectAdmin"
  member = "serviceAccount:${google_service_account.analysis.email}"
}

resource "google_bigquery_dataset_iam_member" "analysis_bq" {
  dataset_id = google_bigquery_dataset.metrics.dataset_id
  role       = "roles/bigquery.dataEditor"
  member     = "serviceAccount:${google_service_account.analysis.email}"
}

resource "google_project_iam_member" "analysis_bq_jobs" {
  project = google_project.rbp.project_id
  role    = "roles/bigquery.jobUser"
  member  = "serviceAccount:${google_service_account.analysis.email}"
}

# ---------------------------------------------------------------------------------------
# Pulling images
# ---------------------------------------------------------------------------------------

# Every job needs to pull from Artifact Registry. Reader only -- a training job has no
# business pushing an image.
# for_each keys are the ACCOUNT IDS, not the resource-computed emails.
#
# Keying on google_service_account.*.email made the whole configuration unplannable from an
# empty state: for_each keys must be known at plan time, and those emails are "known only
# after apply". Terraform then refuses to plan at all, which forces a two-phase -target
# apply and, worse, blocks `terraform import` of anything else in the config.
#
# The email is fully determined by the account id and the project, both known at plan time,
# so building it from var.project_id removes the dependency without changing a single
# resulting grant.
resource "google_artifact_registry_repository_iam_member" "pullers" {
  for_each = toset(["rbp-prep", "rbp-ingest", "rbp-train", "rbp-analysis"])

  depends_on = [
    google_service_account.prep,
    google_service_account.ingest,
    google_service_account.train,
    google_service_account.analysis,
  ]

  location   = google_artifact_registry_repository.images.location
  repository = google_artifact_registry_repository.images.name
  role       = "roles/artifactregistry.reader"
  member     = "serviceAccount:${each.value}@${var.project_id}.iam.gserviceaccount.com"
}

# ---------------------------------------------------------------------------------------
# Off-GCP compute (Modal)
# ---------------------------------------------------------------------------------------

# WHY A SEPARATE IDENTITY FROM rbp-train.
#
# GPUS_ALL_REGIONS is 0 on this project and cannot be raised, and AWS and Azure gate GPU on
# young accounts the same way. Modal does not, so the three large models run there. The
# pipeline needs no code change -- a training task's only cloud dependency is
# google-cloud-storage -- but it does need credentials, and off GCP there is no metadata
# server to mint them.
#
# So this identity exists to hold a real, exportable key. rbp-train deliberately does not:
# it runs on GCE and gets short-lived tokens from the metadata server, so no secret for it
# exists anywhere. Giving it a key would throw that property away for every GCP job at once.
#
# A separate account means the key can be revoked, rotated or deleted without touching
# anything running on GCP, and the audit log distinguishes "a Modal container did this" from
# "a Batch task did this".
resource "google_service_account" "modal" {
  project      = google_project.rbp.project_id
  account_id   = "rbp-modal"
  display_name = "Off-GCP training workers (Modal)"
  depends_on   = [google_project_service.apis]
}

# Reads datasets, panels and the manifest.
resource "google_storage_bucket_iam_member" "modal_derived_read" {
  bucket = google_storage_bucket.derived.name
  role   = "roles/storage.objectViewer"
  member = "serviceAccount:${google_service_account.modal.email}"
}

# Writes ONLY under runs/ and ckpt/, exactly as rbp-train does. The condition matters more
# here than there: this is the one identity in the project whose credential leaves Google's
# network, so it is the one whose blast radius has to be smallest. Even fully compromised it
# cannot alter a dataset, touch the raw bucket, or delete a result outside those prefixes.
resource "google_storage_bucket_iam_member" "modal_derived_write" {
  bucket = google_storage_bucket.derived.name
  role   = "roles/storage.objectAdmin"
  member = "serviceAccount:${google_service_account.modal.email}"

  condition {
    title       = "runs-checkpoints-and-variants-only"
    description = "Sweep outputs, in-flight checkpoints and variant scores, nothing else"
    # variants/ was added when the ClinVar arm moved to Modal. The first probe scored its
    # dataset correctly and then took a 403 on the upload, which is the condition doing
    # exactly its job: a new write path stays denied until someone widens it deliberately.
    # Widened here rather than with a console click so state and reality stay in agreement;
    # a hand-made binding would be invisible to Terraform and reverted by the next apply.
    # rbp-train keeps the narrower condition, because it never writes variant scores.
    expression = <<-EOT
      resource.name.startsWith("projects/_/buckets/${google_storage_bucket.derived.name}/objects/runs/") ||
      resource.name.startsWith("projects/_/buckets/${google_storage_bucket.derived.name}/objects/ckpt/") ||
      resource.name.startsWith("projects/_/buckets/${google_storage_bucket.derived.name}/objects/variants/")
    EOT
  }
}

# NO key resource here, and that is deliberate. google_service_account_key would put the
# private key in Terraform state, and state lives in a GCS bucket -- turning one secret into
# two copies of a secret, one of them in a place nobody thinks of as a secret store. The key
# is minted with gcloud, piped straight into a Modal secret, and the local file deleted.

output "modal_service_account" {
  value       = google_service_account.modal.email
  description = "Identity for off-GCP workers. Key is minted out of band; see docs/49."
}
