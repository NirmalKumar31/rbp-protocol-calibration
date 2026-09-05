output "project_id" {
  value       = google_project.rbp.project_id
  description = "Use with: gcloud config set project <this>"
}

output "buckets" {
  value = {
    raw       = google_storage_bucket.raw.name
    derived   = google_storage_bucket.derived.name
    artifacts = google_storage_bucket.artifacts.name
    tfstate   = google_storage_bucket.tfstate.name
  }
}

output "image_repo" {
  value       = "${var.region}-docker.pkg.dev/${google_project.rbp.project_id}/${google_artifact_registry_repository.images.repository_id}"
  description = "Prefix for image tags. Use a DIGEST for anything whose output is recorded: <this>/cpu@sha256:... . The :latest alias is for interactive convenience and must never be the recorded producer of a published number."
}

# EVERY PRIVILEGED IDENTITY, not the four that run jobs. This listed ingest, prep, train and
# analysis and omitted the two that matter most in a security review: `modal`, which is the only
# one that gets a downloadable private key, and `killswitch`, which holds billing.admin at the
# billing-account scope and can detach the project. An operations summary that omits the two
# broadest identities is the one you would want to be complete.
output "service_accounts" {
  description = "All service accounts this configuration creates, with what each is for."
  value = {
    ingest     = google_service_account.ingest.email
    prep       = google_service_account.prep.email
    train      = google_service_account.train.email
    analysis   = google_service_account.analysis.email
    modal      = google_service_account.modal.email
    killswitch = google_service_account.killswitch.email
  }
}

output "privileged_identities" {
  description = "The two identities that are not job runners. Review these first."
  value = {
    modal = {
      email = google_service_account.modal.email
      note  = "The only identity with a downloadable key. Terraform creates no key resource on purpose, because google_service_account_key stores the private key in state in plaintext; the key is minted out of band and should be revoked when the run ends. See docs/REPRODUCE.md."
    }
    killswitch = {
      email = google_service_account.killswitch.email
      note  = "Holds roles/billing.admin at the BILLING ACCOUNT scope, which is far broader than detaching one project, plus roles/viewer on the project. cloud/killswitch/main.py checks the two permissions that actually authorise an unlink -- resourcemanager.projects.deleteBillingAssignment on the project, or billing.resourceAssociations.delete on the billing account -- because a dry run proves only the read, and because updateBillingInfo is an API method name rather than a permission. Neither check proves an org policy or link lock will permit the write; only a rehearsal in a disposable project does."
    }
  }
}

output "bigquery_runs_table" {
  value = "${google_project.rbp.project_id}.${google_bigquery_dataset.metrics.dataset_id}.${google_bigquery_table.runs.table_id}"
}
