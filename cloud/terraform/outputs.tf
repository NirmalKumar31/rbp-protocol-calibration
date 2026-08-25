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
  description = "Prefix for image tags, e.g. <this>/cpu:latest"
}

output "service_accounts" {
  value = {
    ingest   = google_service_account.ingest.email
    prep     = google_service_account.prep.email
    train    = google_service_account.train.email
    analysis = google_service_account.analysis.email
  }
}

output "bigquery_runs_table" {
  value = "${google_project.rbp.project_id}.${google_bigquery_dataset.metrics.dataset_id}.${google_bigquery_table.runs.table_id}"
}
