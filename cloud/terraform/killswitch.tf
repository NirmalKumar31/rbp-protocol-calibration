# The hard stop.
#
# main.tf already declares a $100 budget with alerts at 25/50/80/100%. Those alerts send
# EMAIL and nothing else -- GCP has no built-in "stop spending" control. If a job runs away
# overnight the alert arrives and the spend continues.
#
# This file wires the only mechanism that actually halts charges: budget -> Pub/Sub ->
# Cloud Function -> detach the project from its billing account. Every VM stops, Batch
# fails, Cloud Run stops serving. Buckets and their contents survive; storage becomes
# inaccessible until billing is re-attached, so the worst case is an interrupted run and a
# manual re-enable, never lost data.
#
# UI EQUIVALENT, for the record: none. The console can create the budget, the topic and the
# function separately, but there is no console feature that stops spend. This has to be
# assembled, which is exactly why it belongs in Terraform rather than in someone's memory.

# ---------------------------------------------------------------------------------------
# The channel the budget shouts down
# ---------------------------------------------------------------------------------------

resource "google_pubsub_topic" "billing" {
  project = google_project.rbp.project_id
  name    = "billing-alerts"

  depends_on = [google_project_service.apis]
}

# NO EXPLICIT PUBLISHER BINDING, AND THAT NEEDED CHECKING RATHER THAN ASSUMING.
#
# Google's own documentation names `billing-budgets-pubsub.iam.gserviceaccount.com` as the
# publishing identity, so that binding was written first. The IAM API rejects it:
# "Error 400: Invalid service account". It is not a real principal you can grant to.
#
# The Billing Budgets service provisions its own access when the budget's all_updates_rule
# accepts a topic, which it did. That is convenient and also exactly the situation this
# comment exists to warn about: a guardrail that LOOKS installed and publishes nothing is
# worse than no guardrail, because it buys false confidence. So it is verified by
# observation, not by inference -- see the pull-and-confirm step in
# docs/45, and the `killswitch-probe` subscription which exists solely so that
# "did a real budget message arrive?" is answerable at any time.

# ---------------------------------------------------------------------------------------
# Identity allowed to pull the plug
# ---------------------------------------------------------------------------------------

resource "google_service_account" "killswitch" {
  project      = google_project.rbp.project_id
  account_id   = "rbp-killswitch"
  display_name = "Disables billing when the budget is exceeded"
  depends_on   = [google_project_service.apis]
}

# Detaching a project from its billing account is a permission on the BILLING ACCOUNT, not
# on the project. roles/billing.admin is the narrowest role that includes it.
#
# This is the most privileged identity in the project, and it is worth being explicit about
# why that is acceptable: it can do exactly one thing that matters, it is not attached to
# any VM, only the Pub/Sub topic can invoke it, and the alternative is having no hard stop
# at all. It cannot spend money -- only stop money.
resource "google_billing_account_iam_member" "killswitch" {
  billing_account_id = var.billing_account
  role               = "roles/billing.admin"
  member             = "serviceAccount:${google_service_account.killswitch.email}"
}

# Reading the project's own billing state before changing it.
resource "google_project_iam_member" "killswitch_viewer" {
  project = google_project.rbp.project_id
  role    = "roles/viewer"
  member  = "serviceAccount:${google_service_account.killswitch.email}"
}

# ---------------------------------------------------------------------------------------
# The function
# ---------------------------------------------------------------------------------------

# Zipped from the working tree so the deployed code is the code in this repo, not something
# uploaded by hand once and forgotten. The hash in the object name means a change to main.py
# forces a redeploy; without it Terraform would see the same object name and do nothing.
data "archive_file" "killswitch" {
  type        = "zip"
  source_dir  = "${path.module}/../killswitch"
  output_path = "${path.module}/.killswitch.zip"
}

# The Cloud Functions service agent copies the source zip into its own staging bucket, so
# it needs read access to ours. Without it: "Error 403: Could not clone object ... does not
# have storage.objects.get access". The agent only exists once cloudfunctions.googleapis.com
# is enabled, which is why this cannot be granted before the API.
resource "google_storage_bucket_iam_member" "gcf_source_read" {
  bucket     = google_storage_bucket.artifacts.name
  role       = "roles/storage.objectViewer"
  member     = "serviceAccount:service-${google_project.rbp.number}@gcf-admin-robot.iam.gserviceaccount.com"
  depends_on = [google_project_service.apis]
}

resource "google_storage_bucket_object" "killswitch" {
  bucket = google_storage_bucket.artifacts.name
  name   = "functions/killswitch-${data.archive_file.killswitch.output_md5}.zip"
  source = data.archive_file.killswitch.output_path
}

resource "google_cloudfunctions2_function" "killswitch" {
  project  = google_project.rbp.project_id
  name     = "billing-killswitch"
  location = var.region

  build_config {
    runtime     = "python312"
    entry_point = "handle"
    source {
      storage_source {
        bucket = google_storage_bucket.artifacts.name
        object = google_storage_bucket_object.killswitch.name
      }
    }
  }

  service_config {
    # Smallest possible. It runs for a second every half hour and does one API call.
    available_memory   = "256M"
    timeout_seconds    = 60
    max_instance_count = 1
    # No concurrency: two copies racing to disable billing would both succeed, which is
    # harmless, but one is simpler to reason about.
    service_account_email = google_service_account.killswitch.email

    environment_variables = {
      TARGET_PROJECT = google_project.rbp.project_id
      # Deliberately BELOW the $100 budget. The budget's job is to warn; this one's job is
      # to stop, and it should stop while there is still credit left to recover with. The
      # whole project is estimated under $11, so $40 is far above any legitimate run and
      # far below the $300 credit.
      KILL_THRESHOLD_USD = "40"
      # Starts in dry-run: it logs what it WOULD do. Flipped to false only after the path
      # has been proven end to end, because the alternative way to test it is to take the
      # project down.
      DRY_RUN = var.killswitch_armed ? "false" : "true"
    }
  }

  event_trigger {
    trigger_region = var.region
    event_type     = "google.cloud.pubsub.topic.v1.messagePublished"
    pubsub_topic   = google_pubsub_topic.billing.id
    retry_policy   = "RETRY_POLICY_RETRY"
  }

  depends_on = [
    google_project_service.apis,
    google_billing_account_iam_member.killswitch,
    google_storage_bucket_iam_member.gcf_source_read,
  ]
}

output "killswitch" {
  value = {
    function  = google_cloudfunctions2_function.killswitch.name
    topic     = google_pubsub_topic.billing.name
    armed     = var.killswitch_armed
    kills_at  = "40 USD"
    re_enable = "gcloud billing projects link ${google_project.rbp.project_id} --billing-account=${var.billing_account}"
  }
  description = "The hard stop. `armed = false` means it logs instead of acting."
}

# Proof that the budget actually publishes.
#
# The function consumes events through Eventarc, which does not leave messages behind to
# inspect. This subscription receives the same messages independently, so at any point
# `gcloud pubsub subscriptions pull killswitch-probe` answers "has a real budget message
# ever arrived?" -- the one question that distinguishes an installed guardrail from an
# inert one. Costs nothing; Pub/Sub bills on message volume and this is a few per hour.
resource "google_pubsub_subscription" "probe" {
  project = google_project.rbp.project_id
  name    = "killswitch-probe"
  topic   = google_pubsub_topic.billing.id

  # Keep a week, so a message is still there on Monday if something is checked on Friday.
  message_retention_duration = "604800s"
  retain_acked_messages      = false
  ack_deadline_seconds       = 20
}
