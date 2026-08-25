# The whole GCP environment for this project, as one file tree.
#
# WHY INFRASTRUCTURE AS CODE HERE, rather than clicking in the console. Three reasons that
# are specific to this project rather than general good practice:
#
#   1. The research claim is reproducibility. "Run these two commands and you regenerate
#      every number in the paper" is only true if the environment is also reproducible.
#   2. The $300 credit expires around 2026-11-18. `terraform destroy` guarantees nothing is
#      left running and quietly billing after the work is done. Console-created resources
#      get forgotten; a state file does not.
#   3. Everything here can be reviewed before it exists. A budget alert written in a file is
#      visible; one you meant to set up in the console is not.
#
# ORDER MATTERS IN THIS FILE. The budget and its alerts are created with the project, before
# any resource capable of spending money exists. That is deliberate.

terraform {
  required_version = ">= 1.5"
  required_providers {
    google = {
      source  = "hashicorp/google"
      version = "~> 6.0"
    }
  }
  # State lives in GCS, not on the laptop.
  #
  # WHY THIS MATTERS MORE THAN IT LOOKS. Terraform state is the only record of which cloud
  # resources exist and are billing. Lose the local file and Terraform no longer knows about
  # the project's resources -- they keep running, keep charging, and `destroy` cannot reach
  # them. GCS also gives object-level locking, so two concurrent applies cannot corrupt it,
  # and the bucket is versioned, so a bad apply can be rolled back.
  #
  # Chicken-and-egg: this bucket is declared in storage.tf, so the FIRST apply necessarily
  # ran with local state. `terraform init -migrate-state` moved it here afterwards.
  backend "gcs" {
    bucket = "rbp-composition-2026-tfstate"
    prefix = "terraform/state"
  }
}

provider "google" {
  project = google_project.rbp.project_id
  region  = var.region
}

# A provider alias without a project, for the calls that create the project itself.
provider "google" {
  alias  = "bootstrap"
  region = var.region
}

# A third alias, purely for the Billing Budgets API.
#
# WHY THIS EXISTS. billingbudgets.googleapis.com refuses requests made with user
# Application Default Credentials unless a "quota project" is attached -- the project that
# gets billed for the API call itself, as distinct from the project the budget is about.
# Without it the call is attributed to Google's own default client project (764086051850)
# where the API is disabled, and returns SERVICE_DISABLED.
#
# `gcloud auth application-default set-quota-project` fixes this for client libraries but
# the Terraform provider does not pick it up; it needs `billing_project` plus
# `user_project_override = true` on the provider itself.
#
# It cannot be folded into the bootstrap provider, because that one has to create the
# project and cannot reference a project that does not exist yet.
provider "google" {
  alias                 = "billing"
  region                = var.region
  billing_project       = var.project_id
  user_project_override = true
}

# ---------------------------------------------------------------------------------------
# The project
# ---------------------------------------------------------------------------------------

resource "google_project" "rbp" {
  provider = google.bootstrap

  name            = var.project_name
  project_id      = var.project_id
  billing_account = var.billing_account

  # No organisation: this is a personal account, so the project sits directly under the
  # billing account. Stated explicitly because the alternative (an org) changes how IAM and
  # policy inheritance work, and someone reading this will wonder.
  labels = {
    purpose = "research"
    project = "rbp-composition-confound"
  }
}

# ---------------------------------------------------------------------------------------
# Budget guardrails -- FIRST, before anything can spend
# ---------------------------------------------------------------------------------------

resource "google_billing_budget" "guardrail" {
  provider = google.billing

  billing_account = var.billing_account
  display_name    = "rbp-guardrail"

  budget_filter {
    projects = ["projects/${google_project.rbp.number}"]

    # MEASURE GROSS SPEND, NOT SPEND AFTER CREDITS.
    #
    # The default is INCLUDE_ALL_CREDITS, which reports cost NET of credits. With a $300
    # trial credit covering everything, that number is $0.00 and stays $0.00 until the
    # credit is exhausted -- so every threshold below would have been unreachable, the
    # alerts would never have emailed, and the killswitch in killswitch.tf would never
    # have fired. A guardrail that cannot trip is worse than none, because it is trusted.
    #
    # EXCLUDE_ALL_CREDITS reports what the resources actually cost. That is the number
    # worth defending: the point is not to avoid a bill (there isn't one while credit
    # lasts), it is to avoid burning through a fixed, non-renewable $300.
    credit_types_treatment = "EXCLUDE_ALL_CREDITS"
  }

  depends_on = [google_project_service.apis]

  amount {
    specified_amount {
      currency_code = "USD"
      units         = tostring(var.budget_usd)
    }
  }

  # Alerts at 25/50/80/100% of budget. The estimate for the whole pipeline is ~$34 with a
  # ~$51 margin (docs/26), so a $100 budget means the first alert fires at $25 -- below the
  # point where anything has gone seriously wrong, and well below the credit.
  dynamic "threshold_rules" {
    for_each = [0.25, 0.5, 0.8, 1.0]
    content {
      threshold_percent = threshold_rules.value
      spend_basis       = "CURRENT_SPEND"
    }
  }

  # Forecast alert: fires when GCP predicts the month will exceed budget, which catches a
  # runaway job hours before actual spend would.
  threshold_rules {
    threshold_percent = 1.0
    spend_basis       = "FORECASTED_SPEND"
  }

  # THE LINE THAT TURNS A WARNING INTO A STOP. Without this the budget only emails. With
  # it, every re-evaluation (roughly every 20-30 minutes, not only when a threshold trips)
  # publishes the current spend to a topic the killswitch function listens on. See
  # killswitch.tf.
  all_updates_rule {
    pubsub_topic                   = google_pubsub_topic.billing.id
    schema_version                 = "1.0"
    disable_default_iam_recipients = false
  }
}

# ---------------------------------------------------------------------------------------
# APIs
# ---------------------------------------------------------------------------------------

# Enabled explicitly rather than on first use, so the set of capabilities this project has
# is reviewable in one place. `disable_on_destroy = false` because disabling an API on
# teardown can fail if another resource is still draining, and a failed destroy is worse
# than a left-enabled API (an enabled API costs nothing).
resource "google_project_service" "apis" {
  for_each = toset([
    # Required by the budget resource below. Its absence is why the FIRST apply created 47
    # of 48 resources and failed on exactly the one that was supposed to come first -- the
    # project existed with billing attached and no spend alerting. Ordering a resource first
    # does not help if the API it needs is not enabled.
    "billingbudgets.googleapis.com",
    # Scans images on push for known CVEs. Low stakes here -- no inbound services and no
    # secrets in the image -- but it is free at our volume and its absence is the kind of
    # thing that is noticed in a review.
    "containerscanning.googleapis.com",
    "cloudbuild.googleapis.com",       # container builds (no local Docker)
    "artifactregistry.googleapis.com", # image storage
    "storage.googleapis.com",          # data and artifacts
    "batch.googleapis.com",            # the 374 CPU preprocessing tasks
    "aiplatform.googleapis.com",       # Vertex AI pipelines and custom training jobs
    "bigquery.googleapis.com",         # metrics warehouse
    "run.googleapis.com",              # results API
    "logging.googleapis.com",
    "monitoring.googleapis.com",
    "compute.googleapis.com", # underlies Batch and Vertex
    "cloudresourcemanager.googleapis.com",
    "iam.googleapis.com",
    "pubsub.googleapis.com",         # budget -> killswitch channel
    "cloudfunctions.googleapis.com", # the killswitch itself
    "eventarc.googleapis.com",       # how gen2 functions receive Pub/Sub events
    "cloudbilling.googleapis.com",   # the API the killswitch calls to stop spend
  ])

  project            = google_project.rbp.project_id
  service            = each.value
  disable_on_destroy = false
}
