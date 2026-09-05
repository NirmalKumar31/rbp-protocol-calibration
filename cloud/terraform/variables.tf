variable "project_id" {
  description = "GCP project ID. Must be globally unique, 6-30 chars, lowercase."
  type        = string
  default     = "rbp-composition-2026"

  validation {
    condition     = can(regex("^[a-z][a-z0-9-]{5,29}$", var.project_id))
    error_message = "Project IDs are 6-30 chars, lowercase letters, digits and hyphens, starting with a letter."
  }
}

variable "project_name" {
  description = "Human-readable project name shown in the console."
  type        = string
  default     = "RBP composition confound"
}

variable "billing_account" {
  description = "Billing account ID to attach. Find it with: gcloud billing accounts list"
  type        = string
}

variable "region" {
  description = "Primary region. us-central1 has the widest GPU availability and lowest spot prices."
  type        = string
  default     = "us-central1"
}

variable "gpu_regions" {
  description = <<-EOT
    Every region that can actually run the sweep. Two independent things had to line up
    and only the second is visible in the quota console:

      1. PREEMPTIBLE_NVIDIA_V100_GPUS quota is 1 in a great many regions.
      2. V100 HARDWARE exists in only five of them. Checked with
         `gcloud compute accelerator-types list --filter=name=nvidia-tesla-v100`.

    Quota in a region with no V100s buys nothing, which is why the earlier "nine regions"
    plan was wrong. Five regions x 1 GPU = five concurrent workers, using 10 of the 12
    vCPU that CPUS_ALL_REGIONS allows globally.

    Each needs its OWN subnet: a subnetwork is a regional resource and a VM cannot attach
    one from another region.
  EOT
  type        = list(string)
  default     = ["us-central1", "us-east1", "us-west1", "europe-west4", "asia-east1"]
}

variable "bq_location" {
  description = "BigQuery location. Multi-region US, so queries are not tied to one zone."
  type        = string
  default     = "US"
}

variable "budget_usd" {
  description = <<-EOT
    Hard budget for alerting, in USD. The whole pipeline is estimated at ~$34 realistic and
    ~$51 with margin (docs/COST.md), so 100 leaves headroom while the first alert (25%) still
    fires at $25 -- before anything has gone badly wrong.
  EOT
  type        = number
  default     = 100
}

variable "killswitch_armed" {
  description = <<-EOT
    false: the billing killswitch logs what it would do and does nothing.
    true:  it detaches the project from its billing account when spend crosses its
           threshold, which stops every VM in the project.

    Starts false on purpose. The only other way to test a killswitch is to fire it, and
    firing it takes the project down to prove it can take the project down. Dry run proves
    Pub/Sub delivery, message parsing, permissions and the decision; arming it after that
    changes exactly one line of behaviour.
  EOT
  type        = bool
  default     = false
}
