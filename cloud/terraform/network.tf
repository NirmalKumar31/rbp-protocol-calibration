# Network for the compute workloads.
#
# WHY THIS FILE EXISTS AT ALL. Everything ran on the auto-created "default" VPC until a
# 488-task job would only start 3 of its 7 nodes. The cause was IN_USE_ADDRESSES, quota 4:
# a Batch VM takes an EXTERNAL IP by default, and the fifth one cannot be created. Batch
# does not fail on that -- it logs the rejections as OPERATIONAL_INFO, keeps the job
# RUNNING, and quietly delivers a fraction of the requested throughput.
#
# The fix is not "ask for more addresses". These workers have no business having public
# addresses: nothing connects TO them, and everything they connect to (Cloud Storage,
# Artifact Registry, Logging, the Batch control plane) is a Google API. So they get no
# external IP, and the subnet gets Private Google Access, which routes traffic to Google
# API ranges over Google's internal network instead of the public internet.
#
# Consequences, all in the right direction:
#   * address quota stops being a ceiling on parallelism entirely
#   * no public attack surface on the workers
#   * no egress charges for reaching Google APIs
#
# WHY A DEDICATED VPC RATHER THAN TURNING PGA ON FOR THE DEFAULT ONE. The default network
# is created by GCP, not by us, so Terraform does not own it. Flipping a setting on it by
# hand would be exactly the kind of out-of-band change that makes the step 9 destroy-and-
# rebuild test lie: a fresh project would come up with PGA off again and nobody would know
# until throughput was wrong. A VPC we declare is a VPC that gets recreated correctly.

resource "google_compute_network" "main" {
  project = google_project.rbp.project_id
  name    = "rbp-net"

  # Without this, GCP creates a subnet in EVERY region -- about 40 of them, all unused.
  # We want exactly one, in the region we actually run in.
  auto_create_subnetworks = false

  depends_on = [google_project_service.apis]
}

resource "google_compute_subnetwork" "workers" {
  project       = google_project.rbp.project_id
  name          = "rbp-workers"
  region        = var.region
  network       = google_compute_network.main.id
  ip_cidr_range = "10.10.0.0/20" # 4,094 usable addresses, far more than 8 instances need

  # THE LINE THIS FILE IS FOR. Lets a VM with no external IP reach *.googleapis.com and
  # *.pkg.dev over Google's internal network. Without it, a worker with no external IP
  # cannot pull its container image and the task hangs rather than erroring.
  private_ip_google_access = true

  log_config {
    # Sampled flow logs. Off entirely would leave no way to answer "did this VM actually
    # try to reach anything?", which is the first question when a private-IP job hangs.
    # 10% sampling keeps the cost negligible on a workload with no steady traffic.
    aggregation_interval = "INTERVAL_10_MIN"
    flow_sampling        = 0.1
    metadata             = "INCLUDE_ALL_METADATA"
  }
}

# NO INGRESS RULES ARE DEFINED, AND THAT IS DELIBERATE.
#
# A custom VPC starts with an implied deny-all for ingress and an implied allow-all for
# egress. Our workers need to make outbound calls and accept nothing, so the defaults are
# already exactly right. The default VPC ships with allow-ssh, allow-rdp and allow-icmp
# from anywhere; none of those are wanted here, and not copying them across is the point.
#
# Debugging is done through Cloud Logging, which is how every failure so far was diagnosed,
# so there is no need for SSH. If it is ever needed, add an IAP-ranged rule (35.235.240.0/20)
# rather than opening 0.0.0.0/0.

# There is deliberately NO Cloud NAT either.
#
# NAT would let these VMs reach the public internet, and they do not need to. The one job
# that does -- cloud_ingest.py, which downloads from ENCODE, GENCODE and NCBI -- keeps an
# external IP, because it is a single VM that runs rarely. A NAT gateway is billed per
# VM-hour attached plus per GB processed, so adding one to serve an occasional download
# would cost more than the download and would give the other workers internet access they
# have no reason to have.

# ---------------------------------------------------------------------------------------
# GPU subnets, one per region that can actually run the sweep.
#
# WHY THIS IS NOT ONE SUBNET. A VPC is global; a SUBNET is regional. `rbp-workers` above
# lives in us-central1, and a VM in us-west1 cannot attach it -- Batch rejects the job at
# submit time with an invalid-argument error naming the subnetwork. Preprocessing never hit
# this because it is CPU work and CPU quota is regional, so everything ran in one place.
# The GPU sweep cannot: V100 quota is one per region, so five regions is the ONLY way to
# get five GPUs, and five regions means five subnets.
#
# Same posture as rbp-workers: Private Google Access on, no external IPs, no NAT, no
# ingress. The workers reach Cloud Storage, Artifact Registry, Logging and the Batch
# control plane, all of which are Google APIs, and nothing else. Model weights are baked
# into the GPU image precisely so that "nothing else" stays true -- see docker/Dockerfile.gpu.
#
# /20 each out of 10.20.0.0/16, non-overlapping. They do not strictly need to be disjoint
# across regions, but keeping them so means a flow log or firewall rule can be written
# against a single range later without ambiguity.
resource "google_compute_subnetwork" "gpu" {
  for_each = { for i, r in var.gpu_regions : r => i }

  project       = google_project.rbp.project_id
  name          = "rbp-gpu-${each.key}"
  region        = each.key
  network       = google_compute_network.main.id
  ip_cidr_range = cidrsubnet("10.20.0.0/16", 4, each.value)

  private_ip_google_access = true

  log_config {
    aggregation_interval = "INTERVAL_10_MIN"
    flow_sampling        = 0.1
    metadata             = "INCLUDE_ALL_METADATA"
  }
}

output "network" {
  value       = google_compute_network.main.name
  description = "VPC for compute workloads"
}

output "gpu_subnetworks" {
  value       = { for r, s in google_compute_subnetwork.gpu : r => s.name }
  description = "region -> subnet name, for the sweep submitter"
}

output "subnetwork" {
  value       = google_compute_subnetwork.workers.name
  description = "Regional subnet with Private Google Access enabled"
}
