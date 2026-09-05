# Security

This is a research repository, not a service. There is nothing here that accepts input from a
network, and no deployed endpoint. The realistic problems are therefore about credentials and
about cost, not about remote code execution.

## Reporting

Open an issue at
<https://github.com/NirmalKumar31/rbp-protocol-calibration/issues>, or email
thirupallikrishnan.n@northeastern.edu if the issue is a leaked credential and should not be
public first.

## What is guarded, and by what

| risk | control |
|---|---|
| a project or bucket id committed in source | `tests/unit/test_no_hardcoded_project.py`, matching the id's shape rather than one historical name, over every tracked `.py`, `.sh`, `.yaml` and `.json` |
| a billing account id committed anywhere | the same test, over every tracked file including docs |
| a private key in Terraform state | no `google_service_account_key` resource exists; the one key is minted out of band. See `docs/REPRODUCE.md` for its lifecycle and revocation |
| runaway GPU spend | `cloud/modal/guard.py`, which discovers every sweep app by pattern and **fails closed** if it cannot read the Modal CLI |
| runaway GCP spend | `cloud/killswitch/`, which detaches billing. It reports at startup whether it actually holds the permission to do so, because a dry run proves only the read |
| a cost report that cannot distinguish zero from unobservable | `cloud/cost.sh` counts and names failed queries and exits non-zero |

## What is not guarded

- **No credential-history audit has been run.** A pattern scan over the current tree found
  nothing; that is not the same as the history being clean.
- **`roles/billing.admin` at the billing-account scope** is broader than the one action the
  killswitch needs. GCP offers no narrower role that permits `projects.updateBillingInfo`. If
  you deploy this, scope it to a billing account you are willing to expose.
- **Base images are pinned by tag, not digest.** A rebuild is not bit-reproducible. The digests
  of the images that produced the published run are recorded in `docs/COST.md` terms only, not
  in a manifest.
- **Dependencies have lower bounds, not a lock file**, outside `docker/requirements-*.txt`,
  which are pinned. `pip install -e .` can resolve a future incompatible version.

## Scope of support

The published results are frozen. The code is provided so they can be checked and extended;
it is not maintained as a library, and there is no supported API surface.
