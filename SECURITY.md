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

## The history scan, and what it found

Run 2026-09-05 over all 202 commits on every ref (`git log --all -p`), looking for private-key
blocks, service-account JSON fields, AWS keys, GitHub and Slack tokens, and billing account IDs.

**Clean:** no private key, service-account key material, or API token appears anywhere in the
history.

**One finding, and it is real.** A live GCP **billing account ID** appears in 11 places in the
history. It was scrubbed from the working tree by commit `f3fab95` ("Submission packaging:
scrub a live billing ID") and `tests/unit/test_no_hardcoded_project.py` has forbidden it in
tracked files ever since, so it is absent from every current file. Git history is not the
working tree: anyone who clones this public repository can recover it.

What that does and does not mean. A billing account ID is not a credential and cannot be used
to authenticate; it is an identifier. It is credential-adjacent, which is why this repository
forbids committing one: it names a real account for anyone constructing a targeted request, and
it appears in support and console URLs.

The options, in order of cost:

1. Accept it, having judged the exposure acceptable for an identifier that grants nothing.
2. Rewrite the history (`git filter-repo`), force-push, and ask collaborators to re-clone. This
   removes it from the canonical repository. Forks and any existing clone keep it.
3. Move the work to a fresh billing account, which is the only action that makes the leaked
   identifier refer to nothing.

**This has not been done, because rewriting published history and force-pushing a public
repository is the repository owner's decision and not an automatic remediation.** It is recorded
here rather than quietly fixed or quietly ignored.

## What is not guarded

- **Forks and existing clones** keep whatever the history contained at the time they were made,
  whatever is done to this repository afterwards.
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
