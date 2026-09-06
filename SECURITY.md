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
| runaway GPU spend | `cloud/modal/guard.py`, which discovers every sweep app by pattern, sums the bound over all of them, and on losing sight of the Modal CLI keeps charging every app it last saw alive rather than the one it started with |
| runaway GCP spend | `cloud/killswitch/`, which detaches billing. It reports at startup whether it actually holds the permission to do so, because a dry run proves only the read |
| a cost report that cannot distinguish zero from unobservable | `cloud/cost.sh` counts and names failed queries and exits non-zero |

## The history scan, and what it found

`scripts/history_scan.py` runs it, over every commit on every ref (`git log --all -p`), looking
for private-key blocks, service-account JSON fields, AWS keys, GitHub, Slack and Google API
tokens, OAuth tokens, and billing account IDs. **The counts are in
`results/tables/history_scan.csv`, generated, not typed here.** They were typed here, as "206
commits at the time of writing"; by the time an external audit read it the repository was larger,
the audit counted a third number, and there was no way to tell which described the release. The
scan is now a gate, and `--check` fails if any FINDING changes, which is the property that
matters rather than the size of the history.

**Clean:** no private key, service-account key material, or API token appears anywhere in the
history, and the scan exits non-zero if that ever stops being true.

**One finding, and it is real.** A live GCP **billing account ID** is in the history. Three
counts, because this document previously gave one of them without saying which: it was
introduced in **4 distinct commits**, touches **9 commit-and-file pairs**, and matches **11 diff
lines** counting both sides, the last being the "11 places" this section used to quote. It was
scrubbed from the working tree by commit `f3fab95` ("Submission packaging:
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

**Status as of the current release candidate: deliberately deferred to the release phase, not
resolved and not dropped.** The scan is a gate now, so the finding cannot fade out of the record
between here and the tag: `scripts/history_scan.py --check` fails if it changes, and it will
still be reported on the commit that is archived. Whoever cuts that release makes the call from
options 1 to 3 above and records which, on the release commit, before the DOI is minted.

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
