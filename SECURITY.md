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
commits at the time of writing"; by the time anyone read it the repository was larger, a fresh
count gave a third number, and there was no way to tell which described the release. The
scan is now a gate, and `--check` fails if any FINDING changes, which is the property that
matters rather than the size of the history.

**Clean:** no private key, service-account key material, or API token appears anywhere in the
history, and the scan exits non-zero if that ever stops being true.

**One finding, and it is real.** A live GCP **billing account ID** is in the history. Three
counts, because this document previously gave one of them without saying which: it was
introduced in **4 distinct commits**, touches **9 commit-and-file pairs**, and matches **11 diff
lines** counting both sides. These totals include the archived working-notes history preserved
by tag `archive/working-notes-2026-09-06`; the published `main` branch contains three of those
commits. Nothing was scrubbed to reduce the totals when the two working branches were removed.

The ID was scrubbed from the working tree by commit `f3fab95` ("Submission packaging:
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

**DECISION FOR v1.0.0: option 1, accept and record.** The owner judged the exposure
acceptable, and the reasoning is recorded here rather than left implicit.

A billing account ID authenticates nothing. It names an account, which is why this repository
forbids committing one, but possessing it grants no access. Set against that, this repository
has been public for weeks, so its history is already in clones, in forks and in GitHub's own
API. Rewriting history would remove the string from the canonical branch and from nowhere else,
while breaking every existing clone and invalidating every commit SHA any reader has cited,
including the ones this paper cites. That is a real cost paid for an incomplete remedy.

Option 3, moving to a fresh billing account, is the only action that would make the identifier
refer to nothing. It is not taken because the account is still in use and the identifier grants
nothing to begin with.

The finding stays gated rather than closed: `scripts/history_scan.py --check` still fails if the
count changes, so accepting the exposure does not stop it being reported. It is visible in
`results/tables/history_scan.csv` in the release snapshot, which is the point of recording a
decision instead of quietly resolving one.

## AI co-author trailers in the commit history

Most commits up to `202953d` carry a `Co-Authored-By: Claude` trailer; the count is in
`results/tables/history_scan.csv`, generated rather than typed, because a number written here is
wrong by the next commit.

The trailer was dropped after `202953d` and then reinstated, so a stretch of commits in between
carries none. Both the gap and its ends are deliberate and neither is a change of disclosure:
the manuscript's Use of AI tools section is the disclosure of record throughout, and no history
has been rewritten. It is written down here because a signal that stops and restarts without
explanation reads worse than either consistent choice.

## What is not guarded

- **Forks and existing clones** keep whatever the history contained at the time they were made,
  whatever is done to this repository afterwards.
- **`roles/billing.admin` at the billing-account scope** is broader than the one action the
  killswitch needs. GCP offers no narrower role that permits `projects.updateBillingInfo`. If
  you deploy this, scope it to a billing account you are willing to expose.
- **Base images are pinned by tag, not digest.** A rebuild is not bit-reproducible, and the
  digests of the images that actually produced the published run were not recorded anywhere.
  `cloud/terraform/outputs.tf` says a digest is what should be recorded for anything whose
  output is published, and that rule was not followed for this run.
- **Dependencies have lower bounds, not a lock file**, outside `docker/requirements-*.txt`,
  which are pinned. `pip install -e .` can resolve a future incompatible version.

## Scope of support

The published results are frozen. The code is provided so they can be checked and extended;
it is not maintained as a library, and there is no supported API surface.
