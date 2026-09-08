# Minting the Zenodo DOI, step by step

A GitHub URL is not a citable reference: repositories get renamed, made private or deleted, and
the code behind a published number has to stay retrievable. Zenodo (run by CERN, free) takes a
snapshot of a GitHub release and mints a DOI that resolves to that exact snapshot forever. That
DOI is what the paper's Code availability section cites.

**Order matters.** Zenodo only archives releases created *after* you switch the repository on.
Creating the release first does nothing, which is the mistake almost everyone makes.

---

## 1. Connect Zenodo to GitHub (once)

1. Go to <https://zenodo.org>.
2. Click **Sign up**, then **Sign up with GitHub** (or **Log in with GitHub** if you have an
   account). Use the same GitHub account that owns the repository, `NirmalKumar31`.
3. Authorise Zenodo when GitHub asks. It requests read access to your repositories and
   permission to create webhooks.

## 2. Switch this repository on

1. In Zenodo, open the account menu (top right) and choose **GitHub**.
   Direct link: <https://zenodo.org/account/settings/github/>
2. Find **NirmalKumar31/rbp-protocol-calibration** in the list of repositories.
   If it is not listed, click **Sync now** at the top right and reload.
3. Flip its toggle to **On**.

That installs a webhook. Nothing is archived yet.

## 3. Add author metadata on Zenodo (do this now, not after)

Still on the Zenodo GitHub settings page, click the repository name, then **edit** the metadata
Zenodo will use. Set at least:

- **Authors:** Thirupallikrishnan Kesavan, Nirmalkumar. Add your ORCID if you have one; if you
  do not, get one first at <https://orcid.org/register>, it takes two minutes and it is what
  makes the deposit findable under your name.
- **Title:** the paper's title, so the software record and the preprint match.
- **Licence:** **both**, not MIT alone. Zenodo supports mixed-license uploads and documents
  the case directly: "files you upload are under different licenses (e.g. software under the
  MIT license, but documentation under CC-BY license). In such cases you can declare all the
  licenses". So record **MIT for the code** and **CC BY 4.0 for `results/` and
  `data/evidence/`**, which is what the repository's own `LICENSE`, `results/LICENSE` and
  `data/evidence/LICENSE` say. This file used to instruct MIT for the whole deposit on the
  ground that Zenodo asks for one licence, which flattened a distinction the paper states.
  See <https://help.zenodo.org/docs/deposit/describe-records/licenses/>.
- **Type:** Software.

If you skip this, Zenodo will guess the author from your GitHub profile name, and the deposit
will be attributed to a username rather than to you.

## 4. Cut the release on GitHub

In the repository on GitHub:

1. **Releases** (right-hand sidebar) then **Draft a new release**.
2. **Choose a tag** then type `v1.0.0` and select **Create new tag: v1.0.0 on publish**.
3. Target: `main`.
4. Release title: `v1.0.0` (or the paper title).
5. Description: one or two lines is enough. For example:
   > Code, committed evidence and manuscript accompanying the preprint. All 1156 verification
   > assertions pass on a clean clone.
6. **Publish release**.

Equivalently, from a terminal in the repository:

```
git tag -a v1.0.0 -m "Preprint release"
git push origin v1.0.0
gh release create v1.0.0 --title "v1.0.0" --notes "Code, committed evidence and manuscript accompanying the preprint."
```

## 5. Collect the DOI

Within a few minutes the repository appears under **Upload** on Zenodo with a DOI badge.
Reload <https://zenodo.org/account/settings/github/> if it does not appear immediately.

Zenodo gives you **two** DOIs and the difference matters:

| DOI | resolves to | use it for |
|---|---|---|
| **concept DOI** | always the latest version | give for discovery, alongside the version DOI |
| version DOI | exactly `v1.0.0` | cite when a specific version matters |

**CITE THE VERSION DOI FOR THE ANALYSED SNAPSHOT.** An audit was right that earlier guidance
here had this backwards. A concept DOI resolves to whatever the latest version happens to be, so
it cannot identify the exact research object the paper's numbers came from; that is the whole
property a reproducibility citation needs. Give both: the version DOI and its tag for the
snapshot the analysis ran on, and the concept DOI beside it so a reader finds later versions.

**THE TAG MUST MATCH THE DECLARED VERSION.** `pyproject.toml` and `CITATION.cff` are the source
of truth for the version string; tag whatever they say. They read 1.0.0 as of the release
commit, raised together, and the decision behind that is recorded in `CHANGELOG.md`: the
version identifies the archived snapshot and does not promise support, which `SECURITY.md`
still declines to offer. Both files must always move together, and a test enforces it.

The concept DOI is the one shown as "Cite all versions" on the record page. Give it in the
manuscript **beside** the version DOI, not instead of it. This paragraph used to say "use it in
the manuscript", which contradicts the table above it by nine lines and was one of three
different workflows the release described. The version DOI names the snapshot the numbers came
from and cannot be superseded; the concept DOI is how a reader finds later ones. Both, always.

**The circularity, and how to resolve it.** The DOI does not exist until the release is
archived, and the manuscript in that release should cite the DOI. There is no way around that
with a single release. Either reserve a DOI in a manual Zenodo draft before archiving, insert
it, and then archive; or archive once, insert the minted version DOI, and cut a second release
whose only change is that insertion. Say in the release notes which of the two was done, because
otherwise the first archive contains a paper citing nothing and looks like an oversight.

## 6. Wire it into the manuscript

Two places, then one rebuild:

- `manuscript/paper.tex`, the Code availability section.
- `manuscript/sections/data-availability.tex`, if the DOI is mentioned there too.

```
cd manuscript && ./build.sh
```

`build.sh` fails on an undefined reference, so a broken edit will not silently produce a PDF.

## 7. The wrinkle, and it is resolved in section 5 rather than here

An earlier version of this section said flatly that the archived snapshot **cannot** contain its
own DOI. That contradicted section 5 above, which correctly offers DOI reservation in a draft as
one of two workflows, and a reader following this file could reasonably have done either.

**Pick one, and section 5's first option is the recommended one:** reserve the DOI in a manual
Zenodo draft, insert it into the manuscript and `CITATION.cff`, rebuild, commit, tag, and then
archive that tag. Zenodo documents reservation directly, and the reserved DOI survives as long
as the draft does. The snapshot then contains its own version DOI and there is no second
release.

The fallback remains valid if reservation is not used: archive once, insert the minted version
DOI, and cut a second release whose only change is that insertion. Whichever is chosen, say so
in the release notes.

## 8. Changing the repo after the DOI exists

Nothing here is one-shot. Zenodo's GitHub integration mints a **new version** each time you cut
a new GitHub release, and every version sits under the same concept DOI. So the loop is:

```
commit -> git tag v1.0.1 -> push the tag -> cut a GitHub release -> Zenodo archives it
```

The concept DOI printed in the paper keeps resolving, now to the newer snapshot. You do not
need to touch the manuscript again.

Two limits. **You cannot swap the files of a version already published**: file changes require a
new version, though metadata (title, description, authors) is editable in place on an existing
record. And **you cannot delete a published record yourself**; withdrawal is a support request
to Zenodo, and the DOI stays registered as a tombstone. So the thing to get right first time is
the *authorship and licence metadata*, not the code, because the code you can always supersede.

## What this is not

Zenodo archives the **code and data**. The **preprint** goes to bioRxiv separately and gets its
own DOI from them. The two cite each other: the paper's Code availability section gives the
Zenodo DOI, and the Zenodo record's metadata can carry the bioRxiv DOI as a related identifier
once you have it.
