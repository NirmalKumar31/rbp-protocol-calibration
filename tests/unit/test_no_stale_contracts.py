"""Present-tense claims in prose and comments must agree with the released evidence.

WHY. An external audit found six places where the repository asserts, in the present tense,
something a committed table or another file in the same repository contradicts. None of them
changed a number. All of them would be found by a referee, and each one costs more credibility
than the defect it describes, because a reader cannot tell a stale comment from a live one.

The shape recurs: a disclosure is written carefully, the defect is then fixed, and the
disclosure stays. `deep_model_contrast.py` disclosed fold leakage in 20 datasets and ended
with "do not restore the same-folds sentence without rerunning the sweep". The sweep was rerun
in 2276eea, `fold_integrity.csv` went from 20 to 0, and the disclosure sat there for three
release candidates telling the reader that a fifth of one arm was contaminated.

Each check below pins one claim to the artefact that decides it, so the claim cannot outlive
the state it describes. These are string checks, which is crude, but the failure they catch is
also crude and nothing subtler was catching it.
"""

import csv
import re
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[2]
SECTIONS = ROOT / "manuscript" / "sections"


def _read(rel):
    p = ROOT / rel
    if not p.exists():
        pytest.skip(f"{rel} not in this checkout")
    return p.read_text()


def test_the_discussion_does_not_deny_an_experiment_the_results_report():
    """P0.1. Both statements were in the same PDF, 850 lines apart."""
    results = _read("manuscript/sections/results.tex")
    disc = _read("manuscript/sections/discussion.tex")
    if "\\label{sec:transport}" not in results:
        pytest.skip("no transport section to be consistent with")
    assert "performs that transfer and we do not" not in disc, (
        "results.tex reports the full train-arm by evaluation-arm experiment under "
        "\\label{sec:transport}, and discussion.tex says the paper did not perform it")
    assert "\\ref{sec:transport}" in disc, (
        "the Discussion's scope paragraph must point at the transport section rather than "
        "describe the design as if it were absent")


def test_the_model_was_refit_per_arm_so_nothing_says_the_model_was_held_fixed():
    """P0.2. The model class was fixed. A separate model was fitted in every arm."""
    bad = re.compile(r"[Hh]olding the model,|\bModel, peak set\b")
    for rel in ("manuscript/sections/introduction.tex", "manuscript/sections/discussion.tex",
                "CITATION.cff"):
        t = _read(rel)
        assert not bad.search(t), (
            f"{rel} says the MODEL was held fixed across arms. It was refit in each arm; what "
            "was held fixed is the model class and its hyperparameters")


def test_the_gc_and_dinucleotide_positives_are_not_called_identical():
    """P0.2. Only 10 of 94 datasets have identical positive sets."""
    t = _read("manuscript/sections/results.tex")
    assert "identical positives" not in t, (
        "the GC and dinucleotide arms retain near-identical, not identical, positives: each "
        "arm keeps only the positives its own matching could pair")
    overlap = ROOT / "results" / "tables" / "positive_set_overlap.csv"
    if not overlap.exists():
        return
    j = [float(r["jaccard"]) for r in csv.DictReader(overlap.open())]
    assert min(j) < 1.0, (
        "every dataset now has identical GC and dinucleotide positives, so the caption may "
        "say identical again and this test should be retired")
    assert sum(x == 1.0 for x in j) < len(j), "same"
    # The caption quotes these three. Pin them so a regenerated table cannot silently
    # contradict the text.
    assert abs(sorted(j)[len(j) // 2] - 0.9972) < 5e-4, "median Jaccard moved"
    assert abs(min(j) - 0.9164) < 5e-4, "minimum Jaccard moved"
    assert sum(x == 1.0 for x in j) == 10, "the count of identical datasets moved"


def test_the_fold_leakage_disclosure_matches_the_fold_integrity_table():
    """P1.13. The disclosure outlived the defect by three release candidates."""
    fi = ROOT / "results" / "tables" / "fold_integrity.csv"
    if not fi.exists():
        pytest.skip("fold_integrity.csv not in this checkout")
    rows = {r["check"]: float(r["value"]) for r in csv.DictReader(fi.open())}
    dirty = rows.get("datasets NOT chromosome-grouped, dn arm")
    src = _read("scripts/deep_model_contrast.py")
    present = re.search(r"For \*\*20 of the 94 dinucleotide-arm\s*\n?datasets\*\* the\s*\n"
                        r"committed", src)
    if dirty == 0:
        assert present is None, (
            "fold_integrity.csv reports 0 non-chromosome-grouped dn-arm datasets, and "
            "deep_model_contrast.py still asserts in the present tense that 20 of them carry "
            "scores from a stratified random partition")
    else:
        assert present is not None, (
            f"fold_integrity.csv reports {dirty} non-chromosome-grouped dn-arm datasets and "
            "deep_model_contrast.py no longer discloses it")


def test_the_baseline_is_not_called_the_whole_story():
    """P1.13. The same script prints that a protocol-specific residual exists."""
    src = _read("scripts/protocol_or_baseline.py")
    assert "A protocol-specific residual DOES exist" in src, (
        "this test pins the docstring to that finding; if the finding moved, move the test")
    for phrase in ("it is the whole story", "carries essentially no information beyond it"):
        assert phrase not in src, (
            f"protocol_or_baseline.py's header says {phrase!r} while the script itself reports "
            "a residual whose interval excludes zero under two designs")


def test_the_killswitch_does_not_promise_that_nothing_is_lost():
    """P1.11. Terraform promised what the function's own docstring retracts."""
    tf = _read("cloud/terraform/killswitch.tf")
    assert "never lost data" not in tf, (
        "killswitch.tf guarantees no data loss on billing detach. Google documents that "
        "disabling billing may delete resources non-recoverably, and cloud/killswitch/main.py "
        "says so. An unsafe guarantee is what persuades an operator to skip the backup")


def test_the_venue_is_not_described_as_arxiv():
    """P1.12. The target has been bioRxiv since the venue was chosen."""
    assert "arXiv PDF" not in _read("manuscript/build.sh")


def test_exit_zero_is_not_described_as_raw_reproduction():
    """P1.7/P0.4. Most released tables are verified, not reconstructed."""
    t = _read("docs/REPRODUCE.md")
    assert "Exit 0 means the science reproduced." not in t, (
        "run.sh regression-verifies most tables against committed evidence rather than "
        "rebuilding them from raw inputs; PROVENANCE.csv says which is which")


def test_the_image_gate_resolves_an_absolute_interpreter():
    """P0.4. `env PATH=/nonexistent python3` can never find python3.

    run.sh gates the paid image build on this script, so a gate that cannot start is a rebuild
    path that is documented and does not run. It failed for everyone, every time, with
    `env: python3: No such file or directory`, and nothing noticed because nothing ran it.
    """
    t = _read("scripts/check_image_tree.sh")
    assert 'PY="${PY:-python3}"' not in t, (
        "PY must be resolved to an absolute path with command -v BEFORE PATH is emptied")
    assert "command -v" in t, "the interpreter is not being resolved at all"


def test_the_image_gate_runs_both_build_selections():
    """P0.4. One simulation cannot stand in for two different images."""
    t = _read("scripts/check_image_tree.sh")
    assert "cloudbuild.cpu.yaml" in t, (
        "the CPU build's --ignore list must be read from the build file, not copied here, "
        "or the simulation drifts from the thing it simulates")
    cpu = _read("docker/cloudbuild.cpu.yaml")
    ignores = re.findall(r"--ignore=tests/[A-Za-z0-9_/.]+", cpu)
    assert ignores, (
        "the parser in check_image_tree.sh finds the CPU ignore list with this pattern; if "
        "cloudbuild.cpu.yaml stops matching it, the gate fails closed but for the wrong reason")


def test_the_neural_stack_is_an_extra_not_a_base_dependency():
    """P0.4. `pip install -e . -c constraints.txt` used to install torch regardless."""
    tomllib = pytest.importorskip("tomllib")
    d = tomllib.loads(_read("pyproject.toml"))
    base = " ".join(d["project"]["dependencies"])
    for pkg in ("torch", "transformers", "multimolecule", "peft"):
        assert pkg not in base, (
            f"{pkg} is a base dependency, so the offline verification install is not minimal. "
            "constraints.txt has always claimed it is")
    assert "neural" in d["project"]["optional-dependencies"]
    rp = d["project"]["requires-python"]
    assert "<" in rp, f"requires-python {rp!r} has no upper bound"


# Entries that legitimately have no DOI. Both are the canonical software citations their
# projects ask to be cited by, and neither publisher minted one.
NO_DOI_BY_DESIGN = {"sklearn", "pytorch"}


def test_every_bibliography_entry_carries_a_resolvable_identifier():
    """P1-1. Five of the seven surveyed methods were discussed with no citation at all.

    Checked OFFLINE, on purpose: a CI job that queries Crossref fails when Crossref is slow,
    and a gate that fails for a reason nobody should act on gets disabled. What this asserts is
    that every entry carries a DOI or is a named exception, so a new entry cannot be added
    without a resolvable identifier and nobody noticing.

    All 29 DOIs were resolved against Crossref by hand on 2026-09-07 and every normalised title
    matched at 0.92 or better. That is a dated manual check, recorded here, and it is not what
    this test does.
    """
    if not (ROOT / "manuscript" / "sections" / "bibliography.tex").exists():
        pytest.skip("no manuscript/ here: this is the container file set")
    text = _read("manuscript/sections/bibliography.tex")
    entries = re.split(r"\\bibitem", text)[1:]
    assert len(entries) >= 25, (
        f"only {len(entries)} bibitems parsed; the pattern stopped matching rather than the "
        "bibliography having shrunk")
    missing = []
    for e in entries:
        key = re.search(r"\]\{([^}]+)\}", e)
        assert key, "a bibitem has no citation key"
        if r"\doi{" not in e and key.group(1) not in NO_DOI_BY_DESIGN:
            missing.append(key.group(1))
    assert not missing, (
        f"bibliography entries with neither a DOI nor an exemption: {missing}. Add the DOI, or "
        "add the key to NO_DOI_BY_DESIGN with the reason")


def test_every_survey_row_is_cited():
    """The seven-method table discussed methods the bibliography did not contain."""
    results = _read("manuscript/sections/results.tex")
    block = re.search(r"\\label\{tab:survey\}(.*?)\\end\{tabular\}", results, re.S)
    assert block, "the survey table no longer matches; the citation check covers nothing"
    rows = [ln for ln in block.group(1).splitlines()
            if "&" in ln and "\\toprule" not in ln and "source &" not in ln]
    assert len(rows) >= 7, f"only {len(rows)} survey rows found"
    uncited = [ln.split("&")[0].strip() for ln in rows if r"\citep{" not in ln]
    assert not uncited, f"survey rows with no citation: {uncited}"


def test_an_escaped_dollar_is_not_a_maths_delimiter():
    """release_consistency.abstract_words() read `\\$115` as an opening maths delimiter.

    It then closed the span at the next real `$` and swallowed 165 words of prose as one token,
    so the abstract counted 300 instead of 465 and the release check reported a stale claim in
    SUBMISSION.md that was not stale. A measurement tool returning a confident wrong number is
    worse than one that fails, because the wrong number gets acted on.
    """
    import sys
    if not (ROOT / "manuscript" / "paper.tex").exists():
        pytest.skip("no manuscript/ here: the container file set carries src, scripts, config "
                    "and tests only, and this test round-trips through paper.tex")
    sys.path.insert(0, str(ROOT / "scripts"))
    rc = pytest.importorskip("release_consistency")

    plain = r"\begin{abstract} one two three four five. \end{abstract}"
    with_cost = r"\begin{abstract} one two costs \$115 three four five. \end{abstract}"
    with_maths = r"\begin{abstract} one two $x + y$ three four five. \end{abstract}"

    def count(tex, tmp=ROOT / "manuscript" / "paper.tex"):
        keep = tmp.read_text()
        try:
            tmp.write_text(tex)
            return rc.abstract_words()
        finally:
            tmp.write_text(keep)

    assert count(plain) == 5
    # "costs", "$", "115" and the five: the dollar amount adds words, it does not delete a span
    assert count(with_cost) >= 7, "an escaped dollar swallowed the rest of the abstract"
    assert count(with_maths) == 6, "a real maths span should count as one token"


def test_the_declared_version_is_one_number_everywhere():
    """P0-5. pyproject.toml, CITATION.cff and any git tag must agree, or the release has no
    single identity and a version DOI cannot point at a definite snapshot.

    They agree today at 0.9.0 and nothing enforced it. The tag does not exist yet; when it
    does, this fails unless it matches, which is the point: the tag is the thing a DOI is
    minted against.
    """
    import subprocess

    tomllib = pytest.importorskip("tomllib")
    pyproject = tomllib.loads(_read("pyproject.toml"))["project"]["version"]
    cff = None
    for line in _read("CITATION.cff").splitlines():
        if line.startswith("version:"):
            cff = line.split(":", 1)[1].strip().strip('"\'')
    assert cff is not None, "CITATION.cff has no version: key"
    assert pyproject == cff, (
        f"pyproject.toml says {pyproject}, CITATION.cff says {cff}")

    r = subprocess.run(["git", "tag", "--points-at", "HEAD"],
                       capture_output=True, text=True, cwd=str(ROOT))
    if r.returncode != 0:
        return                      # not a git repo: an unpacked archive, nothing to check
    tags = [t for t in r.stdout.split() if t.startswith("v")]
    for t in tags:
        assert t.lstrip("v") == pyproject, (
            f"tag {t} does not match the declared version {pyproject}. A version DOI minted "
            "against this tag would name a different release from the one the code declares")


def test_every_repository_path_the_manuscript_names_exists():
    """A paper that tells a reader to look at a file must name a file that is there.

    Added because writing the Methods paragraph for the negative-set survey, I pointed at
    `results/tables/negative_set_survey_sources.csv`, which does not exist; the committed table
    is `negative_set_survey_per_method.csv`. Nothing would have caught it. `verify.py` checks
    values in tables, `audit_manuscript.py` checks numbers, and neither reads a path.

    A reader who follows a dead path in the Methods concludes the evidence is missing.
    """
    if not (ROOT / "manuscript").exists():
        pytest.skip("no manuscript/ here: this is the container file set")
    paths = []
    for f in sorted(ROOT.glob("manuscript/*.tex")) + sorted(ROOT.glob("manuscript/sections/*.tex")):
        for m in re.finditer(
                r"\\texttt\{((?:scripts|src|tests|config|docs|cloud|docker|results|data|"
                r"manuscript)/[A-Za-z0-9_.\\/-]+)\}", f.read_text()):
            rel = m.group(1).replace("\\_", "_").replace("\\", "")
            paths.append((f.relative_to(ROOT), rel))

    assert len(paths) >= 10, (
        f"only {len(paths)} repository paths found in the manuscript; the pattern stopped "
        "matching rather than the paper having stopped citing files")
    missing = sorted({f"{src}: {rel}" for src, rel in paths if not (ROOT / rel).exists()})
    assert not missing, "the manuscript names paths that do not exist:\n  " + "\n  ".join(missing)


def test_every_licence_cross_reference_resolves_from_its_own_directory():
    """P1-2. `data/evidence/LICENSE` pointed at `../LICENSE`, which is `data/LICENSE`.

    The sentence is correct in `results/LICENSE`, where `../LICENSE` IS the repository root, and
    it was copied verbatim one level deeper. So the CC BY 4.0 file that tells a reader where to
    find the code licence pointed at nothing, in the directory whose licensing the paper makes a
    point of separating. Relative paths are only correct relative to something, and nothing here
    had ever resolved one.
    """
    root = ROOT
    if not (root / "LICENSE").exists():
        pytest.skip("no LICENSE here: this is the container file set")
    seen = 0
    for rel in ("LICENSE", "results/LICENSE", "data/evidence/LICENSE"):
        f = root / rel
        if not f.exists():
            continue
        # NOTICE as well as LICENSE. The third-party source list moved out of the root
        # LICENSE into NOTICE, because appending it made GitHub report the repository as
        # NOASSERTION, and the two CC BY files were repointed at it. A pattern that only
        # knew about LICENSE would have stopped covering the reference it was repointed to.
        for ref in re.findall(r"\.\./[./]*(?:LICENSE|NOTICE)", f.read_text()):
            seen += 1
            target = (f.parent / ref).resolve()
            assert target.exists(), f"{rel} points at {ref}, which resolves to a missing {target}"
    assert seen >= 4, (
        f"only {seen} licence cross-references found; the pattern stopped matching rather than "
        "the files having stopped cross-referencing")


def test_the_readers_pins_and_the_ci_pins_are_the_same_environment():
    """Two hand-kept pin files, and nothing compared them.

    README tells a reader `pip install -e . -c constraints.txt`. CI installs
    docker/requirements-cpu.txt instead, so the environment every published number is
    verified in is NOT the one the README hands out, and the only thing that had ever made
    the two agree was that the same person edited both. They do agree today, exactly, on
    all twelve shared packages. A single bump to one file would end that silently, and the
    symptom would be a reader getting a numeric difference nobody can reproduce.

    requirements-cpu.txt is allowed to pin MORE, because CI reads the rendered PDF and a
    reader verifying tables does not need pypdf. It is not allowed to pin the same package
    differently.
    """
    def pins(text):
        out = {}
        for line in text.splitlines():
            m = re.match(r"^([A-Za-z0-9_.-]+)==([^\s;#]+)", line.split("#")[0].strip())
            if m:
                out[m.group(1).lower().replace("_", "-")] = m.group(2)
        return out

    reader = pins(_read("constraints.txt"))
    ci = pins(_read("docker/requirements-cpu.txt"))
    assert reader, "constraints.txt stopped parsing as pins"
    assert ci, "docker/requirements-cpu.txt stopped parsing as pins"
    disagree = {k: (reader[k], ci[k]) for k in reader.keys() & ci.keys() if reader[k] != ci[k]}
    assert not disagree, f"the reader's environment differs from the verified one: {disagree}"
    missing = sorted(reader.keys() - ci.keys())
    assert not missing, (
        f"constraints.txt pins {missing}, which CI never installs, so those pins are "
        "unexercised")


def test_the_column_dictionary_is_a_dictionary_and_is_current():
    """P0.3. It carried four fields, no definitions, and a dtype from the first row alone."""
    cols = ROOT / "results" / "tables" / "COLUMNS.csv"
    if not cols.exists():
        pytest.skip("COLUMNS.csv not in this checkout")
    rows = list(csv.DictReader(cols.open()))
    for field in ("dtype", "unit", "key", "n_rows", "n_missing", "n_distinct", "min", "max",
                  "definition", "producing_script"):
        assert field in rows[0], f"COLUMNS.csv has no {field} column"
    # THE DTYPE BUG. A column typed empty must actually be empty everywhere, which is what
    # reading only the first data row failed to establish for 87 of the 98 it so typed.
    tables = ROOT / "results" / "tables"
    for r in rows:
        if r["dtype"] != "empty":
            continue
        f = tables / r["table"]
        if not f.exists():
            continue
        d = list(csv.DictReader(f.open(newline=""),
                                delimiter="\t" if f.suffix == ".tsv" else ","))
        assert not any((x.get(r["column"]) or "").strip() for x in d), (
            f"{r['table']}:{r['column']} is typed empty but has values further down the file")
    assert int(rows[0]["n_rows"]) >= 0


def test_schema_md_does_not_hand_maintain_counts():
    """P0.3. It carried three drifted figures in the paragraph saying counts are generated."""
    t = _read("results/tables/SCHEMA.md")
    # Ignore the fenced/inline-code and table-cell history; what must not appear is a
    # count stated as fact in the running prose.
    prose = "\n".join(ln for ln in t.splitlines() if not ln.startswith("|"))
    bad = re.findall(r"\b(\d{2,4}) (?:tables|columns|distinct schemas)\b", prose)
    assert not bad, (
        f"SCHEMA.md hand-maintains counts {bad}; they belong in COLUMNS_SUMMARY.csv, which is "
        "generated by scripts/column_dictionary.py")


def test_the_raw_input_manifest_is_complete_and_carries_no_secrets():
    """P1.8. Accessions identify resources; only checksums identify bytes."""
    f = ROOT / "results" / "tables" / "raw_inputs.csv"
    if not f.exists():
        pytest.skip("raw_inputs.csv not in this checkout")
    rows = list(csv.DictReader(f.open()))
    assert len(rows) > 200, f"only {len(rows)} raw inputs; the panel alone needs 94 peak files"
    for r in rows:
        assert r["md5_base64"], f"{r['path']} has no checksum"
        assert r["size_bytes"].isdigit() and int(r["size_bytes"]) > 0
        assert not re.match(r"^/|\.\.|~|\\\\", r["path"]), f"unsafe path: {r['path']}"
    blob = f.read_text()
    for pat in (r"\b\d{6}-[0-9A-F]{6}-[0-9A-F]{6}\b", r"BEGIN [A-Z ]*PRIVATE KEY",
                r"AIza[0-9A-Za-z_-]{35}", r"@[a-z0-9-]+\.iam\.gserviceaccount\.com"):
        assert not re.search(pat, blob), f"credential-shaped string matching {pat}"
    # COVERAGE against the panel of record, not against a count.
    s1 = ROOT / "results" / "tables" / "supplementary_table_s1.csv"
    if s1.exists():
        want = {r["accession"] for r in csv.DictReader(s1.open())
                if r["in_three_arm_panel"] == "True"}
        have = {r["accession"] for r in rows if r["accession"]}
        assert not (want - have), f"panel accessions missing from the manifest: {want - have}"


def test_the_manifest_does_not_claim_fetch_time_provenance():
    """P1.8. The timestamps are bucket writes. No provider checksum was ever recorded."""
    src = _read("scripts/raw_inputs.py")
    assert "uploaded_utc" in src, (
        "the timestamp column must be named for what it is; `retrieval_date` or `fetched_utc` "
        "would assert a fetch-time record that does not exist")
    f = ROOT / "results" / "tables" / "raw_inputs.csv"
    if f.exists():
        head = f.read_text().splitlines()[0]
        for forbidden in ("retrieval_date", "fetched_utc", "download_time"):
            assert forbidden not in head, (
                f"{forbidden} claims a fetch-time record; only the bucket write time exists")
        assert "url_provenance" in head, (
            "recorded and reconstructed URLs are different evidence and must be distinguished")


def test_both_model_loading_paths_honour_a_pinned_revision():
    """P1.9. The bake path honoured it; the runtime path silently did not.

    A revision pinned in config would have pinned the baked weights and left every sweep free
    to fetch whatever the hub's main pointed at, which is worse than no pin because it reads
    as one.
    """
    for rel in ("docker/bake_weights.py", "src/rbp/models/lm.py"):
        src = _read(rel)
        assert "revision" in src, f"{rel} does not pass a revision to from_pretrained"
        assert re.search(r'"revision":\s*(?:rev|spec\["revision"\])', src), (
            f"{rel} mentions revision but does not put it in the from_pretrained kwargs")


def test_the_open_release_decisions_are_still_recorded():
    """P0.6 and P1.14. Both are the author's to make, and both are easy to lose quietly.

    They are deferred, not resolved. The failure mode for a deferred item is that it stops
    being mentioned and then stops being remembered, which is indistinguishable from having
    decided it. This fails if either disappears from the submission checklist.
    """
    sub = _read("SUBMISSION.md")
    assert "AI-use disclosure needs an author decision" in sub, (
        "the AI-use disclosure decision has dropped out of the pre-submission checklist")
    assert "billing-account ID needs an explicit decision" in sub, (
        "the historical billing-account decision has dropped out of the checklist")
    paper = ROOT / "manuscript" / "paper.tex"
    if paper.exists():
        assert "OPEN ITEM, TO BE SETTLED BEFORE SUBMISSION" in paper.read_text(), (
            "the marker above the AI-use paragraph is gone; either the decision was made and "
            "this test should be retired, or it was lost")


def test_the_history_scan_counts_are_generated_not_typed():
    """P1.14 and the standing rule: no hand-maintained count anywhere."""
    sec = _read("SECURITY.md")
    assert not re.search(r"\((\d{3}) at the time of writing\)", sec), (
        "SECURITY.md is hand-maintaining a commit count again; it belongs in "
        "results/tables/history_scan.csv")
    f = ROOT / "results" / "tables" / "history_scan.csv"
    if not f.exists():
        pytest.skip("history_scan.csv not in this checkout")
    rows = {r["check"]: int(r["value"]) for r in csv.DictReader(f.open())}
    assert rows["commits containing credential material of any kind"] == 0, (
        "credential material in the history: stop and rotate")
    # The three counts SECURITY.md quotes must be the ones the scan produced.
    for k in ("commits containing a GCP billing account ID",
              "commit-and-file pairs containing a GCP billing account ID",
              "diff lines containing a GCP billing account ID"):
        assert k in rows, f"{k} is missing from the scan"
        assert f"**{rows[k]} " in sec or f"**{rows[k]}**" in sec, (
            f"SECURITY.md does not state the scanned value for {k} ({rows[k]})")


def test_the_history_scan_degrades_honestly_without_git():
    """The archival case. It crashed with a traceback there, and run.sh gates on it.

    A git export, a Zenodo deposit or a tarball has the files and no .git, so `git rev-list`
    exits 128. The scan must neither crash nor silently pass: there is no history present to be
    clean or dirty, so it reports the committed finding and says it could not rescan. The third
    time this repository has been bitten at the boundary between "cannot look" and "found
    nothing", and the first time in code written to warn about that boundary.
    """
    src = _read("scripts/history_scan.py")
    assert 'if not (ROOT / ".git").exists():' in src, (
        "history_scan.py does not handle a tree without .git, so it raises inside any "
        "unpacked archive and takes run.sh with it")
    assert "no .git and no committed history_scan.csv" in src, (
        "with neither git nor the committed table, nothing establishes what the history "
        "contains, and that must fail rather than pass quietly")


def test_the_shouty_header_habit_does_not_come_back():
    """605 docstring sections opened with the same ALL-CAPS device across 121 files.

    Used once it is emphasis. Used 605 times in one uniform pattern it is a tic, and a
    repository where every module shouts in the same voice reads as though one process wrote
    all of it rather than as though people worked on it. The content was worth keeping and is
    unchanged; only the shouting went.

    A ceiling rather than zero, because genuine emphasis is legitimate and this should not
    become a rule that forbids it. Counted over comments and docstrings only.
    """
    import ast
    import io
    import tokenize
    head = re.compile(r"^\s*#?\s*[A-Z][A-Z0-9 ,'()./\"-]{14,}[.,:]")
    hits = []
    for f in sorted((ROOT / "scripts").glob("*.py")) + sorted((ROOT / "src").rglob("*.py")) \
            + sorted((ROOT / "cloud").rglob("*.py")) + sorted((ROOT / "tests").rglob("*.py")):
        text = f.read_text(errors="ignore")
        try:
            tree = ast.parse(text)
        except SyntaxError:
            continue
        ok = set()
        for node in ast.walk(tree):
            if isinstance(node, (ast.Module, ast.FunctionDef, ast.AsyncFunctionDef,
                                 ast.ClassDef)):
                b = getattr(node, "body", [])
                if (b and isinstance(b[0], ast.Expr) and isinstance(b[0].value, ast.Constant)
                        and isinstance(b[0].value.value, str)):
                    ok.update(range(b[0].lineno, b[0].end_lineno + 1))
        for tok in tokenize.generate_tokens(io.StringIO(text).readline):
            if tok.type == tokenize.COMMENT:
                ok.add(tok.start[0])
        for i, ln in enumerate(text.split("\n"), 1):
            if i in ok and head.match(ln):
                hits.append(f"{f.relative_to(ROOT)}:{i}")
    assert len(hits) <= 60, (
        f"{len(hits)} ALL-CAPS docstring headers, up from 21. The device is back:\n  "
        + "\n  ".join(hits[:15]))


def test_the_trailer_change_is_recorded_as_granularity_not_a_narrowing():
    """206 of the first 228 commits carry an AI co-author trailer; later ones do not.

    A signal that simply stops looks like concealment. A signal that stops with the reason
    written down, the historical ones left in place and the manuscript named as the disclosure
    of record, is a change of granularity. The difference is entirely in whether it is written
    down, so this fails if it stops being.
    """
    sub = _read("SUBMISSION.md")
    assert "Per-commit AI co-author trailers stop after" in sub, (
        "the trailer change has dropped out of the submission checklist, which makes it look "
        "like a signal that quietly stopped rather than one that was deliberately relocated")
    assert "change of granularity and not of disclosure" in sub
    paper = ROOT / "manuscript" / "paper.tex"
    if paper.exists():
        assert "trailers stop after" in paper.read_text(), (
            "the AI-use paragraph must know it is now the disclosure of record, because "
            "narrowing it later would leave nothing in its place")
