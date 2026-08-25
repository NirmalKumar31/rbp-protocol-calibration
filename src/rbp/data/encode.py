"""Query the live ENCODE portal for eCLIP experiments and their peak files.

Nothing here is hard-coded from a previous run: the protein panel and every file
accession are discovered from the API, so `config/proteins.tsv` is reproducible
and auditable.
"""

import requests

SEARCH = "{api}/search/"
FIELDS = [
    "accession",
    "target.label",
    "biosample_ontology.term_name",
    "files.accession",
    "files.output_type",
    "files.assembly",
    "files.file_format",
    "files.status",
    "files.preferred_default",
    "files.biological_replicates",
]


def _get(url, params, timeout=90):
    r = requests.get(url, params=params, headers={"accept": "application/json"}, timeout=timeout)
    r.raise_for_status()
    return r.json()


def experiments(cfg):
    """Released eCLIP experiments for the configured cell line."""
    e = cfg["encode"]
    params = [
        ("type", "Experiment"),
        ("assay_title", e["assay"]),
        ("biosample_ontology.term_name", e["cell_line"]),
        ("status", "released"),
        ("limit", "all"),
        ("format", "json"),
    ] + [("field", f) for f in FIELDS]
    return _get(SEARCH.format(api=e["api"]), params)["@graph"]


def peak_file(exp, assembly, output_type):
    """Accession of the replicate-reproducible peak BED, if present.

    ENCODE's eCLIP pipeline emits one peak file per biological replicate plus one
    merged across replicates, all labelled `output_type: peaks`. We want the merged
    file, so selection is by **most biological replicates**.

    Do not use `preferred_default` for this: it is applied inconsistently to eCLIP.
    For FUS (ENCSR069EVH) and RBFOX2 (ENCSR756CKJ) the flag sits on a single-replicate
    file while the merged [1,2] file is unflagged, so trusting it silently yields
    non-reproducible peaks.
    """
    best, best_n = None, -1
    for f in exp.get("files", []):
        if (f.get("output_type") == output_type
                and f.get("assembly") == assembly
                and f.get("status") == "released"
                and f.get("file_format") == "bed"):
            n = len(f.get("biological_replicates") or [])
            # more replicates wins; the preferred_default flag only breaks ties
            if n > best_n or (n == best_n and f.get("preferred_default") is True):
                best, best_n = f["accession"], n
    return best, max(best_n, 0)


def build_panel(cfg):
    """One row per protein: protein, peak file accession, experiment, cell line.

    Proteins with several eligible experiments keep the first encountered; the
    alternatives are returned separately so the choice is visible rather than silent.
    """
    e = cfg["encode"]
    rows, dupes, skipped = {}, [], []
    for exp in experiments(cfg):
        target = (exp.get("target") or {}).get("label")
        if not target:
            skipped.append((exp.get("accession"), "no target label"))
            continue
        acc, n_reps = peak_file(exp, e["assembly"], e["output_type"])
        if not acc:
            skipped.append((exp.get("accession"),
                            f"no preferred_default bed {e['output_type']} for {e['assembly']}"))
            continue
        name = target.replace("eGFP-", "").strip()
        row = {
            "protein": name,
            "accession": acc,
            "cell_line": e["cell_line"],
            "experiment": exp["accession"],
            "n_replicates": str(n_reps),
        }
        if name in rows:
            dupes.append(row)
        else:
            rows[name] = row
    return sorted(rows.values(), key=lambda r: r["protein"]), dupes, skipped


def peak_url(accession, api):
    return f"{api}/files/{accession}/@@download/{accession}.bed.gz"


def peak_path(root, protein, cell_line):
    """Locate a downloaded peak BED.

    Peaks live under data/raw/peaks/<cell_line>/ rather than one flat directory. The
    same protein has a different accession per cell line, so a flat layout plus a
    `{protein}.*.bed.gz` glob would resolve to whichever file the filesystem returned
    first and silently train on the wrong cell line.
    """
    d = root / "data" / "raw" / "peaks" / cell_line
    hits = sorted(d.glob(f"{protein}.*.bed.gz"))
    if not hits:
        raise FileNotFoundError(f"no peak file for {protein} in {d}")
    if len(hits) > 1:
        raise RuntimeError(f"{protein} in {cell_line} matches {len(hits)} files: "
                           f"{[h.name for h in hits]}")
    return hits[0]
