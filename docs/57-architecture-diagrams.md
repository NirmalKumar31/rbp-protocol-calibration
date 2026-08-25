# 57. Architecture diagrams

Seven Mermaid diagrams. Render in GitHub, VS Code (Markdown Preview Mermaid Support), or
mermaid.live.

---

## 1. The pipeline: fifteen stages, and where each one runs

Colour = execution environment. Note that **stage 6 sits after stage 5**, which is the one
non-obvious edge in the whole graph.

```mermaid
flowchart TD
    subgraph LOCAL["LAPTOP — submits and reads, never computes"]
        S0["stage 0<br/>preflight<br/>$0"]
        S6["stage 6<br/>select panel<br/>$0"]
        S14["stage 14<br/>verify<br/>$0"]
    end

    subgraph GCPCTL["GCP control plane"]
        S1["stage 1<br/>terraform<br/>$0"]
        S2["stage 2<br/>build images<br/>~$0.50"]
    end

    subgraph BATCH["GCP BATCH — CPU fan-out, quota-capped at 8 concurrent"]
        S3["stage 3<br/>ingest<br/>PUBLIC IP"]
        S4["stage 4<br/>panel<br/>PUBLIC IP"]
        S5["stage 5<br/>prep 488 tasks<br/>~$2"]
        S7["stage 7<br/>rehearsal<br/>~$0.60"]
        S8["stage 8<br/>CNN sweep<br/>~$3"]
        S11["stage 11<br/>variants + phyloP<br/>PUBLIC IP"]
        S13["stage 13<br/>analysis<br/>~$0.10"]
    end

    subgraph MODAL["MODAL — GPU, no quota gate"]
        S9["stage 9<br/>SpliceBERT<br/>~$31 REAL MONEY"]
        S10["stage 10<br/>locality probe<br/>~$0.30"]
        S12["stage 12<br/>ClinVar x3 arms<br/>~$0.60"]
    end

    R1(["R1<br/>protocol effect"])
    R2(["R2<br/>four models"])
    R3(["R3<br/>locality"])
    R4(["R4<br/>variant ladder"])

    S0 --> S1 --> S2 --> S3 --> S4 --> S5 --> S6
    S6 --> S7 --> R1
    S6 --> S8 --> R2
    S6 --> S11
    S6 --> S9 --> R2
    S9 --> S10 --> R3
    S9 --> S12
    S11 --> S12 --> R4
    R1 --> S13
    R2 --> S13
    R3 --> S13
    R4 --> S13
    S13 --> S14

    classDef local fill:#e8eef7,stroke:#4878a8,color:#000
    classDef gcp fill:#e6f0ea,stroke:#2b6a4d,color:#000
    classDef modal fill:#fdf0e3,stroke:#e08214,color:#000
    classDef money fill:#fbe4e6,stroke:#b2182b,color:#000,stroke-width:3px
    classDef result fill:#f2f2f2,stroke:#555,color:#000

    class S0,S6,S14 local
    class S1,S2,S3,S4,S5,S7,S8,S11,S13 gcp
    class S10,S12 modal
    class S9 money
    class R1,R2,R3,R4 result
```

**Why stage 6 is after stage 5.** `pairs` — the number the panel is size-ranked on — counts the
positives that could *actually be matched* to a negative. Matching is a search, so `pairs` is a
**result** of preprocessing, not an input. A size-ranked panel cannot be chosen before prep has
produced the counts.

**The three-way fan-out** from stage 6 into stages 7, 8/9 and 11 is real parallelism: those
branches share no inputs and draw on **different quota pools** (GCP vCPU vs Modal containers).
`./run.sh parallel` exploits it.

---

## 2. The two-cloud split, and why

```mermaid
flowchart LR
    subgraph WALL["The quota wall that forced this"]
        direction TB
        Q1["GCP: GPUS_ALL_REGIONS = 0<br/>increase auto-denied<br/>NOT_ENOUGH_USAGE_HISTORY"]
        Q2["AWS: 0 on all four GPU families"]
        Q3["Azure: GPU forbidden on free trial"]
        Q4["Modal: no quota gate"]
    end

    subgraph GCS["GCS — the ONLY shared state"]
        direction TB
        B1[("raw/<br/>immutable inputs")]
        B2[("processed/<br/>matched datasets")]
        B3[("manifest/<br/>task lists + study_panel.tsv")]
        B4[("runs/<br/>weights + fold scores")]
        B5[("variants/<br/>ClinVar, 3 arms")]
        B6[("results/<br/>tables + figures")]
    end

    BATCH["GCP BATCH<br/>CPU fan-out<br/>12 vCPU cap"]
    MODALC["MODAL<br/>A10G / T4<br/>max_containers = 10"]

    Q1 --> MODALC
    Q4 --> MODALC
    BATCH -->|"writes"| B1
    BATCH -->|"writes"| B2
    BATCH -->|"writes"| B3
    B2 -->|"reads"| BATCH
    B2 -->|"reads"| MODALC
    B3 -->|"reads"| MODALC
    MODALC -->|"writes"| B4
    MODALC -->|"writes"| B5
    B4 -->|"reads"| MODALC
    B5 --> B6
    B4 --> B6

    classDef wall fill:#fbe4e6,stroke:#b2182b,color:#000
    classDef store fill:#f2f2f2,stroke:#555,color:#000
    classDef compute fill:#e6f0ea,stroke:#2b6a4d,color:#000
    class Q1,Q2,Q3 wall
    class Q4 compute
    class B1,B2,B3,B4,B5,B6 store
    class BATCH,MODALC compute
```

**What made the split cheap.** A training task's only cloud dependency is
`google-cloud-storage`. It reads its manifest and dataset from GCS, writes results back, and
takes its index from an environment variable. No Batch API calls, no metadata-server
assumptions. So `cloud_train.py` runs on Modal **unchanged**, and `aggregate` cannot tell
afterwards which platform produced a row.

---

## 3. Object layout, and the marker discipline

```mermaid
flowchart TD
    RAW[("BUCKET: PROJECT-raw<br/>immutable")]
    DER[("BUCKET: PROJECT-derived<br/>everything computed")]
    ART[("BUCKET: PROJECT-artifacts<br/>image digests")]
    TFS[("BUCKET: PROJECT-tfstate<br/>terraform state, versioned")]

    RAW --> R1["GRCh38.primary_assembly.genome.fa  3.1 GB"]
    RAW --> R2["gencode.v45.annotation.gtf.gz"]
    RAW --> R3["clinvar.vcf.gz"]
    RAW --> R4["peaks/K562/ + peaks/HepG2/  244 files"]
    RAW --> R5["config/panel_full.tsv  139 K562<br/>config/panel_full_HepG2.tsv  105"]

    DER --> D1["processed/{arm}/{cell}/{protein}/dataset.tsv"]
    DER --> D2["panel/{arm}/panel_final_{cell}_{arm}.tsv"]
    DER --> D3["manifest/prep_tasks.tsv<br/>manifest/rehearsal_tasks.tsv<br/>manifest/sweep_tasks{tag}.tsv<br/>manifest/study_panel.tsv  ← THE panel"]
    DER --> D4["runs/{arm}/{cell}/{protein}/{model}/fold{f}/<br/>best.pt + scores.tsv.gz + metrics.json"]
    DER --> D5["variants/tables/  windows, no FASTA needed<br/>variants/scores_sb/  matched head<br/>variants/scores_mm/  MISMATCHED head"]
    DER --> D6["results/tables/*.csv + results/figures/*"]

    ART --> A1["images/cpu_digest.txt<br/>images/gpu_digest.txt"]

    classDef bucket fill:#e8eef7,stroke:#4878a8,color:#000
    classDef key fill:#e6f0ea,stroke:#2b6a4d,color:#000
    class RAW,DER,ART,TFS bucket
    class D3,D4 key
```

**GCS has no directories.** `processed/dinuc/K562/QKI/dataset.tsv` is one flat object name
containing slashes. Consequence: a wrong path is indistinguishable from missing data, because
there is no schema to violate. Always list the parent prefix before concluding data is absent.

### The write order inside one task

```mermaid
sequenceDiagram
    participant T as task
    participant G as GCS
    Note over T,G: There are NO cross-object transactions.<br/>So order the writes deliberately.
    T->>G: 1. best.pt (75 MB, slow)
    T->>G: 2. scores.tsv.gz
    T->>G: 3. metrics.json  ← the MARKER, written LAST
    Note over T,G: Crash before 3: marker absent,<br/>next run REDOES the work. Costly, correct.
    Note over T,G: Marker first would mean a crash leaves<br/>a marker with no payload: skipped forever,<br/>silently, and the gap is permanent.
```

---

## 4. Identity and blast radius

```mermaid
flowchart LR
    subgraph SA["Five service accounts, one job each"]
        I["rbp-ingest"]
        P["rbp-prep"]
        T["rbp-train"]
        A["rbp-analysis"]
        M["rbp-modal<br/>the ONLY key that<br/>leaves Google's network"]
    end

    RAWB[("raw/")]
    PROC[("processed/")]
    RUNS[("runs/ ckpt/ variants/")]
    RES[("results/")]

    I -->|"write"| RAWB
    RAWB -->|"read"| P
    P -->|"write"| PROC
    PROC -->|"read"| T
    T -->|"write"| RUNS
    PROC -->|"read"| A
    RUNS -->|"read"| A
    A -->|"write"| RES
    M -->|"write, IAM CONDITION"| RUNS

    classDef risky fill:#fbe4e6,stroke:#b2182b,color:#000,stroke-width:3px
    classDef sa fill:#e8eef7,stroke:#4878a8,color:#000
    classDef store fill:#f2f2f2,stroke:#555,color:#000
    class M risky
    class I,P,T,A sa
    class RAWB,PROC,RUNS,RES store
```

`rbp-modal`'s key sits in a Modal secret on a third-party platform, so its binding carries a
CEL condition:

```
resource.name.startsWith(".../objects/runs/")   ||
resource.name.startsWith(".../objects/ckpt/")   ||
resource.name.startsWith(".../objects/variants/")
```

Fully compromised, it cannot alter a dataset under `processed/`, cannot touch the raw bucket,
and cannot delete anything outside those three prefixes.

**This fired in practice.** The first Modal ClinVar probe scored its dataset correctly and then
took a 403 writing `variants/scores_sb/K562_AATF.csv` — because `variants/` was not yet in the
list. That is the guardrail working: a new write path stays denied until somebody widens it
deliberately, in Terraform.

---

## 5. Networking: deliberately not uniform

```mermaid
flowchart TD
    subgraph NET["rbp-net — custom VPC, no auto subnets"]
        SUB["subnet rbp-workers  10.0.0.0/20<br/>Private Google Access ON<br/>noExternalIpAddress: true"]
        GPUSUB["subnet rbp-gpu-us-central1<br/>present, unused: GPU quota is 0"]
    end

    subgraph DEF["default VPC — external IP"]
        EXT["one short-lived VM"]
    end

    W["prep / rehearsal / sweep / analysis<br/>488 workers"] --> SUB
    SUB -->|"PGA route"| GAPI["*.googleapis.com<br/>ONLY"]
    SUB -.->|"NO ROUTE"| WEB["huggingface.co<br/>hgdownload.soe.ucsc.edu"]

    E["ingest / panel / variants"] --> EXT
    EXT --> ENCODE["ENCODE, GENCODE,<br/>NCBI, UCSC phyloP"]

    classDef sealed fill:#e6f0ea,stroke:#2b6a4d,color:#000
    classDef open fill:#fdf0e3,stroke:#e08214,color:#000
    classDef blocked fill:#fbe4e6,stroke:#b2182b,color:#000
    class SUB,W,GAPI sealed
    class DEF,EXT,E,ENCODE open
    class WEB blocked
```

**Why no Cloud NAT.** It bills per gateway-hour *plus* per GB, which for occasional downloads
costs more than the downloads; it hands internet access to 488 workers that have no reason to
have it; and it makes every run depend on third-party sites being up.

**The consequence you must design around.** A sealed worker cannot
`from_pretrained("multimolecule/splicebert")` — that is why model weights are **baked into the
image**, so weights and code share one digest.

**The debugging signature.** A sealed worker attempting a non-Google host does not fail fast; it
**hangs** until `maxRunDuration`, because there is no route and the SYN goes nowhere. If a task
times out doing something that should be quick, suspect the network before the code.

---

## 6. The guards, and what each one catches

```mermaid
flowchart TD
    START(["run a stage"]) --> G0{"preflight<br/>passed?"}
    G0 -->|no| X0["STOP<br/>APIs, billing, quota,<br/>auth, bucket names"]
    G0 -->|yes| G1{"paid stage?"}
    G1 -->|yes| C["confirm<br/>type literal YES<br/>estimate printed"]
    G1 -->|no| G2
    C --> G2{"Modal stage?"}
    G2 -->|yes| GM{"modal gate<br/>secret + credit"}
    GM -->|no| XM["STOP"]
    GM -->|yes| G3
    G2 -->|no| G3{"terraform?"}
    G3 -->|yes| GD{"destroys<br/>in plan?"}
    GD -->|"any"| XD["STOP<br/>state describes<br/>another project"]
    GD -->|zero| RUN
    G3 -->|no| RUN["submit job<br/>count read from MANIFEST"]
    RUN --> W{"WAIT for<br/>completion"}
    W --> MK{"marker<br/>present?"}
    MK -->|no| REDO["redo the work"]
    MK -->|yes| SKIP["skip, cheap rerun"]
    SKIP --> V{"stage 14<br/>verify"}
    REDO --> V
    V -->|"outside tolerance"| XV["STOP<br/>ran, but did not reproduce.<br/>Do not write it up."]
    V -->|"all pass"| DONE(["reproduction verified"])

    classDef stop fill:#fbe4e6,stroke:#b2182b,color:#000
    classDef gate fill:#fdf0e3,stroke:#e08214,color:#000
    classDef ok fill:#e6f0ea,stroke:#2b6a4d,color:#000
    class X0,XM,XD,XV stop
    class G0,G1,G2,G3,GD,GM,W,MK,V gate
    class RUN,SKIP,REDO,DONE,C ok
```

Every one of these exists because of a specific incident:

| guard | the incident |
|---|---|
| preflight | every failure of the first build was a free-to-detect environment problem found *after* spending |
| confirm + printed estimate | a 44%-too-high cost estimate reached $6 out of pocket before anyone compared it to the bill |
| Modal gate separate from GCP gate | the secret needs a service account that stage 1 creates: a single gate is unsatisfiable on a fresh project |
| **destroy guard** | a plan proposing **63 destroys** of the original study's buckets, with `-auto-approve` |
| count from manifest | `COUNT=189` against a 187-row manifest failed a job whose every real task succeeded |
| **wait for completion** | `finalize` would have run seconds after submit and written a **zero-dataset panel** |
| marker written last | a marker-first crash skips the task forever, silently |
| verify | a pipeline that completes and produces different science is worse than one that crashes |

---

## 7. One Batch task, end to end

```mermaid
sequenceDiagram
    autonumber
    participant U as run.sh
    participant S as submit.sh
    participant API as Batch API
    participant MIG as Managed Instance Group
    participant VM as VM + agent
    participant C as container
    participant G as GCS

    U->>S: submit.sh prep
    S->>G: read manifest/prep_tasks.tsv
    G-->>S: 488 rows
    Note over S: COUNT=488, never typed
    S->>G: read images/cpu_digest.txt
    S->>API: jobs.create (spec pinned by digest)
    API-->>S: "successfully submitted"
    Note over S,API: CONTROL PLANE ONLY.<br/>Says nothing about the work.
    API->>MIG: create 2 VMs (3rd refused: QUOTA_EXCEEDED)
    MIG->>VM: boot COS + Batch agent
    VM->>API: register, poll for tasks
    API-->>VM: BATCH_TASK_INDEX = 149, 167, 32, 78
    Note over API,VM: Batch PARTITIONS the index space.<br/>Manifest order is NOT honoured.
    VM->>C: docker run --env BATCH_TASK_INDEX=149
    C->>G: read manifest row 149
    C->>G: read raw/peaks/... + genome
    Note over C: OMP_NUM_THREADS=1<br/>else 16 threads over 4 cores
    C->>G: write processed/.../dataset.tsv
    C->>G: write report (MARKER, last)
    C-->>VM: exit 0
    VM->>API: task SUCCEEDED
    loop every 60s
        S->>API: describe job
        API-->>S: RUNNING, counts
    end
    API-->>S: SUCCEEDED
    S-->>U: exit 0
    U->>U: NOW finalize, safely
```

**The two lines to remember.** Step 6, `"successfully submitted"`, is the control plane
accepting your JSON — nothing more. And step 12: because Batch partitions the index space, you
cannot express "run model X first" by ordering the manifest. Scope is expressed by *which*
manifest, which is why `MANIFEST_TAG` exists.
