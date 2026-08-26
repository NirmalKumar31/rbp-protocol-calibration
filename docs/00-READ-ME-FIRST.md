# The documentation set for the reproducible rebuild

Written 2026-08-25, covering everything after
`rna-binding-proteins/docs/51-chronicle-10-*.md`, which ended with the ClinVar arm complete.

Read in this order. Each file assumes the ones before it.

| file | what it teaches | read it if |
|---|---|---|
| [52-why-a-rebuild.md](52-why-a-rebuild.md) | why the working study was not reproducible, the architecture, every design decision and its alternative | you want the *why* before the *how* |
| [53-cloud-from-zero.md](53-cloud-from-zero.md) | **the big one.** GCP from an empty project, taught three ways for every single thing: CLI/Terraform, the UI click by click, and what the backend actually does | you want to become fluent in this, not just able to run it |
| [54-every-file-every-line.md](54-every-file-every-line.md) | every file in the repo, every function, every non-obvious line, and why it is written that way | you are editing the code |
| [55-the-bug-chronicle.md](55-the-bug-chronicle.md) | all 16 bugs in the order they happened: symptom, mechanism, what it would have cost, the fix, and what class of thinking finds it | you want to develop the instinct rather than memorise the fixes |
| [57-architecture-diagrams.md](57-architecture-diagrams.md) | seven Mermaid diagrams: the stage graph, the two-cloud split, object layout, IAM blast radius, networking, the guard chain, and one task end to end | you think in pictures |
| [56-operating-and-monitoring.md](56-operating-and-monitoring.md) | how to watch a run, read a failure, estimate cost, and decide when to stop | you are running it |
| [58-the-run-chronicle.md](58-the-run-chronicle.md) | **the narrative record of the build-and-run day.** Twenty-four bugs in the order they happened, each with symptom, mechanism, cost, fix and the kind of thinking that finds it | you want to learn from the sequence, not just the list |
| [AGENT-CONTEXT.md](AGENT-CONTEXT.md) | compressed state for an assistant resuming after context loss | you are an AI picking this up mid-flight |

## The one-paragraph summary

The study was scientifically finished and operationally unreproducible: eighteen files
hardcoded one GCP project, the dataset panel existed only as a command-line flag somebody
typed once, four analysis stages ran only on a laptop, and nothing verified that a rerun
produced the same science. `rbp-repro/` is the study rebuilt as fifteen cloud stages with a
preflight that refuses to spend on a broken environment and a final stage that asserts every
published number against a tolerance. Running it on a genuinely new project surfaced sixteen
bugs, one of which would have deleted the original study's results.
