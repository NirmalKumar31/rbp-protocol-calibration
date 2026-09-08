# Documentation

| file | contents |
|---|---|
| [REPRODUCE.md](REPRODUCE.md) | how to check every published number offline, and what the full pipeline needs if you want to rebuild from raw ENCODE files |
| [PANELS.md](PANELS.md) | how the dataset panel was defined, and why the counts 189, 95, 94 and 74 differ |
| [architecture.md](architecture.md) | stage graph; **the provider split as a measured decision** (live quota readings, measured $/Mpair per model, and the accelerator benchmark that picks the middle of the range); object layout; one task end to end |
| [cloud-setup.md](cloud-setup.md) | provisioning from an empty project: Terraform, service accounts, IAM, budget guards and the killswitch |
| [operating.md](operating.md) | running a sweep, reading a failure, estimating cost, and when to stop |

The manuscript, its figure legends and the built PDFs are in
[`../manuscript/`](../manuscript/).

## The pre-specified protocols

These four are cited by name in the paper, because the claims that rest on them are only
checkable if you can read them.

| file | contents |
|---|---|
| [EXTERNAL_BENCHMARK_PROTOCOL.md](EXTERNAL_BENCHMARK_PROTOCOL.md) | the eligibility rule and decision thresholds for the external replication, fixed before the held-out subset was scored |
| [EXTERNAL_BENCHMARK_AMENDMENT.md](EXTERNAL_BENCHMARK_AMENDMENT.md) | the directional estimand and the two fold repairs, committed before either was run |
| [EXTERNAL_CORRECTION_1.md](EXTERNAL_CORRECTION_1.md) | implementation defects found after the results entered the manuscript, with the withdrawn values preserved |
| [SENSITIVITY_SPEC.md](SENSITIVITY_SPEC.md) | the class-ratio and capacity-ladder analyses, specified before running |
