# Documentation

| file | contents |
|---|---|
| [REPRODUCE.md](REPRODUCE.md) | how to check every published number offline, and what the full pipeline needs if you want to rebuild from raw ENCODE files |
| [PANELS.md](PANELS.md) | how the dataset panel was defined, and why the counts 189, 95, 94 and 74 differ |
| [architecture.md](architecture.md) | stage graph; **the provider split as a measured decision** (live quota readings, measured $/Mpair per model, and the accelerator benchmark that picks the middle of the range); object layout; one task end to end |
| [cloud-setup.md](cloud-setup.md) | provisioning from an empty project: Terraform, service accounts, IAM, budget guards and the killswitch |
| [operating.md](operating.md) | running a sweep, reading a failure, estimating cost, and when to stop |
| [ZENODO.md](ZENODO.md) | minting the archival DOI for a release |

The manuscript and its figure legends are in [`../manuscript/`](../manuscript/); the submission
index is [`../SUBMISSION.md`](../SUBMISSION.md).

## Working notes

Development history, the internal review records and the drafting notes live on the
`working-notes` branch rather than here. They document how the study arrived at its current
claim set, including analyses that were withdrawn and the reasons, and are kept because that
record is worth having. They are not part of the release and are not needed to reproduce
anything.

```
git checkout working-notes
```
