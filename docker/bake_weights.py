"""Download every pretrained backbone into the image at BUILD time, then prove it worked.

WHY THIS FILE EXISTS. The sweep workers run with no external IP address, no Cloud NAT, and
Private Google Access. PGA routes to Google APIs only. `huggingface.co` is not a Google
API, so `from_pretrained("multimolecule/splicebert")` on a worker cannot resolve and the
run dies -- or worse, hangs on a retry loop until maxRunDuration. Four of the five models
are pretrained backbones, so that is 80% of the sweep.

Three ways out were available:

  * add Cloud NAT      billed per gateway-hour in five regions, gives the workers general
                       internet access they have no reason to have, and makes every run
                       depend on huggingface.co being up
  * mirror to GCS      works, but adds a download step and a second copy to keep in sync
  * bake into image    the weights become part of the artefact that the digest identifies

The third is the only one where "which weights produced this result?" has the same answer
as "which image produced this result?". Weights are pinned by the same digest as the code.

THE LIST COMES FROM CONFIG, NOT FROM HERE. Adding a sixth model to params.yaml must not
require remembering to edit a Dockerfile. Anything with kind: lm is fetched.

Run in two modes:
    python docker/bake_weights.py fetch    # build time, with network
    python docker/bake_weights.py verify   # build time, with HF_HUB_OFFLINE=1
"""

import os
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from rbp.utils import config as cfgmod  # noqa: E402


def lm_specs(cfg):
    return [(name, s) for name, s in cfg["models"].items() if s.get("kind") == "lm"]


def fetch(cfg):
    """Pull tokenizer and weights for every language-model backbone."""
    import multimolecule as mm
    from multimolecule import RnaTokenizer

    for name, s in lm_specs(cfg):
        repo = s["repo"]
        RnaTokenizer.from_pretrained(repo)
        getattr(mm, s["cls"]).from_pretrained(repo)
        print(f"  baked {name:12} {repo}", flush=True)


def verify(cfg):
    """Build every model with the network unavailable.

    This is the gate, not the fetch above. A fetch that quietly failed would leave an image
    that looks fine and dies on the first GPU task at 2am. Building offline here means a
    missing weight fails the BUILD, which costs a few minutes instead of a night.
    """
    os.environ["HF_HUB_OFFLINE"] = "1"
    os.environ["TRANSFORMERS_OFFLINE"] = "1"
    import torch

    from rbp.models import registry

    # The GPU this image is built for is a V100, compute capability sm_70. CUDA 13 dropped
    # Volta, so a wheel built against it produces "no kernel image is available for
    # execution on the device" -- at run time, on a real GPU, after everything else passed.
    # The arch list is compiled into the wheel and readable on a machine with no GPU at
    # all, which is what makes this checkable here rather than only in the sweep.
    arches = torch.cuda.get_arch_list()
    if "sm_70" not in arches:
        sys.exit(f"torch was built without sm_70 (Volta); a V100 cannot run it. "
                 f"Compiled for: {arches}")
    print(f"  torch {torch.__version__}, cuda {torch.version.cuda}, arches include sm_70",
          flush=True)

    names = registry.names(cfg)
    if not names:
        sys.exit("config defines no models")
    for name in names:
        h = registry.build(name, cfg)
        sz = h.sizes()
        print(f"  {name:12} {h.mode:8} {sz['params_total']:>12,} params "
              f"({sz['trainable_frac']:.2%} trainable)", flush=True)
    n_lm = len(lm_specs(cfg))
    print(f"all {len(names)} models build offline ({n_lm} pretrained backbones)")


if __name__ == "__main__":
    mode = sys.argv[1] if len(sys.argv) > 1 else "fetch"
    c = cfgmod.load()
    {"fetch": fetch, "verify": verify}[mode](c)
