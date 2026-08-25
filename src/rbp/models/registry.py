"""One place that knows about every model, so the trainer stays model-agnostic.

Model definitions, learning rates and LoRA settings all come from config/params.yaml.
Nothing about a model is written twice: adding a backbone is a config edit plus, if it
is a new architecture family, a builder here.

Each handle exposes the same surface -- build, batch, forward -- so the trainer never
branches on which model it is holding.
"""

import torch

from ..utils import config as cfgmod
from . import cnn as cnn_mod
from . import lm as lm_mod

CNN_KIND = "cnn"
LM_KIND = "lm"


def spec(name, cfg=None):
    cfg = cfg or cfgmod.load()
    models = cfg["models"]
    if name not in models:
        raise ValueError(f"unknown model {name!r}; config defines {sorted(models)}")
    return dict(models[name])


def names(cfg=None):
    cfg = cfg or cfgmod.load()
    return tuple(cfg["models"].keys())


def label(name, cfg=None):
    return spec(name, cfg).get("label", name)


def mode(name, cfg=None):
    """How the backbone is adapted: scratch, full or lora.

    `lora_threshold_m` in the config decides for pretrained models, so the rule is
    stated once as a size threshold rather than repeated per model.
    """
    cfg = cfg or cfgmod.load()
    s = spec(name, cfg)
    if s.get("mode"):
        return s["mode"]
    if s["kind"] == CNN_KIND:
        return "scratch"
    return "lora" if s.get("params_m", 0) > cfg["lora_threshold_m"] else "full"


class Handle:
    def __init__(self, name, model, spec_, tcfg):
        self.name = name
        self.model = model
        self.spec = spec_
        self.tcfg = tcfg
        self.kind = spec_["kind"]
        self.mode = spec_["_mode"]
        self.label = spec_.get("label", name)

    def to(self, device):
        self.model.to(device)
        return self

    def batch(self, seqs, device):
        if self.kind == CNN_KIND:
            return (cnn_mod.one_hot_batch(seqs).to(device),)
        return lm_mod.encode(self.model.tok, seqs, device)

    def forward(self, inputs):
        return self.model(*inputs)

    def param_groups(self):
        t = self.tcfg
        if self.kind == CNN_KIND:
            return [{"params": self.model.parameters(), "lr": float(t["lr_head"])}]
        lr_body = float(t["lr_encoder_lora"] if self.mode == "lora"
                        else t["lr_encoder_full"])
        body = [p for p in self.model.encoder.parameters() if p.requires_grad]
        return [{"params": body, "lr": lr_body},
                {"params": self.model.head.parameters(), "lr": float(t["lr_head"])}]

    def state(self):
        if self.kind == CNN_KIND:
            return self.model.state_dict()
        return lm_mod.trainable_state(self.model, self.mode)

    def load(self, state):
        missing, _ = self.model.load_state_dict(state, strict=False)
        return missing

    def sizes(self):
        m = self.model
        total = sum(p.numel() for p in m.parameters())
        train = sum(p.numel() for p in m.parameters() if p.requires_grad)
        return {"params_total": total, "params_trainable": train,
                "trainable_frac": round(train / total, 5) if total else 0.0}


def build(name, cfg=None):
    cfg = cfg or cfgmod.load()
    s = spec(name, cfg)
    s["_mode"] = mode(name, cfg)
    t = cfg["train"]

    if s["kind"] == CNN_KIND:
        model = cnn_mod.DeepBindCNN(
            channels=(s["conv1"]["out"], s["conv2"]["out"]),
            kernels=(s["conv1"]["kernel"], s["conv2"]["kernel"]),
            hidden=s["hidden"], pool=s["pool"], dropout=s["dropout"])
    elif s["kind"] == LM_KIND:
        model = lm_mod.build_from_spec(s, t)
    else:
        raise ValueError(f"unknown model kind {s['kind']!r} for {name!r}")
    return Handle(name, model, s, t)


def device_of(prefer=None):
    """cuda if present, then Apple MPS, else cpu. An explicit choice always wins."""
    if prefer:
        return torch.device(prefer)
    if torch.cuda.is_available():
        return torch.device("cuda")
    if torch.backends.mps.is_available():
        return torch.device("mps")
    return torch.device("cpu")
