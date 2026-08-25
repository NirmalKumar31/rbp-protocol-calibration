"""RNA language-model backbones with a shared classification head.

All four pretrained models get the same pooling and the same head, so a difference in
score is a difference in the backbone rather than in how it was read out.

Pooling is a mean over real nucleotide tokens only. Special tokens are masked out
because cls/eos/pad carry no sequence content and including them would dilute short
windows differently across tokenizers.

Large backbones are adapted with LoRA rather than fully fine-tuned: with a few thousand
windows per protein, updating ~100M weights overfits. Small backbones are tuned end to
end, where full fine-tuning is both affordable and stronger.
"""

import torch
import torch.nn as nn



def masked_mean_pool(hidden, ids, tok):
    keep = torch.ones_like(ids, dtype=torch.bool)
    for sid in (tok.cls_token_id, tok.eos_token_id, tok.pad_token_id):
        if sid is not None:
            keep &= ids != sid
    m = keep.unsqueeze(-1).to(hidden.dtype)
    return (hidden * m).sum(1) / m.sum(1).clamp(min=1.0)


class LMClassifier(nn.Module):
    def __init__(self, tokenizer, encoder, hidden=128, dropout=0.3):
        super().__init__()
        self.tok = tokenizer
        self.encoder = encoder
        self.head = nn.Sequential(
            nn.Linear(encoder.config.hidden_size, hidden), nn.ReLU(),
            nn.Dropout(dropout), nn.Linear(hidden, 1))

    def forward(self, ids, attention_mask):
        out = self.encoder(input_ids=ids, attention_mask=attention_mask)
        pooled = masked_mean_pool(out.last_hidden_state, ids, self.tok)
        return self.head(pooled).squeeze(-1)

    def apply_lora(self, r=8, alpha=16, dropout=0.05, targets=("query", "value")):
        from peft import LoraConfig, get_peft_model
        self.encoder = get_peft_model(self.encoder, LoraConfig(
            r=r, lora_alpha=alpha, lora_dropout=dropout,
            target_modules=list(targets), bias="none"))

    @property
    def n_trainable(self):
        return sum(p.numel() for p in self.parameters() if p.requires_grad)

    @property
    def n_total(self):
        return sum(p.numel() for p in self.parameters())


def build_from_spec(spec, train_cfg):
    """Tokenizer plus classifier for one `models:` entry, LoRA attached if the spec says so.

    Everything -- repo, class, head shape, LoRA rank -- comes from config, so this
    function holds no model-specific knowledge.
    """
    import multimolecule as mm
    from multimolecule import RnaTokenizer
    tok = RnaTokenizer.from_pretrained(spec["repo"])
    encoder = getattr(mm, spec["cls"]).from_pretrained(spec["repo"])
    model = LMClassifier(tok, encoder,
                         hidden=train_cfg.get("head_hidden", 128),
                         dropout=train_cfg.get("head_dropout", 0.3))
    if spec["_mode"] == "lora":
        lc = train_cfg.get("lora", {})
        # a model may override the target names: BERT-style backbones expose
        # query/value, ESM and MSA-style ones expose q_proj/v_proj
        targets = spec.get("lora_targets") or lc.get("targets", ("query", "value"))
        model.apply_lora(r=lc.get("r", 8), alpha=lc.get("alpha", 16),
                         dropout=lc.get("dropout", 0.05), targets=tuple(targets))
    return model


def encode(tok, seqs, device="cpu"):
    enc = tok(list(seqs), return_tensors="pt", padding=True)
    return enc["input_ids"].to(device), enc["attention_mask"].to(device)


def trainable_state(model, mode):
    """What to save. LoRA runs store only trainable tensors, keeping checkpoints small."""
    state = model.state_dict()
    if mode != "lora":
        return state
    keep = {n for n, p in model.named_parameters() if p.requires_grad}
    return {k: v for k, v in state.items() if k in keep}
