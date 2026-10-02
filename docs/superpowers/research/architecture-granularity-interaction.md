# Does architecture interact with tokenization granularity?

Research notes. Every claim carries a URL and year. Quotes are short and marked. Where a paper
reports a comparison, the actual numbers are given. Findings only — no recommendations.

---

## SECTION 1 — Does architecture interact with tokenization granularity?

### 1.1 MambaByte (Wang, Gangavarapu, Yan, Rush)

arXiv 2401.13660, v1 Jan 2024, v2 Apr 2024, v3 Aug 2024.
https://arxiv.org/abs/2401.13660 · https://arxiv.org/html/2401.13660v2 · OpenReview https://openreview.net/forum?id=X1xNsuKssb (2024)

**The exact claim about why SSMs suit byte-level modelling.** The argument is memory, not
inductive bias:

> "unlike Transformers, whose memory scales linearly in sequence length, Mamba maintains a large
> fixed-size memory state, which makes it suitable for direct byte-level modeling"

Formalised in the paper: an *m*-layer Mamba with hidden state h(t) ∈ R^(n_state × d) maintains
`m × n_state × d` floats, independent of L_ctx. The paper's key inequality is that
`m × n_state × d ≫ L_ctx` "in all but extreme cases", so the state has room to encode L_ctx
bytes regardless of whether the input is bytes or subwords. Therefore "if Mamba can be used for
tokenized models, MambaByte should enable modeling byte-level sequences without the need for
length-compression trade-offs".

Secondary argument, also explicitly length-framed: training cost is O(L_ctx) for Mamba versus
O(L_ctx²/p² + L_ctx·p) for MegaByte at patch size p — O(L_ctx^4/3) even at the favourable patch
size L_ctx^(1/3).

**What was held fixed.** Two separate settings, not one:

- *Compute-matched* (FLOPs per training byte). Table 1 of the paper:
  MegaByte-758M+262M : MambaByte-353M = **1.02 : 1**;
  MegaByte-1.3B+350M : MambaByte-972M = **0.54 : 1** and MegaByte-1.3B+218M : MambaByte-972M = **0.40 : 1**.
  MegaByte uses patch size 8.
- *Parameter-matched* (Figure 1, PG19, 8,192 consecutive bytes): standard Transformer, MegaByte,
  gated diagonalised S4, MambaByte at a fixed parameter budget.
- Common recipe across all models: contiguous 8,192-byte sequences, one per document from a random
  position, BF16, same optimizer/LR schedule. Datasets: PG19, Stories, Books, ArXiv, Code.
- **Training bytes were NOT matched** in the headline tables — this is the weak point (see 1.2).
  MambaByte-353M saw 30B bytes versus 80B for the baselines; MambaByte-972M saw 150B versus 400B.

**Measured numbers — byte-level, medium scale (Table 2, test bits-per-byte, lower better).**
MegaByte-758M+262M and MambaByte-353M use the same FLOPs per byte.

| Model | Context | Bytes trained | PG19 | Stories | Books | ArXiv | Code |
|---|---|---|---|---|---|---|---|
| Transformer-320M | 1,024 | 80B | 1.057 | 1.064 | 1.097 | 0.816 | 0.575 |
| PerceiverAR-248M | 8,192 | 80B | 1.104 | 1.070 | 1.104 | 0.791 | 0.546 |
| MegaByte-758M+262M (patch 8) | 8,192 | 80B | 1.000 | 0.978 | 1.007 | 0.678 | 0.411 |
| **MambaByte-353M** | 8,192 | **30B** | **0.930** | **0.908** | **0.966** | 0.663 | **0.396** |

So MambaByte beats the compute-matched MegaByte on all five datasets with **0.63× less compute and
training data**.

**Byte vs subword, large scale (Table 3, PG-19, word-level perplexity).**

| Model | Vocab | Effective ctx (bytes) | Effective bytes trained | Val PPL | Test PPL |
|---|---|---|---|---|---|
| Transformer-XL (36L) | 32K | 2,048/4,096 | 400B | 45.5 | 36.3 |
| Compressive (36L) | 32K | 2,048/2×2,048 | 400B | 43.4 | 33.6 |
| Routing-490M | 82K | 32,768 | 330B | – | 33.2 |
| PerceiverAR-974.6M | 32K | 8,192 | 1.68T | 45.9 | 28.9 |
| Block-Recurrent-1.3B | 32K | 4,096/recurrence | – | – | **26.5** |
| **Mamba-1.03B (subword)** | 32K | 8,192 | 150B | 40.7 | 34.6 |
| Transformer-320M (byte) | 256 | 8,192 | 400B | 81.6 | 69.4 |
| PerceiverAR-248M (byte) | 256 | 8,192 | 400B | 119.1 | 88.8 |
| MegaByte-1.3B+350M (byte) | 256 | 8,192/patch 8 | 400B | 42.8 | 36.4 |
| **MambaByte-972M (byte)** | 256 | 8,192 | 150B | **39.6** | 33.0 |

**The single most relevant number for a granularity study.** MambaByte explicitly retrained a
subword Mamba-1.03B "to control for the benefits of the Mamba architecture", and found:

> "at the same parameter size and amount of training bytes, the (subword) Mamba and MambaByte
> perform similarly"

with the paper's own explanation being that "these models effectively have the same memory capacity
despite significant differences in the input sequence length". In other words: **inside the Mamba
architecture, moving from subword to byte granularity is roughly performance-neutral on language
modelling loss** (39.6/33.0 vs 40.7/34.6 PPL). The byte advantage shows up elsewhere — in noise
robustness and length extrapolation, not in loss.

**Where the granularity effect actually appears (Table 6, PG-19 word-PPL degradation under noise;
lower is better).** Mamba-1.03B subword vs MambaByte-972M:

| Noise | p | Mamba (subword) | MambaByte |
|---|---|---|---|
| Drop | 0.05 / 0.3 | +16.9 / +213.2 | **+8.5 / +31.7** |
| Repeat | 0.05 / 0.3 | +6.3 / +28.4 | **+6.2 / +26.6** |
| Antspeak | – | +58,300.0 | **+28.3** |
| Uppercase | 0.05 / 0.3 | +5.4 / +18.3 | **+1.6 / +5.5** |
| Random case | – | +20.8 | **+7.7** |
| Swap | 0.05 / 0.3 | +29.0 / +630.6 | **+9.3 / +28.7** |

**Length extrapolation.** Trained at 8,192 bytes, MambaByte extrapolates at least 4× (Figure 3) and
up to 64× with a sliding window (Figure 5) without degradation, whereas "limited by the position
embeddings, Transformer models don't extrapolate beyond the training length". MambaByte extrapolates
slightly better than subword Mamba, attributed to MambaByte seeing 4× longer sequences at training
for the same byte budget.

**Generation.** Generating 8,192 bytes on one A100-80GB: MambaByte-972M **29 s** vs
MegaByte-1.3B+218M **265 s** (same hardware, open-source MegaByte) — 2.6× faster parameter-matched
vs the 93 s figure quoted from Yu et al. With speculative subword drafting (Mamba-110M drafter,
top-3 acceptance, 3 subwords per iteration): 2.6× relative speedup, log-odds ratio 0.89 vs 0.10 for
substituting subword Mamba-1.03B outright.

### 1.2 Independent critique of MambaByte's control — SpaceByte

Slagle, arXiv 2404.14408 (2024). https://arxiv.org/abs/2404.14408

SpaceByte states that MambaByte's one subword control was matched on parameters and training data
(14 epochs of PG-19) **but not on compute**, and that "MambaByte was trained using roughly four
times as much compute than Mamba". Under SpaceByte's own compute- and inference-FLOP-controlled
protocol, a Transformer-family byte model matches MambaByte on PG-19:
**SpaceByte-793M+184M 0.918 BPB vs MambaByte-353M 0.930 BPB** at comparable inference FLOPs/byte
(~728M), with arXiv 0.663 vs 0.663 and Github 0.411 vs 0.396. SpaceByte also reports that
"byte-level Transformer and MegaByte models can require roughly 10 times more training FLOPs to
achieve the same performance as a subword-level Transformer" — i.e. the granularity penalty for
Transformers is large and measurable, which is the other half of the interaction.

### 1.3 Studies that cross architecture with MORE THAN ONE granularity

**(a) Genomics — the clearest reported reversal.**
Lindsey, Pershing, Habib, Stephens, Blaschke, Sundar, "A Comparison of Tokenization Impact in
Attention Based and State Space Genomic Language Models", bioRxiv 10.1101/2024.09.09.612081 (2024).
https://www.biorxiv.org/content/10.1101/2024.09.09.612081v1.full-text

- Design: 3 attention gLMs (Nucleotide Transformer, DNABERT1, DNABERT2) + 3 state-space gLMs
  (HyenaDNA, Mamba, Caduceus) + 2 baselines (3-layer CNN, pretrained GPT-2), on 45 downstream tasks
  from the Genomic Benchmark, Nucleotide Transformer Tasks, and GUE. Each fine-tuning task run 10×.
- Controlled ablation: 4-layer and 8-layer Mamba trained with **both character and BPE tokenization**,
  "leaving all other model parameters fixed", vocabulary capped at 4,096.
- **Reported interaction, verbatim from the discussion:**
  > "When using attention-based models, tokenization methods that compress the input … are
  > preferred. In state-space models, where a limited context window is not a concern, our
  > experiments indicate that character-based tokenization are the best choice for all genomic
  > language tasks except epigenetic mark prediction."
- Direction of the exception: on the Nucleotide Transformer tasks, Mamba-char wins splice site
  prediction, Mamba-BPE wins epigenetic marks. On GUE, Mamba-char wins transcription factors and
  promoters and "there are no tasks where the Mamba-BPE model has consistently better results".
  On the Genomic Benchmark "there is minimal difference".
- They also note the prior single-architecture conclusions disagreed with each other: DNABERT2
  (transformer) found BPE ≫ k-mer; HyenaDNA (SSM) found k-mer "degraded the results". Neither
  controlled input segment length.
- Caveat on strength of evidence: the architecture axis uses **published checkpoints with different
  pretraining data**, so only the tokenizer axis is a clean within-architecture control. The abstract's
  claim is conditioned on controlling input sequence length.

**(b) Chemistry — a genuine full crossed matrix.**
Wang, Sultan, Volkamer, Klakow, "Chemical Language Models for Natural Products: A State-Space Model
Approach", arXiv 2602.13958 (2026). https://arxiv.org/html/2602.13958v1

- Design: **3 architectures (Mamba, Mamba-2, GPT-2) × 8 tokenizers × 2 data splits = 48 models
  pre-trained from scratch** on 1,030,273 natural-product SMILES. Tokenizers span the exact
  granularity ladder of interest: character-level (vocab 48), Atom-in-SMILES / atom-level
  (vocab 1,023), a general PubChem BPE (7,924), and NP-specific BPE at vocab 60, 100, 1,000,
  7,924 and 30k. Each model generates 100k SMILES; property prediction via repeated 5×5-fold CV.
- Architecture marginals (validity / uniqueness / novelty, averaged over tokenizers and splits):
  **Mamba 76.27 / 74.77 / 62.23 · GPT 75.23 / 73.50 / 63.47 · Mamba-2 74.60 / 73.05 / 60.35**.
- Tokenizer marginals (same metrics, averaged over architectures): **AIS 81.26 / 79.99 / 69.63** ·
  NPBPE1000 80.23 / 78.25 / 62.68 · NPBPE100 80.81 / 79.45 / 67.46 · NPBPE60 79.69 / 78.65 / 68.06 ·
  BPE 78.44 / 77.34 / 65.99 · NPBPE7924 68.77 / 66.32 / 52.16 · **NPBPE30k 53.77 / 51.19 / 41.58**.
- **The headline comparison of effect sizes**, from the appendix MCSim analysis:
  > "differences between tokenizers are more pronounced than between models"
  The tokenizer span is ~28 points of validity (81.26 → 53.77); the architecture span is ~1.7
  points (76.27 → 74.60).
- **Orderings that do flip**, as reported:
  - By *metric*: Mamba wins validity and uniqueness by 1.4%; GPT wins novelty by ~2% and produces
    1–2k more novel scaffolds.
  - By *target*: AIS is best for molecules, but "tokenizers with more fine-grained tokens, like
    character-level and NPBPE60, however, yielded more unique and novel scaffolds … contrasting
    with the earlier observation that the AIS tokenizer is more effective for molecules."
    Scaffold marginals: character-level 38,821 unique / 29,446 novel; NPBPE60 39,480 / 30,204;
    AIS 36,160 / 25,781; NPBPE30k 25,577 / 19,472.
  - By *split*: Mamba/Mamba-2 beat GPT by 0.02–0.04 MCC under random split; under scaffold split
    "all models perform comparably".
- **Mechanistic link to granularity.** Mamba makes "6.15% and 3.73% fewer long-range syntax errors
  than GPT and Mamba-2"; syntax errors are ~80% of rejected molecules, roughly half of them
  unclosed rings. Larger-vocabulary tokenizers "reduce the average tokenized sequence length by
  approximately 3- to 6-fold (from around 70 tokens with a character-level tokenizer to 10–20
  tokens)" and raise training speed 2–4× — so granularity and sequence length are being varied
  together, and the SSM's advantage is specifically on the long-range error class.
- **Architecture-dependent hyperparameters.** "hyperparameter search also consistently leads to
  Mamba and Mamba-2 models adopting deeper architectures (with more stacked Mamba blocks) than GPT
  models". This is an interaction at the level of the search space itself.
- The paper's own framing: "Specific model-tokenizer pairs, like M1+NPBPE100 and GPT+NPBPE100,
  balance outstanding overall performance and speed, with over 50% faster training and up to 12
  times faster inference than the worst-performing ones … These results underscore the importance
  of aligning model-tokenizer choices with specific goals for optimal outcomes."
- **What it does not do:** it reports marginals and selected cells, not a formal interaction test,
  and it does not state whether the architecture ranking changes *within* a fixed tokenizer. The
  48-cell table exists (Appendix Figures 19–20) but the per-cell architecture ordering is not
  summarised in the text.

**(c) Near-miss — OCR.** Agbeti-Messan et al., "A Benchmark of State-Space Models vs. Transformers
and BiLSTM-based Models for Historical Newspaper OCR", arXiv 2604.00725 (2026).
https://arxiv.org/html/2604.00725v2 — advertises line vs paragraph granularity and char vs BPE, but
only one of six models (DANIEL) uses BPE, so tokenization is confounded with architecture. Mamba-AR
1.83% CER / 53.3 ms vs DAN 1.83% CER / 156.5 ms on Antiqua lines; DANIEL (BPE) 2.37% Antiqua but
6.18% Fraktur, attributed to "BPE tokenizer vocabulary mismatch". Useful as an example of exactly
the confound the researcher is worried about, in reverse.

### 1.4 The negative finding, stated plainly

**No study was found that crosses Transformer against an SSM at more than two tokenization
granularities in natural language under a single controlled protocol, and reports whether the
architecture ordering changes.** The two studies that genuinely cross the two axes are both outside
NLP — genomics (Lindsey et al. 2024, 2 granularities, architecture axis via public checkpoints) and
chemistry (Wang et al. 2026, 8 granularities, 3 architectures, full matrix but no interaction test).
MambaByte itself is a single-architecture-family study with a byte-vs-subword contrast inside Mamba
plus Transformer baselines at byte level only; it does not run a subword Transformer and a subword
Mamba and a byte Transformer and a byte MambaByte under one budget. The nearest thing to a 2×2 is
assembled across two papers (MambaByte 2024 + SpaceByte 2024) with incompatible compute controls.

---

## SECTION 2 — SSM weaknesses that are granularity-dependent

### 2.1 Jelassi, Brandfonbrener, Kakade, Malach — "Repeat After Me"

ICML 2024, arXiv 2402.01032. https://arxiv.org/abs/2402.01032 · https://arxiv.org/html/2402.01032v2

**Theory.** For a dictionary of size D and a length-L uniform copy distribution:

- *Theorem 2.3*: for every n there is a depth-2 Transformer of dimension O(n log D) whose copy error
  is below p_n-gram(D_L) for all 2n ≤ L ≤ D^n.
- *Lemma 2.4*: p_n-gram(D_L) < L² D^(−n), so error decays exponentially in n.
- *Corollary 2.5*: a depth-2 Transformer of dimension O(log(L/ε)·log D) copies with error < ε —
  **parameters logarithmic in sequence length**.
- *Theorem 2.7*: **any** GSSM with state space S has copy error > 1 − |S|/D^L.
- *Corollary 2.8*: any GSSM with mem(S) < L log D − 1 has error > 1/2.

The length dependence is the entire content of the lower bound: for a **fixed** state, error → 1 as
L grows. This is the formal statement of "the degradation is worse for long sequences".

**Synthetics (Mamba and Transformers both ≈160M params; LSTM ≈40M, "the largest LSTM we managed to
train").**

- *Learning efficiency*, train on strings ≤300, evaluate string-level accuracy at 300:
  "the transformers need **100x less samples** than the best GSSMs to learn the copy task".
  The LSTM never learns at length 300; the authors show it learns at shorter lengths and conclude
  "string length is the bottleneck".
- *Length generalization*, train ≤50, test to 1000: all models solve in-distribution; the GSSMs'
  (LSTM and Mamba) "performance … drops to zero almost immediately when increasing the input
  length", Transformer performance "decays much more gradually", and Hard-ALiBi achieves "almost
  perfect length generalization up to sequences of length 1000".
- *Mechanism*: accuracy holds for duplicated n-grams up to n ≤ 4 and drops from n = 5, matching a
  5-gram retrieval algorithm.

**The important nuance — the deficit is task-shaped, not blanket.** Two variants of n-gram lookup,
both trained on strings ≤30:

- *Suffix key* (key given after the sequence; requires storing the whole context): Transformers
  degrade only slightly out to length 100; "GSSMs, however, perform poorly beyond their training
  distribution".
- *Prefix key* (key given first; can be matched on the fly): "GSSMs achieve perfect
  length-generalization on this variant" and **outperform the NoPE and ALiBi Transformers**
  (though not Hard-ALiBi).

The authors' own summary: "GSSMs seem to be memory limited, but can be effective when the tasks only
require a summary of the inputs rather than storing the entire context."

**Pretrained models (Pythia 410M/1.4B/2.8B vs Mamba 360M/1.4B/2.8B, both on the Pile, same
tokenizer; Mamba has the *lower* perplexity).**

- Copying C4 strings: "Even the smallest transformer model dramatically outperforms the largest GSSM."
- Shuffled word order (less compressible): both drop, GSSMs more; "the largest GSSM now gets zero
  accuracy on strings of length 300".
- Phone-book lookup: "even the smallest transformer (410M parameters) outperforms the largest GSSMs
  (2.8B parameters) when the phone-book size is long enough (L ≥ 70)".
- SQuAD QA (F1, binned by paragraph length, 2.8B models, 50 questions per bin): comparable on short
  paragraphs; "the performance of Mamba degrades more quickly with the paragraph length".

### 2.2 Zoology — associative recall, and the width-vs-length theorem

Arora, Eyuboglu, Timalsina, Johnson, Poli, Zou, Rudra, Ré, arXiv 2312.04927 (2023), ICLR 2024.
https://arxiv.org/abs/2312.04927

- 17 models pretrained across 4 scales (70M–1.4B) and 5 architectures on identical data and
  infrastructure. SoTA gated-convolution architectures "still underperform attention by up to
  **2.1 perplexity points** on the Pile".
- Errors on "AR hits" account for **82% of the perplexity gap** on average while being only **6.4%
  of tokens**. A **70M attention model outperforms a 1.4B gated-convolution model** (20× larger) on
  associative recall. The gap persists at 7B comparing RWKV to Llama-2.
- **The theorem that makes this granularity-relevant.** Theorem 4.4: for BaseConv — which provably
  simulates architectures built from gating and convolution (H3, Hyena, RWKV, RetNet) — the model
  dimension required to solve multi-query associative recall **grows with the input sequence
  length**, whereas (Proposition 4.3) attention solves MQAR with model dimension **independent of
  sequence length**. A finer tokenization raises L for the same content, so it moves a
  fixed-width recurrent model down this curve while leaving attention unaffected.
- Prior single-architecture evidence was misleading: earlier synthetic AR formulations used one
  query at a fixed position from a small vocabulary (|V| < 50), on which gated convolutions score
  perfectly; MQAR (multiple recalls, varying positions, vocabulary larger than model dimension)
  reproduces the real gap.
- Hybrids with input-dependent sparse attention close **97.4%** of the gap to attention.

### 2.3 Based — the recall–throughput (i.e. recall–state-size) tradeoff

Arora, Eyuboglu, Zhang, Timalsina, Alberti, Zinsley, Zou, Rudra, Ré, arXiv 2402.18668 (2024), ICML 2024.
https://arxiv.org/abs/2402.18668

- "we identify a key tradeoff between a model's **state size** and **recall ability**" — the tradeoff
  is reported to "hold across architecture classes", established by varying hyperparameters that set
  the size of the recurrent state during generation.
- Attention beats Mamba on real-world recall-intensive tasks by **32.2 accuracy points** (Table 1).
- Based (linear attention + 64–128 token sliding window) matches Mamba in perplexity and beats it on
  recall-intensive tasks by **10.36 accuracy points**, at up to 1.3B params.
- Linear attention alone "struggles to solve associative recall"; sliding-window recall range is
  "limited by the width of the windows", and widening the window grows the recurrent state linearly.
- 24× higher generation throughput than FlashAttention-2 at 1,024 generated tokens, 1.3B params.

### 2.4 Just Read Twice — recall hardness is a communication-complexity problem in L

Arora, Timalsina, Singhal, Spector, Eyuboglu, Zhao, Rao, Rudra, Ré, arXiv 2407.05483 (2024).
https://arxiv.org/abs/2407.05483

- "a **2.8Bn** parameter Mamba LM trained on **300Bn** tokens of the Pile underperforms a **1.3Bn**
  param (2.2× smaller) Transformer LM trained on **50Bn** tokens (6× fewer tokens) by **5 points**",
  averaged over recall-intensive ICL tasks. Scale does not buy the gap back.
- The hardness of in-context recall is reduced to **set disjointness**; the recurrent memory needed
  to solve it changes with the order in which the sets appear. Memory requirement is a function of
  how much must be held, i.e. of sequence length and ordering, not of semantic content.
- JRT-Prompt (repeat the context) gives **+11.0 ± 1.3 points** averaged over 16 recurrent LMs and
  6 ICL tasks, with 11.9× higher prefill throughput than FlashAttention-2 at 32k length, batch 16, H100.
- JRT-RNN recovers 99% of Transformer quality at 360M/30B tokens and 96% at 1.3B/50B.

### 2.5 Bounded state is also not a state-tracking advantage

Merrill, Petty, Sabharwal, "The Illusion of State in State-Space Models", ICML 2024,
arXiv 2404.08819. https://arxiv.org/abs/2404.08819

- S4, Mamba and related SSMs are shown to lie in **TC⁰**, "limited very similarly to transformers",
  and provably cannot solve permutation composition (S₅), chess move tracking in source-target
  notation, code evaluation, or long-narrative entity tracking.
- Experiments confirm: Transformers and these SSMs cannot learn to compose permutations with a fixed
  number of layers, "whereas RNNs can compose permutations with just a single layer".
- Relevant because it removes the obvious counter-hypothesis that a fine granularity would favour
  Mamba through superior sequential state tracking. It would not.

### 2.6 Synthesis: the mechanism by which granularity and architecture interact

Every documented SSM deficit is indexed by sequence length, not by content:

1. Jelassi Thm 2.7: error > 1 − |S|/D^L. Fixed |S|, growing L → error → 1.
2. Zoology Thm 4.4: required model dimension for MQAR grows with L for gated-conv/SSM-class models;
   constant in L for attention.
3. Based: state size ↔ recall is a Pareto frontier; the position on it is set by how much of the
   context must be retained.
4. Jelassi §3.5: GSSMs are fine when the task needs a running summary (prefix key) and fail when it
   needs the whole context (suffix key). Longer inputs push more tasks into the second category.

Finer tokenization of the same object (residue → monomer → character → atom) multiplies L for
constant information content. On that reading, fine granularity is a pure tax on a fixed state and
no tax at all on attention — the interaction is predicted, not merely possible.

**Counter-evidence that must be stated.** MambaByte's own argument is that `m × n_state × d ≫ L_ctx`,
so the state is not the binding constraint at 8k bytes; empirically subword Mamba and MambaByte
performed equally at matched parameters and bytes, and MambaByte extrapolated 4×–64× past its
training length. So the theoretical tax is not detectable at 8k-byte context on language-modelling
loss. It is detectable on copying, phone-book lookup and MQAR at those same lengths. Which of these
a molecular generation metric resembles is the open question.

---

## SECTION 3 — Hybrids and what they imply

| Model | Gap it cites | Fix | Scale / numbers |
|---|---|---|---|
| **Jamba** (AI21), arXiv 2403.19887 (2024), https://arxiv.org/abs/2403.19887 | Transformers: KV-cache memory, no summary state. SSMs: "more efficient to train than RNNs and … more capable at handling long distance relationships, but still lag behind the performance of comparably sized Transformer language models". Notes Gu & Dao found interleaving Mamba+attention "only slightly better than pure Mamba in terms of perplexity, with models up to 1.3B" — i.e. the gap is invisible in perplexity at small scale. | Interleave Transformer and Mamba layers at a chosen ratio, MoE on alternate MLP layers (16 experts, top-2). | 12B active / 52B total, fits one 80GB GPU at 8-bit even past 128K context; 256K context; comparable to Mixtral-8x7B and Llama-2 70B; 3× Mixtral's throughput at long context. Ablations up to 7B / 250B tokens. |
| **Zamba** (Zyphra), arXiv 2405.16712 (2024), https://arxiv.org/abs/2405.16712 | Cites Jelassi et al. 2024 directly: SSMs "do not fully match the expressivity and performance of transformers at scale", "with several works highlighting in particular in-context learning (ICL) weaknesses (Park et al., 2024; Grazzi et al., 2024)". | Mamba backbone + a **single shared attention module**, "obtaining the benefits of attention at minimal parameter cost". | 7B, 1T tokens; "best non-transformer model at this scale"; faster inference and much lower generation memory than comparable transformers. Cites Poli et al. 2024: ~¼ of layers as self-attention is optimal. |
| **Samba** (Microsoft), arXiv 2406.07522 (2024/2025), https://arxiv.org/abs/2406.07522 | "SSMs struggle with memory recall due to their recurrent nature (Arora et al., 2023)", and retrieval results "have further shown that SSMs are not as competitive as their attention-based counterparts" (Fu 2023, Wen 2024, Arora 2024). | Layer-wise interleave Mamba + SwiGLU + **Sliding Window Attention**: "Mamba layers capture the time-dependent semantics … while SWA fills in the gap modeling complex, non-recurrent dependencies." | 421M–3.8B, 3.2T tokens. 3.8B post-trained: MMLU **71.9**, HumanEval **62.8**, GSM8K **87.6**. Pretrained at 4K, extrapolates to **1M** zero-shot (256×). After 500 steps of 4K instruction tuning, **perfect** passkey recall at 256K, while "the fine-tuned SWA-based model simply cannot recall memories beyond 4K length". Closes the Phonebook (Jelassi et al. 2024) retrieval gap with full attention. 3.73× throughput at 128K prompts, 3.64× at 64K generation. |
| **Zoology hybrids**, arXiv 2312.04927 (2023) | Associative-recall gap, 82% of a 2.1-ppl deficit. | Minimal input-dependent sparse attention on exact-match repeated bigrams. | Closes **97.4%** of the gap to attention while staying sub-quadratic. |
| **Based**, arXiv 2402.18668 (2024) | Recall–state-size tradeoff; Mamba 32.2 points behind attention on recall-intensive tasks. | Linear attention (2nd-order Taylor softmax feature map) + 64–128 token sliding window. | +10.36 points over Mamba on recall-intensive tasks at ≤1.3B; 24× throughput over FA2. |

**What this implies about how a pure-Mamba result would be read.** The hybrid literature is a
documented, converging consensus that pure SSMs have a *specific* deficit — in-context recall,
copying and retrieval — which is **not** visible in perplexity (Jamba says interleaving is "only
slightly better than pure Mamba in terms of perplexity"; Jelassi's Mamba checkpoints have *lower*
Pile perplexity than the Pythia models they lose to). So: a pure-Mamba result is representative of
pure-SSM behaviour and would be read as such. It would be read as *idiosyncratic* precisely to the
degree that the evaluated metric is recall-shaped and the sequences are long — which is also the
regime where fine-grained tokenization lives.

---

## SECTION 4 — Mamba in chemistry

### 4.1 SATURN (Guo & Schwaller) — the exact architecture comparison

Preprint: arXiv 2405.17066 (2024), https://arxiv.org/abs/2405.17066 · https://arxiv.org/html/2405.17066
NeurIPS 2024 AI4Mat workshop spotlight, https://openreview.net/forum?id=Hbpzrh7JbN ·
Published: *Nature Machine Intelligence*, 2026, doi 10.1038/s42256-026-01200-4,
https://link.springer.com/article/10.1038/s42256-026-01200-4 (carries "Extended Data Table 1
Sample efficiency across architectures").

- **Claim of novelty**: "the first application of the Mamba architecture for generative molecular
  design", and "the first application of Mamba … specifically for goal-directed generation with
  reinforcement learning".
- **Architectures and sizes**: RNN (LSTM) **5.8M**, decoder transformer **6.3M**, Mamba **5.2M**.
  Backbone of the Augmented Memory algorithm (experience replay + SMILES augmentation) on top of
  REINVENT. ">5,000 experiments" across the grid.
- **The replicate-success figure, verbatim from the paper:**
  > "Mamba with 10 augmentation rounds successfully generates 100 molecules above the reward
  > threshold (OB 100 metric) in **10/10 replicates**, compared to only **5/10** and **4/10**
  > successful replicates for RNN and transformer, respectively (Table 1)."
- Other reported findings: "Across the Yield and OB metrics, Mamba consistently outperforms both the
  RNN and transformer backbones." "Increasing augmentation rounds decreases diversity and
  inconsistently improves Yield and OB for RNN and transformer. Mamba more consistently benefits
  from increasing augmentation rounds." Mamba at 10 rounds: IntDiv1 0.714 ± 0.035.
- **Why Mamba wins, per the authors — this is an interaction, not a main effect.** "Mamba (5.2M) and
  RNN (5.8M) have similar parameter counts but during pre-training, the former converges to a lower
  loss …, indicating a better match to the data distribution." Under RL, the average max conditional
  token probability during generation "approaches 1, and near collapses to a Dirac delta function
  (less so for RNN)", producing repeated SMILES. The paper's framing throughout is that Mamba
  "synergistically exploits" the experience-replay-plus-augmentation mechanism and that the three
  architectures respond differently to the same augmentation schedule. That is an
  architecture × training-pipeline interaction measured directly.
- **Scope relevant to a granularity study**: the comparison is run at **one tokenization** —
  SMILES, with SMILES randomisation/enumeration as data augmentation (a data-level, not
  granularity-level, manipulation). SATURN provides **no** evidence about how the RNN/transformer/
  Mamba ordering would change under a different molecular granularity.
- Downstream: Saturn outperforms 22 models (16 in the journal version) on MPO tasks and is
  sample-efficient enough to directly optimise DFT as a high-fidelity oracle.

### 4.2 Other chemical / molecular SSM work

- **Natural-product chemical language models with Mamba, Mamba-2 and GPT** — Wang, Sultan, Volkamer,
  Klakow, arXiv 2602.13958 (2026), https://arxiv.org/html/2602.13958v1. Full details in §1.3(b).
  **This is the only chemistry paper found that compares architectures at more than one tokenization**
  (8 tokenizers × 3 architectures × 2 splits = 48 models). Self-described as "the first extensive
  experimental comparison of selective state-space models (S6) and transformers in NP-focused tasks".
  Its summary sentence: "for information-dense sequences, selective SSMs tend to preserve structural
  validity and long-range consistency, whereas self-attention favors recombination and novelty."
- **S4 for chemical language modelling** — Özçelik et al. (2024), cited in the above as having shown
  S4 "capture[s] long-distance dependencies and reduce[s] SMILES design errors more reliably than
  GPT". Pre-dates the Mamba work; same qualitative direction.
- **Materials** — "A Mamba-based foundation model for materials", *npj Artificial Intelligence* (2025),
  https://www.nature.com/articles/s44387-025-00009-7.
- **Retrosynthesis** — "A Comparison of Selective State Space Models and Transformers for Single-Step
  Retrosynthetic Reaction Prediction", BSc thesis, University of Hamburg (2024),
  https://www.inf.uni-hamburg.de/en/inst/ab/lt/teaching/theses/completed-theses/2024-ba-roth.pdf.
  PDF is image-based; contents not machine-extractable in this session. Title indicates a single
  tokenization.
- **Protein / biology** — Mamba protein language models (Sgarbossa et al. 2024) and Caduceus
  (Schiff et al. 2024, bi-directional Mamba for DNA) are the adjacent precedents; Caduceus is one of
  the state-space gLMs in the Lindsey et al. benchmark.
- **Tokenization-only chemistry studies (one architecture)**: "Comparing SMILES and SELFIES
  tokenization for enhanced chemical language modeling" (2024),
  https://pmc.ncbi.nlm.nih.gov/articles/PMC11499904/ — BPE vs a new Atom Pair Encoding, SMILES vs
  SELFIES, all inside BERT-based models; APE+SMILES beats BPE on HIV, Tox21 and BBBP ROC-AUC. No
  architecture axis.

---

## SECTION 5 — How papers scope a single-architecture result

### 5.1 What limitation statements actually say

The striking pattern: tokenizer papers scope their conclusions to **model size** and **language**
almost universally, and to **architecture** almost never.

- **Ali et al., "Tokenizer Choice For LLM Training: Negligible or Crucial?"**, arXiv 2310.08754
  (2023), NAACL 2024 Findings. https://ar5iv.labs.arxiv.org/html/2310.08754 — Limitations §8 names
  three things: no seed repetition (because "training all 2× models only once required 59,000 GPU
  hours"), "we did not investigate whether the results obtained could be extrapolated to larger
  model sizes", and no few-shot evaluation. **Architecture is not named as a limitation.**
- **"Beyond Text Compression: Evaluating Tokenizers Across Scales"**, arXiv 2506.03101 (2025),
  ACL 2025. https://arxiv.org/html/2506.03101 — "Our study focuses on decoder-only models up to 2.7B
  parameters … we have not verified whether these trends hold for larger architectures." Also names
  language coverage, benchmark variance and seeds. The word "architectures" is used to mean *scales*.
- **"Unpacking Tokenization"**, ACL Findings 2024,
  https://aclanthology.org/2024.findings-acl.134.pdf — limitations are compute and non-English
  coverage only.
- **TokSuite**, arXiv 2512.20757 (2026), https://arxiv.org/pdf/2512.20757v1.pdf ·
  https://openreview.net/pdf?id=iExjy56t3o — the paper's whole premise is the confound the
  researcher is guarding against: "it would be fraught to try to compare the Qwen 3 and Llama 3
  tokenizers by studying the respective models because differences in training data, training
  duration, and **architectural details** make it difficult to attribute performance differences
  specifically to tokenization", hence "reliable comparison can only be made through models that are
  completely identical apart from the tokenizer used". They release 14 LMs with "identical
  initialization, architecture, and training data composition, varying only in the tokenizer used".
  Their own limitations section names languages, domains and scale — **not** the single architecture.
- **TokLens**, ACL SRW 2026, https://aclanthology.org/2026.acl-srw.18.pdf — the one source found
  that names the limit explicitly: controlled tokenizer-only experiments (citing Ali et al. 2024 and
  Altıntaş et al. 2025) "provide stronger causal evidence but require retraining from scratch and
  **are limited to a single architecture each**."
- **Lindsey et al. (2024)** scope to vocabulary: "The results from this study were limited to
  vocabularies of size 4096, but larger token vocabularies may have different properties."
- **Wang et al. (2026)** scope to domain: "Future work should examine whether these insights extend
  beyond the domain of NPs, while also aiming to … investigate model size, architecture, and
  training data volume trade-offs."

So the conventional scope statement is: *size, language, domain, vocabulary, seeds*. The
architecture limit is normally left implicit, and when a paper does cross architectures it tends to
advertise that as a contribution rather than treat the single-architecture case as a defect.

### 5.2 Do reviewers accept it?

Only indirect evidence was obtainable. MambaByte (OpenReview X1xNsuKssb, 2024), Jelassi et al.
(ICML 2024), Zoology (ICLR 2024), Based (ICML 2024), Merrill et al. (ICML 2024) and SATURN
(*Nature Machine Intelligence*, 2026) all published with single-architecture or single-tokenization
framing in their headline claims. OpenReview review text was behind a verification wall in this
session, so **no direct reviewer quotation could be obtained** and no claim is made here about what
reviewers said. What *is* visible is peer critique in subsequent papers (§5.3), which is where the
scoping gets tested in practice.

### 5.3 Cases where a single-architecture (or single-protocol) conclusion was overturned

1. **Long Range Arena — the strongest precedent.** Tay et al. (2020) concluded from randomly
   initialised models that Transformers are poor at long-range dependencies (LRA average **53.66**,
   chance on PathX) and that SSMs are dramatically better. Amos, Berant & Gupta,
   "Never Train from Scratch", **ICLR 2024**, arXiv 2310.02980, https://arxiv.org/abs/2310.02980,
   showed that "random initialization leads to **gross overestimation of the differences between
   architectures**" and that with self-pretraining on the downstream task data alone, "vanilla
   Transformers … match the performance of S4 on Long Range Arena" — a >30% mean absolute improvement
   for Transformers — and the best SSM result on PathX-256 improves by **20 absolute points
   (67 → 87)**. They further show the hand-crafted structured parameterisations of S4 "become mostly
   redundant in the presence of data-driven initialization". An architecture-ranking conclusion was
   reversed by a *training protocol* change, with no new architecture involved.
2. **Synthetic associative recall.** Prior work established that gated convolutions "can perfectly
   solve synthetic tests for AR capability", which was taken as evidence of parity with attention.
   Zoology (arXiv 2312.04927, 2023) showed those synthetics used one query at a fixed position from
   a <50-token vocabulary, and that under multi-query AR a **70M attention model beats a 1.4B
   Hyena** — the parity conclusion did not survive a change of benchmark.
3. **Byte-level patching.** MegaByte (Yu et al. 2023) concluded, from Transformers, that byte-level
   modelling requires fixed-size patch compression. MambaByte's Figure 1 reports that "patching can
   also lower the model performance compared to the standard Transformer" and that an
   un-compressed SSM beats patched Transformers at matched FLOPs. The necessity of patching was an
   architecture-local conclusion.
4. **…and then partially back again.** SpaceByte (arXiv 2404.14408, 2024) showed MambaByte's one
   byte-vs-subword control was not compute-matched ("roughly four times as much compute than
   Mamba"), and that a Transformer with dynamic, word-aligned patching reaches **0.918 BPB on PG-19
   vs MambaByte's 0.930** at comparable inference FLOPs/byte. A single-architecture conclusion about
   granularity was overturned by a different architecture under a tighter control.
5. **Genomic tokenization.** DNABERT2 (transformer) concluded BPE ≫ k-mer; HyenaDNA (SSM) concluded
   k-mer tokenization "degraded the results" for state-space models. Lindsey et al. (2024) found
   character tokenization best for Mamba on every task category except epigenetic marks while
   attention models prefer compressing tokenizers. Neither original conclusion generalised across
   architectures, and both were drawn without controlling input segment length.
6. **Perplexity parity.** Gu & Dao's Mamba was reported as matching or exceeding Transformers.
   Jelassi et al. (2024) then showed Mamba checkpoints with *lower* Pile perplexity than Pythia are
   "dramatically" worse at copying and retrieval. Equivalence on the training metric did not transfer
   to the capability.

---

## SECTION 6 — Cost of a second architecture as a robustness check

### 6.1 Papers that ran the full matrix

- **Wang et al. 2026 (chemistry NP CLMs)**, arXiv 2602.13958 — **full 3 architectures × 8 tokenizers
  × 2 split protocols = 48 pretrained models**, plus 100k generated SMILES per model (4.8M molecules
  total) and repeated 5×5-fold cross-validation with hyperparameter tuning and early stopping for
  each of three downstream tasks. Reported cost envelope: parameters **8.96M–69.7M** per model,
  A100 GPUs, per-model training time "ranges from hours to a week". The full matrix was affordable
  only because the models are tiny (<70M) and the corpus is ~1M SMILES. This is the closest existing
  precedent for a chemistry tokenization study that also varies architecture.
- **Zoology**, arXiv 2312.04927 (2023) — **17** models across 4 scales (70M–1.4B) and 5
  architectures on identical data and infrastructure. Note: 4×5 = 20, so even this is an incomplete
  grid; three cells were dropped.

### 6.2 Reduced designs actually used

| Paper | Full axis | Reduced axis, concretely |
|---|---|---|
| **Lindsey et al. 2024** (genomics), bioRxiv 2024.09.09.612081 | 45 downstream tasks × 8 models × 10 runs each | The **tokenizer ablation was retrained on one architecture only**: Mamba at 4 and 8 layers × {char, BPE} = **4 models**, "leaving all other model parameters fixed". The architecture axis used **published checkpoints** (NT, DNABERT1, DNABERT2, HyenaDNA, Mamba, Caduceus) with different pretraining data. Vocabulary fixed at 4,096 throughout. "A small hyperparameter search was performed for the state space models, since they are extremely sensitive to these parameters." |
| **Jelassi et al. 2024** (ICML), arXiv 2402.01032 | Synthetics: Transformer variants (RoPE/NoPE/ALiBi/Hard-ALiBi) vs Mamba vs LSTM, all trained from scratch | Parameters matched at ~160M for Mamba and Transformers but **LSTM capped at ~40M** ("the largest LSTM we managed to train"). The real-data arm **does not retrain anything** — it uses off-the-shelf Pythia/Mamba pairs chosen because they share the Pile and the tokenizer. Evaluation budget cut: 10 batches × 128 for synthetics, 10 × 64 for pretrained copying, and **50 questions** for SQuAD at **2.8B only**, because "smaller models were unable to achieve reasonable and consistent performance". |
| **MambaByte 2024**, arXiv 2401.13660 | 5 datasets × several byte architectures | **Two matched settings instead of a grid**: compute-matched *or* parameter-matched, never both at once. The 5-dataset BPB table runs at medium scale only (353M). The large-scale byte-vs-subword comparison runs on **PG-19 only**. One extra model (subword Mamba-1.03B) was retrained purely as a control "to control for the benefits of the Mamba architecture". Training bytes were allowed to differ (30B vs 80B; 150B vs 400B) — which is what SpaceByte later attacked. |
| **SATURN 2024/2026**, arXiv 2405.17066 | >5,000 experiments, 3 architectures × batch size × augmentation rounds, 10 replicates per cell | All at **~5–6M parameters** and a **1,000-oracle test budget** for the grid search before the real MPO tasks. The mechanistic analysis is reduced to **Mamba vs RNN only**: "Given Mamba's superior sample efficiency, we focus our analysis on comparing it to the RNN baseline in the remainder of this section (transformer results are provided in Appendix B.3)." One tokenization throughout. |
| **Based 2024**, arXiv 2402.18668 | Architecture sweep on synthetic MQAR with varied state size | Headline LM claims at **up to 1.3B** only; the Pareto-frontier sweep is on synthetics, not pretrained LMs. |
| **Ali et al. 2023**, arXiv 2310.08754 | 24 models varying the tokenizer | **One architecture; one seed each.** Stated cost: "training all 2× models only once required **59,000 GPU hours**", explicitly given as the reason seeds were not repeated. |
| **TokSuite 2026**, arXiv 2512.20757 | 14 LMs, identical architecture / data / budget / initialization, varying only the tokenizer | **One architecture.** The design philosophy is to spend the whole budget on tokenizer coverage. |

### 6.3 What the cost evidence says

- The only hard GPU-hour figure located is Ali et al.: **59,000 GPU hours for 24 single-architecture
  models**. Doubling the architecture axis doubles pretraining cost; nobody in this set absorbed that
  at scale.
- The four observed ways of avoiding the doubling, in order of frequency:
  1. **Substitute public checkpoints for the architecture axis** and retrain only the tokenizer axis
     (Lindsey et al.; Jelassi et al.'s pretrained arm).
  2. **Shrink the models** until the full matrix fits (Wang et al.: 48 models, all <70M).
  3. **Run the second architecture on a subset of conditions** — one dataset, one scale, one metric
     (MambaByte's PG-19-only large-scale run; SATURN's appendix-only transformer analysis).
  4. **Match on one budget axis, not all of them** (MambaByte's compute-matched *or*
     parameter-matched settings; its unmatched training bytes).
- Where a full matrix exists (Wang et al. 2026), the finding it buys is precisely the one a
  single-architecture study cannot assert: that the tokenizer main effect (≈28 points of validity)
  dwarfs the architecture main effect (≈1.7 points), and that orderings flip by metric and by target.

---

## GAPS

- **No controlled Transformer × SSM × ≥3 granularity study exists in any domain.** The two crossed
  studies are genomics (2 granularities, architecture axis via public checkpoints with different
  pretraining data) and chemistry natural products (8 granularities, 3 architectures, full matrix but
  no interaction test and no per-cell ordering reported).
- **No paper reports a formal architecture × granularity interaction test.** Wang et al. 2026 report
  marginals and an MCSim pairwise significance heatmap; nobody reports an ANOVA-style interaction
  term, so "the ordering changes" is currently an eyeball claim everywhere it is made.
- **The SSM recall deficit has never been measured at varying tokenization granularity on the same
  content.** Jelassi, Zoology, Based and Just Read Twice all vary sequence length by adding content,
  not by re-tokenizing fixed content — so the length-indexed bounds have not been tested in the one
  regime that a granularity study inhabits.
- **The one near-decisive MambaByte control is compute-unmatched.** The "subword Mamba ≈ MambaByte at
  matched parameters and training bytes" result — the single best evidence that granularity is
  architecture-neutral inside Mamba — used ~4× the compute on the byte side per SpaceByte (2024). No
  one has rerun it compute-matched.
- **No chemistry work connects the SSM recall literature to molecular metrics.** Wang et al. 2026
  measure long-range *syntax* errors (unclosed rings) and explicitly say "long-context modeling is
  not central to our tasks, as NP SMILES are relatively short (avg. <200 characters)". Whether
  validity, uniqueness, novelty or RL sample efficiency are recall-shaped tasks in the Zoology/Based
  sense is unestablished.
- **SATURN's Mamba advantage is confounded with the training algorithm, not isolated as an
  architecture effect.** Its own explanation is "synergistic" interaction with experience replay plus
  SMILES augmentation, and the three backbones respond differently to the same augmentation
  schedule — so even the one chemistry architecture comparison with clean replicate counts (10/10 vs
  5/10 vs 4/10) is not a pure architecture main effect.
