# Scale sweeps, ranking stability, embedding fraction, weight tying, reduced subsets, seeds

Research notes on what the ML literature treats as sufficient evidence that a ranking of arms is
not an artefact of model capacity. Every claim carries a URL and year. Quotes are short and marked.

Contents: Section 2 (ranking flips) first, because it decides whether a capacity sweep is mandatory;
then 1 (how many sizes), 4 (weight tying), 3 (embedding fraction), 5 (reduced subsets), 6 (seeds), GAPS.

---

## SECTION 2 — DO RANKINGS FLIP WITH SCALE?

**Short answer: yes, documented repeatedly, in several literatures, and at spans as small as one
order of magnitude in parameters. Ranking flips are not exotic.** Below, the cases, strongest first.

### 2.1 A complete ranking inversion in a controlled, seeded, factorial study

**Li 2026, "Do Post-Training Algorithms Actually Differ? A Controlled Study Across Model Scales
Uncovers Scale-Dependent Ranking Inversions"** — https://arxiv.org/abs/2603.19335 (2026)

This is the single best methodological template for the researcher's problem, because it anticipates
and rules out the confound. **Caveat: single-authored arXiv preprint, no citation record, not known
to be peer-reviewed.** Treat the design as a template to imitate, not as a citable authority; the
peer-reviewed load for "rankings flip" should rest on Tay et al. 2023 and Wei et al. 2022.

- Protocol: 8 algorithms × **4 model scales (0.5B, 1.5B, 3B, 7B)**, 3 evaluation domains, plus a
  20-variant DPO taxonomy run only at 1.5B with **5 seeds each (100 runs)**. ~240 training runs total.
- Result: at 1.5B the order is SFT (54.4%) > IPO (52.2%) > KTO (51.2%) > DPO (49.1%) > SimPO (38.7%)
  on GSM8K. At 7B the order is SimPO (85.8%) > DPO (83.9%) > IPO (81.4%) > KTO (80.2%) > SFT (76.4%).
  **SimPO goes from worst to best — a ~21 pp reversal relative to SFT.**
- Crucially for the researcher: they anticipated the objection that the flip was an artefact of a
  nuisance factor (LoRA, used at 7B) and ran a **2×2 factorial (3B/7B × full-FT/LoRA)** to isolate it.
  At 3B LoRA has no detectable effect (SFT differs by 0.31 pp; DPO gap 1.56 pp, inside the seed
  variance σ=2.01). They conclude scale, not LoRA, drives the inversion.
- They also show rankings invert **across tasks and across metrics**, not only across scale: the
  19.3 pp spread on GSM8K collapses to 0.54 pp on MATH (36×) and 0.47 pp on general benchmarks (41×),
  and on MATH the ordering itself changes (SGRPO 1st→4th, SimPO last→2nd).
- Their hierarchy of leverage: scale (~50 pp) ≫ paradigm (~10 pp) ≫ online/offline (~9 pp) ≫
  loss function (~1 pp).

The lesson the researcher should take: **a 2×2 factorial that crosses the arm of interest with the
suspected confound, at two sizes, is what a 2026 paper does to defend "this is not capacity".**

### 2.2 Architecture rankings explicitly flip with scale

**Tay et al. 2023, "Scaling Laws vs Model Architectures: How Does Inductive Bias Influence Scaling?"
(Findings of EMNLP 2023)** — https://arxiv.org/abs/2207.10551 (2022/2023),
https://aclanthology.org/2023.findings-emnlp.825/

- Ten+ architectures: Transformer, Evolved Transformer, Universal Transformer, Switch Transformer,
  Performer, Funnel Transformer, ALBERT, Mixture-of-Softmaxes Transformer, GLU-Transformer,
  Lightweight/Dynamic Convolutions, MLP-Mixer.
- **Five size buckets spanning 15M to ~40B parameters**: Tiny (15–27M), Small (60–96M),
  Base (210–324M), Large (738M–1.2B), XL (2.3–29.6B).
- Headline, in their words: **"the best performing model can fluctuate at different scales"**, and
  models that do well in one compute region are "not necessarily the best in another compute-region".
- Named reversals: **Evolved Transformer** does well against the standard Transformer at tiny-to-small
  scale on downstream tasks, but "this quickly changes when scaling the model up". **MoS-Transformer**
  outperforms vanilla Transformers at some regions and not others.

This is the paper to cite for the general claim that architecture ranking at one size does not
predict ranking at another.

### 2.3 Tokenizer rankings crossing with scale

**Schmidt et al. 2024, "Tokenization Is More Than Compression" (EMNLP 2024)** —
https://arxiv.org/abs/2402.18376 (2024)

- **64 models: 54 at 350M, 6 at 1.3B, 4 at 2.4B** (span ~7×), 18 tokenizer variants, 3 vocabulary
  sizes (32,768 / 40,960 / 49,152); the larger models used only the middle vocabulary.
- Finding: the relative performance of the tokenizers **"do vary by model size"**, with **"the
  prevalence of crossing lines"** in the per-size plots.
- So a tokenizer comparison produced rank crossings over a span of less than one order of magnitude.

**Limisiewicz et al. 2026, "Compute Optimal Tokenization" (FAIR / UW)** —
https://arxiv.org/abs/2605.01188 (2026), https://arxiv.org/html/2605.01188v1, https://co-tok.github.io

The largest tokenization sweep in the literature and a clean scale-dependent reversal:

- **988 BLT (latent-tokenized) models from 50M to 7B parameters**, plus **320 subword-tokenized
  models** spanning 4B to 1.1T training bytes; compute budgets **5×10^18 to 2×10^21 FLOPs**.
- The optimal compression rate T* (bytes per token) is **not fixed — it decreases with compute**:
  T* = 3.69 at C=10^20 FLOPs, T* = 3.33 at C=2×10^21 FLOPs. Fitted T* = T₀/C^δ with T₀=18.2, δ=0.035.
- The reversal: **models with 90% and 75% of the BPE vocabulary masked outperform the original BPE
  tokenizer in the largest runs**, "even though both spend the same compute in the embedding and
  de-embedding layers". At smaller budgets standard BPE is near-optimal. So the *winner changes with
  scale*, and the authors explicitly rule out an embedding-compute explanation for it.
- They also generalise Chinchilla: compute-optimal English models want ≈**60 bytes per parameter**,
  not ≈20 tokens per parameter, and the 20-tokens rule holds only at BPE's specific compression rate.
- Parameter accounting: when fitting the first-stage law they count **only latent-module parameters,
  excluding BLT encoder/decoder and excluding subword embedding parameters**, specifically to make
  cross-tokenizer comparison fair.

### 2.4 Inverse scaling, and inverse scaling reversing again

**McKenzie et al. 2023, "Inverse Scaling: When Bigger Isn't Better"** —
https://arxiv.org/abs/2306.09479 (2023)

- 11 datasets from a public contest. Four identified causes of inverse scaling: preference for
  repeating memorised sequences over following in-context instructions; imitation of undesirable
  training-data patterns; an easy *distractor* task embedded in the hard real task; and misleading
  few-shot demonstrations.
- The paper itself notes these tasks drove discovery of "U-shaped and inverted-U scaling trends,
  where an initial trend reverses".

**Wei et al. 2022, "Inverse scaling can become U-shaped"** — https://arxiv.org/abs/2211.02011 (2022)

The strongest quantified flip rate in the literature:

- Re-ran the 11 prize tasks on models up to **540B** (≈5× the compute of the 280B maximum used in the
  original prize evaluation).
- **Only 4 of 11 tasks remained inverse-scaling. 6 of 11 became U-shaped** (performance falls then
  rises again). 1 became positive-scaling. **7 of 11 trends (64%) reversed sign** once the sweep was
  extended by roughly one further scale step.
- 1-shot prompting and chain-of-thought mitigate the bad trends further — the direction of the trend
  depends on the evaluation protocol as well as the size.

Caveat for the researcher: these tasks were *selected* for an unusual trend, so 64% is not a base
rate. It is, however, a hard bound on how safely a monotone trend can be extrapolated one step out.

### 2.5 Apparent flips that are metric artefacts

**Schaeffer, Miranda, Koyejo 2023, "Are Emergent Abilities of Large Language Models a Mirage?"
(NeurIPS 2023)** — https://arxiv.org/abs/2304.15004 (2023)

- Core claim: nonlinear or discontinuous metrics "produce apparent emergent abilities, whereas linear
  or continuous metrics produce smooth, continuous predictable changes".
- Evidence in four prongs: (1) swapping exact-match/multiple-choice accuracy for a continuous metric
  (token edit distance, Brier score) makes the jump disappear; (2) three predictions about metric
  choice made and confirmed on the InstructGPT/GPT-3 family; (3) two predictions validated by
  meta-analysis of BIG-Bench emergence claims; (4) **they induced never-before-seen "emergent"
  abilities in vision models across diverse architectures purely by choosing a discontinuous metric**
  — a constructive demonstration that the phenomenon can be manufactured.
- They also note alleged emergence "evaporates ... with better statistics" — i.e. more seeds and more
  eval items remove some claimed discontinuities. This links directly to Section 6.

Consequence: if arm A and arm B are compared on a thresholded or exact-match metric, a crossover may
be a metric artefact. The defence is to report a continuous metric (loss, per-token likelihood, edit
distance, FCD) alongside any thresholded one.

### 2.6 Gap compression — the quieter failure mode

**Goldman et al. 2024, "Unpacking Tokenization"** — https://arxiv.org/abs/2403.06265 (2024)

- **Three sizes: 10M, 128M, 1B** (embedding layers excluded from the counts), six tokenizers (support
  data from 1M documents down to 1 document, plus character-level), 200k steps at batch size 512,
  span corruption on C4; a representative subset replicated on Turkish.
- Pearson correlation between compression and downstream *generation* performance:

  | Model size | Correlation range |
  |---|---|
  | 1B | −0.994 to −0.976 |
  | 128M | −0.988 to −0.949 |
  | 10M | −0.933 to −0.870 |

- Their abstract: correlations are "more pronounced for generation tasks (over classification) or for
  **smaller models** (over large ones)"; smaller models are "especially vulnerable to poor
  tokenizations", with the 10M model dropping more than the 1B under an inferior tokenizer.

**Over-Tokenized Transformer 2025** — https://arxiv.org/abs/2501.16975 (2025) shows the same shape in
the opposite direction: the loss gain from a 12.8M-entry input vocabulary is **0.082 at OLMoE-1.3B
and 0.076 at OLMoE-7B** — rank preserved, magnitude slightly compressed.

So the dominant observed pattern for a representation/tokenizer change is: **the better arm stays
better, but the margin narrows as capacity grows**, because a larger model compensates for a worse
input representation. For a claim "A beats B", a shrinking gap threatens external validity almost as
much as a sign reversal, since the gap can shrink into the seed noise.

### 2.7 A case where nothing flips

**NovoMolGen 2025** — https://arxiv.org/abs/2508.13408 (2025), full text
https://arxiv.org/html/2508.13408v1

- Llama-architecture molecular generators, 1.5B molecules from ZINC-22. Representations: **SMILES,
  SELFIES, SAFE, DeepSMILES**. Tokenizers: **Atomwise** and **BPE (vocab 500)**.
- **Three sizes: 32M, 157M, 300M** (9.4× span). Three random seeds per setting are stated for the PMO
  hyperparameter search; no uniform per-cell seed count is stated for the pretraining grid.
- Conclusion: "only modest performance differences overall"; SMILES+BPE gives "the most consistent
  (albeit marginal) gains" across FCD, PMO and docking, offered as a practical default.
- On scale: "scaling from 32M to 300M parameters yields only modest and inconsistent gains". No flip —
  but note that **the scale axis itself barely moved the metric**, which weakens the sweep as evidence
  either way. A sweep over a range where capacity does not bind cannot rule out a capacity confound.

### 2.8 The counterweight — evidence that small-scale conclusions often do transfer

**Khaddaj, Engstrom, Mądry 2025, "Small-to-Large Generalization: Data Influences Models Consistently
Across Scale"** — https://arxiv.org/abs/2505.16260 (2025). Finding: small- and large-scale model
predictions "(generally) do highly correlate across choice of training data", while also conceding
"changes in data do not influence smaller and larger models identically". Used to justify proxy
models for data attribution and dataset selection.

**"Most Transformer Modifications Still Do Not Transfer at 1–3B"** — https://arxiv.org/abs/2605.20798
(2026), a 2026 update to **Narang et al. 2021** (https://arxiv.org/abs/2102.11972). Of 20 post-2021
modifications, only two survive Bonferroni correction at 1.2B and **one of those two diverges at 3B**;
six significantly degrade the baseline and eight fall inside noise. The reproduction of Narang's
null result at larger scale is itself a scale-transfer finding.

### 2.9 How often, and under what conditions — plain summary

Ranking flips are **not rare, and they are concentrated in identifiable conditions**:

1. **Thresholded / discontinuous metrics** (Schaeffer 2023). Highest-risk condition.
2. **Tasks with a distractor sub-task or a memorisation-vs-instruction conflict** (McKenzie 2023).
3. **Arms close in effect size relative to seed noise** — then any "flip" may be sampling (Li 2026's
   DPO variants cluster inside 3.2 pp and none beats vanilla DPO after correction).
4. **When the arm changes the parameter budget.** A vocabulary/representation change moves embedding
   parameters, which are a large share of a small model and a negligible share of a large one
   (Section 3). This is the single most likely mechanical cause of a capacity-confounded ranking in a
   tokenizer study, and the one the researcher must address.
5. **When the small-scale regime is capability-limited** and the large-scale regime is not (Li 2026:
   the inversion "requires sufficient model capacity (≥7B)").

Conversely, the modal outcome for input-representation changes specifically (Goldman 2024,
Over-Tokenized 2025, SuperBPE 2025, NovoMolGen 2025) is **rank preserved, gap compressed**.

**Bottom line.** The literature does not support a default assumption of scale-invariant ranking. It
supports the narrower position that a ranking measured at one size is credible when (a) the swept
span crosses the regime where the suspected confound changes materially, (b) the metric is continuous
or reported alongside a continuous one, and (c) the gap at each size is large relative to the
measured seed noise.

---

## SECTION 1 — HOW MANY MODEL SIZES DO ABLATION PAPERS USE?

| Paper (year, URL) | # sizes | Smallest | Largest | Span | Seeds / cell | What was concluded |
|---|---|---|---|---|---|---|
| **Kaplan et al. 2020** https://arxiv.org/abs/2001.08361 | many (shape grid) | shape (2,128), ~10^3 non-emb params | >1B non-emb | "more than seven orders of magnitude" in N | not reported | Loss is a power law in non-embedding params N, data, compute; width/depth matter little within a wide range |
| **Hoffmann et al. 2022 (Chinchilla)** https://arxiv.org/abs/2203.15556 | 400+ models, 9 IsoFLOP budgets (6×10^18–3×10^21) | **70M** | **>16B** | ~230× | not reported; 80%-subsample bootstrap ×100 for CIs on exponents | Params and tokens should scale ~equally (~20 tokens/param); validated by training Chinchilla-70B |
| **Tay et al. 2023 (archs)** https://arxiv.org/abs/2207.10551 | 5 buckets × 10+ archs | **15M** | **~40B** | ~2700× | not reported | Architecture ranking **fluctuates** across scale; small-scale arch comparison does not predict large-scale |
| **Limisiewicz et al. 2026 (Compute Optimal Tokenization)** https://arxiv.org/abs/2605.01188 | 988 BLT models + 320 subword | **50M** | **7B** | 140× | not reported per cell; multi-restart BFGS + Hessian CIs on the fitted law | ≈60 bytes/param compute-optimal; optimal compression rate decreases with compute; 90%-masked BPE beats BPE at the largest budget |
| **Schmidt et al. 2024 (tokenization)** https://arxiv.org/abs/2402.18376 | 3 | **350M** (54 models) | **2.4B** (4 models) | ~7× | no seed repeats; averaged over 3 vocab sizes as a noise proxy; one-sided Wilcoxon signed-rank over 30 paired scores at 350M | Compression alone does not explain downstream performance; pre-tokenization matters; **rankings vary by model size** |
| **Goldman et al. 2024** https://arxiv.org/abs/2403.06265 | 3 | **10M** | **1B** | 100× | not reported | Compression correlates with downstream performance; correlation **stronger for smaller models** and for generation |
| **SuperBPE 2025** https://arxiv.org/abs/2503.13423 | ~6 configs | **680M** | **8B** (plus an 11B compute-matched variant) | ~12–16× | not reported | +4.0% absolute mean over 30 tasks, +8.2% MMLU, wins 25/30, 27% less inference compute; advantage holds "at every model size and every training budget tested" |
| **Li 2026 (post-training algos)** https://arxiv.org/abs/2603.19335 | 4 | **0.5B** | **7B** | 14× | **5 seeds** for the 20-variant DPO grid at 1.5B; 3 seeds for several headline cells | **Complete ranking inversion** between 1.5B and 7B; confirmed not-LoRA by a 2×2 factorial |
| **Ali et al. 2024 (tokenizer choice)** https://arxiv.org/abs/2310.08754 | **1** | 2.6B | 2.6B | — | **none** — explicitly "we did not investigate the effect of different random seeds ... due to the additional computational costs" | Tokenizer choice materially affects downstream performance and cost (up to +68% training cost for English-centric multilingual tokenizers); fertility/parity are not always predictive |
| **Dagan et al. 2024** https://arxiv.org/abs/2402.01035 | 2 | **1.5B** (GPT-2 XL, 2T tokens) | **7B** (Llama, 500B-token fine-tunes) | 4.7× | not reported | Vocab sizes 32k/64k/128k/256k; larger base models tolerate larger vocabularies because the softmax cost is amortised; >50B fine-tuning tokens allow retokenising a pretrained LLM |
| **MobileLLM 2024** https://arxiv.org/abs/2402.14905 | 2 (main sweep) | **125M** (9 architectures) | **350M** (10 architectures) | 2.8× | not reported | Deep-and-thin wins at sub-billion scale; embedding sharing + layer sharing pay for themselves |
| **NovoMolGen 2025** https://arxiv.org/abs/2508.13408 | 3 | **32M** | **300M** | 9.4× | 3 seeds for the PMO hyperparameter search; not uniform across the pretraining grid | SMILES+BPE a marginal, consistent default; 32M→300M gives "modest and inconsistent gains" |
| **Skinnider 2024 (NMI, SELFIES)** https://doi.org/10.1038/s42256-024-00821-x | **1 capacity** (3-layer LSTM, hidden 1024, embedding 128); robustness to "architecture" in Extended Data only | — | — | — | **10 seeds per condition; 180 models total** (90 SMILES, 90 SELFIES) | SMILES robustly beats SELFIES; invalid outputs act as a likelihood filter; causally confirmed by removing SELFIES valency constraints |
| **Narang et al. 2021 / 2026 update** https://arxiv.org/abs/2102.11972, https://arxiv.org/abs/2605.20798 | 1 main + 1 check | **1.2B** (20 modifications) | **3B** (robustness subset) | 2.5× | 3 seeds on the baseline only (for the noise floor); single seed elsewhere | Most Transformer modifications do not transfer; 2/20 survive Bonferroni at 1.2B and 1 of those diverges at 3B |
| **Press & Wolf 2017** https://arxiv.org/abs/1608.05859 | several configs | **2.65M** (small NNLM, tied) | **66M** (large NNLM) | ~25× | not reported | Weight tying reduces perplexity and halves NMT model size |
| **Inan et al. 2017** https://arxiv.org/abs/1611.01462 | **3** | 200 hidden units | 1500 hidden units | — | not reported | Tying + augmented loss improves perplexity at all three sizes; the relative benefit of the two components **changes with size** |

### What span is treated as credible

Reading the table as a norm rather than a rule:

- **Scaling-law papers** (Kaplan, Chinchilla, Tay, Limisiewicz) use **≥4–5 sizes spanning 100×–1000×**,
  with 400–1000 runs. Nobody reviews a scaling-law claim built on fewer.
- **Ablation papers that claim "across scale"** converge on **3 sizes spanning roughly one order of
  magnitude** (Schmidt ~7×, Goldman 100×, NovoMolGen 9.4×, SuperBPE ~12×, Li 14×). Three sizes over
  ~10× appears to be the de facto minimum for the phrase to pass.
- **Two sizes** (Dagan 1.5B+7B; Narang-update 1.2B+3B) is presented as *confirmation at a second
  point*, not as a sweep, and the papers word their claims accordingly.
- **One size** is publishable (Ali 2024, Skinnider 2024) **only with an explicit disclaimer**. Ali
  et al. state plainly: "we did not investigate whether the results obtained could be extrapolated to
  larger model sizes", and their only supporting argument is that the tokenizer in question has been
  used up to 65B — which they concede is indicative, not empirical.

**Minimum that appears to be treated as credible: three sizes spanning about one order of magnitude,
with the largest size at or above the regime where the suspected confound stops binding.** Below
that, the paper must either disclaim (Ali) or substitute seeds for sizes (Skinnider: 10 seeds × 180
models at one capacity, with capacity robustness relegated to Extended Data).

---

## SECTION 4 — WEIGHT TYING ABLATIONS

Short answer: **tying is almost never swept across sizes as a factor in its own right.** It is run
once at a reference configuration, or — most commonly in modern work — simply asserted as a default
with no ablation at all. The two exceptions below are both from 2017.

### Press & Wolf 2017 — tied vs untied at several configurations, one task family

https://arxiv.org/abs/1608.05859 (2016/2017). Datasets: PTB, text8, IMDB, BBC; NMT on WMT'14 EN→FR
and WMT'15 EN→DE. Reported tied-vs-untied pairs:

| Configuration | Params untied → tied | Test PPL untied → tied | Δ |
|---|---|---|---|
| Large NNLM + dropout (PTB) | 66M → 51M | 78.4 → 74.3 | **−4.1** |
| Large + Bayesian dropout | 51M | 75.2 → 73.2 | **−2.0** |
| Small NNLM, no dropout (PTB) | 4.65M → 2.65M | 114.5 → 112.4 | **−2.1** |
| Small + projection regularisation | 2.69M | 111.7 → 100.9 | **−10.8** |
| RHN + Bayesian dropout | 32M → 24M | 68.5 → 66.0 | **−2.5** |

NMT: tying reduces model size to **less than half** (~52% parameter reduction) without harming
performance. Note the ablation varies configuration (dropout regime, architecture) more than it
varies *scale*, and the effect size is not monotone in size (−10.8 at 2.69M vs −2.0 at 51M).

### Inan, Khosravi & Socher 2017 — tying ablated at three sizes

https://arxiv.org/abs/1611.01462 (2016/2017). PTB (10k vocab) and WikiText-2 (33,278 vocab); 2-layer
LSTMs at **200 / 650 / 1500 hidden units**. PTB test perplexity:

| Size | Baseline | +AL (augmented loss) | +RE (reused embeddings = tying) | +REAL (both) |
|---|---|---|---|---|
| Small (200) | 87.3 | 82.9 | 85.1 | 82.7 |
| Medium (650) | 77.7 | 74.7 | 73.9 | 73.2 |
| Large (1500) | 72.6 | 71.2 | 69.0 | 68.5 |

Note the scale interaction: **AL beats RE at the small size (82.9 vs 85.1) and RE beats AL at medium
and large** (73.9 vs 74.7; 69.0 vs 71.2). The paper states that RE "significantly outperforms AL for
larger networks" while augmented loss helps smaller models more. **This is itself a documented rank
flip between two components, at hidden size 200→650.** It is also the only paper in this set that
runs tying as a crossed factor at more than one size.

### MobileLLM 2024 — tying as a parameter-reallocation decision, not a quality decision

https://arxiv.org/abs/2402.14905 (2024). Vocab 32k.

- Sharing input/output embeddings at 125M saves ~16M parameters = **11.8% of total parameters**.
- Immediate accuracy cost: **−0.2 pp** on zero-shot reasoning.
- Reinvesting the saved parameters as two extra layers (30 → 32) yields **+0.4 pp net gain** while
  still ending with **10M fewer parameters** than the original 135M model.
- Framing: tying is worth it *because the freed parameters buy depth*, which is the right framing for
  a small model where embeddings are a large share of the budget.

### Over-Tokenized Transformer 2025 — tying partially retained, decoupled from the main lever

https://arxiv.org/abs/2501.16975 (2025). The standard 1-gram input embedding **remains tied to the
output layer**; the added 2-gram and 3-gram embeddings are computed separately via modulo hashing on
a tiled matrix. Input vocabulary is scaled from ~0.1M to **1.2M and 12.8M** entries (12× and 128×
baseline). Scaling law: **L = 2.6754 − 0.0256·log₁₀(m)**, i.e. **every 4× increase in input vocabulary
m reduces training loss by 0.015**. Models: OLMo2-151M / 400M / 1B dense (400B–1T tokens), OLMoE-1.3B
(260M active), OLMoE-7B (1.3B active), plus an in-house 400M-active/4B-total MoE. Throughput cost at
m=10^7 on FSDP is **under 5%**. There is **no tied-vs-untied ablation at matched parameter count** —
tying is a structural given, not a factor.

### Limisiewicz 2026 — a tying-adjacent control

https://arxiv.org/abs/2605.01188 (2026). Their masked-vocabulary comparison holds **the same compute
in the embedding and de-embedding layers** across arms, which is the strongest "embedding budget is
not the explanation" control found in this survey. It is a compute control, not a parameter control.

### Summary for Section 4

- **Swept across sizes as a factor: only Inan et al. 2017** (3 sizes), and there the tying-vs-loss
  ordering flips between the smallest and the middle size.
- **Run as a factor at one configuration: Press & Wolf 2017, MobileLLM 2024.**
- **Asserted, not ablated: Over-Tokenized 2025**, and effectively all modern LLM pretraining work
  (tying is a default below ~1B, untied above).
- **Effect sizes where reported:** PTB perplexity −2 to −11 (Press & Wolf, scale-dependent and
  non-monotone); −0.2 pp accuracy at 125M before reallocation, +0.4 pp after (MobileLLM);
  PTB 85.1→73.9→69.0 across three sizes with tying (Inan).
- **Nobody reports a tied-vs-untied comparison at matched total parameter count** in the modern
  LLM-scale literature. MobileLLM comes closest by reallocating the saved parameters into depth.

---

## SECTION 3 — EMBEDDING FRACTION AS THE RELEVANT VARIABLE

### The one paper that quantifies it directly

**MobileLLM 2024** — https://arxiv.org/abs/2402.14905 (2024). Vocabulary **32k**; embedding dimension
512 at 125M and 960 at 350M.

- "with an embedding dimension of 512 and a vocabulary size of 32k, the input and output embedding
  layers each comprise **16 million parameters**. Together, these embedding layers account for
  **more than 20% of the total parameters of a 125M-parameter model**."
- Contrast given in the same passage: embeddings are **3.7% of LLaMA-7B** and **0.7% of LLaMA-70B**.
- Sharing (tying) the two tables removes 16M = **11.8% of the total**.

That one sentence is the cleanest statement in the literature of the fact the researcher needs: the
embedding share falls by roughly **30× between 125M and 70B**, so any effect mediated by the embedding
budget must behave differently at the two ends of that range.

### The convention split on whether embeddings even count

- **Kaplan et al. 2020** (https://arxiv.org/abs/2001.08361) **excludes** them, and says why explicitly:
  "To observe these trends it is crucial to study performance as a function of N; if we instead use
  the total parameter count (including the embedding parameters) the trend is somewhat obscured."
  And: "When we include embedding parameters, performance appears to depend strongly on the number of
  layers in addition to the number of parameters." Excluding embeddings makes models of different
  depths "converge to a single trend". They add: "This suggests that the embedding matrix can be made
  smaller without impacting performance." Vocabulary 50,257.
- **Hoffmann et al. 2022 (Chinchilla)** (https://arxiv.org/abs/2203.15556) does the **opposite**:
  "We include all training FLOPs, including those contributed to by the embedding matrices... we also
  count embeddings matrices in the total parameter count." The NeurIPS 2024 reconciliation paper
  (https://proceedings.neurips.cc/paper_files/paper/2024/hash/b6341525cd84f3be0ef203e4d7cd8556-Abstract-Conference.html)
  finds this accounting difference is one of the causes of the Kaplan/Hoffmann discrepancy:
  accounting for the decoding layer's compute "shifts compute-optimal scaling toward a more constant
  token-to-parameter ratio".
- **Goldman et al. 2024** reports its 10M / 128M / 1B sizes with **embedding layers excluded**.
- **Limisiewicz et al. 2026** fits its first-stage law counting **only latent-module parameters**,
  excluding BLT encoder/decoder and subword embeddings, to compare tokenizers fairly.

So the field has two incompatible conventions, and at least one high-profile discrepancy traced
partly to the difference. Any tokenizer paper has to declare which it uses.

### Arithmetic the researcher can use

Embedding share for a decoder-only Transformer, untied, is

    f_emb = 2·V·d / (2·V·d + 12·L·d²)        (tied: drop one factor of 2 on V·d)

With V = 32k that gives, for typical shapes:

Computed (V = 32,768; non-embedding = 12·L·d², the standard approximation):

| Non-embedding params | d | L | f_emb untied | f_emb tied |
|---|---|---|---|---|
| 25M | 512 | 8 | **57.1%** | 40.0% |
| 85M | 768 | 12 | **37.2%** | 22.9% |
| 302M | 1024 | 24 | **18.2%** | 10.0% |
| 1.21B | 2048 | 24 | **10.0%** | 5.3% |
| 6.4B | 4096 | 32 | **4.0%** | 2.0% |
| 64B | 8192 | 80 | **0.8%** | 0.4% |

The 4.0% and 0.8% rows bracket MobileLLM's reported 3.7% (LLaMA-7B) and 0.7% (LLaMA-70B), which is a
sanity check on the formula — the small residual is LLaMA's V=32,000 and SwiGLU FFN.

Vocabulary scales the numerator linearly. At the 302M-non-embedding shape above, the untied embedding
share is **0.34% at V=500** (NovoMolGen's BPE vocab), **18.2% at V=32,768**, and **57.6% at V=200,000**
(SuperBPE's fixed vocab). The practical implication for the researcher: with a small chemical
vocabulary the embedding confound is negligible at every size and the sweep argument is weak; with a
word-level or super-word vocabulary it dominates the budget below ~1B and the sweep argument is the
whole ballgame.

### Is "sweep until the embedding fraction changes materially" an established criterion?

**No. Nobody in this survey states it as a criterion.** What exists is the raw material for it:

1. MobileLLM quantifies the fraction and shows it drives architecture choice at small scale (2024).
2. Kaplan gives the mechanistic reason it distorts trends and excludes embeddings on that basis (2020).
3. Chinchilla includes them and the NeurIPS 2024 reconciliation shows the accounting choice
   measurably moves the fitted exponent.
4. Limisiewicz 2026 is the only tokenizer paper that explicitly **holds embedding compute equal across
   arms** and then reports a ranking change anyway — the closest thing to a direct test that the
   ranking change is *not* an embedding-budget artefact.

So the criterion is defensible and well-supported by component facts, but it is a **synthesis the
researcher would be proposing, not a norm they can cite as established practice.** The honest framing
is: "we choose our smallest and largest sizes so that the untied embedding share falls from X% to Y%,
following the observation of Yin et al. (MobileLLM, 2024) that this share moves from >20% at 125M to
0.7% at 70B, and Kaplan et al.'s (2020) finding that embedding parameters obscure scaling trends."

---

## SECTION 5 — REDUCED-SUBSET SWEEPS

**Yes, running the full matrix at one reference configuration and a reduced subset elsewhere is
standard and generally accepted. Every multi-scale ablation paper in this survey does it.** What
varies is how explicitly it is justified and how the subset is chosen.

### Documented instances

| Paper | Full matrix at | Reduced subset at | How the subset was chosen |
|---|---|---|---|
| **Schmidt 2024** https://arxiv.org/abs/2402.18376 | 350M: 54 models, 18 variants × 3 vocab sizes | 1.3B: 6 models; 2.4B: 4 models, **middle vocab size only** | Not stated explicitly; they keep the vocab size that averaged best and carry the leading variants up |
| **Li 2026** https://arxiv.org/abs/2603.19335 | 1.5B: 20 DPO variants × 5 seeds (100 runs) | 0.5B / 3B / 7B: 8 headline algorithms only | The 20-variant taxonomy is a *within-family* question; the 8 algorithms are the *cross-family* question that scale could plausibly affect |
| **MobileLLM 2024** https://arxiv.org/abs/2402.14905 | 125M (9 archs) and 350M (10 archs) depth/width grid | larger sizes carry only the chosen configuration | Grid is at the sizes where the design question (deep vs thin) is live |
| **Narang-update 2026** https://arxiv.org/abs/2605.20798 | 1.2B: 20 modifications | 3B: "robustness check" on a subset | The subset is the modifications that survived significance at 1.2B — and one of them diverged at 3B |
| **Dagan 2024** https://arxiv.org/abs/2402.01035 | 1.5B: tokenizer design ablations across 32k/64k/128k/256k | 7B: three Llama fine-tunes | The 7B runs test the one practical question (can you retokenise a pretrained model) |
| **SuperBPE 2025** https://arxiv.org/abs/2503.13423 | 8B: full tokenizer comparison on 30 tasks | 680M / 1.9B (and compute-matched 910M / 2.5B) | Smaller sizes carry only the scaling question, not the full benchmark suite |
| **Resolving Discrepancies in Compute-Optimal Scaling (NeurIPS 2024)** https://proceedings.neurips.cc/paper_files/paper/2024/hash/b6341525cd84f3be0ef203e4d7cd8556-Abstract-Conference.html | Hyperparameter sweep (LR, batch size, β₂) at **5M–108M** | **Validation sweeps only at 220M and 901M** | Fit power laws for optimal LR and batch size on the small grid, extrapolate, then validate at two larger points |

### The justifications papers actually give

The NeurIPS 2024 compute-optimal-scaling paper gives the most explicit defence found, and it is worth
copying as a template because it names the threat and argues against it rather than waving at compute:

> "Due to limited compute budgets, our hyperparameter sweep only targeted the smaller models in our
> grid, and furthermore trained each model for only 20N steps... This raises the concern that the
> hyperparameters we chose unfairly favor models trained for that particular token-to-parameter
> ratio... We believe this is unlikely: at small scales (where hyperparameter tuning is crucial) our
> original set of hyperparameters favored higher token-to-parameter ratios..."

Note the structure: (1) state the reduction, (2) name the specific bias it could introduce, (3) give a
directional argument that the bias, if present, runs *against* the paper's conclusion. That third step
is what makes the reduction acceptable rather than merely confessed.

Ali et al. 2024's justification for not sweeping at all is weaker and reads as a limitation rather
than a defence: "the fact that this tokenizer has been used for training state-of-the-art models up to
65B might indicate that our results also transfer to larger model sizes" — appeal to usage, not
evidence (https://arxiv.org/abs/2310.08754, 2024).

### The validity argument for proxy-at-small-scale

Two 2025–2026 papers argue directly that a reduced/proxy design is sound when the quantity of interest
is a *region* rather than a point:

- **Khaddaj et al. 2025** https://arxiv.org/abs/2505.16260 — data influence correlates highly across
  scale, justifying proxy models for data attribution and dataset selection, while conceding the
  influences are not identical.
- **"Scaling Near-Optimal SFT–RL Annotation Budget Allocation"** https://arxiv.org/abs/2609.01573
  (2026) — introduces the **near-optimal region** (allocations within 2–10% of peak) rather than the
  optimum, shows it is wide, **widens with model scale**, and transfers from small proxy models to
  large targets, so "small proxy-model experiments suffice". Notably this paper cites Li 2026
  precisely because "post-training algorithm rankings can flip across scale" — it is arguing for proxy
  designs *in spite of* the flip literature, by weakening the claim from "the best" to "a near-best
  region".

That weakening is the clean reviewer-facing move available to the researcher: if the full matrix only
runs at one size, state the conclusion as a near-optimal set rather than a strict total order.

### Criticism

Direct reviewer-facing guidance on scale sweeps in ablations is thin (see GAPS). What exists is
case-level: the 2026 Narang update (https://arxiv.org/abs/2605.20798) shows a modification that passed
at 1.2B **diverging at 3B**, which is the concrete argument that carrying only the winners upward can
mislead. And the Pith referee reports surveyed (e.g. https://pith.science/paper/2510.22228, 2025;
https://pith.science/paper/2607.22577, 2026) consistently flag two things: **matched-compute claims
contradicted by the paper's own FLOP accounting**, and **three seeds reported with no error bars**.
Both are more commonly the stated objection than the number of sizes.

---

## SECTION 6 — SEEDS AND VARIANCE AT EACH CELL

### What papers actually do

| Paper | Seeds per cell | Statistics reported |
|---|---|---|
| **Skinnider 2024 (NMI)** https://doi.org/10.1038/s42256-024-00821-x | **10** independent training datasets per parameter set; **180 models total** (90 SMILES, 90 SELFIES), crossing seed × training-set size (30k/100k/300k) × diversity (Tc ≥ 0.0/0.05/0.10/0.15) × augmentation (canonical/10×/30×) × database (ChEMBL vs GDB-13) | Distribution-learning metrics (FCD) with significance tests across replicates; correlations between invalid-output rate and performance gain |
| **Li 2026** https://arxiv.org/abs/2603.19335 | **5** for the 20-variant DPO grid; 3 for several headline cells | Mean ± SE; **Bonferroni correction over m=19/20 comparisons**; reports seed σ (e.g. DPO σ=2.01 at 3B) and uses it to judge whether a 1.56 pp gap is real |
| **Narang-update 2026** https://arxiv.org/abs/2605.20798 | **3 on the baseline only**; single seed for the other 19 methods | **Noise floor σ_baseline = 0.00208 CLIMB-avg units** from the 3 baseline seeds; **bootstrap resampling N=10,000**; flag |z| > 2; **Bonferroni at m=19** |
| **Schmidt 2024** https://arxiv.org/abs/2402.18376 | **0** (no repeats) | Averages over the 3 vocabulary sizes as a noise proxy; **one-sided Wilcoxon signed-rank test** on 30 paired accuracy scores (3 vocab sizes × 10 tasks) at 350M. Limitation stated: further experiments would be needed "to gain a better estimate of any potential noise" |
| **Ali 2024** https://arxiv.org/abs/2310.08754 | **0**, stated explicitly as a limitation | None |
| **Chinchilla 2022** https://arxiv.org/abs/2203.15556 | not reported | 80% subsample bootstrapped 100×, reporting 10th/90th percentiles on fitted exponents |
| **NovoMolGen 2025** https://arxiv.org/abs/2508.13408 | 3 for the PMO hyperparameter search | not uniform across the grid |

**Modal practice for a pretraining ablation: 1 seed per cell, with variance estimated once from a
small baseline replicate set and used as a noise floor for all cells.** That is the Narang-update
design and it is the cheapest defensible protocol. **3–5 seeds per cell is what a paper does when the
ranking itself is the claim** (Li 2026). **10 seeds per cell is an outlier** and comes from a paper
that spent its budget on seeds instead of sizes (Skinnider 2024).

### Does Agarwal et al. apply?

**Agarwal et al. 2021, "Deep Reinforcement Learning at the Edge of the Statistical Precipice"
(NeurIPS 2021, outstanding paper)** — https://arxiv.org/abs/2108.13264 (2021), library `rliable`.

Their recommendations transfer directly, because the structural situation is identical: few runs per
cell, many tasks, aggregate comparison of arms.

- Diagnosis: "5 or less runs are common in the field". Achieving statistically defensible *median*
  scores on Atari 100k would need roughly **50–100 runs**.
- Their alternative, workable at **3–10 runs per task**:
  1. **Interquartile mean (IQM)** — "discards the bottom and top 25% of the runs and calculates the
     mean score of the remaining 50% runs". Robust to outliers, "considerably less bias than median",
     and reaches small confidence intervals with far fewer runs than the median.
  2. **Stratified bootstrap confidence intervals** rather than point estimates.
  3. **Performance profiles** (score distributions across tasks) instead of a single aggregate number.
- Demonstrated consequence: a published ranking reversed under proper uncertainty treatment — "der may
  in fact be better than otr, unlike what the reported point estimates suggest".

**Where the analogy holds for the researcher:** N tokenizations × M downstream tasks is exactly the
arms × tasks structure IQM and performance profiles were built for. Reporting a per-arm IQM over the
task suite with stratified bootstrap CIs, plus a performance profile, is a stronger and cheaper claim
than per-task means with 3 seeds.

**Where it does not hold:** in RL the run-to-run variance comes from environment stochasticity and is
large; in LM pretraining the dominant variance at a fixed data order is smaller and the per-run cost is
far higher, so the trade is tilted toward fewer seeds and more cells. The relevant NLP evidence for
seed variance is **Dodge et al. 2020, "Fine-Tuning Pretrained Language Models: Weight Initializations,
Data Orders, and Early Stopping"** — https://arxiv.org/abs/2002.06305 (2020) — which fine-tuned BERT
hundreds of times per task over 4 GLUE datasets (2,100 released trials) varying only the seed, and
found weight initialisation and data order contribute **comparably** to out-of-sample variance, with
substantial spread. That is fine-tuning rather than pretraining, so it bounds the downstream-eval
component of the researcher's noise, not the pretraining component.

### What a defensible claim looks like, assembled from the above

For "arm A beats arm B and this is not capacity", the literature's assembled standard is:

1. **≥3 model sizes spanning ≥10×**, chosen so the suspected confound (here, embedding share) changes
   materially between the ends (Section 1 norm + MobileLLM's numbers).
2. **A noise floor** measured from ≥3 seeds at one reference cell, applied as a threshold to all cells
   (Narang-update 2026), **or** 3–5 seeds per cell if the ranking itself is the claim (Li 2026).
3. **Multiple-comparison correction** across arms (Bonferroni in both Li 2026 and Narang-update 2026)
   or a paired non-parametric test across tasks (Wilcoxon signed-rank, Schmidt 2024).
4. **A continuous metric reported alongside any thresholded one** (Schaeffer 2023).
5. **An explicit factorial against the suspected confound** at ≥2 sizes (Li 2026's 2×2) or an explicit
   equalisation of the confound across arms (Limisiewicz 2026's matched embedding compute).
6. **IQM + stratified bootstrap CIs + performance profiles** when aggregating over a task suite
   (Agarwal 2021).

---

## GAPS

- **No paper states a rule for how wide a scale sweep must be.** The ~3 sizes over ~10× norm is
  inferred from what gets published, not from any methodological guidance. There is no NeurIPS/ICLR
  checklist item, no survey, and no reviewer guide on the question.
- **"Sweep until the embedding fraction changes materially" is nowhere stated as a criterion.** The
  component facts exist (MobileLLM's >20% at 125M vs 0.7% at 70B; Kaplan's rationale for excluding
  embeddings; the Kaplan/Hoffmann accounting discrepancy), but the synthesis is the researcher's, and
  must be presented as a proposed criterion rather than cited as practice.
- **No tied-vs-untied ablation at matched total parameter count exists in the modern LLM-scale
  literature.** The only multi-size tying ablation is Inan et al. 2017 at 200/650/1500 hidden units.
  Everything since either asserts tying (Over-Tokenized 2025) or treats it as a budget-reallocation
  decision (MobileLLM 2024).
- **Base rates for ranking flips are unknown.** Wei et al.'s 6-of-11 U-shaped figure comes from tasks
  *selected* for an anomalous trend; Li 2026 and Tay 2023 report single clear inversions; Schmidt 2024
  reports "crossing lines" without quantifying how many pairs cross. Nobody has measured, over a
  representative sample of ablations, what fraction of pairwise orderings survive a decade of scale.
- **Seed-variance for pretraining (as opposed to fine-tuning) at a fixed configuration is barely
  measured.** Dodge et al. 2020 covers fine-tuning; the Narang-update's σ = 0.00208 CLIMB-avg at 1.2B
  is one of the only published pretraining noise floors, from three seeds at one size. Nobody reports
  how that noise floor scales with model size, which is exactly what a cross-size ranking claim needs.
- **How to choose the reduced subset is never justified methodologically.** Papers carry the winners
  upward (Narang-update, Schmidt) — which is the selection rule most likely to miss a flip, since a
  flip requires carrying up a loser too. No paper argues for its subset rule; the 2026 Narang update's
  finding that a 1.2B winner diverged at 3B is the strongest evidence that the usual rule is unsafe.
