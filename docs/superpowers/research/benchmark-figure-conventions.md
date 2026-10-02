# Figure conventions in molecular / peptide generative-design benchmark papers

Research note, compiled 2026-10-02. Every claim carries a URL and year. Figure
inventories were extracted from open-access full text (arXiv PDF/HTML, PMC, JMLR,
publisher OA); anything not directly read is flagged inline.

Sections are ordered 1–6 as requested; they were written 5, 4, 2, 6, 1, 3.

---

## SECTION 1 — WHAT BENCHMARK PAPERS ACTUALLY SHOW

Figure inventories below were extracted from the open-access PDFs (arXiv) by text
extraction of every `Figure N:` / `Table N:` caption, in document order.

### 1.1 PMO — Gao, Fu, Sun, Coley, "Sample Efficiency Matters", NeurIPS 2022 D&B

https://arxiv.org/abs/2206.12411 (2022) · code https://github.com/wenhao-gao/mol_opt ·
read from arXiv PDF, 55 pp. 25 methods × 23 oracles, 10,000-oracle budget, metric = AUC of
the top-10 average vs oracle calls.

**Direct answer to the variance question: yes, PMO reports variance, and more of it than
any other paper in this section.** Every cell of every results table is **mean ± standard
deviation over 5 independent runs** with different random seeds ("Reported results are
from 5 independent runs with various random seeds"). It is the only paper here that also
**discloses and visualises its hyperparameter tuning** (Appendix D.2, 17 figures, tuned on
3 runs of two GuacaMol tasks). It reports no confidence intervals, no bootstrap, no
significance test, and its main-text *figures* carry no variance bands in their captions —
variance lives entirely in the tables.

**Its main results object is Table 2, not a figure.** Only three figures appear in the
main text, and none of them is a results leaderboard:

| # | What it is | Claim it supports |
|---|---|---|
| Table 1 | Taxonomy grid: rows = optimization algorithm (GA, MCTS, BO, VAE, GAN, HC, RL), columns = molecular assembly strategy (SMILES / SELFIES / graph-atom / graph-fragment / synthesis); cells = method names. | Frames the two factors the paper later ablates. A *design-space* table, not results. |
| **Table 2** | **The main results object.** Rows = 23 tasks, columns = the 10 best methods. Cells = **mean ± std of AUC Top-10 over 5 runs**; best per task bolded. | "Older algorithms still win" — REINVENT and Graph GA top the table despite being years old. |
| Table 3 | Rank table: rows = all 25 methods, columns = 6 metrics (AUC Top-1/10/100 and Top-1/10/100) plus mean rank. Integer ranks only. | Ranking is largely metric-invariant *except* where it is not — SMILES-LSTM-HC drops under AUC vs final-score metrics, i.e. it is sample-*inefficient*. The rank-stability argument. |
| Figure 1 | Two optimization curves, x = oracle calls, y = top-10 average: `isomers_c9h10n2o2pf2cl` (isomer-type) and `celecoxib_rediscovery` (similarity-type). Only 8 of 25 methods shown "for clarity"; full set in Appendix A. | The *shape* of the budget curve, not the endpoint, separates methods. |
| **Figure 2** | **Two parity scatter panels.** (a) x = AUC Top-10 of the SMILES variant, y = AUC Top-10 of the SELFIES variant **of the same method**; one point per task; colour = optimization algorithm; **the fraction of tasks above the parity line printed in the legend parentheses**. (b) the same construction for model-free vs model-based methods (with a caption note that the GA pair "is not a head-to-head comparison"). | **The controlled ablation, as a main-text headline figure.** One factor isolated, everything else held fixed, 23 tasks aggregated into a single "fraction above parity" number per algorithm family. |
| Figure 3 | Clustered heatmap: rows = oracles, columns = methods, colour = AUC Top-10 **min-max normalised within each task**; UPGMA dendrogram over oracles. | Oracles fall into a similarity cluster, four isomer-based oracles, and unclustered outliers. "Different types of landscape are more suitable for different kinds of methods" — a method × task-family interaction. |

Appendix: Figures 4–8 are per-oracle-family optimization curves for all 25 methods;
Figure 9 SA_Score of top-100 by method; Figure 10 internal diversity, showing "the stronger
a model is in optimization, the less diverse the results are"; **Figures 14–30 are one
hyperparameter-sensitivity plot per algorithm** (W&B, endpoint = sum of AUC Top-10 on
`zaleplon_mpo` + `perindopril_mpo` over 3 runs — e.g. "sigma (σ) has large impact… and the
optimal value is much larger than the default setting in the original paper"); Figure 31
the distribution of Zaleplon MPO values in ZINC-250k; Tables 4–22 the full 25 × 23
leaderboards for all six metrics, every cell mean ± std over 5 runs.

- **Leaderboard tables vs plots:** both, tables dominate; the three figures do interpretive
  work rather than reporting scores.
- **Per-task breakdown:** yes, exhaustively — every table is tasks × methods.
- **Variance:** ± std over 5 runs in every table cell; none advertised in figure captions.
- **Controlled ablation as headline figure: yes — Figure 2.** PMO is the only classical
  molecular benchmark in this set that puts an isolated-factor comparison in the main text.

### 1.2 GuacaMol — Brown, Fiscato, Segler, Vaucher, JCIM 59(3):1096–1108 (2019)

https://arxiv.org/abs/1811.09621 (2018) · DOI 10.1021/acs.jcim.8b00839 (2019) ·
https://github.com/BenevolentAI/guacamol. Inventory read from the arXiv PDF; the published
JCIM version has the same figure/table structure (flagged in GAPS as not independently
re-verified panel-by-panel).

| # | What it is | Claim |
|---|---|---|
| Figure 1 | Chemical structure drawings of compounds 1–3, generated by a SMILES-LSTM and experimentally active against RXR/PPAR. | Motivation; not a results figure. |
| Table 1 | Distribution-learning leaderboard: rows = 6 models (Random sampler, SMILES LSTM, Graph MCTS, AAE, ORGAN, VAE), columns = Validity, Uniqueness, Novelty, KL divergence, FCD. **Single numbers.** | The **random sampler scores KL 0.998 and FCD 0.929** — the trivial baseline nearly saturates the distribution metrics, and its novelty is zero by construction. |
| Table 2 | Goal-directed leaderboard: 20 task rows × 5 models + "Best of dataset", plus a **Total row** (sum over tasks, e.g. Best-of-dataset 12.144). Single numbers. | SMILES LSTM scores 1.000 on every rediscovery/similarity task; Graph GA is competitive; "Best of Dataset" is a virtual-screening baseline. |
| Figure 2 | "Results of the Goal-Directed Benchmarks" — a grouped bar chart over the 20 tasks, one bar per model per task. Visual restatement of Table 2. | Per-task breakdown made scannable. |
| Figure 3 | Bar chart: x = model (Best of ChEMBL, Graph GA, Graph MCTS, SMILES GA, SMILES LSTM), y = "Ratio of acceptable compounds" (fraction of the top-100 passing quality filters), 0–0.8. | **The key negative result**: Graph GA 40%, SMILES GA 36%, Graph MCTS 22%, versus **77% for both Best-of-ChEMBL and the LSTM**. High benchmark scores ≠ chemically acceptable molecules. |
| Figure 4 | Plots of the score-modifier functions (MaxGaussian, MinGaussian, Gaussian, Thresholded…). | Methods figure. |

Appendix: Table 3 scoring-function specifications; **Table 4 a list of deliberately
*trivial* objectives** (logP targets −1.0 and 8.0, TPSA 150, CNS MPO, QED, C7H8N2O2,
Pioglitazone MPO); **Table 5 their results, where almost every model scores 1.000** —
supporting the explicit conclusion that "such trivial benchmarks are not suitable for the
assessment of generative models for de novo molecular design."

- **Leaderboard tables with bar charts as visual duplicates.** No optimization curves, no
  oracle-budget axis — exactly what PMO later attacked.
- **Per-task breakdown:** yes (Table 2, Figure 2).
- **Variance: none.** Single numbers throughout; no repeated runs, no standard deviations
  anywhere in the baseline results. The largest statistical gap in the paper and the norm
  that MOSES (3 inits) and PMO (5 runs) reacted against.
- **Controlled ablation:** no. Figure 3 is a post-hoc diagnostic of the winners, not a
  one-variable-isolated study.

### 1.3 MOSES — Polykovskiy et al., Front. Pharmacol. 11:565644 (2020)

https://arxiv.org/abs/1811.12823 (2018) · https://doi.org/10.3389/fphar.2020.565644 (2020) ·
https://github.com/molecularsets/moses. Inventory read from the arXiv PDF (flagged in GAPS:
the published Frontiers version may renumber or add figures).

| # | What it is | Claim |
|---|---|---|
| Figure 1 | Pipeline schematic: dataset → baseline models → evaluation metrics. | Framing. |
| Figure 2 | Grid of drawn example molecules from the MOSES dataset. | Dataset characterisation. |
| Tables 1–4 | Four stacked metric tables, rows = 8 baselines + a **"Train" reference row**, **mean ± std over three independent model initialisations**. T1: Valid, Unique@1k, Unique@10k. T2: Filters, Novelty, IntDiv1, IntDiv2. T3: FCD and SNN, each split into **Test and TestSF (scaffold split)**. T4: Frag and Scaf similarity, Test and TestSF. | The Train row is the reference oracle; the Combinatorial (BRICS-fragment) and NGram non-neural baselines sit close to the neural models on most metrics. |
| Figure 3 | Overlaid kernel-density curves for four physico-chemical properties (MW, logP, QED, SA), one panel per property, one curve per model plus the test set, **with the Wasserstein-1 distance to the MOSES test set printed in brackets in the legend**. | Distribution learning shown as *distributions*, not a scalar — with a scalar distance annotated so it stays rankable. The one results plot in the paper. |
| Figures 4–8 | Methods/architecture schematics (vanillin as graph/string/image; autoencoder; AAE; combinatorial generator) and a diverse-subset molecule gallery. | Methods. |

- **Pure leaderboard-table paper.** No per-task axis exists — there are no tasks, only
  metrics. Figure 3 is the only results plot.
- **Variance:** yes — **mean ± std over 3 model initialisations** in every cell. Three is
  small by Section 5 standards, but MOSES is the earliest paper here to report it at all.
- **Controlled ablation:** no. The nearest thing is the **Test vs TestSF column pairing**
  inside Tables 3–4, which isolates one factor (the split) and shows every model degrades
  on it. A column-pair ablation, never lifted into a figure.

### 1.4 Tartarus — Nigam, Pollice, Friederich, Aspuru-Guzik et al., NeurIPS 2023 D&B

https://arxiv.org/abs/2209.12487 (2022/2023) · https://github.com/aspuru-guzik-group/Tartarus

The cleanest example of a benchmark paper whose **main-text figures are all schematics and
whose entire evidentiary content is tables**.

| # | What it is |
|---|---|
| Figure 1 | Two-panel framework schematic (task definition + simulation workflows; evaluation pipeline). |
| Figure 2 | Simulation-workflow flowchart for organic photovoltaics: SMILES → Open Babel → crest/xtb conformer search → GFN2-xTB single point → HOMO/LUMO, gap, dipole → Scharber-model PCE. |
| Table 1 | OPV leaderboard: 8 generative models + a **"Dataset"** row (best molecule in the training subset of the Harvard Clean Energy Project DB). Cells = **mean ± std of the best objective value over 5 independent runs**. |
| Figure 3 | Organic-emitter workflow flowchart (xtb → TD-DFT with pyscf → singlet–triplet gap, oscillator strength, vertical excitation). |
| Table 2 | Emitter leaderboard, same format, trained on a GDB-13 subset. |
| Table 3 | Protein-ligand leaderboard: 3 targets (1SYH, 6Y2F, 4LDE), each with docking score ΔE **and a structural-filter success rate SR column**; includes both a "Dataset" row and a **"Native Docking"** row (the original crystal-structure ligands). |
| Figure 4 | Docking workflow flowchart (Open Babel → QuickVina2 → smina re-scoring). |
| Figure 5 | Reaction-substrate workflow flowchart (reactant/product optimization → SEAM transition-state guess → constrained conformational sampling). |
| Table 4 | Reaction-substrate leaderboard; baselines = best-in-dataset **and the unsubstituted parent substrate**. |
| Figure 6 | **Three-panel bar chart of model timing** (A training time, B sample time, C single-epoch time with and without GPU), mean ± std over 5 runs. The only quantitative plot in the main text. |

SI: Figures S1–S4 drawn best-molecule galleries per benchmark family; Table S5 raw timing
numbers; Figure S5 estimated CPU training time.

- **Tables only** for performance. Six main figures: five workflow diagrams and one cost
  bar chart. No optimization curves, no per-task plot, no parity plot.
- **Variance:** mean ± std over 5 runs in every table cell.
- **Controlled ablation:** none. The SMILES-VAE/SELFIES-VAE and SMILES-LSTM-HC/
  SELFIES-LSTM-HC row pairs *are* matched representation pairs, but they stay adjacent
  table rows and are never lifted into a figure the way PMO Figure 2 does.
- The distinctive move is **baselines in the table**: "Dataset", "Native Docking" and
  "Parent Substrate" give the reader three reference lines a model must beat.

### 1.5 Therapeutics Data Commons — Huang, Fu, Gao, … Zitnik, NeurIPS 2021 D&B

https://arxiv.org/abs/2102.09548 (2021) · https://tdcommons.ai. 66 datasets, 22 learning
tasks, 23 evaluation strategies, 17 oracles, 33 data processors, **29 leaderboards**.

| # | What it is |
|---|---|
| Figure 1 | Platform overview schematic. |
| Figure 2 | "Therapeutics Machine Learning" conceptual diagram across the discovery→development pipeline. |
| Figure 3 | Tiered-design diagram: 3 problems → learning tasks → datasets. |
| Table 1 | The 22 learning tasks × therapeutic-product type (small molecule / macromolecule / cell & gene therapy) checkbox matrix. |
| Table 2 | The 66-dataset catalogue: size, feature type, task type, **suggested metric, recommended split**. |
| Table 3 | ADMET benchmark-group leaderboard: 22 endpoints × featurisation/pretraining strategies. GIN + context prediction best on 8 endpoints, GIN + attribute masking on 5, **RDKit2D expert descriptors on 5**, SMILES CNN on 1. "the ML SOTA models do not work well consistently for these novel realistic endpoints." |
| **Figure 4** | **Heatmap of domain-generalization results** for the DTI task: in-distribution (2013–2018) vs out-of-distribution (2019–2021) Pearson correlation across methods. In-distribution sits at ~0.7 PCC and is stable across years; **OOD degrades by 33.9% to 43.6%**. Best methods are MMD and CORAL, but "the standard training strategy has similar performances as current ML SOTA domain generalization algorithms". |
| Table 4 | DTI-DG leaderboard: In-Distribution and Out-of-Distribution PCC, **average ± standard deviation across five random runs**, best bolded and second underlined. |
| Table 5 | Molecule-generation leaderboard under oracle budgets of **100 / 500 / 1,000 / 5,000** calls. Only at 5,000 do Graph-GA (−14.811) and SMILES-LSTM (−13.017) surpass the best-in-data docking score. "Graph-GA dominates the leaderboard with 0 learnable parameters." |

- **Leaderboard tables are the product** (29 of them, hosted live); the paper's four figures
  are three schematics plus one heatmap.
- **Per-task breakdown:** yes — Table 3 is 22 endpoints wide.
- **Variance:** mean ± std over 5 random runs in Table 4; sparser elsewhere.
- **Controlled ablation:** closest is Figure 4's in-distribution vs out-of-distribution
  heatmap — one factor (temporal domain shift) isolated and shown as a two-block colour
  matrix. TDC Table 5 also pre-announces PMO's central result a year early.

### 1.6 Ehrlich functions / holo-bench — Stanton, Alberstein, Frey, Watkins, Cho (2024)

https://arxiv.org/abs/2407.00236 (2024) · https://github.com/prescient-design/holo-bench.
Plot mechanics and metric definitions are in §4.1; here, the figure strategy.

**Only six figures, four of them non-quantitative, and two controlled ablations.** This is
the paper whose figure strategy is closest to what a NeurIPS evaluation-track paper should do.

| # | What it is | Claim |
|---|---|---|
| Fig. 1 | A render of the **Ackley function**. | The framing argument: a good test function earns its place through *geometric* similarity (many local minima, changing curvature), not semantic realism. Opening with someone else's test function is a deliberate move. |
| Fig. 2 | Four-panel biology schematic: (a) arginine–glutamate salt bridge; (b–c) antibody–epitope binding; (d) two antibodies binding homologous epitopes as shared motifs in sequence space (PDB 3gbn, 4fqi). | Justifies the **construct**: why spaced motifs with quantised partial credit are the right abstraction for affinity maturation. |
| Fig. 3 | Schematic of an epistatic second-order interaction. | Justifies non-additivity. |
| **Fig. 4** | **Four-panel controlled difficulty sweep**, x = Function Evaluations (M), y = Simple Regret, **optimizer held fixed** at the GA baseline, one Ehrlich parameter varied per panel off a fixed base config (L=256, c=4, k=4, q=k). 10/50/90 quantile bands over 32 trials. | "The benchmark has a difficulty dial, and here is the dial working." Quantization q dominates. |
| **Fig. 5** | **Three-panel reverse ablation**: function fixed (k=8, q=4), two GA hyperparameter configurations; panels = Simple Regret, Cumulative Regret, **Feasible Particles** (a diagnostic, not a score). Same quantile bands. | "Configuration A outperformed B on most random seeds, but some seeds show essentially no improvement" — the honest-variance statement. The tuning was itself a 512-configuration W&B Bayesian search consuming >200B function evaluations. |
| Fig. 6 | Python code listing (minimal usage example, BoTorch + PyTorch optim APIs). | Usability. |

- **No leaderboard table at all.** No baselines other than the authors' own GA. The
  contribution is the test function, and the paper refuses to pretend it is a ranking.
- **Variance: 10/50/90 quantiles over 32 trials**, stated once in §4.2 and applied to every
  plot. Not mean ± std, not a CI on a mean. The single most transferable convention here.
- **Controlled ablation as headline figure: yes, twice, and in both directions** — fix the
  optimizer and sweep the problem (Fig. 4), then fix the problem and sweep the optimizer
  (Fig. 5). The two-sided design is what makes the difficulty claim credible rather than
  asserted.
- Known optimum by construction (`f* = 1` for feasible sequences), so **simple regret has a
  true zero** and "solved" is well defined.

### 1.7 Cross-cutting pattern

| Paper | Results in | Per-task | Variance | Controlled ablation in main text |
|---|---|---|---|---|
| GuacaMol 2019 | tables + duplicate bar chart | yes | **none** | no |
| MOSES 2020 | tables | n/a | ± std, 3 inits | no (TestSF column pair only) |
| TDC 2021 | tables (29 live leaderboards) | yes | ± std, 5 runs (partial) | Fig. 4 heatmap, in-dist vs OOD |
| PMO 2022 | tables | yes | ± std, 5 runs | **yes — Fig. 2 parity scatter** |
| Tartarus 2023 | tables only | yes | ± std, 5 runs | no |
| Ehrlich 2024 | **plots only** | n/a | **10/50/90 quantiles, 32 trials** | **yes — Figs. 4 and 5** |

The trend over six years is monotone: toward budget-aware curves, toward reported variance,
and (only in the two most recent) toward one-factor-isolated figures. **Nobody in this
corpus uses a critical-difference diagram, a performance profile, a bootstrap confidence
interval or any significance test** — the Section 5 conventions have not reached this
literature. That is simultaneously the gap and the opportunity.

## SECTION 2 — CRITIQUE PAPERS AND WHAT THEY SHOW

**Author-list correction up front:** the failure-modes paper is Renz, **Van Rompaey,
Wegner, Hochreiter, Klambauer** (JKU–Janssen), not Van Deursen/Reymond. The rebuttal is
**Langevin, Vuilleumier, Bianciotto** (J Cheminform 2022); the Langevin–Grebner–Matter
applicability-domain paper is a *different*, constructive paper in ACS Omega 2023.

### 2.1 Renz, Van Rompaey, Wegner, Hochreiter, Klambauer, "On failure modes in molecule generation and optimization", Drug Discovery Today: Technologies 32–33:55–63 (2019/2020)

DOI https://doi.org/10.1016/j.ddtec.2020.09.003 · open-access PDF
https://epub.jku.at/download/pdf/5687408.pdf (2020) · SI https://epub.jku.at/obvulioa/download/res/9/5687411

**Three main figures, one main table.** The split-scoring-function design: ChEMBL data
for JAK2/EGFR/DRD2 is binarised and split in half *before* any model is fitted. Three
random forests on ECFP4-1024: **OS** (split 1, drives optimization), **MCS** (split 1,
*different random seed* → isolates model-specific bias), **DCS** (split 2 → isolates
data-specific bias). Table S1 shows all three have statistically indistinguishable
held-out ROC AUC, which is what makes any later divergence damning. Three optimizers
(Graph GA, SMILES-LSTM hill-climb, Particle Swarm in CDDD latent space), 10 runs, 150
iterations.

| # | What it shows | Device |
|---|---|---|
| Table 1 | GuacaMol distribution-learning leaderboard **adapted from the benchmark being criticised**, with one extra column: **AddCarbon** (insert a random `C` into a random training SMILES; keep if valid and novel). Validity 1.000, Uniqueness 0.999, Novelty 1.000, KL 0.982, FCD 0.871 — beats every baseline but the LSTM. | **Trivial baseline in the last column of the opponent's own table.** The FCD row is the honest exception and they say so. |
| Figure 1 | Drawn structures: (a) 2 Graph-GA compounds with a reactive diene and an N–F bond; (b) 2 SMILES-LSTM compounds with –O–O–, S–S, N–O chains; (c) 3 real DRD2 actives from the training set. Table S2 gives each unusual substructure as SMARTS with its frequency in 1,591,378 ChEMBL compounds. | **Gallery of embarrassing winners**, with a quantitative backstop so it is not taste. |
| **Figure 2** | **3 × 4 grid.** Rows = optimizer, columns = iteration 0 / 50 / 100 / 150. Every panel: x = Optimization score (0–1), y = Data Control score (0–1), **solid y = x diagonal**, grey background points = random ChEMBL, coloured points = generated (one colour per run), **grey contours marking where split-1 actives lie**. Clouds start on the diagonal near the origin and end far right-and-below it (OS ≈ 0.7–0.9, DCS ≈ 0.1–0.3), inside the split-1 contour. | **The identity line as implicit contract.** Departure from y = x is the violation, and four time steps make it kinetic rather than a single number. The contour overlay does the memorization argument simultaneously. |
| **Figure 3** | **3 × 3 grid** (rows = optimizer, columns = DRD2/EGFR/JAK2). x = Iteration 0–150, y = Score. Three curves: **OS blue solid**, **MCS green dashed**, **DCS red dash-dot**; bold line = median of per-run means, shaded = IQR over 10 runs. OS climbs; MCS lags; DCS is flat and far below. | **The train/test overfitting idiom transplanted into generative chemistry.** The blue–green gap and green–red gap are two separately named mechanisms you can measure with a ruler, and the IQR bands establish the gap exceeds seed noise. |

SI worth noting: Fig. S2 overlaid predicted-score histograms for train/test
actives/inactives (the bias mechanism in one panel); Fig. S3 the Fig. 2 grid for DRD2 and
JAK2, whose caption adds a neat symmetry argument ("Due to symmetry the region of split 2
actives could be obtained by mirroring the contour lines across the diagonal"); Fig. S8
t-SNE showing each seed finding a private corner; Fig. S10 paired nearest-neighbour
similarity distributions to split 1 vs split 2.

### 2.2 Langevin, Vuilleumier, Bianciotto, "Explaining and avoiding failure modes in goal-directed generation of small molecules", J Cheminform 14:20 (2022)

DOI https://doi.org/10.1186/s13321-022-00601-y · https://pmc.ncbi.nlm.nih.gov/articles/PMC8973583/ (2022) ·
code https://github.com/Sanofi-Public/IDD-papers-avoiding_failure_modes

Thesis: the Sopt/Smc/Sdc divergence is real but is a property of the *predictive models*,
not evidence of generator bias-exploitation. On an untouched held-out sample of the
original data, optimization and control models **already disagree by about the same
amount** at high Sopt. **Eight figures, no main tables.**

| # | What it shows | Device |
|---|---|---|
| Fig. 1 | Schematic redrawing of Renz's experimental setup (dataset → split 1/2 → three models → generation). | Restate the opponent's design faithfully first. |
| Fig. 2 | **Exact reproduction** of Renz Fig. 3 on all 9 target × optimizer cells, 10 runs, median of per-run means with 97.5% envelope. | **Replicate before disputing.** Removes "you ran it wrong" as a counter. |
| **Fig. 3** | **3 columns (targets) × 3 rows.** Row 1: **hexbin density, log colour**, x = Sopt, y = Sdc, **dashed y = x**; the dense blob already sits well below the diagonal (EGFR's high-Sopt cluster at Sopt ≈ 0.3–0.5, Sdc ≈ 0.1). Row 2: mean absolute difference Sdc−Sopt as a function of a *threshold* on Sopt (at abscissa x, MAD is over molecules with Sopt > x) — monotone, reaching 0.3–0.4. Row 3: boxplots of Sdc binned by Sopt with 95% CI. Held-out set 10× augmented with Topliss-tree analogues to de-noise the high-Sopt tail (Fig. S9 = un-augmented). **No generative model anywhere in the figure.** | **The mirror attack.** Same axes, same diagonal, real held-out molecules instead of generated ones. Same picture with no generator ⇒ the picture was never about the generator. Scatter → hexbin also neutralises the "look how far the coloured points moved" reading. |
| **Figs. 4–5** | **Renz Figure 3's exact layout** — 3 × 3, same colours, same line styles, same IQR envelopes — plus **red box-and-whisker tolerance intervals every 30 iterations**, computed as the distribution of control score expected given the observed Sopt distribution, using P[Sdc\|Sopt] estimated from held-out data. The red Sdc curve lies **inside the whiskers everywhere**. | **Overlay a null expectation on the opponent's own figure.** The cheapest rebuttal graphic in the corpus and the most devastating: concede every data point, change only the reference frame. |
| Fig. 6 | Paired similarity distributions (Tanimoto, Morgan r=2, 1024 bits, **log y-axis**) of the top-5%-by-Sopt molecules in an *external test set* to split 1 vs split 2. The split-1 skew is already present in untouched data. | Neutralises Renz Fig. S10 by the same method. |
| Figs. 7–8 | The Fig. 3 composite, then median OS/MCS/DCS trajectories, for **ALDH1 and a modified-architecture JAK2** — two constructed cases where optimization and control models agree a priori. The generator then scores highly on *both*. | **The constructive close / positive control.** Converts "your evidence doesn't support your claim" into "here is what the uncontaminated experiment looks like." |

Honest concession worth copying: they state Renz's two problems — score divergence and
poor chemical quality — are **independent**, and that even in their clean tasks "the
molecules generated can be irrelevant from a drug-discovery perspective" (Figs. S2, S5
show non-drug-like LSTM output: repeated tetrazoles, S–S chains). **Renz's Figure 1
survives the rebuttal intact; only Figures 2–3 are reinterpreted.** They also attack the
raw data, not just the statistics: the datasets are small (842 / 667 / 842 molecules; 40 /
140 / 59 actives) and Table S1 lists questionable DRD2 "actives". Figs. S10–S11 repeat the
Fig. 3 composite with physico-chemical descriptors and Atom-Pair fingerprints, showing the
disagreement is descriptor-dependent.

Related but distinct: Langevin, Grebner, Güssregen, Sauer, Li, Matter, Bianciotto,
"Impact of Applicability Domains to Generative Artificial Intelligence", *ACS Omega*
8(25):23148–23167 (2023), https://doi.org/10.1021/acsomega.3c00883 (PMC10308412).
Constructive, 14 figures and 12 tables. Three transferable devices: **Fig. 3**
binary-vs-count fingerprint comparison showing the mechanism by which a bit-vector
applicability domain waves through a molecule with ten repeats of a fragment; **Fig. 10**
a *molecular Turing test* with 15 human participants (black bar = mean, box = 90%
interval) used as the ground truth the metrics are validated against; **Table 10** a
gallery of "typical problematic structures generated with the three ADs that are able to
optimize the scoring function while producing high scores on evaluation metrics". Their
critique content: classic QSAR applicability domains have ~100% rejection rates for
drug-likeness, and **QED can be optimized in unintended ways**, so it is not a valid AD.

### 2.3 Tripp & Hernández-Lobato, "Genetic algorithms are strong baselines for molecule generation", arXiv:2310.09267 (2023)

https://arxiv.org/abs/2310.09267 (2023). Six pages.

**Zero figures. Two tables.** The absence is the rhetoric: the claim is that the field's
*numbers*, not its concepts, are the problem, so only numbers are presented.

- **Table 1** — ZINC-250k unconditional generation. Columns Validity / Novelty@10k /
  Uniqueness; eleven rows (JT-VAE, GCPN, MolecularRNN, GraphNVP, GraphAF, MoFlow, GraphCNF,
  GraphDF, ModFlow, GraphEBM) plus **AddCarbon (Renz et al.) 100 / 99.94 / 99.86** and
  **MOL_GA 99.76 / 99.94 / 98.60**. A footnote marks which rows are copied from which
  source paper's own table; another notes "the lack of error bars is because most papers do
  not report them."
  **Device: the saturated leaderboard.** Eleven rows of near-identical numbers; the reader
  does the "wait, these are all the same" work. Importing AddCarbon explicitly chains the
  two critiques. The error-bar footnote is a quiet second indictment.
- **Table 2** — PMO AUC Top-10, 23 task rows, three columns: REINVENT and Graph GA (*both
  copied from Gao et al.*) and MOL_GA (their runs), each ± std. Then three summary rows:
  **Sum** (14.196 / 13.751 / **14.708**), **Old Rank** (1 / 2 / N/A), **New Rank** (2 / 3 /
  **1**). One hyperparameter changed — offspring per generation ≈100 → 5, giving ≈2000
  improvement steps inside a 10,000-oracle budget instead of ≈100.
  **Device: print the leaderboard being rewritten.** The Old Rank / New Rank row pair is
  the entire instrument. The authors then disown the win — "We believe this result is
  likely an artifact of the tuning of the baselines in PMO, rather than MOL_GA being an
  especially good method (given that MOL_GA is essentially Graph GA)" — converting a SOTA
  claim into an indictment of tuning effort as a confounder. Appendix A (the whole method)
  is one paragraph of quantile-based sampling; its triviality is part of the argument.
  Prescription: the **"GA criterion"** — referees should require a new method to show an
  empirical or conceptual advantage over a GA.

### 2.4 Tripp & Hernández-Lobato, "Diagnosing and fixing common problems in Bayesian optimization for molecule design", arXiv:2406.07709 (2024)

https://arxiv.org/abs/2406.07709 (2024) · ICML 2024 **AI for Science workshop** (not the
main track) · code https://github.com/AustinT/basic-mol-bo-workshop2024. Four figures, one
table. The same move aimed at themselves: "poor BO performance in prior works may
essentially be due to poor tuning of hyperparameters."

- **Fig. 1** — a 1D toy objective; black dots (known data) clustered near a local optimum,
  red dashed line marking a global optimum nowhere near the data. Deliberately
  "qualitatively similar to molecular design."
- **Figs. 2–3 — the hyperparameter dial figures.** Two column-pairs each; **top row = GP
  posterior, bottom row = PI acquisition function**, for prior width σ = 1.0 vs 0.1
  (Fig. 2) and lengthscale ℓ = 0.05 vs 50.0 (Fig. 3). At σ = 0.1 the acquisition is flat on
  the unexplored half — the optimizer will never go there.
  **Device: same model, one hyperparameter changed, opposite conclusion — with the causal
  chain stacked vertically** (hyperparameter → belief → where the algorithm looks). The
  cleanest way to argue that a reported benchmark number measures a *setting*, not a method.
- **Fig. 4** — two pairs of molecules with **identical binary Morgan r=2 fingerprints**:
  two alkanes of different length, and celecoxib vs a much larger analogue with repeated
  substructures (SMILES in Appendix C). The Tanimoto-on-binary-fingerprints kernel used by
  PMO's own GP BO baseline literally cannot tell them apart.
  **Device: the counterexample pair** — one picture falsifying a representation assumption.
- **Table 1 (Appendix B)** — PMO AUC Top-10, 23 rows, four columns: REINVENT* (Gao et al.),
  MolGA** (their 2023 paper), Genetic GFN† (Kim et al. 2024), **Our GP BO**. Sums 14.196 /
  14.708 / 16.213 / **16.303**. Their GP BO beats Gao's own GP BO implementation by **>3.0
  points — "about the same as the score difference between the best and 10th best methods
  from Gao et al."**
  **Device: the magnitude comparison.** "Retuning one baseline moves it as far as the entire
  top-10 spread" is a quantitative statement that rank is dominated by tuning effort, not
  method identity. The discussion then volunteers against itself: no ablation, no
  acquisition-function sweep, "a very limited pilot study."

### 2.5 Thomas, O'Boyle, Bender, de Graaf, "Re-evaluating sample efficiency in de novo molecule generation", arXiv:2212.01385 (2022)

https://arxiv.org/abs/2212.01385 (2022), NeurIPS 2022 AI4Science workshop. This is the
paper that attacks PMO's **metric** rather than its tuning, and the closest match to
"budget sensitivity".

- **Fig. 1** — per-replicate distributions of MW, logP, and **ratio of de-novo fingerprint
  bits** (0–10% of bits never seen in ZINC-250k) for the top-10 molecules of each of 5
  REINVENT replicates on JNK3, against the ZINC-250k reference distribution. 4 of 5
  replicates drift far outside.
  **Device: "same method, different seed", split out per replicate** — shows the drift is
  not a fluke *and* that replicates disagree, both hidden by a single aggregate AUC.
- **Fig. 2** — 2 × 5 grid of drawn structures, labelled Run 1 … Run 5; oversized molecules
  with repeating substructures. **Laid out by replicate, which pre-empts "you picked the
  worst one".**
- **Fig. 3 — the purest "this comparison is invalid" graphic in the whole corpus.** Grouped
  bar chart. x = ~25 generative models (SMILES-AHC … MolDQN); **y = Rank, axis inverted
  (1 at top, 25 at bottom)**; **four bars per model** = AUC Top-10, AUC Top-10 (Filtered),
  AUC Top-10 (Diverse), AUC Top-10 (Combined). Evolutionary/rule-based methods (Graph GA,
  GP BO, SMILES GA, STONED) drop sharply under the filtered/diverse variants because they
  never learned the ZINC-250k distribution.
  **Device: the dependent variable *is the rank itself*, and within-group bar-height
  disagreement is the message.** If a model swings from 4th to 18th across defensible
  variants of the same metric, the published rank is a property of the metric, not the
  model. The inverted y-axis keeps "tall = good" so the eye reads ragged groups as
  instability before reading the caption.
- **Fig. 4** — per-oracle grid/heat-map of AUC Top-10 (Combined) by model × objective.
- **Fig. 5 (appendix)** — hyperparameter optimization of Augmented Hill-Climb on two test
  objectives, yielding batch_size = 256, **σ = 120**, K = 0.25. PMO tuned REINVENT to
  **σ = 500**, and σ ≥ 240 is already known to push generation outside the training property
  space. **The benchmark's own tuning choice produced the pathological chemistry.**
- **Figs. 6–7 (appendix)** — property space of the benchmark reference molecules vs
  ZINC-250k: six PMO objectives have reference molecules in the **lowest 0.01%** of
  ZINC-250k MW and logP space, i.e. unreachable by construction for distribution-based models.
- Self-flagged limitation: the re-ranked models were not re-tuned against the new metrics.

### 2.6 Confounded-protocol papers

#### Gao & Coley, "The Synthesizability of Molecules Proposed by Generative Models", JCIM 60(12):5714–5723 (2020)

https://arxiv.org/abs/2002.07007 (2020). **Four main figures, one main table**, plus
Figs. S1–S14. The protocol criticism: benchmarks score with heuristic proxies (SA Score,
SCScore) that do not measure what they claim, validated here against a real CASP oracle
(ASKCOS).

- **Fig. 1** — 7-panel taxonomy schematic (a)–(g) of ways to handle synthesizability, from
  "ignore it" to "construct from buyable building blocks".
  **Device: a map, not a complaint** — it pre-allocates every remedy a slot, so the rest of
  the paper reads as filling in which slots work.
- **Fig. 2 — five panels, and the layout is the argument.**
  (a) **stacked bar**, y = fraction synthesizable, x = dataset (**MOSES 0.898, ChEMBL 0.683,
  ZINC 0.613, Sheridan 0.567, GDB 0.034**), stack coloured by ASKCOS step count 0–10. GDB —
  widely used in generative ML — is **3.4% synthesizable**.
  (b) stacked bar within MOSES (Random **0.898** > VAE 0.870 > LSTM 0.803 > AAE 0.790) and
  within ChEMBL (Random **0.683** > AAE 0.643 > LSTM 0.603 > VAE 0.600). The **"Random" bar
  is leftmost and tallest in each group** — the trivial baseline again, in bar form.
  (c)–(e) line plots with a grey histogram behind: x = **SA_Score**, **SCScore**, and — the
  knife — **raw SMILES string length**; left y = fraction CASP-solvable, right y = molecule
  count; one coloured line per dataset.
  **Device: dataset-dependence as the attack.** If the heuristic measured synthesizability
  the lines would superimpose; they do not. And SMILES string length produces a monotone
  curve of comparable quality to both purpose-built chemistry heuristics. **Putting the
  absurd proxy in the same visual grammar as the respectable ones, as a third panel of the
  same figure, is the whole argument.** Fig. S3 converts the visual trend into ROC AUCs so
  it is not left to the eye.
- **Fig. 3 — a 4 × 10 trellis of two-point slope graphs.** Rows = Best-from-Data / SMILES
  LSTM / SMILES GA / Graph GA; columns = 10 MPO objectives. Each tiny panel has **exactly
  two x-ticks: "No Bias" and "SA_Score Biased"**, y = 0.0–1.0. Green solid = fraction of
  top-100 synthesizable (green dashed = ChEMBL reference); red solid = objective value of
  the top *synthesizable* molecule; red dashed = top molecule regardless. **A missing solid
  red line means no synthesizable structure existed in the top 100.**
  **Device: paired before/after slope charts en masse.** Forty two-point charts in one
  figure; green up, red down, trade-off legible at a glance across every method × objective
  cell. The *absence* of a line is itself a data point — blank cells read as total failure
  without needing a zero.
- **Fig. 4** — four worked examples (a)–(d); (c,d) are cases where the top-1 objective value
  *decreased* while the top-1 **synthesizable** objective value *increased*. "Your headline
  metric moved the wrong way for the right reason."
- **Table 1** — fraction synthesizable in the top-100, by training database × task
  difficulty × bias. ChEMBL-hard **30.2% → 80.2% (SA) → 55.4% (SC)**; MOSES-hard 32.7% →
  77.2% → 58.0%; trivial tasks ≈60% → ≈91% / ≈78%. On hard objectives **~70% of a
  generative model's best proposals are not makeable.**
- **Fig. S14** correlates ASKCOS step count against expert-provided scores — **validating
  the oracle before using it as ground truth.** That move is what makes the whole critique land.

#### Cieplinski, Danel, Podlewska, Jastrzębski, "We Should at Least Be Able to Design Molecules That Dock Well", arXiv:2006.16955 (2020, v5 2023; JCIM 2023)

https://arxiv.org/abs/2006.16955 · benchmark https://github.com/cieplinski-tobiasz/smina-docking-benchmark

- **Fig. 1** — pipeline schematic: proposed molecule → SMINA docking → score from the pose.
- **Table 1** — dataset sizes (1,193–4,199 molecules for 5-HT1B, 5-HT2B, ACM2, CYP2D6).
  Establishes the realistic-small-data premise *before* any result, so "it would work with
  more data" is answered in advance.
- **Fig. 2** — t-SNE of ECFP fingerprints, four panels. Training set **red**; REINVENT output
  coloured by **which objective was optimized** — blue = docking score, orange = hydrogen
  bonding, green = repulsion. The generated clouds are small, tight, and sit on the same
  side of the map **regardless of objective**.
  **Device: colour by objective, not by method.** If the colours separated, the model would
  be responding to the objective. They do not, so the location is a property of the ChEMBL
  prior.
- **Table 2** — three sub-tables, one per objective (docking ↓, repulsion ↓, H-bonding ↑);
  9 rows × 4 target columns; every cell = mean over 250 generated molecules **with internal
  diversity in parentheses**. Rows: CVAE, GVAE, REINVENT, then **Train (50/10/1%)** and
  **ZINC (50/10/1%)** percentile baselines. For 5-HT1B docking: REINVENT −9.774 (**0.506**)
  vs **ZINC-10% −9.894 (0.862)** and ZINC-1% −10.496 (0.861); training-set diversity is
  0.787. "—" marks models that could not produce 250 drug-like molecules at all.
  **Device: the parenthetical disqualifier inside the same cell.** It becomes impossible to
  read the score alone; the one good number sits adjacent to its own refutation. The
  percentile baseline rows turn the table into a **dominance test with a pre-registered
  pass mark** ("task solved = beat ZINC top-1% at training-set diversity") rather than a
  relative leaderboard. The authors concede the threshold "is necessarily arbitrary".
- **Figs. 3–5** — grids of best-scoring molecules per task for REINVENT, CVAE, GVAE.
- **Figs. 6–7 — the confounder scatter.** Four panels each (one per target) with marginal
  histograms; x = docking score, y = **molecular weight** (Fig. 6) / **number of rotatable
  bonds** (Fig. 7); training set red, generated coloured by objective; **Pearson r printed
  in each panel: −0.79 / −0.68 / −0.79 / −0.01** and −0.68 / −0.64 / −0.64 / −0.34.
  **Device: plot the metric against something nobody wants to optimize, and print r.** The
  marginal histograms show distribution shift and correlation in one object. The near-zero
  CYP2D6 r is left in and explained (larger binding site), which buys credibility for the
  other three. They also note generated compounds are shifted toward better docking at
  *equal* molecular weight, crediting REINVENT with a genuine improvement — the concession
  is what makes the attack stick.

#### Wallach & Heifets, "Most Ligand-Based Classification Benchmarks Reward Memorization Rather than Generalization", JCIM 58(5):916–932 (2018)

https://arxiv.org/abs/1706.06619 (2017/2018). The foundational "your split is leaking"
paper for molecular ML and the direct ancestor of the data-leakage literature.
**19 figures, 1 main table.** **AVE bias = (AA − AI) + (II − IA)**, from the distributions
of nearest-neighbour similarity between validation actives/inactives and training
actives/inactives.

- **Fig. 1** — a worked toy counterexample: two four-molecule sets ('+'/'−' labels, blue =
  training, red = validation, arrows annotated with pairwise similarities) with very
  different absolute train–val similarity but **identical AVE bias = 0.8**
  ((0.9−0.5)+(0.8−0.4) = (0.6−0.2)+(0.5−0.1)). "relative distance analysis does not simply
  consider the domain of applicability."
  **Device: pre-empt the reviewer's first objection as Figure 1, before any result.**
- **Figs. 2, 3, 4, 6–15, 17, 18 — one template, fifteen repetitions.** Top-left: AVE bias
  per target, sorted. Bottom-left: AUC for RF / LR / SVM / 1-NN **plotted in AVE-bias
  order**, so the AUC trace visibly climbs left-to-right. Right: a 2 × 2 grid of
  **AUC-vs-AVE-bias scatters**, one per classifier, each annotated with **r², Spearman ρ and
  Kendall τ** (MUV: RF r² = 0.88, ρ = 0.94, τ = 0.79). Applied to J&J (560 targets, already
  single-linkage clustered), MUV (already de-biased by its own procedure), ChEMBL
  v19→v20 time-split, PCBA/ToxCast/SIDER under **random *and* Murcko-scaffold** splits, and
  proprietary Merck sets. Fig. 5 adds a t-SNE of 4 MUV datasets showing the spatial clumping
  that produces the bias.
  **Device: argument by induction over benchmarks, including the ones already de-biased.**
  Sorting the AUC panel by bias converts a correlation into a monotone staircase the reader
  sees before reading any r².
- **Figs. 15, 17, 18 — the shuffled-label experiment, the rhetorical climax.** Train on
  target X's actives (everything >0.7 Dice-similar to any Y active removed) plus target Y's
  inactives; evaluate on Y. **Every training point could legitimately have been used as a Y
  inactive.** Expected AUC 0.5. Observed over 90 pairs: **RF 0.59, SVM 0.60, LR 0.56, 1-NN
  0.52 — and where AVE bias was highest, every algorithm except 1-NN exceeded AUC 0.9.**
  Fig. 18 (kinase vs non-kinase pairs) caption: "When the bias is high, a model trained for
  unrelated proteins can still miraculously separate active and inactive molecules."
  **Device: a model that has never seen a single true positive scoring well above chance.**
  No plot can argue around it.
- **Fig. 16** — a Voronoi schematic showing why fixed-distance exclusion filters (the
  standard remedy) cannot work: fixed-radius circles either swallow neighbouring cells
  (discarding useful hard examples) or leave most of the cell outside.
- **Table 1** — r² between AVE bias and each algorithm's AUC, with a column headed
  **"Previous Unbiasing"** naming the de-biasing the benchmark's own authors already applied.
  **MUV, de-biased by the MUV procedure itself: 0.88 / 0.84 / 0.86 / 0.73.** PCBA Murcko
  0.51/0.66/0.60/0.47; ToxCast Murcko 0.56/0.42/0.29/0.46; SIDER Murcko 0.51/0.63/0.50/0.64;
  Merck shuffled 0.73/0.64/0.65/0.66.
  **That "Previous Unbiasing" column is the entire argument of the table** — and including
  Murcko-scaffold splits is what makes it apply to essentially every published molecular
  property-prediction result.

### 2.7 The shared visual grammar of "this comparison is invalid"

Six reusable devices:

1. **Identity line with a drifting cloud.** Renz Fig. 2; Langevin Fig. 3 (same axes, same
   diagonal, real data, hexbin). The diagonal encodes the implicit contract; departure is
   the violation.
2. **Train/test-style diverging curves.** Renz Fig. 3. Borrows the overfitting idiom so the
   reader convicts the method before reading the text.
3. **Overlay a null expectation on the opponent's own figure.** Langevin Figs. 4–5. Concede
   every data point; change only the reference frame.
4. **Tables where the trivial baseline wins, with the rank printed.** Renz Table 1
   (AddCarbon); Tripp & H-L Tables 1–2 (**Old Rank / New Rank** rows); Gao & Coley Fig. 2b
   ("Random" tallest); Cieplinski Table 2 (percentile baselines as a pre-registered pass
   mark). Strongest form: make the **rank itself the dependent variable** (Thomas Fig. 3).
5. **Confounder scatter with the correlation printed.** Cieplinski Figs. 6–7 (docking score
   vs MW, r = −0.79); Wallach & Heifets' AUC-vs-AVE-bias scatters with r², ρ, τ in every
   panel. Plot the metric against something nobody wants to optimize.
6. **Gallery of embarrassing winners, with a quantitative backstop.** Renz Fig. 1 + Table S2
   (SMARTS counts over 1.59M ChEMBL); Thomas Fig. 2 (per replicate, so cherry-picking
   cannot be alleged); Gao & Coley Fig. 4; Langevin ACS Omega Table 10.

Four specialist moves worth naming separately: **mass paired two-point slope charts where a
missing line is a data point** (Gao & Coley Fig. 3); **colour by objective rather than by
method** to show the objective had no effect (Cieplinski Fig. 2); **posterior and
acquisition function stacked with one hyperparameter changed** (Tripp & H-L 2024 Figs. 2–3);
**one template repeated across fifteen benchmarks, including pre-de-biased ones** (Wallach
& Heifets).

**The asymmetry worth planning around:** the critique papers use **far more varied and more
statistically serious figures than the benchmark papers they attack**. Section 1's benchmark
papers produce leaderboards; Section 2's critiques produce diagnostics. An
evaluation-track paper that wants to survive review has to look like Section 2, not
Section 1.

## SECTION 3 — CONTROLLED ABLATION PRESENTATION

The question is: what visual form does a **representation × task interaction** take in
print? There are four answers in use, in increasing order of how much they commit to.

### 3.1 NovoMolGen — Chitsaz, Balaji, Fournier, Bhatt, Chandar (2025)

https://arxiv.org/abs/2508.13408 (2025) · https://github.com/chandar-lab/NovoMolGen ·
weights https://huggingface.co/chandar-lab/NovoMolGen. Read from the arXiv HTML.

Claimed scope, verbatim: "the largest systematic study (**>30,000 experiments**) to date on
Mol-LLMs by evaluating the effects of molecular representation, tokenization, model
scaling, and dataset size on de novo generation." The factor grid is 4 representations
(SMILES, SELFIES, SAFE, DeepSMILES) × 2 tokenizers (Atomwise, BPE) × 3 model sizes (32M,
157M, 300M) × intermediate pretraining checkpoints × RL hyperparameter sweeps, pretrained
on 1.5B molecules.

**How >30,000 experiments are compressed into print — five devices:**

1. **The heatmap + bar chart composite. This is the load-bearing form.** Figure 2 (and its
   three appendix clones I23, I24, I25) is one figure in two halves: *left*, a heatmap with
   **rows = model variant, columns = the 23 PMO tasks**, cell colour = AUC Top-10
   (normalised in the appendix versions); *right*, a **horizontal bar chart of the
   per-variant total score**. REINVENT and f-RAG baselines are drawn in as reference
   rows/bars. The heatmap carries the interaction; the bar carries the ranking. Reading
   left you see *where* a variant wins, reading right *whether* it wins.
   - Figure 2: rows = model sizes (SMILES only).
   - Figure I23: rows = sizes × tokenizers.
   - Figure I24: rows = intermediate pretraining checkpoints.
   - Figure I25: rows = all four molecule types at 32M.
   **The trick: one factor per heatmap, rows = factor levels, columns = tasks, same
   template instantiated four times.**
2. **The scatter that kills a metric.** Figure 4: x = aggregated PMO benchmark score
   (higher better), y = FCD (lower better), one point per (representation, tokenization)
   model. Text reports **r = 0.376, p = 0.358**. "models that achieve the best FCD scores
   (e.g., DeepSMILES (BPE)) often exhibit only modest downstream performance, whereas
   top-performing models on the PMO benchmark (e.g., SAFE) have sub-optimal FCD scores. The
   low correlation… challenges the utility of distribution-based metrics like FCD as
   reliable predictors of a model's functional capabilities in goal-directed generation."
   This is the paper's most striking figure and it is a *negative* claim — the same device
   the critique papers use (§2.7 item 5).
3. **Saturation curves.** Figure 3: key metrics for 32M–300M at intermediate checkpoints vs
   number of molecules seen; "performance saturates early, with diminishing returns for
   extended training." Figures D6–D9: train/val loss vs molecules seen, one figure per
   tokenizer and per molecule type, with **solid lines = random split and dashed = scaffold
   split on the same axes** — a split × representation interaction shown without a second figure.
4. **Parallel-coordinates plots for the hyperparameter sweep.** Figures E10–E21, **twelve of
   them**, one per (representation × tokenizer × size) cell, from the W&B sweep: axes =
   hyperparameters, final axis = Aggregated Score, lines coloured by score, "showing the
   importance and effects of hyperparameters on the Aggregated Score." **This is where the
   30,000 experiments physically live, and it is entirely in the appendix** — no
   parallel-coordinates plot appears in the main text.
5. **Tables with significance colouring.** Table 1 (MOSES-style distribution metrics),
   Table 2 (docking, 3 runs), Tables G4–G6, I7–I8, J9: **mean (std) over three independent
   model initialisations**, with *blue* = best and *pink* = second-best, **gated on
   p < 0.05**. Colour-as-significance inside a table is a cheap substitute for a
   critical-difference diagram, and it is the only significance machinery in the paper.

Supporting: Figure B5 (batch diversity — length distribution + pairwise Tanimoto);
Figure H22 (8-panel property-distribution overlay: QED, SA, logP, MW, TPSA, Bertz
complexity, rotatable bonds, rings — generated vs ZINC-Random vs ZINC-Scaffold);
Figures I26/J28 (top-k reward vs oracle calls); Figures I27/J29 (drawn top molecules with
scores); Figure K30 (BPE substructure inventory).

**Honest compression ratio: >30,000 runs → 4 main-text figures and 2 main-text tables.**
The main text shows one heatmap+bar composite, one saturation plot, one
metric-invalidation scatter, and a pipeline schematic. Everything factorial is appendix.

### 3.2 Skinnider, "Invalid SMILES are beneficial rather than detrimental to chemical language models", Nature Machine Intelligence 6:437–448 (2024)

https://www.nature.com/articles/s42256-024-00821-x · DOI 10.1038/s42256-024-00821-x (2024).
Read from the publisher PDF.

The methodological gold standard in this corpus for isolating one variable, and it uses
**none** of the Section 1 conventions. Design: **n = 10 paired models per condition**, same
training molecules, representation swapped — so paired tests are legitimate and the plots
can show *differences* rather than two clouds.

| Fig | Panels and chart types |
|---|---|
| Fig. 1 | a–c schematics (chemical-space exploration; a single-character substitution breaking caffeine's SMILES; the benchmark framework). **d** paired dot plot of % valid molecules, SMILES vs SELFIES, n = 10 each, **paired t-test P = 2.0 × 10⁻¹⁰**. **e** same for Fréchet ChemNet distance, P = 1.1 × 10⁻⁹. **f** scatter: x = proportion of valid SMILES generated, y = **difference in FCD between each model and its matched SELFIES twin**, with linear regression + **95% CI band** and Pearson r and P inset. |
| Fig. 2 | **a** loss distributions of valid vs invalid SMILES from a representative model (n = 10⁷ SMILES, two-sided t-test P < 10⁻¹⁵). **b** **Cohen's d dot plot across the n = 10 models** — one effect size per model, one-sample t-test on the effects (P = 1.5 × 10⁻¹³). **c** losses split into **six RDKit error categories**. **d** Cohen's d per error category across 10 models (all P ≤ 1.4 × 10⁻¹⁰). **e** frequencies of each error type as mean proportion across ten models. **f** % valid SMILES **per decile of loss** in samples of 500,000 strings, with a **Jonckheere–Terpstra trend test** (P < 10⁻¹⁵). |
| Fig. 3 | **a** FCD for SMILES vs default SELFIES vs **"Texas SELFIES"** (valency constraints modified to allow pentavalent carbon) vs **unconstrained SELFIES**, n = 10 each, paired t-tests, both P ≤ 3.0 × 10⁻⁵ vs default. **This is the causal manipulation: the validity constraint is the dial, turned in three positions — a dose–response design, not an A/B.** **b–c** loss distributions and Cohen's d for valid vs invalid SELFIES. **d** scatter of % valid SELFIES vs ΔFCD with regression + 95% CI. |
| Fig. 4 | **a/c** overlaid count distributions of aromatic (a) and aliphatic (c) rings for SMILES-generated, SELFIES-generated and training molecules. **b/d** Cohen's d dot plots for the same (P = 1.2 × 10⁻¹¹, 1.2 × 10⁻⁷). **e** **volcano plot**: x = mean difference in effect size between SMILES- and SELFIES-generated molecules, y = statistical significance (paired t-test), dotted line at P = 0.05, across many structural descriptors. **f/g** Cohen's d for valid vs invalid SELFIES (P = 7.6 × 10⁻¹¹, 2.6 × 10⁻¹²). **h** scatter: x = SMILES-vs-SELFIES effect sizes, y = valid-vs-invalid-SELFIES effect sizes, regression + 95% CI — **the mediation argument made visual**: the representation effect *is* the validity-filtering effect. |
| Fig. 5 | **a** valid novel molecules per 100M samples from models trained on 1M GDB-13 molecules, SMILES vs SELFIES (P = 1.9 × 10⁻¹⁰). **b** molecules sampled **outside** GDB-13 chemical space (P = 8.5 × 10⁻⁷). **c** fraction of the full ~975M-molecule GDB-13 reproduced (P = 1.1 × 10⁻⁷). **d** **saturation curve**: x = number of valid molecules sampled, y = proportion of GDB-13 reproduced, one curve per representation. |

**What to steal:**
- **Paired design**, not independent arms. Every SELFIES model has a SMILES twin on the same
  training molecules.
- **Cohen's d dot plots across replicates** instead of a bar with an error whisker. The
  reader sees 10 effect sizes and whether they all point the same way — this is the direct
  replacement for the mean ± std bar in a representation comparison.
- **The volcano plot** (Fig. 4e) for "which of many structural properties shifted" — effect
  size on x, significance on y. Borrowed from genomics, essentially unused in
  molecular-generation benchmarking, and exactly the right object when you have many
  metrics and one factor.
- **Scatter-with-CI-band as the mechanism figure** (1f, 3d, 4h): if the causal story is
  "X works *because of* Y", plot X's effect against Y's effect and show they lie on a line.
- **The three-position dial** (default → Texas → unconstrained SELFIES) rather than a binary
  comparison.

### 3.3 PepFoundry, JCIM 66(2):1264 (2026)

https://doi.org/10.1021/acs.jcim.5c02629 (2026, CC-BY) · preprint
https://chemrxiv.org/doi/full/10.26434/chemrxiv-2025-6tr6g (2025) ·
https://github.com/ — repo linked from the ChemRxiv record.

**I could not read the figure inventory** — pubs.acs.org and chemrxiv.org both served bot
challenges, and I did not attempt to bypass them. Flagged in GAPS. What is established
from the abstract (read):

- The claim: "atomic-level representations of peptides containing non-canonical amino acids
  **consistently outperform sequence-level representations, regardless of model type**."
  The phrase *regardless of model type* is the structural giveaway that the evidence is a
  **representation × model-type grid**, and that the claim is one of a **main effect with
  no interaction** — the opposite of what PMO Figure 2 and NovoMolGen Figure 2 are built to
  surface. A grid showing no interaction is a legitimate result, but it needs the grid shown.
- "We additionally explore the representation of non-canonical peptides through **latent
  space visualization** and show that models with atomic-level information can effectively
  learn relationships between analogous sequences of L-peptides, D-peptides, and peptoids."
  A latent-space projection used as the **mechanism figure** — the role Skinnider gives to
  his scatter-with-CI plots, and Cieplinski to his t-SNE (§2.6).
- Mechanism: peptides handled as SMILES in CHUCKLES format → atom-mapped RDKit molecule
  objects → atom-level features (Morgan fingerprints, graph representations).

### 3.4 HELM-BERT, JCIM 66(14):7900–7916 (2026)

https://arxiv.org/abs/2512.23175 (2025, preprint, read) · DOI 10.1021/acs.jcim.6c00451 ·
PMC13417886 · https://github.com/clinfo/HELM-BERT. Encoder-only DeBERTa-style transformer
(6 layers, H = 768, A = 12) pretrained with MLM on 39,079 peptides in HELM notation.

| # | What it is |
|---|---|
| Figure 1 | Architecture schematic: HELM tokenisation → span masking → hybrid first layer (disentangled self-attention + nGiE) → five transformer blocks → Enhanced Mask Decoder with two weight-tied refinement steps → MLM head. |
| **Table 1** | **The controls table.** For each encoder (HELM-BERT, MoLFormer, PeptideCLM) × each protocol (Full FT / Head FT / Linear probe): head architecture, encoder params, head params, learning rates. Explicit note: **"All encoders are of comparable scale (43–54M parameters)."** |
| Table 2 | Same for the PPI task, adding ESM-2 at 35M / 150M / 650M and a peptide-descriptor baseline, with concat dim and head params. |
| Figure 2 | Pretraining train/val MLM loss curves, 127 epochs, early stopping patience 20. |
| Table 3 | CycPeptMPDB permeability leaderboard: R², Pearson r, RMSE, MAE; **mean ± std**; **a † marker on every value statistically significantly different from HELM-BERT**; reported separately for all three protocols. |
| Figure 3 | **Ablation loss curves**: validation MLM loss vs step for HELM-BERT and five architectural variants (w/o Span Masking, w/o EMD, w/o nGiE, w/o Disentangled Attention, Vanilla-BERT), overlaid on one axis. |
| Table 4 | **Architecture ablation on the downstream task**: the same six variants × four metrics, mean ± std, † for significance. Disentangled attention and Vanilla-BERT are the only significant drops. |
| Table 5 | **Pretraining-data-composition ablation**: full / w/o ChEMBL / w/o Propedia / w/o CycPeptMPDB / from scratch. Reported with an **effect size**: from-scratch R² = 0.664, "p = 0.002, d = 2.48". Data composition otherwise had limited impact. |

**The transferable move is Table 1 / Table 2** — a *controls table* published alongside the
results, stating parameter counts, head architectures and learning rates for every arm,
with an explicit sentence that the encoders are matched in scale. That pre-empts "your
baseline was undertuned", which is the most common reviewer objection to a representation
comparison (§6.4 item 2).

Second transferable move: **the ablation is split across one figure and two tables** — the
figure shows its effect on the *pretraining objective*, the tables on the *downstream
metric*. Keeping those separate is deliberate, and it is the same pretraining-metric vs
downstream-metric disconnect that NovoMolGen Figure 4 makes its headline.

### 3.5 What a representation × task interaction looks like in print — the four forms

| Form | Example | Shows | Cost |
|---|---|---|---|
| **Parity scatter**, one point per task, x = factor level A, y = factor level B, fraction above the line annotated | PMO Fig. 2 (2022) | whether the factor flips *per task*, and by how much | handles only 2 levels of 1 factor |
| **Heatmap (rows = factor levels, cols = tasks) + total bar chart** | NovoMolGen Fig. 2, I23–I25 (2025) | many factor levels at once; where each wins | colour hides magnitude; the normalisation choice must be stated |
| **Paired dot plot + Cohen's d per replicate + volcano plot** | Skinnider Figs. 1d–e, 2b, 4e (2024) | effect size, its consistency across replicates, and which of many metrics moved | needs a paired replicate design (n ≈ 10) from the start |
| **Adjacent table rows / column pairs** | Tartarus Tables 1–4; MOSES Test vs TestSF | nothing visually — the reader does the subtraction | free, and ignored |

If one headline figure must carry a representation × task interaction: the **parity
scatter** is the most honest (it cannot hide per-task reversals), the **heatmap + bar** is
the most scalable (>2 levels), and the **Skinnider apparatus** is the most credible but has
to be designed in before any model is trained. Add the **controls table** (HELM-BERT
Table 1) whichever you pick — it is cheap and it closes the undertuned-baseline objection.

## SECTION 4 — REGRET AND SUCCESS-RATE VISUALISATION

### 4.1 Ehrlich functions / holo-bench — Stanton, Alberstein, Frey, Watkins, Cho (2024)

Exact title: **"Closed-Form Test Functions for Biophysical Sequence Optimization
*Algorithms*"**. https://arxiv.org/abs/2407.00236 (2024, v1 only) · ICML 2024 Machine
Learning for Life and Material Sciences workshop · https://github.com/prescient-design/holo-bench

**Six figures; only two are results plots.** Metric definitions (§4.2, verbatim): simple
regret `r_t = f* − f(x̂*_t)`; cumulative regret `R_t = Σ_{j=1..t} r_j`. **`f* = 1` by
construction** for feasible sequences (`−∞` for infeasible), so regret has a true zero.
Feasibility = fraction of the GA population satisfying the discrete-Markov-process
constraint at iteration t. And the one sentence that governs every plot: **"In all plots
we show the 10%, 50%, and 90% quantiles of the reported performance metrics, estimated
from 32 trials."**

- **Fig. 4** — four side-by-side panels. **y = Simple Regret** (linear, 0.0–1.0);
  **x = Function Evaluations (M)**, millions, linear, 0–60 / 0–300 / 0–500 / 0–500
  depending on panel. One colour per parameter value; legends *Sequence Length*
  (128/256/512), *Motif Count* (4/8/16), *Motif Length* (4/8/10), *Objective Quantization*
  (2/4/8). Curves are monotone-decreasing **staircases** — regret only drops at discrete
  quantization levels. Median line plus **shaded 10–90% quantile band over 32 trials**
  (not a CI, not standard error). **Optimizer held fixed** at the GA baseline.
- **Fig. 5** — three panels sharing **x = Function Evaluations (M)**, 0–500M:
  (1) **y = Simple Regret** (0.3–1.0), (2) **y = Cumulative Regret** (0–300),
  (3) **y = Feasible Particles** (0–1.0). Two series (configuration A vs B), median +
  10–90% quantile ribbon over 32 trials. **Function held fixed** (k=8, q=4).
- Figs. 1–3 and 6 are non-quantitative: Ackley surface render, antibody/epitope cartoon,
  epistasis cartoon, Python code listing.

Repo-level confirmation: `scripts/benchmark_optimizer.py` logs `simple_regret_best`,
`simple_regret_last`, `cumulative_regret`, `frac_particles_feasible`, `timestep` to W&B
and **stops early on `simple_regret_best == 0`**. There is **no plotting module in
holo-bench** — the published figures came from the W&B dashboard. So "success" is encoded
as a *solve time* (regret hits exactly 0), not as a plotted success-rate curve.

**Convention: simple-regret-vs-budget trace with quantile ribbons** ("known-optimum
anytime regret curve"). Answers: how fast does the optimizer close the gap to a known
`f*`, and how reliably across seeds.
**Failure mode:** with `f*` known and bounded at 1, curves saturate and go uninformative
once a method stalls (Fig. 4's hardest line is flat at 1.0 for 500M evaluations — all you
learn is "failed"). Quantile ribbons also hide the bimodality the text itself admits:
"Configuration A outperformed B on most random seeds, but some seeds show essentially no
improvement."

### 4.2 FLEXS / AdaLead — Sinai, Wang, Whatley, Slocum, Locane, Kelsic (2020)

https://github.com/samsinai/FLEXS · paper https://arxiv.org/abs/2010.02141 (2020).
FLEXS has no separate paper; the README directs you to cite AdaLead.

`flexs/evaluate.py` defines **three sweep functions and no plotting at all**:
`robustness(signal_strengths=[0, 0.5, 0.75, 0.9, 1])` against degraded-fidelity surrogate
oracles; `efficiency(budgets=[(100,500),(100,5000),(1000,5000),(1000,10000)])` over
`(sequences_batch_size, model_queries_per_batch)`; `adaptivity(num_rounds=[1,10,100])` at
a fixed total budget. It returns `(setting, dataframe)` tuples; the user plots.

**AdaLead Figure 2 — four panels, four different conventions**, and the de-facto FLEXS
figure. Caption: "We record the cumulative maximum over all sequences generated by each
algorithm when run with 10 batches of size 100… Scores are normalized to top known or top
imputed."

- **A — box + strip plot.** y = **"Cumulative max y"** (0.94–1.00); x = categorical
  algorithm. Box = quartiles over 13 initializations × 5 TF landscapes, individual points
  overplotted.
- **B — the distinctive FLEXS plot: a robustness curve.** x = **α (surrogate model
  quality), 0.0 → 1.0**; y = **mean cumulative max y**. One line per algorithm.
  **No bands or error bars.** Tests "does the method actually benefit from a better model,
  and does it survive a useless one (α = 0)".
- **C — anytime curve.** x = **# samples (ground-truth oracle calls), 1→1000**;
  y = **cumulative max y** (0–1.0); **shaded band over 5 initializations**. Two stacked
  sub-panels: perfect oracle (α=1) above, 3×CNN ensemble below.
- **D — grouped bar chart with black error bars.** y = cumulative max y; x = landscape
  class, grouped by algorithm, split into ENS-3xCNN and α=1 blocks.
- **Table 1** — number of distinct local optima found above thresholds y_τ ∈ {0.75, 0.9,
  1.0}, with `# peaks` (brute-forced total) as the denominator. This is FLEXS's
  success/diversity count, not a success-rate curve.

**Key difference from Ehrlich:** FLEXS plots **cumulative-max fitness normalised to the
known/imputed optimum (y* = 1)**, *not* regret.
**Failure mode:** cumulative-max is monotone and ceiling-bounded, so near-optimal methods
become visually indistinguishable (panel A spans only 0.94–1.00); and "normalized to top
known **or top imputed**" means the denominator is sometimes a guess, silently changing
the y-scale across landscapes.

### 4.3 COCO / BBOB — the bootstrapped ECDF of runtimes

Origin: Hansen, Auger, Ros, Mersmann, Tušar, Brockhoff, **"COCO: A Platform for Comparing
Continuous Optimizers in a Black-Box Setting"**, *Optimization Methods and Software*
36(1):114–144, https://arxiv.org/abs/1603.08785 (arXiv 2016, journal 2021). Measurement
conventions are in the companion **"COCO: Performance Assessment"**,
https://arxiv.org/abs/1605.03560 / http://numbbo.github.io/coco-doc/perf-assessment/ (2016).

- **Runtime** = the number of function evaluations needed to reach a quality-indicator
  target for the first time. CPU time is deliberately avoided.
- **Targets** = `I^target,θ = I^ref,θ + ΔI`, where `I^ref` is the **known optimal function
  value**. The standard grid is **51 targets uniform on a log scale between 10^+2 and
  10^−8** (ratio 10^0.2 between neighbours). This is literally a distance-to-optimum
  ladder re-expressed as "which precision levels have you reached".
- **aRT (average runtime)** = `(Σ RT^s_i + Σ RT^us_j) / n_s = #FEs / n_s` — total function
  evaluations over **all** trials divided by the number of **successful** ones. Estimates
  the expected runtime of an idealised restart algorithm,
  `E(RT) = E(RT^s) + ((1−p_s)/p_s)·E(RT^us)`. aRT is the renamed successor of **ERT**
  (expected running time) in older BBOB papers. Caveat from the source: averaging is only
  meaningful if instances have similar, non-heavy-tailed distributions.
- **ECDF-of-runtimes plot:** **x = budget in function evaluations (log scale, usually
  normalised by dimension, "evaluations / dimension"); y = fraction of (function,
  instance, target) problems solved within that budget**, 0→1. "for any budget we see the
  fraction of problems solved within the budget as y-value."
- **Aggregation** over the 51 target precisions **and** 15 instances **and** (at suite
  level) all 24 functions — but **never over dimension**, deliberately, because dimension
  is a legitimate algorithm-selection variable.
- **Bootstrapped / simulated-restarts ECDF:** repeatedly draw uniformly with replacement
  from the K trials until a success is encountered, summing evaluation counts; repeat
  hundreds-to-thousands of times. This is how unsuccessful runs enter the plot instead of
  being dropped. A **cross marker** on the curve indicates at least one unsuccessful run.
- **Variance** is not drawn as a band — the spread *is* the ECDF; uncertainty enters via
  the bootstrap over trials. Aggregation is over **both** seeds/instances and
  tasks/targets.

**Failure mode**, from COCO's own documentation: when runs from several instances are
aggregated "the association to the single run is lost, as is the association to the
function when aggregating over several functions" — a method excellent on 6 functions and
hopeless on 18 can trace the same ECDF as one that is mediocre everywhere. Easy targets
(10^+2) are reached immediately by anything, so the left half is near-uninformative, and
the choice of target grid directly determines the curve's shape. See also Hansen, Auger,
Brockhoff, Tušar, "Anytime Performance Assessment in Blackbox Optimization Benchmarking",
*IEEE TEVC* 26(6):1293–1305 (2022).

### 4.4 Bayesian-optimization / HPO benchmarks

**BayesMark** (uber/bayesmark, Apache-2.0, 2019/2020) —
https://bayesmark.readthedocs.io/en/latest/scoring.html
- Normalised against **both** the known optimum and random search:
  `norm-mean-perf = (mean-perf − opt_p) / (clip_p − opt_p)`, where `opt_p` is the
  estimated global minimum and `clip_p` the median score after a single evaluation; scores
  clipped to [−1, 1]. Random search then "performs as a straight line at 1 for all t";
  **0 = optimum, 1 = random-search level, lower is better.**
- Axes (from `notebooks/plot_mean_score.ipynb`): `xlabel("evaluation")`,
  `ylabel("normalized median score")`, plus a second figure with `ylabel("mean score")`.
- **Variance:** `fill_between(LB, UB)` bands. From `bayesmark/experiment_analysis.py`, the
  **median** band is a distribution-free order-statistic quantile CI (`quantile_and_CI`)
  and the **mean** band is a **t-distribution CI** (`t_EB`, α = 0.05). Per-problem bands
  are computed **across trials/seeds** (`axis=1`); the aggregate band is a t-interval
  **across test cases/tasks**.
- Failure mode: t-intervals assume approximate normality of a mean over a small,
  heterogeneous task set; and the random-search denominator moves the y-scale if that
  baseline is re-estimated.

**HPO-B** — Pineda Arango, Jomaa, Wistuba, Grabocka, https://arxiv.org/abs/2106.06257
(2021, NeurIPS D&B)
- **Normalized regret** = `min_{x ∈ X_e}(f(x) − y*_min) / (y*_max − y*_min)`. Because the
  benchmark is tabular, the optimum is *known exactly* as the best cell.
- Figures: **x = trial number (1–100, excluding 5 seed configurations); y = normalized
  regret, or mean rank**, linear scale. Aggregated over 5 seeds per task, then tasks, then
  search spaces.
- Notably, it **also reports critical-difference diagrams at trials 25, 50 and 100** — CD
  diagrams used as the significance layer on top of regret-vs-trial curves. This is the
  cleanest precedent for combining §4 and §5 conventions in one paper.

**HPOBench** — Eggensperger, Müller, Mallik, Feurer, Sass, Klein, Awad, Lindauer, Hutter,
https://arxiv.org/abs/2109.06716 (2021, NeurIPS D&B). The earlier tabular-benchmark paper
is Klein & Hutter, "Tabular Benchmarks for Joint Architecture and Hyperparameter
Optimization", https://arxiv.org/abs/1905.04970 (2019), which introduced exhaustively
enumerated grids where the global optimum is known and therefore plotted as *regret*.
- Convention: **x = wall-clock time or function evaluations (log); y = validation regret
  `f(x̂) − f*` (log)**, variance as **quantile bands over repeated runs**. *(Medium
  confidence on the exact band definition — the PDF resisted text extraction.)*

**YAHPO Gym** — Pfisterer, Schneider, Moosbauer, Binder, Bischl, PMLR 188 (AutoML 2022),
https://proceedings.mlr.press/v188/pfisterer22a/pfisterer22a.pdf · https://arxiv.org/abs/2109.03670
- Figure 2 caption, verbatim: "Mean normalized regret (top) and mean ranks (bottom) of
  different HPO methods on different benchmarks. **Ribbons represent standard errors.**
  The gray vertical line indicates the cumulative budget used for the initial design of BO
  methods… **30 replications.**"
- **Two stacked rows — y = mean normalized regret and y = mean rank; x = cumulative
  budget.** The **vertical rule marking the end of BO's random initial design** is a detail
  most papers omit and then over-interpret; worth copying.
- Explicit recommendation (§E): "We encourage reporting mean normalized regret and mean
  ranks for the anytime performance of an optimizer… In order to assess variance, we
  encourage reporting averages and standard errors across **30 replications** with
  differing random seeds." Multi-objective: normalized hypervolume indicator.
- Their own finding, worth flagging for anyone building a surrogate benchmark: **tabular
  benchmarks distort the ranking** — surrogate-vs-real consensus rankings differ by
  permutation order 2, tabular-vs-real by order 5.

### 4.5 Dolan–Moré performance profile — the origin of "performance profile"

Dolan & Moré, **"Benchmarking optimization software with performance profiles"**,
*Mathematical Programming* 91(2):201–213 (2002), https://arxiv.org/abs/cs/0102001

- Performance ratio `r_{p,s} = t_{p,s} / min{t_{p,s} : s ∈ S}` (∞ if solver s fails on p).
- Profile `ρ_s(τ) = (1/n_p)·|{p : r_{p,s} ≤ τ}|`.
- **x = τ, the factor-of-the-best threshold, conventionally log₂; y = ρ_s(τ) ∈ [0,1], the
  fraction of problems on which solver s is within a factor τ of the best solver.**
- `ρ_s(1)` = fraction of problems where s is the outright winner; `ρ*_s = lim_{τ→∞} ρ_s(τ)`
  = fraction it solves at all (its reliability).
- **Failure mode — Gould & Scott, "A note on performance profiles for benchmarking
  software", *ACM TOMS* 43(2), article 15, https://doi.org/10.1145/2950048 (2016).** With
  more than two solvers you can read off who is most likely within a factor τ of the best,
  but **you cannot read off the relative performance of two non-best solvers**. In their
  artificial 5-problem/3-solver example, Solver B wins on *nothing* when A is in the pool,
  yet becomes best on 60% of the set once A is removed. Recommended workaround: a
  *sequence* of profiles with the leading solver removed each time. Verbatim conclusion:
  "performance profiles must be used with care."

### 4.6 Reinforcement-learning conventions

**Human-normalized score (Atari).** Mnih et al., *Nature* 518:529–533 (2015),
https://www.nature.com/articles/nature14236, on Bellemare et al.'s ALE, *JAIR* (2013),
https://arxiv.org/abs/1207.4708.
- `100 × (score_agent − score_random) / (score_human − score_random)`. 0 = random policy,
  100 = professional human tester.
- Learning-curve convention: **x = environment frames/steps; y = median (or mean)
  human-normalized score across 57 or 26 games.** Mean over episodes within a run, then
  median across games. Historically a **point estimate with no interval at all**, from 1–5
  runs.
- This is a *known-reference* rather than known-optimum normalization — "human" is not the
  optimum, which is why scores above 100% are routine.
- Failure mode: the mean is dominated by a handful of games with astronomical normalized
  scores (Atlantis, Video Pinball); the median is insensitive — zero on nearly half the
  games does not move it. Exactly what rliable was written to fix.

**rliable's three plot forms** (Agarwal et al. 2021, https://arxiv.org/abs/2108.13264) —
the run-score distribution, the interval-estimate bar, and the IQM sample-efficiency curve
— are specified in §5.2 above. The one point to carry here: rliable's **sample-efficiency
curve** is the direct replacement for the Atari convention: **x = training frames, y = IQM
human-normalized score, shaded = pointwise 95% percentile stratified-bootstrap CI**,
aggregating over **tasks and runs simultaneously**. And its performance profile is the
**tail** distribution (higher is better), the mirror of Dolan–Moré's CDF, and thresholds
against an **absolute normalized score** rather than against the best competitor.

### 4.7 Critical difference diagrams

Origin Demšar, *JMLR* 7:1–30 (2006), https://jmlr.org/papers/v7/demsar06a.html —
construction, anatomy and known defects are in §5.1. For Section 4 purposes: **x = average
rank (no y-axis); aggregation over tasks/datasets only; variance is not drawn at all**,
being folded into a single CD bar. Modern usage: time-series classification (Ismail Fawaz
et al., *DMKD* 33:917–963, 2019, https://arxiv.org/abs/1809.04356 — their widely-copied
`cd-diagram` code swaps Nemenyi for **Wilcoxon signed-rank with Holm correction**), HPO-B
at trials 25/50/100, and the `autorank` (Herbold, *JOSS* 2020,
https://doi.org/10.21105/joss.02173), `scikit-posthocs`, `aeon.visualisation.plot_critical_difference`
and `mlr` implementations.

### 4.8 Summary table — convention, origin, axes, when to use

| Convention | Origin | x | y | Use when |
|---|---|---|---|---|
| **Simple-regret trace with quantile ribbons** | Stanton et al. 2024 (arXiv 2407.00236) | function evaluations | simple regret, f*−f(x̂) | the optimum is known in closed form |
| **Cumulative-regret trace** | classical bandit/BO; here Stanton et al. 2024 | function evaluations | Σ regret | you care about the cost of the whole search, not just the best point |
| **Cumulative-max-fitness curve, normalised to y\*=1** | Sinai et al. 2020 (arXiv 2010.02141) | oracle calls | best-so-far, normalised | the optimum is estimated rather than known |
| **Robustness-vs-oracle-fidelity curve** | Sinai et al. 2020, FLEXS `evaluate.robustness` | surrogate quality α | mean cumulative max | the method depends on a learned surrogate |
| **Bootstrapped ECDF of runtimes** | Hansen et al., COCO (arXiv 1603.08785 / 1605.03560) | budget (log, FEs/dim) | fraction of (function, instance, target) solved | many tasks, many precision levels, known f* |
| **aRT / ERT** | COCO | — | expected restart runtime | you need one number per (function, target) |
| **Normalized-regret-vs-trial + mean-rank panel** | HPO-B (arXiv 2106.06257), YAHPO Gym (PMLR 188) | trial / cumulative budget | normalized regret; mean rank | tabular or surrogate HPO where f* is the best table cell |
| **Dolan–Moré performance profile** | Dolan & Moré 2002 (arXiv cs/0102001) | τ, factor of the best (log₂) | fraction of problems within τ of best | comparing solvers on runtime over a problem suite |
| **rliable run-score distribution** | Agarwal et al. 2021 (arXiv 2108.13264) | τ, normalized score | fraction of **runs** above τ | you have multiple runs on multiple tasks and want the whole distribution |
| **IQM interval estimate / sample-efficiency curve** | Agarwal et al. 2021 | — / training budget | IQM with 95% stratified bootstrap CI | any aggregate claim across tasks and runs |
| **Critical difference diagram** | Demšar 2006 (JMLR) | average rank | — | many methods × many tasks, one score each, significance wanted |

### 4.9 The structural pattern worth planning around

Every "aggregate many tasks into one curve by comparing against the best" convention here
— Dolan–Moré profiles, CD diagrams, and to a lesser degree COCO ECDFs — has **the same
documented defect: the result for a pair of methods depends on which *other* methods or
problems are in the pool.** Gould & Scott 2016 for profiles; Benavoli et al. 2016 and
Ismail-Fawaz et al. 2023 for CD diagrams; COCO's own docs for ECDF aggregation.

rliable's run-score distribution is the only one in this list that is **pool-independent by
construction**, because it never ranks or ratios against the best competitor — it
thresholds against an absolute normalized score. **Benchmarks that know `f*` (Ehrlich,
tabular HPO, COCO) can do the same thing: threshold against distance-to-optimum rather than
against the best competing method.** For a benchmark paper whose selling point is a known
optimum, this is the argument that justifies the whole design.

### 4.10 Variance: bands, bars, or bootstrap? Over seeds or tasks?

| Convention | Variance shown as | Over what |
|---|---|---|
| Ehrlich / holo-bench Figs. 4–5 | shaded **10/50/90% quantile ribbon**, no CI | seeds only (32 trials), one task per curve |
| FLEXS/AdaLead Fig. 2C | shaded band | initializations/seeds (5 or 13) |
| FLEXS/AdaLead Figs. 2A, 2D | boxplot / black error bars | initializations × landscapes (both) |
| FLEXS/AdaLead Fig. 2B | **nothing** | — |
| COCO/BBOB ECDF | no band; the spread *is* the curve; bootstrap over trials for simulated restarts | both — instances (seeds) and functions × targets (tasks) |
| COCO aRT | point estimate (sometimes bootstrapped percentile ticks) | instances, per function × target |
| BayesMark | `fill_between`: **t-CI** for mean, order-statistic quantile CI for median | per-problem = seeds; aggregate = **tasks** |
| HPO-B | mean normalized regret / mean rank curves + CD diagrams | 5 seeds × tasks × search spaces |
| HPOBench | quantile bands | repeated runs (seeds) |
| YAHPO Gym | ribbons = **standard errors**, 30 replications | seeds; the rank panel folds in tasks |
| Atari HNS curves (classic) | usually none | 1–5 runs |
| rliable | **pointwise 95% percentile stratified-bootstrap CI bands** and horizontal CI bars | **both**, by construction |
| Dolan–Moré profile | none | tasks only |
| CD diagram | none (one CD bar) | tasks/datasets only |

## SECTION 5 — STATISTICAL PRESENTATION NORMS

**Headline: no methodology source in this literature endorses a mean ± stddev bar
chart.** Every one replaces it with either (a) rank-based significance diagrams
(Demšar lineage) or (b) interval estimates plus whole distributions (rliable lineage).

### 5.1 Demšar, "Statistical Comparisons of Classifiers over Multiple Data Sets", JMLR 7:1–30 (2006)

https://jmlr.org/papers/v7/demsar06a.html · PDF https://www.jmlr.org/papers/volume7/demsar06a/demsar06a.pdf (2006)

**Seven figures.** Figures 1 and 2 are the CD diagrams; Figures 3–7 are an empirical
power/replicability study, not presentation advice.

- **Figure 1** — "Visualization of post-hoc tests for data from Table 6", two panels.
  (a) Nemenyi, all-pairs: "Groups of classifiers that are not significantly different
  (at p = 0.10) are connected." (b) Bonferroni–Dunn, one control vs the rest: "All
  classifiers with ranks outside the marked interval are significantly different
  (p < 0.05) from the control."
- **Figure 2** — a six-method CD diagram for feature-selection measures.
- **Figures 3–7** — power and replicability curves: t-test vs Wilcoxon (3),
  replicability R(p) vs R(e) (4), ANOVA vs Friedman (5), Tukey vs Nemenyi (6),
  Bonferroni-Dunn/Holm/Hochberg vs Dunnett (7).

**How a CD diagram is built**

1. Rank algorithms **within each dataset** (best = 1), average ranks for ties;
   R_j = (1/N) Σ_i r_i^j.
2. **Friedman test** χ²_F = [12N/k(k+1)]·[Σ_j R_j² − k(k+1)²/4], k−1 df. Rule of thumb
   N > 10 and k > 5.
3. **Iman–Davenport correction** (the paper insists on it): Friedman's χ²_F is
   "undesirably conservative"; use F_F = (N−1)χ²_F / [N(k−1) − χ²_F], F-distributed with
   k−1 and (k−1)(N−1) df.
4. If H₀ rejected: **Nemenyi** for all-pairs, **CD = q_α·√(k(k+1)/6N)** with q_α "based
   on the Studentized range statistic divided by √2" (q_0.05 = 1.960, 2.343, 2.569,
   2.728, 2.850, 2.949, 3.031, 3.102, 3.164 for k = 2…10). Or **Bonferroni–Dunn** for one
   control vs the rest (same formula, critical values at α/(k−1): 1.960, 2.241, 2.394,
   2.498, 2.576, 2.638, 2.690, 2.724, 2.773). Explicit: "We thus should not make pairwise
   comparisons when we in fact only test whether a newly proposed method is better than
   the existing ones."
5. **Drawing it:** "The top line in the diagram is the axis on which we plot the average
   ranks of methods. The axis is turned so that the lowest (best) ranks are to the right
   since we perceive the methods on the right side as better." There is **no y-axis**.
   The **horizontal bar connects groups that are *not* significantly different** — a
   clique bar, **not an error bar**, carrying no variance information. The CD itself is
   drawn as a labelled segment above the axis. For Bonferroni–Dunn, mark one CD left and
   right of the control's rank; anything outside differs.
6. Demšar declines to draw per-classifier adjusted intervals for Holm/Hochberg — it
   "would need to plot a different adjusted critical interval for each classifier…
   which could easily become confusing."

**What he argues against (verbatim)**

- **Averaging accuracies across datasets.** Quoting Webb (2000): "it is debatable whether
  error rates in different domains are commensurable". "If the results on different data
  sets are not comparable, their averages are meaningless." "Averages are also
  susceptible to outliers. They allow classifier's excellent performance on one data set
  to compensate for the overall bad performance."
- **Paired t-test across datasets.** Three named defects: commensurability ("using the
  paired t-test for comparing a pair of classifiers makes as little sense as computing
  the averages over data sets"); normality (needs ~30 datasets, and with few you are
  "unlikely to detect abnormalities"); outlier sensitivity, which "skew[s] the test
  statistics and decrease[s] the test's power".
- Framing: "We shall warn against the widely used t-test as usually conceptually
  inappropriate and statistically unsafe."

**What he recommends instead**

- Two methods over many datasets → **Wilcoxon signed-ranks test** ("more sensible than
  the t-test… outliers have less effect"); sign test as a weaker fallback.
- Multiple methods → **Friedman + Iman–Davenport, then Nemenyi or Bonferroni–Dunn,
  presented as a CD diagram.** "CD diagrams are 'space-friendly'… yet they present the
  order of the algorithms, the magnitude of differences between them (in terms of ranks)
  and the significance of the observed differences much more clearly than it can be done
  in textual or in a pure numerical form."
- "the actual experiments should be conducted on as many data sets as possible."
- Non-dogmatic closer: "statistical tests should not be the deciding factor for or
  against publishing the work."

**Known defect of the CD diagram.** Benavoli, Corani & Mangili, "Should we really use
post-hoc tests based on mean-ranks?", *JMLR* 17(5), 2016,
https://www.jmlr.org/papers/v17/benavoli16a.html — "the outcome of the mean-ranks test
depends on the pool of algorithms originally included in the experiment… the difference
between A and B could be declared significant if the pool comprises algorithms C,D,E and
not significant if the pool comprises algorithms F,G,H. To overcome these issues, we
suggest instead to perform the multiple comparison using a test whose outcome only
depends on the two algorithms being compared, such as the **sign-test or the Wilcoxon
signed-rank test**." Practically: **adding or removing a weak baseline can flip a
significance claim on your CD diagram.** See also Ismail-Fawaz, Dempster, Tan, Webb et
al., https://arxiv.org/abs/2305.11921 (2023) — CD diagrams "are open to both inadvertent
and intentional manipulation"; they propose the **Multiple Comparison Matrix (MCM)**, a
pairwise matrix of mean-score difference, win/draw/loss counts and Wilcoxon p-values.
Also García & Herrera, *JMLR* 9:2677–2694 (2008),
https://jmlr.org/papers/v9/garcia08a.html, for more powerful n × n procedures and
"adjusted and comparable p-values".

Clique bars are additionally **non-transitive in display**: A can share a bar with B,
B with C, yet A be significantly better than C. The diagram cannot render that, and
readers habitually misread the bars as equivalence classes.

### 5.2 Agarwal, Schwarzer, Castro, Courville, Bellemare, "Deep RL at the Edge of the Statistical Precipice", NeurIPS 2021 (Outstanding Paper)

https://arxiv.org/abs/2108.13264 (2021) · library https://github.com/google-research/rliable

**Twelve main figures**, in order:

1. Number of runs used in RL papers over time — ≤5 runs has been common since DQN.
2. (L) sampling distribution of median normalized score across subsampled runs;
   (R) 95% CI *width* for median vs IQM as a function of run count.
3. Expected sample median of task means — substantial bias, varying by algorithm and N.
4. 95% CIs for detecting score improvements, median (L) vs IQM (R).
5. Normalized DER scores — non-standard evaluation protocols explain gains attributed to
   SUNRISE and CURL.
6. **Coverage check** of 95% stratified bootstrap CIs for median and IQM across run counts.
7. Performance profiles on Atari 100k: **score distributions (left, recommended) vs
   average-score distributions (right)**, with 95% bands.
8. The "anatomy" figure: **IQM = red shaded region, optimality gap = orange shaded region
   of a performance profile.**
9. Atari 200M: **interval estimates with 95% CIs for Mean / Median / IQM**, 7 algorithms,
   55 games.
10. (L) Atari 200M score distributions; (R) **IQM sample-efficiency curve** vs training
    frames with 95% bands.
11. DM Control: (a) interval estimates, (b) score distributions, (c) rank /
    probability-of-improvement comparisons, 9 algorithms.
12. Procgen: (L) score distributions; (R) probability of improvement, 6 algorithms, 16 tasks.

**Definitions, as stated**

- **IQM** — "discards the bottom and top 25% of the runs and calculates the mean score of
  the remaining 50% runs." Computed over **all runs × tasks pooled**, not per-task means.
  It is a 25%-trimmed mean.
- **Stratified bootstrap CI** — "re-samples runs with replacement **independently for each
  task**… from which we calculate a statistic and repeat this process many times."
  Figure 6 validates coverage: "percentile CIs provide good interval estimates for as few
  as N = 10 runs for both median and IQM scores."
- **Performance profile / run-score distribution** — F̂_X(τ) = (1/M)·Σ_m (1/N)·Σ_n
  1[x_{m,n} > τ]: the **empirical tail distribution** of a random score (higher curve is
  better), with pointwise percentile-bootstrap confidence bands. Unbiased, a step function
  in 1/(MN) rather than 1/M, and a single outlier run can move it by at most 1/(MN) —
  unlike the older **average-score distribution** (fraction of *tasks* whose mean exceeds
  τ), which is biased. Readable facts: "The τ value where the profiles intersect y = 0.5
  shows the median while for a non-negative random variable, area under the performance
  profile corresponds to the mean."
- **Probability of improvement** — "how likely it is for X to outperform Y on a randomly
  selected task", averaged over tasks, with CIs. Explicitly does *not* account for effect size.
- **Optimality gap** — "the amount by which the algorithm fails to meet a minimum score of
  γ = 1.0", i.e. mean of max(0, γ − x). The RL analogue of regret against a declared
  "good enough" target.

**What it says about point estimates, metric choice and run counts**

- Point estimates "ignore the statistical uncertainty implied by the use of a finite
  number of training runs" and "evade the question: Would similar findings be obtained
  with new independent runs?"
- **Median**: "only depends on the performance ordering across tasks and not on the
  magnitude except at most 2 tasks"; "zero performance on nearly half of the tasks does
  not change it"; statistically inefficient. **Mean**: "often dominated by performance on
  a few outlier tasks." **IQM**: "considerably less bias than median" while "robust to
  outliers."
- Contradicts "folk wisdom in experimental RL suggesting that 20 or 30 runs are enough":
  in their case study, "statistically defensible improvements with median scores is only
  achieved for 25 runs." The Atari-100k literature they audit used 3, 5, 10 or 20.
- "we emphasize the importance of published papers providing results for **all runs** to
  allow for future statistical analyses."
- Overlapping CIs are **not** a no-difference verdict: compute CIs for score *differences*
  (their Figure A.16).

**What rliable produces by default** — the README is an explicit three-row replacement
table for the bar-chart convention:

| Desideratum | Criticised current approach | Recommendation |
|---|---|---|
| Uncertainty in aggregate performance | "**Point estimates**: Ignore statistical uncertainty. Hinder results reproducibility" | "Interval estimates using **stratified bootstrap confidence intervals**" |
| Variability across tasks and runs | "**Tables with task mean scores**: Overwhelming beyond a few tasks. Standard deviations frequently omitted. Incomplete picture for multimodal and heavy-tailed distributions" | "**Score distributions** (performance profiles)… Easily read any score percentile" |
| Aggregate metric | "**Mean**: often dominated by outlier tasks. **Median**: statistically inefficient… 0 scores on nearly half the tasks doesn't change it" | "**Interquartile Mean (IQM)** across all runs… report *Probability of improvement* and *Optimality gap*" |

API defaults: `rly.get_interval_estimates(..., reps=50000)` (2000 for probability of
improvement); `plot_utils.plot_interval_estimates` renders Median / IQM / Mean /
Optimality Gap as **horizontal interval bars with 95% stratified bootstrap CIs**;
`create_performance_profile` + `plot_performance_profiles`; `plot_sample_efficiency_curve`;
`plot_probability_of_improvement`.

Self-criticism worth knowing: IQM discards 50% of runs, so a method excellent on a quarter
of tasks and catastrophic on a quarter looks identical to a uniformly mediocre one; the
optimality gap depends entirely on an arbitrary γ; and the performance profile inherits
Dolan–Moré's crossing problem (§4.5).

### 5.3 NeurIPS checklist and Evaluations & Datasets track guidance

**Paper Checklist (version live as of this writing: NeurIPS 2026)** —
https://neurips.cc/public/guides/PaperChecklist

- **Item 1, Claims** — "Do the main claims made in the abstract and introduction
  accurately reflect the paper's contributions and scope?" Contributions must be stated
  with important assumptions and limitations; aspirational goals are acceptable only if
  clearly marked as unattained.
- **Item 4, Experimental Result Reproducibility** — "If the contribution is a dataset or
  model, what steps did you take to make your results reproducible or verifiable?"
  NeurIPS "does not require releasing code" but requires "some reasonable avenue for
  reproducibility."
- **Item 6, Experimental Setting/Details** — "did you specify all the training details
  (e.g., data splits, hyperparameters, how they were chosen)?"
- **Item 7, Experiment Statistical Significance — the load-bearing one.** "Does the paper
  report **error bars suitably and correctly defined** or other appropriate information
  about the statistical significance of the experiments?" The guidance bullets require you
  to: include error bars, CIs or statistical tests for main claims; **state what factor of
  variability the error bars capture** (splits, initialization, overall run); **explain how
  they were calculated** (closed form, library call, bootstrap); **state assumptions** such
  as normality; **distinguish standard deviation from standard error**; one-sigma
  acceptable only if stated, **two-sigma preferred when normality is not verified**; and
  **avoid symmetric error bars that produce impossible values** (accuracy > 1 or < 0).

  *This is the hook.* A mean ± std bar chart with no stated variability factor, no
  calculation method, and symmetric whiskers crossing a bound is a **direct item-7
  failure**, independent of anything else in the paper.

**Evaluations & Datasets track (the renamed D&B track), NeurIPS 2026** —
https://neurips.cc/Conferences/2026/EvaluationsDatasetsReviewerGuidelines ·
CFP https://neurips.cc/Conferences/2026/CallForEvaluationsDatasets ·
announcement https://blog.neurips.cc/2026/03/23/introducing-the-evaluations-datasets-track-at-neurips-2026/

- Four criteria: **Quality** ("technical soundness, well-supported claims, completeness,
  honest evaluation of strengths and weaknesses"), **Clarity** (incl. "sufficient detail
  for reproducibility"), **Significance**, **Originality** (explicitly not requiring a new
  method).
- **Double-blind is now the default**; single-blind needs a scientific/ethical justification.
- Per-contribution-type language that bears directly on evidence: a dataset must clarify
  "its evaluative role: **what specific claims it supports, under what assumptions, and
  what limitations apply**"; benchmark design is judged on whether it "address[es] a
  meaningful gap or enable[s] **reliable comparison across methods**";
  **evaluation-methodology papers are judged on statistical soundness**; reproducibility
  and stress-testing papers must be "grounded in rigorous, systematic analysis" with
  negative results "rigorously supported".

### 5.4 Other methodology sources

**Henderson, Islam, Bachman, Pineau, Precup, Meger, "Deep Reinforcement Learning That
Matters", AAAI 2018** — https://arxiv.org/abs/1709.06560 (2017/2018)

- The seed experiment (**Figure 5**): 10 trials, *same* hyperparameters, varying only the
  random seed; split into two groups of 5 and averaged. Caption: "TRPO on HalfCheetah-v1
  using the same hyperparameter configurations averaged over two sets of 5 different
  random seeds each. The average 2-sample t-test across entire training distribution
  resulted in **t = −9.0916, p = 0.0016**." I.e. **two groups of the same algorithm come
  out "significantly different" from each other.**
- "even averaging several learning results together across totally different random seeds
  can lead to the reporting of misleading results"; "it is not uncommon for the top-N
  trials to be selected from among several trials or averaged over only small number of
  trials (N<5)… Our experiment with random seeds shows that this can be potentially
  misleading."
- Their replacement: **Table 3, "Bootstrap mean and 95% confidence bounds… 10k bootstrap
  iterations and the pivotal method"** — e.g. HalfCheetah-v1 DDPG 5037 (3664, 6574) vs PPO
  3043 (1920, 4165). Overlapping intervals a bar chart would have hidden. Note they report
  standard *error*, not std, in their own curves, and show it is still insufficient.

**Bouthillier et al., "Accounting for Variance in Machine Learning Benchmarks", MLSys
2021** — https://arxiv.org/abs/2103.03098

- Models data sampling, augmentation, initialization **and the whole hyperparameter-
  optimization process** as variance sources. "there are other, larger, sources of
  uncontrolled variation [than weight init] and the risk is that conclusions are driven by
  differences due to arbitrary factors, such as data order, rather than model improvements."
- Three recommendations, verbatim: **(1)** "As many sources of variation as possible should
  be randomized whenever possible… including weight initialization, data sampling, random
  data augmentation and the whole hyperparameter optimization." **(2)** "Deciding of whether
  the benchmarks give evidence that one algorithm outperforms another should not build
  solely on comparing average performance but account for variance. We propose a simple
  decision criterion based on requiring a high-enough probability that in one run an
  algorithm outperforms another." **(3)** "Resampling techniques such as out-of-bootstrap
  should be favored instead of fixed held-out test sets."
- The criterion: **ℙ(A > B) ≥ γ with γ = 0.75** (statistical significance:
  ℙ(A>B) − CI_min > 0.5; meaningfulness: ℙ(A>B) + CI_max > γ).
- Counter-intuitive and directly useful: a cheap estimator that randomizes *many* sources
  (one HOpt run, then k splits randomizing everything else, O(k+T) trainings) approaches
  the ideal O(k·T) estimator — ~21 h vs ~1,070 h in their experiments. **Adding more
  sources of randomization to a cheap estimator beats adding more seeds to a narrow one.**

**Colas, Sigaud, Oudeyer, "How Many Random Seeds?", 2018** — https://arxiv.org/abs/1806.08295

- **Figure 1 is exactly the artifact in question**: two algorithms' mean ± 95% CI over 5
  seeds, barely overlapping — "We might consider that Algo1 outperforms Algo2 because there
  is not much overlap… But is it sufficient evidence…? Below, we show that the performances
  of these algorithms is actually the same."
- Use **Welch's t-test**, not Student's — the equal-variance assumption "rarely holds when
  comparing two different algorithms". Rules out the KS test: it "is unable to prove any
  order relation."
- Sample sizes: the bootstrap CI test "should not be used with less than N = 20 samples";
  run "at least n = 20 samples in the pilot study"; "use larger sample size N than the one
  prescribed by the power analysis." Under non-normality the empirical false-positive rate
  exceeds nominal, so "tighter significance level should be used (<0.05)."
- Multiple comparisons: with N_E experiments the false-positive rate grows linearly; apply
  Bonferroni α/N_E.

### 5.5 Concretely, what to show instead of mean ± stddev bars

1. **Release per-run, per-task raw scores.** (Agarwal 2021; NeurIPS checklist item 4.)
2. **Aggregate with IQM, plus optimality gap — not mean, not median — drawn as horizontal
   interval estimates with 95% stratified bootstrap CIs.** (Agarwal 2021 / rliable.)
3. **A performance profile / score distribution with CI bands**, so the whole distribution
   rather than one moment is visible. (Agarwal 2021; echoes Demšar's outlier objection.)
4. **For pairwise claims: probability of improvement with CIs** (Agarwal 2021) or
   **ℙ(A > B) ≥ 0.75** (Bouthillier 2021). Never "the bars don't overlap".
5. **Across many tasks: Friedman + Iman–Davenport, then a CD diagram** (Demšar 2006) — with
   the Benavoli caveat that Nemenyi bars depend on the method pool, so prefer per-pair
   Wilcoxon/sign tests with adjusted p-values, or an MCM matrix.
6. **Randomize more than seeds** — splits, init, augmentation, HPO (Bouthillier 2021);
   budget N by power analysis; never ≤5 seeds (Colas 2018, Henderson 2018).
7. **In every caption: state what the interval is, what variability it captures, how it was
   computed, and whether it is 1σ or 2σ.** (NeurIPS checklist item 7.)
8. Where the optimum is known, **threshold against distance-to-optimum rather than against
   the best competing method** — the only aggregation in this corpus that is
   pool-independent by construction (see §4.8).

## SECTION 6 — WHAT GETS CRITICISED

### 6.1 Reviewer-facing checks that get a benchmark paper flagged

From the **NeurIPS 2026 Evaluations & Datasets reviewer guidelines**
(https://neurips.cc/Conferences/2026/EvaluationsDatasetsReviewerGuidelines) and the 2025
CFP and hosting pages (https://neurips.cc/Conferences/2025/CallForDatasetsBenchmarks,
https://neurips.cc/Conferences/2025/DataHostingGuidelines), the hard failure modes are:

- **Missing or invalid Croissant metadata.** Reviewers are explicitly told to flag it, and
  an **automated compliance report** is attached to the submission as a checklist — though
  the guidelines add that "expert judgment remains essential."
- **Data not reachable at submission** without contacting the PI, or not on a persistent
  ML-dataset host. "New datasets should be hosted at one of the hosting sites dedicated to
  ML datasets (Dataverse, Kaggle, Hugging Face, or OpenML) or at a bespoke hosting site";
  "Datasets and code should be available and accessible to all reviewers, ACs and SACs at
  the time of submission."
- **Missing or broken `code_URL`**, or code buried in supplementary rather than hosted,
  anonymized and documented. Code is **mandatory for tool/platform contributions** and
  "highly encouraged" elsewhere, with a code-submission justification demanded when absent.
- **Missing Responsible AI metadata** documenting "biases, limitations, intended uses, and
  sensitive information."
- **Anonymity**: double-blind is now default; single-blind requires a scientific or ethical
  justification on the submission form.
- **Ethics**: reviewers flag for ethics review when Code-of-Ethics issues arise and they
  cannot judge severity. Human-centred evaluation must document protocols and **"fair
  compensation for all participants."**
- **Evidence quality.** Quality = "technical soundness, well-supported claims… honest
  evaluation of strengths and weaknesses." A dataset must state "what specific claims it
  supports, under what assumptions, and what limitations apply." **Evaluation-methodology
  contributions are judged on statistical soundness.** Reproducibility/stress-testing work
  must be "grounded in rigorous, systematic analysis" with negative results "rigorously
  supported." Empirical analyses must have claims that are "independently verifiable."
- Reviewers may not use LLMs in reviewing; authors must disclose non-routine LLM/agent use.

Scale and compliance data from the **D&B chairs' retrospectives**
(https://blog.neurips.cc/2025/09/30/reflecting-on-the-2025-review-process-from-the-datasets-and-benchmarks-chairs/,
https://blog.neurips.cc/2025/12/05/neurips-datasets-benchmarks-track-from-art-to-science-in-ai-evaluations/):
**1,995 submissions in 2025** (1,820 in 2024). On initial submission, Croissant gaps were
**license 11.9%**, dataset description 4.9%, URLs 3.5% missing. These are the cheapest
possible desk-reject reasons and they still catch roughly one paper in eight.

### 6.2 Retrospectives on benchmark quality

**Raji, Bender, Paullada, Denton, Hanna, "AI and the Everything in the Whole Wide World
Benchmark", NeurIPS 2021 D&B** — https://arxiv.org/abs/2111.15366 (2021)
"There is a tendency across different subfields in AI to valorize a small collection of
influential benchmarks… frequently framed as foundational milestones on the path towards
flexible and generalizable AI systems… we explore the limits of such benchmarks in order to
reveal the **construct validity** issues in their framing as the functionally 'general'
broad measures of progress they are set up to be." Three criticisms: construct-validity
failure (the benchmark does not measure the construct claimed), overstated generality from
narrow task scope, and unjustified community valorization of a few benchmarks.
*Direct implication for a molecular-design benchmark: say what construct the oracle is a
proxy for, and show the proxy is valid — which is exactly what Gao & Coley's Fig. S14 does
and what GuacaMol/Tartarus do not.*

**Liao, Taori, Raji, Schmidt, "Are We Learning Yet? A Meta-Review of Evaluation Failures
Across Machine Learning", NeurIPS 2021 D&B** —
https://datasets-benchmarks-proceedings.neurips.cc/paper/2021/file/757b505cfd34c64c85ca5b5690ee5293-Paper-round2.pdf ·
https://openreview.net/forum?id=mPducS1MsEK (2021)
Meta-review of **107 survey papers** across NLP, recommender systems, CV, RL, computational
biology and graph learning. Taxonomy split into **internal validity** (improper baseline
comparisons, **insufficiently tuned baselines**, overfitting from test-set reuse) and
**external validity** (whether progress on one learning problem transfers to seemingly
related tasks). *"Insufficiently tuned baselines" is precisely the Tripp & Hernández-Lobato
finding in §2.3–2.4, independently catalogued as a field-wide failure mode.*

**Koch, Denton, Hanna, Foster, "Reduced, Reused and Recycled: The Life of a Dataset in
Machine Learning Research", NeurIPS 2021 D&B** — https://arxiv.org/abs/2112.01716 (2021)
Usage patterns 2015–2020: "increasing concentration on fewer and fewer datasets within task
communities, significant adoption of datasets from other tasks, and concentration across
the field on datasets that have been introduced by researchers situated within a small
number of elite institutions." Implications flagged for "scientific evaluation, AI ethics,
and equity/access."

### 6.3 Published criticism of specific benchmarks and of leaderboard practice

**Reuel, Hardy, Smith, Lamparth, Hardy, Kochenderfer, "BetterBench", NeurIPS 2024 D&B
(Spotlight)** — https://arxiv.org/abs/2411.12990 (2024)
"we develop an assessment framework considering **46 best practices across an AI
benchmark's lifecycle and evaluate 24 AI benchmarks** against it. We find that there exist
large quality differences and that commonly used benchmarks suffer from significant
issues. We further find that **most benchmarks do not report statistical significance of
their results nor allow for their results to be easily replicated.**" They ship a
minimum-quality-assurance checklist and a living repository of benchmark assessments.
*This is the single most citable sentence for justifying a statistics-forward benchmark
paper: not reporting statistical significance is the modal failure, measured.* (The arXiv
abstract page does not enumerate the 46 criteria; those need the PDF appendix.)

**Dehghani, Tay, Gritsenko, Zhao, Houlsby, Diaz, Metzler, Vinyals, "The Benchmark
Lottery"** — https://arxiv.org/abs/2107.07002 (2021)
"The benchmark lottery postulates that **many factors, other than fundamental algorithmic
superiority, may lead to a method being perceived as superior**… we show that the relative
performance of algorithms may be altered significantly simply by choosing different
benchmark tasks." Named lotteries: the **benchmark lottery**, **task-selection bias**
(rankings flip with the task subset), **community bias**, **statefulness** (benchmarks
accumulate tricks over time, so *when* you enter matters), **annotation bias**, and the
related **hardware lottery** (Hooker 2020). Their recommendations overlap heavily with
Section 5: treat a benchmark metric "as a single sample from the distribution describing
the model's performance" and compare distributions with appropriate tests; use multiple
fixed dataset splits because "dataset split contributes the most to model variance";
replace plain averaging with **macro-averaging by domain, geometric mean, or average rank /
robust average rank**; report model size, energy and latency alongside accuracy; publish
data statements; use benchmarking checklists in review; build living benchmarks with
test-query limits to prevent "creeping overfitting"; and **fix the hyperparameter-tuning
budget so improvements are not purely compute**.
*That last item is the direct remedy for the Tripp & Hernández-Lobato critique, and a
benchmark that states a tuning budget pre-empts it.*

**Singh, Nan, Wang, D'Souza, Kapoor, Üstün, Koyejo, Deng, Longpre, Smith, Ermis, Fadaee,
Hooker, "The Leaderboard Illusion"** — https://arxiv.org/abs/2504.20879 (2025)
On Chatbot Arena: "undisclosed private testing practices benefit a handful of providers who
are able to test multiple variants before public release and retract scores if desired…
the ability of these providers to choose the best score leads to **biased Arena scores due
to selective disclosure**." Numbers: **27 private Llama variants** tested before the
Llama-4 release; proprietary closed models sampled at higher battle rates and deprecated
less often; **Google ~19.2% and OpenAI ~20.4% of all arena data** versus **29.7% for 83
open-weight models combined**; extra arena data yields "relative performance gains of up to
112% on the arena distribution." Conclusion: "these dynamics result in **overfitting to
Arena-specific dynamics rather than general model quality**."
*Relevant to a molecular benchmark mainly as the strongest recent statement that submission
and retraction policy is part of a leaderboard's validity — a hosted leaderboard needs a
stated policy on repeated submissions.*

### 6.4 What this means concretely for a molecular/peptide evaluation paper

Synthesising §6.1–6.3 against what Section 1 shows the field actually does, the
distrust-triggers in descending order of how likely they are to apply:

1. **No statistical significance reported** — BetterBench 2024 measured this as the modal
   failure across 24 benchmarks, and NeurIPS checklist item 7 asks for it directly. GuacaMol
   reports no variance at all; MOSES uses 3 seeds; PMO and Tartarus 5.
2. **Insufficiently tuned baselines** — Liao et al. 2021 catalogue it as a field-wide
   internal-validity failure, and §2.3–2.4 demonstrate it specifically for PMO, where
   retuning one baseline moved it as far as the entire top-10 spread.
3. **Unstated or unfixed tuning budget** — Dehghani et al. 2021's explicit recommendation;
   PMO is the only paper in Section 1 that discloses its tuning sweep.
4. **Construct validity of the oracle not established** — Raji et al. 2021. Gao & Coley's
   validate-the-oracle-first move (Fig. S14) is the model answer.
5. **Task-selection sensitivity not tested** — Dehghani et al. 2021; Thomas et al. 2022
   Fig. 3 shows ranks swinging from 4th to 18th across defensible metric variants.
6. **Croissant/licensing/hosting/code gaps** — cheapest to fix, still caught ~11.9% of 2025
   D&B submissions on license alone.

## GAPS

- **PepFoundry's figure inventory is unread.** pubs.acs.org (DOI 10.1021/acs.jcim.5c02629)
  and chemrxiv.org both served bot challenges; §3.3 rests on the abstract alone. The
  representation × model-type grid that supports "regardless of model type" is the one
  artifact I most wanted and do not have.
- **Published-version figure numbering unverified for GuacaMol (JCIM 2019) and MOSES
  (Front. Pharmacol. 2020).** Both inventories come from the arXiv preprints. The claims
  that matter — GuacaMol reports no variance, MOSES reports mean ± std over 3
  initialisations — are solid in the preprints, but panel numbering may differ in print.
- **No molecular or peptide generative benchmark in this corpus uses any Section 5
  convention.** No critical-difference diagram, no performance profile, no bootstrap CI, no
  significance test across methods. The closest approximations are NovoMolGen's
  colour-coded p < 0.05 table cells and HELM-BERT's † significance markers. HPO-B
  (arXiv 2106.06257) is the only benchmark found anywhere that combines regret curves with
  CD diagrams — there is no chemistry precedent to cite.
- **Nobody reports enough runs.** 3 (MOSES, NovoMolGen), 5 (PMO, Tartarus, TDC), 10
  (Skinnider, Renz), 32 (Ehrlich). Agarwal et al. 2021 found "statistically defensible
  improvements with median scores is only achieved for 25 runs"; Colas et al. 2018 say the
  bootstrap CI test "should not be used with less than N = 20 samples". Only Ehrlich and
  Skinnider clear either bar, and neither is a leaderboard.
- **The seed-variance question for molecular optimization has not been asked directly.**
  Nothing in this literature is the equivalent of Henderson et al.'s Figure 5 — split one
  algorithm's runs into two groups and show they look "significantly different" from each
  other. Thomas et al. 2022 Fig. 1 is the nearest, and it is a side observation in a
  workshop paper.
- **Ehrlich's quantile convention has no second adopter.** The 10/50/90-over-32-trials band
  is the most defensible variance presentation in the corpus, and as of this writing it
  appears in exactly one paper, from one group, in a workshop.
