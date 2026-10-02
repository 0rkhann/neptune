# Docking as a deterministic oracle for peptide generative design — evidence file

Scope: Vina-family docking (AutoDock Vina, QuickVina2, smina, GNINA, Vina-GPU / QuickVina2-GPU,
AutoDock-GPU) applied to peptide and macrocyclic ligands, 60–390 heavy atoms, 15–25 rotatable
bonds, non-canonical residues; thousands of calls; determinism required.
Facts only. Every claim carries a URL and a year.

---

## QUESTION 1 — IS VINA-FAMILY DOCKING VALID AT THIS LIGAND FLEXIBILITY?

### 1.1 Short answer

**No validated regime exists at 15–25 rotatable bonds.** There is no single published
"collapse threshold" stated as a hard number by the Vina authors, but three independent
lines of evidence converge on **~10 rotatable bonds as the edge of the validated envelope**,
and the Vina developers themselves place the *search feasibility* ceiling at 20–50 rotatable
bonds — which is a statement about whether the global optimum is findable at all, not about
accuracy. For peptides specifically, the honest statement is:
**flexible-ligand Vina has been assessed on peptides, it fails on most of them, and the
assessments cover peptides far smaller than the ligands in this benchmark.**

### 1.2 The original Vina validation set — what it actually covers

Trott & Olson, *J Comput Chem* 31(2):455–461, 2010.
https://pmc.ncbi.nlm.nih.gov/articles/PMC3041641/

- Validation set: 190 protein–ligand complexes, the AutoDock 4 training set. Rigid receptor,
  flexible ligand, **active rotatable bonds ranging 0 to 32**.
- Overall RMSD < 2 Å success: **78% (Vina) vs 49% (AutoDock 4)** on the 190; on the 116
  complexes not in PDBbind, **80% vs 53%**. (2010)
- Accuracy vs. flexibility is shown only as **figures (Fig. 3 and Fig. 5), binned by number of
  active rotatable bonds — no tabulated per-bin success rates are given in the text.** The paper
  states the degradation qualitatively: Vina's lower standard error on binding free energy is
  attributed to "the ligands with many active rotatable bonds, for which AutoDock has difficulty
  finding the correct bound conformation, as can be seen in Figure 3".
- Note the composition: this is a *drug-like* set. The high-Nrot tail (20–32) is a handful of
  complexes, not a validated regime.

**Consequence for citation hygiene:** anyone quoting "Vina is 78% accurate" is quoting a
drug-like-ligand number. That figure carries no information about 15–25 rotatable bonds.

### 1.3 The ~10 rotatable-bond threshold — published statements

| Source | Year | Statement |
|---|---|---|
| Röhrig et al., *Attracting Cavities 2.0*, JCIM 63(12):3925 — https://pubs.acs.org/jcisd8/article/63/12/3925/982687/Attracting-Cavities-2-0-Improving-the-Flexibility | 2023 | "ligands with more than 10 rotatable bonds display significantly lower success rates" (evaluated across docking programs incl. Vina-family) |
| Yang et al., *Dockformer*, arXiv:2411.06740 / IEEE — https://arxiv.org/html/2411.06740v4 | 2024/2026 | "the success rates of most traditional docking programs significantly decrease when the number of rotatable bonds increases, and satisfactory predictions are obtained only for compounds with **fewer than 10 rotatable ligand bonds**" |
| Wang et al., *Comprehensive evaluation of ten docking programs*, PCCP 18(18):12964 — https://pubs.rsc.org/cp/article/18/18/12964/517163 | 2016 | Sampling power degrades with rotatable-bond count; >half (40/72) of the hardest ligands carry >10 rotatable bonds |
| Sarigun et al., *PocketVina*, arXiv:2506.20043 — https://arxiv.org/html/2506.20043 | 2025 | Stratifies PDBbind2020-timesplit and PoseBusters by flexibility: **rigid ≤5 RB, medium 6–10 RB, flexible >10 RB**. Even the best method (QuickVina2-GPU 2.1 over multiple pockets) reaches ~80% physically-valid success on rigid, >60% on medium, and shows "a significant decrease in valid success rate for flexible ligands". The three-bin split is itself the field's operational admission that >10 RB is a separate, degraded regime. |

**Note the ceiling of the stratification: "flexible" means >10 RB. The peptides in question
(15–25 RB) sit beyond the top bin of every published stratified benchmark.**

### 1.4 The developers' own ceiling

AutoDock-Vina GitHub issue #420, "What is the practical limit on ligand size?" (2025)
https://github.com/ccsb-scripps/AutoDock-Vina/issues/420

Diogo Santos-Martins (Forli Lab / AutoDock maintainer), 2025-04-11:
> "For the sake of search efficiency, the limit for what's possible to find the best solution
> should be somewhere between 20 and 50 rotatable bonds."

This is a **search-feasibility** statement, not an accuracy statement: above it, the global
optimum is not reliably findable. 15–25 RB sits at or just under this stated boundary.
The thread also documents hard failure (OOM kill) on a large flexible ligand (>500 torsional
DOF) on a 64 GB / 32-core machine.

### 1.5 Published assessments on peptide ligands — with numbers

**Rentzsch & Renard, *Briefings in Bioinformatics* 16(6):1045–1056, 2015.**
"Docking small peptides remains a great challenge: an assessment using AutoDock Vina"
https://academic.oup.com/bib/article/16/6/1045/225862 · PDF: https://academic.oup.com/bib/article-pdf/16/6/1045/606428/bbv008.pdf

- Dataset: **47 complexes, peptides of 2–5 residues** (i.e. far smaller than this benchmark's ligands).
- Protocol: Vina 1.1.2, semi-blind (30 Å box), exhaustiveness 8→1024, 5 replicate runs each,
  3,760 dockings total.
- **Flexible-ligand docking (FD_512) top-pose success: ~1/4 of cases. Rigid redocking
  (crystal conformation preserved, RD_128): ~3/4 of cases.** Verbatim: "FD is (reproducibly)
  successful in only 1/4 of cases, whereas RD is successful for 3/4 cases."
  (Their success band is top-pose backbone RMSD < 1.5 Å; "bad" is ≥ 3.5 Å.)
- **"There is no overall correlation between pose rank and similarity to the native pose."**
  i.e. the scoring function does not rank the correct pose first.
- Reproducibility finding (relevant to Q2): Vina's own output is seed-dependent and pose ranking
  is unstable across replicate runs; at exhaustiveness 8 only 24/47 complexes passed their
  stability test in flexible mode, rising to 42/47 at exhaustiveness 1024 **without converging**.
- Secondary summary of this work (Agrawal et al. 2019, below): "The overall performance of
  AutoDock Vina was case dependent and **poor for peptides with more than 4 residues**."

**Hauser & Windshügel, LEADS-PEP, *JCIM* 56(1):188–200, 2016.**
https://pubmed.ncbi.nlm.nih.gov/26651532/

- 53 protein–peptide complexes, **peptide length 3–12 residues**.
- Conclusion: "most tested programs were capable to generate native-like binding modes of small
  peptides, **only Surflex-Dock and AutoDock Vina performed reasonably well for peptides
  consisting of more than five residues**." (Vina is the *best* of the small-molecule programs
  here — which is the ceiling, not a pass.) Per-program RMSD<2 Å tables are behind the ACS paywall
  and were not retrieved; the qualitative conclusion above is from the abstract.

**Agrawal et al., *BMC Bioinformatics* 20:33, 2019.** (PPDbench)
https://link.springer.com/article/10.1186/s12859-018-2449-y

- 133 protein–peptide complexes, **peptides 9–15 residues** — the length range closest to this
  benchmark. Six docking programs benchmarked (ZDOCK, FRODOCK, Hex, PatchDock, ATTRACT, pepATTRACT).
- Blind docking, top pose: best method FRODOCK, **average L-RMSD 12.46 Å**. Best-of-20 poses:
  FRODOCK 3.72 Å. Re-docking top pose: best ZDOCK, **8.60 Å**; best-of-20: 2.88 Å.
- Success rate within 2.0 Å L-RMSD, **top pose: ZDOCK 32.33%, FRODOCK 39.09%**.
- **Vina itself was only run on the 40-complex ≤5-residue Rentzsch subset, not on the 9–15-residue
  main set** — where it reached best-pose (not top-pose) average L-RMSD 2.09 Å on re-docking.
- Reported failure cause directly relevant here: several programs **"fail on certain PDB IDs …
  due to the presence of non-natural residues (example ACE, ACY) in the peptide"**.

**Holcomb, Chang, Goodsell & Forli, *QRB Discovery* 3:e18, 2022.** (macrocycles)
https://www.cambridge.org/core/journals/qrb-discovery/article/performance-evaluation-of-flexible-macrocycle-docking-in-autodock/D8417BC284AEE198EC6AF25C7E677249

This is the only systematic macrocycle assessment by the AutoDock authors. AutoDock-GPU v1.5.3,
Meeko ring-breaking + anisotropic closure potential, 90 curated PDB macrocycle complexes,
7- to 33-membered rings, 20 GA runs each.

- **Rigid redocking (crystal ring conformation given): 76% (68/90) within 2 Å.**
- **Flexible redocking (ring conformation unknown — the prospective case): 53% (48/90).**
- Restricted to rings <15 atoms: **69% rigid vs 59% flexible** (N=54).
- 23/90 complexes did not converge (top cluster ≤3 of 20 poses). Rerunning those with
  autostop/heuristics disabled at **100M evaluations** improved only 3 of 23:
  **success on that hard subset went 17% → 29%.**
- Cyclic peptides were **deliberately excluded from the main set** and handled as 6 hand-picked
  case studies (HIV-1 protease inhibitors 1b6j/1b6p/4cpw; arylomycin C 3s04, 24 torsions;
  darobactin 7nrf, 32 torsions; vancomycin 1rrv, 39 torsions), each run at 100M evaluations with
  heuristics off. 5 of 6 top poses <2 Å; 4cpw failed at top pose (correct pose was cluster 3).
- The authors' own caveat: "Peptide ligands tend to be larger than typical organic macrocyclic
  structures, resulting in a very large number of active torsions. **Specialized software with
  ad-hoc protocols such as AutoDock CrankPep … may be better suited to this task**" and
  "with respect to cyclic peptides … **a systematic treatment would be challenging because of the
  large number of active torsions**."

**Earlier macrocycle datapoint:** AutoDock-GPU reached <2 Å on all 19 D3R Grand Challenge 4
macrocycles — but **"using visual inspection to select the best pose"** (Santos-Martins et al.,
*J Comput Aided Mol Des* 33:1071, 2019; recounted in the QRB paper above). Human pose selection
is not available to an oracle.

### 1.6 smina, GNINA, QuickVina2 on peptides specifically

**Nobody has properly assessed them on peptides.** Searches across the peptide-docking
benchmarking literature (LEADS-PEP 2016, PPDbench 2019, Rentzsch 2015, the ADCP/AlphaFold2
peptide comparison of Shanker & Sanner 2023 — https://pmc.ncbi.nlm.nih.gov/articles/PMC13459648/)
turn up **no dedicated peptide or macrocycle pose-prediction benchmark for smina, GNINA,
QuickVina2, or any Vina-GPU variant.** These programs appear in peptide papers only as
*components* (e.g. QuickVina2-GPU 2.1 inside PocketVina) evaluated on drug-like benchmarks
(PDBbind, PoseBusters, Astex, DockGen), all of which are small-molecule sets.

GNINA's CNN rescoring is trained on drug-like protein–ligand complexes (PDBbind/CrossDocked);
no published training or validation set covers peptides or non-canonical residues. Treat any
GNINA score on a 25-rotatable-bond peptide with a non-canonical residue as **extrapolation
outside the training distribution**, not as a measurement.

### 1.7 Verdict on Question 1

For the stated ligand class (60–390 heavy atoms, 15–25 RB, non-canonical residues):

1. **Not validated.** Every stratified benchmark tops out at ">10 rotatable bonds"; the ligands
   here start above that bin. The Vina maintainer's own search-feasibility ceiling is 20–50 RB.
2. **Where peptides *have* been measured, the numbers are bad.** Flexible-ligand Vina on
   2–5-residue peptides: ~25% top-pose success (2015). Flexible macrocycle redocking with the
   authors' own best protocol: 53% (2022), and that set deliberately excluded peptides.
3. **Non-canonical residues are a documented hard failure mode**, not merely an accuracy
   degradation — docking programs error out on ACE/ACY-containing peptides (2019).
4. **The successes that do exist on peptide-like macrocycles required protocol escalation that
   defeats throughput** (heuristics off, 100M evaluations) **and in one published case human
   visual pose selection.**
5. **Ranking is the specific failure.** Rentzsch & Renard (2015) found no correlation between
   pose rank and native-likeness; the QRB macrocycle study found the correct pose sampled but
   mis-ranked (4cpw). An oracle consumes the top-ranked score, so sampling-level partial success
   does not transfer.

### 1.8 Addendum — AutoDock maintainers on torsion count (added while researching Q2)

AutoDock-GPU issue #297, maintainer `atillack`, 2025:
https://github.com/ccsb-scripps/AutoDock-GPU/issues/297
> "This likely means your ligand has too many torsions for AD-GPU to be able to traverse search
> space efficiently. If possible, we recommend to keep the number of torsions **below 20** by
> selectively 'freezing' some or by docking smaller sub-units. While Vina typically does a bit
> better with more torsions as it uses a different search algorithm and force field;
> **even for Vina, too many torsions will inhibit convergence**."

AutoDock-GPU README (https://github.com/ccsb-scripps/AutoDock-GPU/blob/develop/README.md):
> "For molecules with many rotatable bonds (e.g. **about 15 or more**) it may be advisable to
> increase `--heurmax`."

i.e. 15 RB is where the shipped evaluation heuristic is already documented as insufficient, and
20 torsions is the maintainers' recommended practical cap. The target ligand class is 15–25 RB.

---

## QUESTION 2 — IS QUICKVINA2-GPU DETERMINISTIC?

### 2.1 Short answer

**No tool in the Vina GPU family documents or guarantees bit-identical output, and the
AutoDock-GPU maintainer has stated in writing that a fixed seed does not guarantee it.**
CPU AutoDock Vina is the only member with a stated (and historically broken) fixed-seed
reproducibility property.

| Tool | `--seed` exposed? | Documented determinism guarantee? | Evidence of run-to-run variance |
|---|---|---|---|
| AutoDock Vina (CPU) | Yes, documented CLI option | Not a formal guarantee. Docs say the algorithm "is non-deterministic"; reproducibility with a fixed seed is asserted in the manual and demonstrated once in the 2010 paper | **Yes — was a real bug in v1.2.3, fixed in v1.2.4** |
| smina | Yes (`--seed`) | None found | not measured in literature found |
| GNINA | Yes (`--seed`) | None found. GPU CNN scoring (torch/cuDNN) carries no determinism statement | GNINA 1.0 paper reports single- vs double-precision pose differences |
| AutoDock-GPU | Yes (`--seed`, up to 3 ints; default = system time + pid) | **Explicitly disclaimed by a maintainer** (see 2.4) | Yes (issue #262, 2024) |
| Vina-GPU / Vina-GPU 2.0 / 2.1, QuickVina2-GPU (2.0/2.1) | **Not listed in any README argument table.** The binary prints "Using random seed: N" at runtime, so the inherited Vina option exists but is undocumented for the GPU builds | **None anywhere** — no determinism or reproducibility statement in any README or in the Vina-GPU 2.1 paper | Yes — see 2.3 |

### 2.2 What CPU AutoDock Vina actually guarantees with a fixed `--seed`

Official FAQ (https://autodock-vina.readthedocs.io/en/latest/faq.html), undated, current:
> "The docking algorithm is non-deterministic."

Trott & Olson 2010 (https://pmc.ncbi.nlm.nih.gov/articles/PMC3041641/) is the one primary
demonstration: after a single-threaded run, "Vina was rerun with the same random seed and 8-way
multithreading … **This latter Vina run produced identical results**, but executed faster."
So: same seed + different thread count → identical output, as of 2010. The Vina manual's standing
phrasing (quoted widely, e.g.
https://bioinformatics.stackexchange.com/questions/15812/) is that exact reproducibility is assured
by the same seed **"only if all other inputs and parameters are the same as well. Even minor
changes to the input can have an effect similar to a new random seed."**

**This guarantee has been empirically broken.** AutoDock-Vina issue #90 (2022):
https://github.com/ccsb-scripps/AutoDock-Vina/issues/90
A user docked 1000 receptor structures twice **without changing the seed** and got score
correlation ≈0.8 with differences up to ~5 kcal/mol, on cases where the *pose was identical*.
Maintainer `diogomart`, 2022-03-04: **"This is a bug, and it's fixed by … PR #81. I'll merge the PR
soon and release v1.2.4."** He later added (2023-06-08) that the same issue "Might have been
affected by #200" — a separate macrocycle-specific bug
(https://github.com/ccsb-scripps/AutoDock-Vina/issues/200, 2023) where ligands carrying
**macrocycle closure atoms (CG0, G0, CG1 …) produced absurd free energies, e.g. −50 kcal/mol**,
from a mismatch between intramolecular and unbound energies. The verification run for that fix
used exhaustiveness 64 and **the same random seed (42) in all runs** across ~34k ligands.

Two things follow, both directly relevant here: (a) fixed-seed determinism in CPU Vina is a
property of a *particular build*, not a contract — it regressed and was repaired; (b) the macrocycle
glue-atom machinery that this benchmark's ligands would use is precisely where the scoring bug lived.

### 2.3 QuickVina2-GPU / Vina-GPU 2.x

**No documentation statement exists.** The Vina-GPU 2.1 README
(https://github.com/DeltaGroupNJUPT/Vina-GPU-2.1, and
https://raw.githubusercontent.com/DeltaGroupNJUPT/Vina-GPU-2.1/main/README.md) lists the complete
argument table — `--config`, `--receptor`, `--ligand_directory`, `--output_directory`, `--lbfgs` /
`--rilc_bfgs`, `--thread`, `--search_depth`, `--center_x/y/z`, `--size_x/y/z`,
`--opencl_binary_path` — and **`--seed` does not appear in it.** The same is true of the
QuickVina2-GPU 2.1 subdirectory README. The word "deterministic", "reproducib*" or "seed" appears
nowhere in either README. The Vina-GPU 2.1 paper (Tang et al., IEEE/ACM TCBB, 2024,
https://pubmed.ncbi.nlm.nih.gov/39320991/; preprint
https://www.biorxiv.org/content/10.1101/2023.11.04.565429v1.full-text) reports speed and
enrichment, not reproducibility.

The seed nevertheless exists at runtime: binaries print `Using random seed: <N>` (e.g.
https://bioinformatics.stackexchange.com/questions/22394/, 2024, QuickVina2-GPU 2.1 showing
`Using random seed: 433324452`).

**GitHub issue on exactly this question — Vina-GPU #36 (2023), still open:**
https://github.com/DeltaGroupNJUPT/Vina-GPU/issues/36
Title: "run Vina-GPU.exe multiple times with the same example config.txt but get different Vina
dock scores". Developer `Glinttsd`, 2023-03-21, the only reply:
> "the reason lies in the random seed, you should use the same random seed to produce the same results"

That is an assertion by the developer, **not a measurement, and the issue was left open.** No
reproduction, no test, and the argument table still does not document the flag.

**The one published numeric measurement of run-to-run spread** comes from the derivative project
Vina-CUDA (Li et al., *JCIM* 65(10):4751–4759, 2025, DOI 10.1021/acs.jcim.4c01933), whose README
lists under "Limitation" (https://github.com/HPCLab-933/Vina-CUDA):
> "Due to the variability of random seed values, **the docking score error ranges from 0.1 to 0.5
> in each run** (this phenomenon also occurs in **Vina-GPU 2.1**)."

0.1–0.5 kcal/mol of run-to-run spread, stated by the authors as also applying to Vina-GPU 2.1.
Note the attribution is to *seed* variability (the GPU builds seed from the clock by default),
not to a fixed-seed reduction-order effect — the latter has not been separately measured
in any source found.

**Downstream users treat it as stochastic.** Guo et al., arXiv:2407.12186 (2024), reproducing
RGFN's QuickVina2-GPU-2.1 oracle, report in their appendix that they "validated their docking
protocol by re-docking the reference ligand multiple times with varying random seeds and observed
**notable differences in docking scores and poses**", and note that "When executing
QuickVina2-GPU-2.1, if a seed is not specified, a random seed is used" — i.e. the default behaviour
of the exact tool used as a generative-design oracle is non-reproducible unless the caller
intervenes, and the original RGFN paper did not state whether it had.
(Summary: https://www.emergentmind.com/open-problems/rgfn-docking-seed-usage)

### 2.4 AutoDock-GPU — an explicit maintainer disclaimer

AutoDock-GPU issue #126, "Random Seed", maintainer `atillack`
(https://github.com/ccsb-scripps/AutoDock-GPU/issues/126):
> "By default when AD-GPU starts, the first two seeds are set to the system time and the process id …
> **Although the seeds can be used to have identical results as the RNGs we are using are
> deterministic … this is not necessarily the case all the time. When running things in parallel
> the sequence of random numbers — even with the same seeds — can diverge when worker threads using
> them are run out-of-order, which can easily happen on a GPU when running more work units than
> compute units or when using different GPUs. This is not a bug at this point but more a 'potential
> future feature'** :-)"

This is the clearest statement available anywhere in the Vina/AutoDock ecosystem: **a fixed seed on
GPU does not guarantee identical results, by design, and the maintainers do not consider it a bug.**
The stated mechanism is RNG-stream divergence from out-of-order work-unit scheduling, which is
*in addition to* the floating-point summation-order effect.

Corroborating report: AutoDock-GPU issue #262, "AutoDock GPU Batch Docking Result Unreproducible"
(2024, https://github.com/ccsb-scripps/AutoDock-GPU/issues/262) — a user could not reproduce a
−14 kcal/mol result from a first batch in any of 100 later runs of the same ligand, including
"by using the random seed in the first batch dlg file". Maintainer reply: "this is unexpected,
I'll try to reproduce"; the user later attributed part of it to edited map files. Not a clean
measurement, but the thread shows seed reuse did not restore the result.

### 2.5 The underlying mechanism, with references

GPU non-determinism from floating-point non-associativity is well characterised outside docking
and is the correct prior:
- NVIDIA CCCL issue #5550, "[EPIC] Floating-Point Deterministic Algorithms"
  (https://github.com/NVIDIA/cccl/issues/5550) defines the three tiers the field uses:
  **GPU-to-GPU** (bitwise identical on any GPU), **run-to-run** (identical on the same GPU),
  **not guaranteed** — and states "The non-associativity of floating-point arithmetic and the lack
  of order of execution guarantees from parallel execution generally mean that two identical
  invocations of an algorithm on identical floating-point inputs may result in different results."
  Deterministic reductions are an opt-in feature NVIDIA is still adding.
- Oak Ridge et al., arXiv:2408.05148 (2024),
  "Impacts of floating-point non-associativity on reproducibility for HPC and deep learning
  applications": atomicAdd-based reductions are non-deterministic; pairwise/two-pass schemes are
  deterministic by construction but must be chosen deliberately. "The variability from
  non-deterministic reductions can approach the tolerance thresholds used in high-accuracy
  molecular simulation correctness tests."

No Vina-family GPU project documents which reduction scheme it uses or claims either tier.

**The one quantitative docking-specific study of reduction order** is Solis-Vasquez et al.,
"Accelerating Drug Discovery in AutoDock-GPU with Tensor Cores", arXiv:2410.10447 (2024)
(https://arxiv.org/html/2410.10447v1). They replace AutoDock-GPU's reduction with a tensor-core
implementation, run **1000 dockings per complex with the pseudo-random generator initialized with
the same seed for both their code and the baseline**, and report box-and-whisker distributions of
best energy. Findings: for 1ac8 and 3tmn "the best energy values show no significant variance
between runs for both implementations"; for 1stp, 7cpa and 3ce3 the distributions are merely
"similar"; and **"for all tested complexes, the relative difference between the average best scores
for each method is below 0.18%."** They validate by *distributional similarity and average relative
error*, never by bit-identity — which is the operative tell: the people who changed the reduction
order inside AutoDock-GPU did not have, and did not claim, bitwise reproducibility.

### 2.6 GNINA / smina specifics

- Both expose `--seed` (`gnina --seed arg  explicit random seed`,
  https://github.com/gnina/gnina/blob/master/README.md). Neither README states a determinism guarantee.
- GNINA 1.0 (McNutt et al., *J Cheminform* 13:43, 2021,
  https://pmc.ncbi.nlm.nih.gov/articles/PMC8191141/) documents a precision change versus smina:
  "**unlike Smina, Gnina does computation with single (32 bit) precision rather than double (64 bit)
  precision due to the need to shift calculations to the GPU**". Comparing GNINA (CNN off,
  autobox_extend off) against smina on redocking: "**A majority of the output poses are exactly the
  same, with slight differences seen for some output poses.**" A majority — not all.
- GNINA's default `--cnn_scoring rescore` runs a 3-model CNN ensemble on the GPU (Caffe in 1.0,
  torch from 1.3). The CNN pass therefore sits on a deep-learning stack whose own determinism is
  framework- and kernel-dependent; no statement is made about it in GNINA's documentation.

### 2.7 What DOCKSTRING actually froze

García-Ortegón, Simm, Tripp, Hernández-Lobato, Bender & Bacallado, *JCIM* 62(15):3486–3502, 2022.
https://pubs.acs.org/doi/full/10.1021/acs.jcim.1c01334 · open copy
https://pmc.ncbi.nlm.nih.gov/articles/PMC9364321/ · code https://github.com/dockstring/dockstring

The claim: "**Subsequently, we fixed the random seed to obtain a fully deterministic pipeline.**"
Their complaint about predecessors (TDC, Cieplinski et al.) was precisely this: "Neither TDC nor
Cieplinski et al. control sources of randomness during the docking procedure (e.g., random seeds
input into the docking program or **the conformer generation routines**), leading to the potential
for considerable variance between runs on the same molecule."

**The frozen steps, in order:**

1. **Targets frozen as files, not as a procedure.** 58 pre-prepared receptor PDBQT files ship with
   the package. 57 come from DUD-E; PDBs standardized with Open Babel, polar hydrogens added,
   PDBQT conversion by AutoDock Tools. The 58th (DRD2, PDB 6CM4) was prepared once by hand:
   ligand/waters removed in PyMOL, `obminimize` with GAFF, protonation at pH 7.4 with PROPKA,
   polar H + PDBQT via AutoDock Tools.
2. **Search boxes frozen** — fixed per-target box shipped alongside each receptor file.
3. **Ligand input restricted.** Sanity check rejects radicals and multi-fragment molecules.
4. **Protonation frozen**: (de)protonation at pH 7.4 with Open Babel.
5. **Conformer generation frozen**: a **single** 3D conformer from RDKit ETKDG, refined with MMFF94
   — one conformer, not an ensemble, so no selection step to vary.
6. **Stereochemistry frozen**: determined stereocenters preserved; "any undetermined stereocenters
   are assigned randomly (**but consistently across different runs to ensure the reproducibility of
   docking scores**)".
7. **Charges frozen**: Gasteiger charges, PDBQT written by Open Babel.
8. **Docking parameters frozen**: AutoDock Vina with defaults — exhaustiveness 8, num_modes 9,
   energy_range 3; only the lowest (best) of the up-to-9 scores is used.
9. **Random seed fixed** — after first measuring the sensitivity: "we investigated this dependence
   and found **no target–ligand combination for which the docking scores deviated by more than
   0.1 kcal/mol**."
10. **Dependency versions pinned, and verified by a test.** The repo pins
    `python<=3.10, openbabel=3.1.1, rdkit<=2022.03` and ships
    `tests/test_dataset_matching.py`, which re-docks N random molecules from the published dataset
    and checks the scores match. The README's own acceptance language is the honest part:
    "**If 99%+ of scores match then it is probably ok to use dockstring in the benchmarks, but there
    will of course be some error and this should be noted.**"

**Three limits on transferring this to the present problem:**
- DOCKSTRING is **CPU AutoDock Vina only**. Nothing in it addresses GPU reductions.
- Its 0.1 kcal/mol seed-insensitivity finding was measured on **drug-like ZINC-scale ligands**
  (260k molecules × 58 targets). Seed sensitivity grows with the size of the search space;
  Rentzsch & Renard (2015, §1.5) found peptide flexible-ligand docking *does not* converge with
  respect to sampling even at exhaustiveness 1024. The 0.1 kcal/mol figure must not be carried over.
- Determinism was achieved by **removing degrees of freedom from the pipeline** (one conformer,
  one protonation state, fixed box, fixed receptor file, pinned library versions) — not by making
  the docking engine deterministic. The engine was simply given identical bytes each time.

---

## SECTION 3 — CIEPLIŃSKI ET AL.: SCORE–SIZE CORRELATION AND DEGENERATE OPTIMA

Ciepliński, Danel, Podlewska & Jastrzȩbski,
"Generative Models Should at Least Be Able to Design Molecules That Dock Well: A New Benchmark",
*JCIM* 63(11):3238–3247, 2023 (preprint arXiv:2006.16955, v1 2020, v5 2023).
https://pubs.acs.org/doi/full/10.1021/acs.jcim.2c01355 · open copy
https://pmc.ncbi.nlm.nih.gov/articles/PMC10268949 · preprint https://arxiv.org/html/2006.16955
Code: https://github.com/cieplinski-tobiasz/smina-docking-benchmark

Setup: docking with **SMINA at default settings**, 8 targets (5HT1B, 5HT2B, ACM2, CYP2D6, ADRB1,
MOR, A2A, D2), binding-site coordinates chosen manually, receptors prepared in Schrödinger.
Three optimization tasks: the full docking score, the **repulsion** term alone (minimize), and the
**non_dir_h_bond** term alone (maximize). Scores averaged over the top 5 poses. Generated sets are
filtered by Lipinski and MW > 100. Score decomposition as printed in the preprint:
`−0.035579·gauss(o=0,w=0.5) − 0.005156·gauss(o=3,w=2) + 0.840245·repulsion
 − 0.035069·hydrophobic − 0.587439·non_dir_h_bond`
(these are the AutoDock Vina default weights; the JCIM text describes the components as Vinardo
terms — a labelling discrepancy between the two versions, not load-bearing here).

### 3.1 The score–molecular-weight correlation: **no numeric coefficient was published**

This is the honest answer to "the exact correlation they measured." The paper reports it
**qualitatively and graphically only**:

> "We noticed a **moderately strong correlation** between docking scores and the number of
> rotatable bonds or molecular weight."
> "with the increasing number of rotatable bonds or molecular weight, the docking scores improve."
> "the distribution of generated compounds is shifted towards better docking scores and smaller
> molecular weights … molecules achieve better docking scores at the same molecular weight after
> the optimization. **The correlations are weaker for CYP2D6, which may be caused by a bigger
> binding site of this enzyme.**"

The evidence is two scatter-plot figures, neither of which carries a coefficient:
- **Figure 6** — "Correlation between docking score and molecular weight. The training set is
  marked with red dots, and the compounds generated by REINVENT by enhancing different optimization
  targets are colored in blue (Docking Score Function), orange (Hydrogen Bonding), and green
  (Repulsion)."
- **Figure 7** — "Correlation between docking score and the number of rotatable bonds." (same colour key)

**No Pearson or Spearman r is given in the text, in either figure caption, or in the abstract.**
Anyone citing a number for this correlation is citing something other than this paper.

### 3.2 Degenerate optima by scoring-function component — verbatim

> "there are clear patterns visible in the top molecules of CVAE and GVAE (Figures 4 and 5).
> For example, **CVAE generates macrocycles in the task of docking score optimization**, while
> **GVAE generates long chains with no cycles** when optimizing the same objective. These models
> also create **oxygen or nitrogen chains when optimizing Hydrogen Bonding**, and **very small
> molecules (often less than 3 heavy atoms) for the Repulsion task.**"

Mapping to objective:
| Objective optimized | Degenerate structure emitted |
|---|---|
| Full docking score | **macrocycles** (CVAE); **long acyclic chains** (GVAE) |
| `non_dir_h_bond` (hydrogen bonding) alone | oxygen or nitrogen chains |
| `repulsion` term alone | molecules **often under 3 heavy atoms** |

Supporting quantitative context from Table 2/3: on the Repulsion task CVAE reaches 1.148 / 1.001 /
1.132 (5HT1B / 5HT2B / ACM2) against a ZINC-1% reference of 0.613 / 0.625 / 0.612 — i.e. the models
that emit sub-3-heavy-atom molecules still do not beat the ZINC top 1%, and on the full docking
score CVAE and GVAE reach only ≈ −4.2 to −5.4 kcal/mol against ZINC-1% of −8.8 to −11.5. Dashes in
the tables mark cases where "the model failed to generate 250 molecules that satisfy drug-like
filters".

### 3.3 Why this matters for the stated plan

The objective being proposed is a docking score over **peptides and macrocycles at 60–390 heavy
atoms**. Ciepliński et al.'s measured pathology is that an unconstrained docking score rewards
**exactly** mass, rotatable-bond count, and macrocyclicity — and that the degenerate optimum of the
full docking score for one of their two VAEs *was a macrocycle*. A benchmark whose legitimate
chemical space *is* large flexible macrocycles has no size filter available to separate the reward
signal from the artifact: the Lipinski + MW>100 filter that Ciepliński et al. used to keep the
benchmark honest cannot be applied here, because every valid ligand violates it.

---

## VERDICT

**No. Vina-family docking cannot serve as a deterministic oracle for peptide/macrocycle ligands at
15–25 rotatable bonds over thousands of calls. It fails three independent requirements, any one of
which is disqualifying.**

**1. Validity — fails.** The ligand class sits above the top bin of every published stratified
benchmark (">10 rotatable bonds"), above the AutoDock-GPU maintainers' recommended cap
(<20 torsions, with the eval heuristic already flagged as insufficient at ~15 RB), and inside the
Vina maintainer's stated search-feasibility boundary of 20–50 RB where "the best solution" is no
longer reliably findable. Where peptides have actually been measured: **~25% top-pose success for
flexible-ligand Vina on 2–5-residue peptides (2015)**, and **53% for flexible macrocycle redocking
with the AutoDock authors' own best protocol on a set that deliberately excluded peptides (2022)**.
Non-canonical residues are a documented hard-failure mode (2019), not a graceful degradation.
smina, GNINA, QuickVina2 and every GPU variant have **never been assessed on peptides at all**.

**2. Determinism — fails, and fails specifically on GPU.** No Vina-family GPU tool documents
bit-identical output. QuickVina2-GPU 2.1 does not list `--seed` in its argument table at all and
makes no reproducibility statement anywhere; its only determinism "evidence" is one sentence from a
developer on an issue that is still open. The AutoDock-GPU maintainer has stated in writing that
**fixed seeds do not guarantee identical GPU results**, by design, because RNG streams diverge under
out-of-order work-unit scheduling — before floating-point reduction order is even considered. The
one published number for run-to-run spread in this family is **0.1–0.5 kcal/mol per run**, stated by
the Vina-CUDA authors to apply to Vina-GPU 2.1 as well. CPU Vina is the only member with a stated
fixed-seed reproducibility property, and that property was **broken in v1.2.3 and repaired in
v1.2.4** — with a related macrocycle-glue-atom scoring bug (absurd −50 kcal/mol energies) living in
exactly the code path this ligand class would use.

**3. Objective validity — fails.** Ciepliński et al. (2023) showed the docking score is
"moderately strongly" correlated with molecular weight and rotatable-bond count, and that the
degenerate optimum of the full score for a CVAE was **macrocycles** — which here are not an artifact
to be filtered out but the target chemistry itself. The size filter that made their benchmark
honest is unavailable.

**The single hardest requirement to recover is #2 on GPU.** Determinism on CPU Vina is in principle
achievable by the DOCKSTRING method — freeze the receptor and box as files, freeze protonation,
generate exactly one conformer, fix stereochemistry assignment, pin library versions, fix the seed,
and verify with a re-docking regression test. That is a recipe for making the *inputs* identical,
not for making the engine deterministic, and it has never been validated at this flexibility
(DOCKSTRING's 0.1 kcal/mol seed-insensitivity measurement was on drug-like ZINC ligands, and
peptide docking is documented not to converge with respect to sampling even at exhaustiveness 1024).
Moving to GPU to afford thousands of calls reintroduces non-determinism that the maintainers
explicitly decline to guarantee away. The throughput requirement and the determinism requirement
are therefore in direct conflict, and requirement #1 means that even if both were satisfied, the
resulting number would be a reproducible measurement of something that is not the binding pose.
