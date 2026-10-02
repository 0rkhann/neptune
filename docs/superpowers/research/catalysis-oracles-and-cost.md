# Semiempirical QM descriptors as generative-design oracles, and what catalysis data exists

Research notes, compiled 2026-10-02. Every claim carries a URL and year. Quotations are short phrases only.
Numbers labelled **[measured here]** were produced on this machine during this research and are reproducible from the commands given; everything else is from the cited source.

---

## SECTION 1 — xTB cost, measured

### 1.1 Direct measurement on this machine [measured here]

Hardware: Intel Core i7-14700HX (28 logical cores), Linux. Software: `xtb` 6.7.1 (locally compiled, 2025-04-24).
Inputs: PDB Chemical Component Dictionary "ideal" SDF geometries (https://files.rcsb.org/ligands/download/<ID>_ideal.sdf) and PubChem 3D SDF for amoxicillin (CID 33613, https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/33613/SDF?record_type=3d).
Commands: `xtb in.sdf --gfn 2 --sp|--opt --parallel N` and `xtb in.sdf --gfnff --sp|--opt --parallel 1`, gas phase, default accuracy and default (`normal`) optimisation thresholds. Wall time measured around the process.

| molecule | heavy atoms | GFN2 single point, 1 core | GFN2 single point, 4 cores | GFN2 `--opt`, 1 core | GFN2 `--opt`, 4 cores | GFN-FF single point, 1 core | GFN-FF `--opt`, 1 core |
|---|---|---|---|---|---|---|---|
| amoxicillin | 25 | 0.30 s | 0.25 s | 2.97 s | 3.86 s | 0.05 s | 0.21 s |
| ritonavir (RIT) | 50 | 0.65 s | 0.76 s | 34.6 s | 30.9 s | 0.08 s | 0.79 s |
| erythromycin A (ERY) | 51 | 0.81 s | 0.95 s | 51.4 s | 38.9 s | 0.09 s | 0.80 s |
| tacrolimus (FK5) | 57 | 0.96 s | 0.85 s | 39.8 s | 32.5 s | 0.09 s | 0.99 s |
| rapamycin (RAP) | 65 | 1.36 s | 1.11 s | 346 s | 261 s | 0.10 s | 3.40 s |
| cyanocobalamin (B12, contains Co) | 91 | 4.69 s | 4.47 s | 601 s | 835 s | 0.18 s | 5.18 s |

Note on B12: a 180-atom Co(III) corrinoid is outside the organic regime of interest here, so treat it as a scaling anchor rather than a representative case. Two things it shows: the GFN2 **SCF** step grows steeply (≈ O(N³)) from 1.4 s at 65 heavy atoms to 4.7 s at 91; and the GFN2 **optimisation** took 601 s on 1 core and *more* on 4 cores (835 s) — at this size the optimiser trajectory is not thread-deterministic, which is itself a reproducibility datum for §2. GFN-FF stayed at 0.18 s (single point) and 5.2 s (optimisation), i.e. ~26× and ~120× cheaper.

Readings that matter for the oracle-budget decision:

1. **GFN2-xTB single point on one conformer is sub-second to ~1.5 s** for 25–65 heavy atoms on one core, and essentially flat in thread count at this size (threading overhead cancels the gain below ~100 atoms; the same effect is reported by xtb users and developers at https://github.com/grimme-lab/xtb/issues/1195 , 2025).
2. **GFN-FF single point is 0.05–0.10 s**, i.e. of order 0.1 s, not seconds — a factor 7–15 below GFN2 at these sizes. GFN-FF is a force field with an iterative Hückel treatment for π atoms, so it yields charges, bond orders and geometries but **no SCF orbital energies**; HOMO/LUMO-derived descriptors are not available from GFN-FF alone (method paper: https://onlinelibrary.wiley.com/doi/10.1002/anie.202004239 , 2020; a GFN2 single point on a GFN-FF geometry is the standard patch, see CREST composite calculators, https://crest-lab.github.io/crest-docs/page/examples/composite.html ).
3. **Geometry optimisation is where the cost is, and it is the least predictable step.** GFN2 `--opt` ranged from 3 s to 346 s over five organic molecules of 25–65 heavy atoms — two orders of magnitude of spread, only loosely related to size. These runs start from CCD "ideal" / PubChem 3D geometries, which are far from the GFN2 minimum, so the cycle counts are a pessimistic case. The spread, not the median, is the planning risk.
4. GFN-FF optimisation is 0.2–3.4 s, 30–100× cheaper than GFN2 optimisation on the same inputs.

**Follow-up measurements on the two obvious cost levers [measured here]**, 1 core, same inputs:

| lever | ritonavir (50 heavy) | erythromycin (51) | rapamycin (65) |
|---|---|---|---|
| GFN2 `--opt` from the raw CCD geometry (baseline) | 34.6 s | 51.4 s | 346 s |
| GFN2 `--opt crude` from the raw geometry | **8.0 s** | **10.7 s** | not run |
| GFN-FF `--opt` then GFN2 `--opt` (total = FF time + the figure shown) | 57.8 s | 45.5 s | **71.0 s** |

Two findings worth stating because they are counter-intuitive: (a) a **GFN-FF pre-optimisation does not reliably reduce** the subsequent GFN2 optimisation cost — it made ritonavir *worse* (34.6 → 57.8 s) and erythromycin only slightly better — but it removes the pathological tail (rapamycin 346 → 71 s, ~5×), which is exactly the variance that breaks a fixed compute budget; (b) relaxing the convergence threshold to `crude` buys a clean **4–5×** on 50-heavy-atom molecules and is the cheaper lever if the descriptor tolerates a looser geometry.

### 1.2 Published timings, same quantities

- **GFN-FF vs GFN2 speed, general**: GFN-FF reaches O(N²) scaling versus O(N³) for all GFNn-xTB methods, and is reported as "about two orders and three orders of magnitude faster than GFN0-xTB and GFN1-/GFN2-xTB" for single point + gradient on proteins of 300–6000 atoms (Bannwarth et al., *WIREs Comput. Mol. Sci.* 2021, https://wires.onlinelibrary.wiley.com/doi/10.1002/wcms.1493 ; summarised at https://exa.ai/library/publication/cxmvcy7p7gx ). For thermostatistics, GFN-FF is quoted as 2–3 orders of magnitude faster than GFN2-xTB and ~5 orders faster than DFT (Spicher & Grimme, *J. Phys. Chem. Lett.* 2020, https://pubs.acs.org/doi/abs/10.1021/acs.jpclett.0c01930 ).
- **Geometry optimisation, CPU time per molecule in a full pipeline**: a 2025 GFN benchmark on QM9-derived and Harvard CEP molecules reports average CPU times per molecule of GFN-FF 14.9 s, GFN0 41.0 s, GFN2 130.9 s, GFN1 168.1 s (QM9 subset) and GFN-FF 56.3 s, GFN0 302.2 s, GFN2 1358.5 s, GFN1 1608.9 s (CEP subset), against BP86/def2-SVP DFT at ~9931 s; hardware "two 16-core Intel Xeon Gold 6142@2.6 GHz", 8 cores per task (https://arxiv.org/html/2505.09606v2 , 2025). **Important**: their pipeline is force-field pre-opt → xTB pre-opt → **CREST conformer search** → final xTB optimisation with ALPB solvation, so these are pipeline costs, not single-optimisation costs. They are therefore a reasonable proxy for "oracle call that includes a conformer search", and they are 1–3 orders of magnitude above a single-point call.
- **Drug-sized minimisation in a commercial stack**: Cresset report GFN2-xTB minimisation + single point on biomolecules up to erythromycin as "on the order of minutes and not days", on an 8-core CPU (https://cresset-group.com/about/news/minimization-using-semi-empirical-gfn2-xtb/ , 2023).
- **In-loop use in a generative run**: a 2026 inverse-design study using GFN2-xTB IP/EA as the in-loop reward describes "subsecond cost per molecule" for frontier-orbital properties and evaluated 179,316 valid 3D geometries out of 360,846 candidates (https://exa.ai/library/publication/4vwtsbm586x , 2026). This is the existence proof that a single-point xTB oracle runs at a 10⁵-call budget.

### 1.3 CREST conformer-ensemble workflow cost

- CREST's own documentation states a conformer search for alanylglycine, **20 atoms**, needs **~2.8 × 10⁵ energy+gradient calls**, taking **18 s at GFN-FF** (8 threads; total wall time 17.96 s, 2.8463 × 10⁵ calls) — "at a DFT level it would take weeks" (https://crest-lab.github.io/crest-docs/page/examples/composite.html , accessed 2026).
- The CREST 3.0 paper states that iMTD-GC on n-octane, amoxicillin and Tamiflu requires energy+gradient evaluations "on the order of O(10⁵)", with wall times measured on "16 threads of Intel Xeon Gold 6326 CPUs (2.90 GHz)"; CREST 3.0 is 1.6–4.7× faster than 2.12 (Pracht et al., *J. Chem. Phys.* 160, 114110, 2024; https://pubs.aip.org/aip/jcp/article/160/11/114110/3278084 , full text PDF https://publications.rwth-aachen.de/record/986468/files/986468.pdf ). The figure itself is logarithmic and the numeric bar values are not given in text.
- An HPC centre's empirical guidance, the most concretely useful figure found for peptide-sized systems: "for a flexible organic molecule of 65 atoms, conformational search using GFN-FF method and 24 cores took about 15-20 minutes and semiempirical GFN2 needed 5-8 hours" (TalTech HPC docs, https://docs.hpc.taltech.ee/chemistry/crest.html , accessed 2026).
- Kraken built its ligand library with CREST (GFN2-xTB, GBSA toluene) and reports **33,489 unique conformers for 190 DFT-level ligands, 20.8 conformers per ligand on average** (https://descriptor-libraries.molssi.org/kraken/ , accessed 2026) — i.e. the ensemble multiplies every downstream per-conformer cost by ~20.

### 1.4 RDKit ETKDG conformer generation cost

- RDKit's own maintained benchmark page ("How Long Things Take", release 2025.09.1, machine "stoat" = Intel i9-12900K, 64 GB RAM, Ubuntu): generating **one** conformer for **1000 molecules** takes **44.1 s ± 0.17 s with plain DG** and **1 min 5 s ± 0.245 s with ETKDGv3** — i.e. **44 ms and 65 ms per molecule** (https://github.com/rdkit/rdkit/wiki/How-Long-Things-Take and https://greglandrum.github.io/rdkit-blog/posts/2025-10-31-how-long-does-it-take.html , 2025).
- Multi-conformer cost, same author, 100 drug-like molecules (8 physical performance cores): 100 conformers each = 159.0 s on 1 thread → 26.7 s on 8 threads; 400 conformers each = 629.9 s on 1 thread → 98.1 s on 8 threads. Scaling is good to ~6 threads and poor beyond 8 (https://greglandrum.github.io/rdkit-blog/posts/2025-08-30-confgen-scaling.html , 2025). Per molecule: ~1.6 s for 100 conformers, ~6.3 s for 400 conformers, single-threaded.
- Historical calibration: ETKDG costs ~1.7–3× plain DG per conformer but needs roughly a quarter of the conformers for equivalent crystal-structure recovery; DG+UFF optimisation takes ~1.97× the ETKDG runtime (Riniker & Landrum, *JCIM* 2015, discussed at https://www.blopig.com/blog/2016/06/advances-in-conformer-generation-etkdg-and-etdg/ and http://rdkit.blogspot.com/2017/05/looking-at-platinum-dataset.html ).

### 1.5 What one oracle call costs — the arithmetic

For a 60–120-heavy-atom bifunctional peptidomimetic, at a 10,000-call budget on one 8-core node:

| oracle definition | per call | 10,000 calls |
|---|---|---|
| ETKDG 1 conformer + GFN-FF single point | ~0.15 s | ~25 min |
| ETKDG 1 conformer + GFN2 single point (no optimisation) | ~0.7–1.5 s | 2–4 h |
| ETKDG 1 conformer + MMFF + GFN2 optimisation + GFN2 single point | ~10–350 s (measured spread) | 1–40 CPU-days |
| ETKDG 20 conformers + GFN-FF opt each + GFN2 single point on the best | ~5–10 s | 14–28 h |
| CREST (GFN-FF) ensemble + GFN2 re-ranking | ~15–20 min (65-atom reference) | ~3–14 CPU-months |
| CREST (GFN2) ensemble | 5–8 h (65-atom reference) | infeasible |

The decision therefore turns on one line: a **single-point** descriptor on a **single, fixed** conformer is affordable at 10⁴ calls (minutes to hours); anything that re-optimises geometry per call is 1–2 orders of magnitude more expensive and, worse, has 100× variance per molecule; a CREST ensemble per call is not a benchmark oracle at any budget that includes thousands of evaluations per run.

---

## SECTION 2 — Conformer determinism

### 2.1 DOCKSTRING's freezing protocol (the template)

DOCKSTRING (García-Ortegón et al., *J. Chem. Inf. Model.* 62, 3486–3502, 2022; https://pubs.acs.org/doi/full/10.1021/acs.jcim.1c01334 , open text at https://pmc.ncbi.nlm.nih.gov/articles/PMC9364321/ ) fixes every stochastic degree of freedom in the ligand pipeline. Its steps, in order:

1. **Sanity filter** — reject radicals and multi-fragment ligands.
2. **Protonation** — (de)protonate at pH 7.4 with Open Babel, called out as a step that competing wrappers omit.
3. **One conformer only** — "a single three-dimensional (3D) conformation is generated with the Euclidean distance geometry algorithm ETKG as implemented in RDKit".
4. **Force-field refinement** — MMFF94 on that single conformer.
5. **Stereochemistry rule** — defined stereocentres are preserved; undefined ones are "assigned randomly (but consistently across different runs to ensure the reproducibility of docking scores)".
6. **Charges and file conversion** — Gasteiger charges, PDBQT via Open Babel.
7. **Seed audit, then seed freeze** — they first measured the seed dependence: "no target–ligand combination for which the docking scores deviated by more than 0.1 kcal/mol", and then "fixed the random seed to obtain a fully deterministic pipeline".
8. **Fixed docking settings** — AutoDock Vina with default exhaustiveness 8, 9 modes, energy range 3; the lowest of up to 9 scores is the reported value.
9. Cost disclosed: "computing a score with eight CPUs takes around 15s".

The explicit contrast they draw is with TDC and Cieplinski et al., which "[do not] control sources of randomness during the docking procedure (e.g., random seeds input into the docking program or the conformer generation routines), leading to the potential for considerable variance between runs on the same molecule".

### 2.2 Evidence that conformer handling moves computed descriptors

- **Same ensembles, any answer you like.** For transition-state ensembles under Curtin–Hammett treatment, processing the *same* computed TS ensembles in different but defensible ways gives "virtually any selectivity": Boltzmann-weighting the 10 lowest TSs gives ΔΔG‡ = 0.48 kcal/mol (31:69), the 20 lowest gives 0.07 kcal/mol (49:51), and correct filtering of non-interconvertible families gives 1.75 kcal/mol (5:95) — a ~1.7 kcal/mol swing driven purely by the selection rule ("Overcoming the Pitfalls of Computing Reaction Selectivity from Ensembles of Transition States", 2024, https://pmc.ncbi.nlm.nih.gov/articles/PMC11284845/ ; the `marc` tool is their proposed fix).
- **Conformer-search randomness alone changes published descriptors.** The refactored Kraken workflow states: "The conformer search produces slightly different conformers each run, so results vary slightly (around 1%) between runs", and, comparing DFT properties for 60+ monophosphines between workflow versions, "For a small number of properties (typically octant volumes), the difference between the old and new values for a given monophosphine exceeds 75% of the original value. These differences likely arise from randomness in the conformer search" (https://github.com/SigmanGroup/kraken , 2024).
- **Program-version changes silently change descriptors.** Same source: "Several descriptors vary substantially with xTB 6.7.0 or greater (EA/IP descriptors, nucleophilicity) because IPEA-xTB is not used for vertical IP/EA calculations", and "CREST 2.12 produces many more conformers than CREST 2.8".
- **Ensemble-condensation rule is a modelling choice, not a detail.** Kraken publishes five different condensed measures per property — Boltzmann average, max, min, std, the value at the minimum-buried-volume conformer, and max−min delta — precisely because a single number per ligand is underdetermined (https://descriptor-libraries.molssi.org/kraken/library_details ). Its conformer selection rule is explicit: conformers that minimise or maximise any of ~19 xTB steric properties, plus up to 20 conformers within 3 kcal/mol, RMSD-pruned.
- **CREST ensembles over-count flexibility relative to DFT.** For 24 Rh precatalysts, CREST(GFN2//GFN-FF) produced 678 conformers in total (average 23 per ensemble) which collapsed to an average of **2 per ensemble** after PBE0-D3(BJ)/def2-SVPP optimisation; energy-based filtering of the xTB ensemble was found "ineffective", geometry-based (DBSCAN/RMSD) filtering worked better (https://pmc.ncbi.nlm.nih.gov/articles/PMC12120983/ , 2025).
- **Conformer quality propagates into learned representations.** "The impact of conformer quality on learned representations of molecular conformer ensembles" (https://arxiv.org/abs/2502.13220 , 2025) studies exactly this dependence for ensemble-based ML models.
- Related: "Impact of Model Selection and Conformational Effects on the Descriptors for In Silico Screening Campaigns: A Case Study of Rh-Catalyzed Acrylate Hydrogenation" (https://pmc.ncbi.nlm.nih.gov/articles/PMC12025388/ ) and "Bisphosphine Ligand Conformer Selection to Enhance Descriptor Database Representation" (2025, https://exa.ai/library/publication/qvl6r6qtm7l ) both report descriptor values shifting with the conformer-selection protocol.

### 2.3 What must be frozen for a conformer-dependent oracle to reproduce across machines

Pinning, in the order in which each one bites:

1. **The input string normalisation** — canonical SMILES from a pinned toolkit version, explicit protonation state and formal charges (a superbase oracle is meaningless if the base is silently protonated differently), explicit stereo on every stereocentre, and a fixed rule for undefined stereo (DOCKSTRING: random but seeded).
2. **The embedding** — ETKDG variant (ETKDG / v2 / v3 / srETKDGv3 are different algorithms), `randomSeed`, `numThreads = 1` (multithreaded embedding is order-dependent; RDKit's own scaling study runs conformers in parallel and notes they do not pruning-independently), `useRandomCoords`, and the number of conformers requested.
3. **The FF relaxation** — force field (MMFF94 vs MMFF94s vs UFF), convergence thresholds, iteration cap.
4. **The ensemble rule** — energy window, RMSD-pruning threshold, max conformers kept, and the condensation statistic (Boltzmann at which T, min, max, delta). Per §2.2 this choice alone is worth >1 kcal/mol.
5. **The QM step** — xtb version string (6.7.x changed IP/EA behaviour), `--acc`, electronic temperature, solvation model and solvent, optimisation level, and whether IPEA-xTB or GFN2 is used for vertical IP/EA.
6. **Numerics and environment** — thread count (results of parallel optimisation loops can reorder), BLAS/LAPACK implementation, and OS/compiler, since xtb is compiled locally.
7. **Publication practice** — distribute a frozen reference set of (SMILES → descriptor) pairs plus a container image, as DOCKSTRING distributes a precomputed 260,000-molecule × 58-target matrix so that users can check their installation reproduces it.

The cheapest route to full determinism, and the one DOCKSTRING chose, is **one conformer, one seed, no ensemble**. An ensemble oracle can be made reproducible, but only by freezing items 2 and 4 as hard as the seed.

---

## SECTION 3 — Conceptual-DFT descriptors in catalysis ML

### 3.1 Origin papers

- **Electrophilicity index ω = μ²/2η (= χ²/2η)**: Parr, von Szentpály, Liu, "Electrophilicity Index", *J. Am. Chem. Soc.* 121, 1922–1924 (1999), prompted by Maynard et al., *PNAS* 95, 11578 (1998). It "measures the second-order energy change of an electrophile as it is saturated with electrons" (https://exa.ai/library/publication/h150k4pbt8k ; ~7,500 citations).
- **Nucleophilicity index N = E_HOMO(Nu) − E_HOMO(TCE)**: Domingo and co-workers (2008), an empirical index referenced to tetracyanoethylene so that values are positive, built on Koopmans/Kohn–Sham HOMO energies (review: Domingo, Ríos-Gutiérrez, Pérez, *Molecules* 21, 748, 2016, https://pmc.ncbi.nlm.nih.gov/articles/PMC6273244/ ).
- Underlying quantities: electronic chemical potential μ (Parr et al., *J. Chem. Phys.* 68, 3801, 1978) and chemical hardness η (Parr & Pearson, *JACS* 105, 7512, 1983), with μ ≈ −(I+A)/2, η ≈ (I−A)/2, and I ≈ −E_HOMO, A ≈ −E_LUMO.
- Scale thresholds (B3LYP/6-31G(d), the level the scales were defined at): strong electrophiles ω ≥ 1.5 eV, moderate 0.8–1.5, marginal < 0.8, superelectrophiles ≥ 4.0; strong nucleophiles N ≥ 3.0 eV, moderate 2.0–3.0, marginal < 2.0, supernucleophiles ≥ 4.0 (Domingo, *Scientiae Radices* 3(3), 2024, https://sci-rad.com/scirad2024v3i3a02/ ).

### 3.2 What is validated against measured reactivity

Validated, repeatedly, against **small-molecule kinetics and barriers**:

- ω vs computed Diels–Alder activation energies for twelve ethylenes: **R² = 0.92**; ω of benzaldehydes vs barriers for cyanide attack: **R² = 0.95**; N of 5-substituted indoles vs experimental ln k with a benzhydryl cation: **R² = 0.89**; N vs barriers for a 12-member bicyclic diene series: **R² = 0.99** (Domingo review 2016, https://pmc.ncbi.nlm.nih.gov/articles/PMC6273244/ ; Domingo 2024, https://sci-rad.com/wp-content/uploads/2024/09/scirad2024v3i3a02.pdf ).
- Against Mayr's experimental E/N scales (Mayr & Patz, *Angew. Chem. Int. Ed.* 1994, https://onlinelibrary.wiley.com/doi/10.1002/anie.199409381 ): a 2025 reformulation using cubic energy expansions correlates CDFT indices with Mayr parameters for benzhydrylium ions and 15 nucleophiles and reports "the descriptors correlated better with the electrophilicity parameter (r2 = 0.981) than with the nucleophilicity parameter (r2 = 0.827)" (https://pubs.rsc.org/en/content/articlehtml/2025/cp/d5cp00994d ).
- ML on the same target: "Predicting experimental electrophilicities from quantum and topological descriptors" (*J. Comput. Chem.*, 2020, https://onlinelibrary.wiley.com/doi/10.1002/jcc.26376 ).
- Two caveats from the same literature. (i) Correlations degrade across chemical families: a systematic study of theoretical electrophilicity measures against Mayr data concludes correlations at HF/6-31G(d) "are of semi-quantitative value" and that "working with families of compounds with similar functional groups is indispensable" (https://repositorio.uchile.cl/handle/2250/150056 ). (ii) The numbers are **level-of-theory specific**: the reference scales were fixed at B3LYP/6-31G(d) and a 2023 paper exists purely to regress ω and N across 48 other DFT levels so they can be compared at all (https://www.sciencedirect.com/org/science/article/pii/S0894323023001273 , 2023).

### 3.3 Do these descriptors predict *organocatalyst* performance?

Short answer: **no direct evidence was found that global CDFT indices (ω, N, η, μ, gap) predict organocatalyst rate, yield or ee.** What exists is weaker and of three kinds.

1. **They appear as features inside larger descriptor sets.** The 2024 review of ML in enantioselective organocatalysis lists, among typical inputs, "Charton or Sterimol values, NBO charges, NMR chemical shifts, bond distances and angles, HOMO–LUMO gaps, local electro/nucleophilicity, or RDKit descriptors", with outputs "ΔΔG‡, e.e., or yield" (Beilstein J. Org. Chem. 20, 2024, https://www.beilstein-journals.org/bjoc/content/pdf/1860-5397-20-196.pdf ). So they are used; the review does not identify them as the carrying feature.
2. **They are used to map catalyst space, validated against acidity/basicity rather than catalysis.** OSCAR uses the **nucleophilicity N-index** as one of two axes for its NHC chemical-space map and the **LUMO energy** for H-bond donors, justifying them as "an indirect estimate of the catalysts' Brønsted acidity/basicity" and checking against experimental pK_a of 23 NHC precursors (https://pubs.rsc.org/en/content/articlehtml/2022/sc/d2sc04251g , 2022). That is descriptor → thermodynamic property, not descriptor → catalytic outcome.
3. **The models that actually predict ee well use other descriptors.** Denmark's CPA models (MAD 0.161–0.238 kcal/mol on external test sets, q² = 0.748) are built on grid-based average steric occupancy (ASO) and electrostatic-indicator fields, not CDFT globals (*Science* 2019). Wennemers/Denmark peptide models (MAE_test 0.22 kcal/mol for ee, 0.10 for dr) use ASO/AEIF descriptors (*ACS Cent. Sci.* 2024). Miller/Sigman peptide models use conformer-specific geometric and NBO parameters (*JACS* 2018). The best transferable models in 2026 use **transition-state-derived** features explicitly because "simple stereoelectronic parameters may fail to describe mechanistically complex transformations" (*Nature*, 2026, https://www.nature.com/articles/s41586-026-10239-7 ).

For a superbase/H-bond-donor bifunctional catalyst specifically, the physically closest validated quantity is not ω or N but **computed pK_a / proton affinity**, where DFT protocols reach R² > 0.99 and MUE ≤ 1.0 pK_a unit against experiment for phosphorus, nitrogen and carbon bases (https://www.mdpi.com/1422-0067/23/18/10576 , 2022). That is a defensible oracle target with literature backing; "electrophilicity index of the catalyst correlates with ee" is not.

### 3.4 Does xTB track DFT closely enough to substitute?

Mixed, and property-dependent.

| quantity | xTB vs reference | source |
|---|---|---|
| HOMO–LUMO gap, QM9-like small molecules | MAD vs B3LYP: GFN0 1.26 eV, GFN1 1.32 eV, **GFN2 1.54 eV**, GFN-FF 2.03 eV; gaps systematically **underestimated** | https://arxiv.org/html/2505.09606v2 (2025) |
| HOMO–LUMO gap, extended π systems (CEP) | MAD: GFN1 0.091 eV, GFN-FF 0.124 eV, **GFN2 0.143 eV**, GFN0 0.278 eV | same |
| IP, EA, gap, optical gap, dipole (50 polyimide repeat units) | Pearson R vs B3LYP/aug-cc-pVTZ: **EA 0.98, IP 0.93, gap 0.86**, optical gap 0.89, dipole 0.76; but absolute IP off by ~1.6–1.8 eV and EA by ~0.7 eV | https://exa.ai/library/publication/4vwtsbm586x (2026) |
| Redox potentials, ROP313 (313 experimental values) | MAD vs experiment: **GFN2 0.30 V**, GFN1 0.31 V, PM6-D3H4 0.61 V, PM7 0.60 V, B97-3c 0.25 V; organometallics much worse (GFN2 0.74 V) | Neugebauer, Bohle, Bursch, Hansen, Grimme, 2020, https://exa.ai/library/publication/cl9xyrpfxxr |
| Vertical IP specifically | GFN1 IPs "are not sufficiently accurate", which is why the special-purpose **IPEA-xTB** reparameterisation exists; GFN2 IPs are "qualitatively correctly" computed | Koopman & Grimme, *ACS Omega* 4, 15120 (2019), https://pubs.acs.org/doi/pdf/10.1021/acsomega.9b02011 |
| Conformational energies (relevant if the descriptor is Boltzmann-averaged) | CREST 3.0 paper, correlation to ωB97X-V reference: **GFN2 R² = 0.889, r_s = 0.925; GFN-FF R² = 0.687, r_s = 0.787** | https://publications.rwth-aachen.de/record/986468/files/986468.pdf (2024) |
| Stability of xTB's own values | Kraken: "Several descriptors vary substantially with xTB 6.7.0 or greater (EA/IP descriptors, nucleophilicity) because IPEA-xTB is not used for vertical IP/EA calculations" | https://github.com/SigmanGroup/kraken (2024) |

Reading: **ranking** is largely preserved (R ≈ 0.86–0.98 for IP/EA/gap; Spearman 0.93 for conformer energies at GFN2), **absolute values are not** (eV-scale systematic offsets for IP/EA, 1.5 eV MAD on small-molecule gaps). For a benchmark oracle whose job is to rank candidates monotonically, that is usable; for any claim about an absolute ω or N on the Domingo B3LYP/6-31G(d) scale, it is not — and the thresholds in §3.1 are defined on that scale. Note also that newer **g-xTB** reportedly halves GMTKN55 WTMAD-2 (9.3 vs ~18 kcal/mol) and specifically improves gaps at ~30% extra cost; this was found only in a secondary summary (https://docs.mqs.dk/sections/section_005_qc/ ) and should be checked against the primary paper before relying on it.

---

## SECTION 4 — Superbase and peptide catalysis data

### 4.1 Inventory

| resource | size | measured quantities | machine-readable? | licence | download |
|---|---|---|---|---|---|
| **Wennemers tripeptide ML study** (Schnitzer, Schnurr, Zahrt, Sakhaee, Denmark, Wennemers, *ACS Cent. Sci.* 10, 367–373, 2024) | **200 experiments** (50 dPro-containing peptides × 4 aldehyde/nitroolefin combinations); separately a **universal training set of 161 catalysts**; in-silico library of **30,276** tripeptides from 174 commercial amino acids | conversion, dr (58:42–98:2), ee (10–98%) | SI only (synthetic protocols, analytical data, computational details); descriptor profiles reported as study output | ACS (open access on PMC) | https://pmc.ncbi.nlm.nih.gov/articles/PMC10906243/ ; https://doi.org/10.1021/acscentsci.3c01284 ; listed at https://moleculemaker.org/datasets/machine-learning-todeveloppeptide-catalysts-successes-limitations-and-opportunities/ |
| **Wennemers primary catalysis literature** (Pro-Pro-Xaa aldol and 1,4-additions, 2008–2013+) | screens of order 10–15 tripeptides per paper (e.g. "a screening of a collection of 15 tripeptides of the type Pro-Pro-Xaa") | yield, dr, ee | no — SI tables in individual papers | journal | e.g. https://onlinelibrary.wiley.com/doi/10.1002/anie.200704972 (2008); review https://ojs.chimia.ch/chimia/article/download/2013_279/4672/15357 (2013) |
| **Miller-group peptide catalysts, atroposelective bromination** | initial library **24 sequences**, expanded to **54 peptides**; **35** peptides structurally characterised (53 distinct conformational states) | er/ee (up to 97:3 er at 1 mol %), plus X-ray/NMR conformations | no — SI and CSD depositions | journal / CCDC | https://pubs.acs.org/doi/full/10.1021/jacs.6b11348 (2017); https://pubs.acs.org/doi/full/10.1021/jacs.7b11303 (2018) |
| **Miller/Sigman peptide parameterisation** (Crawford, Stone, Metrano, Miller, Sigman, *JACS* 140, 868–871, 2018) | multivariate models from an **initial training set of nine catalysts** varying only at i+2, split by β-turn type (I′ vs II′) | ΔΔG‡ from ee | no — SI only | journal | https://pubs.acs.org/doi/10.1021/jacs.7b11303 ; open copy S-EPMC5817992 |
| **OSCAR** (Gallarati, van Gerwen, Laplaza, Vela, Fabrizio, Corminboeuf, *Chem. Sci.* 13, 13782, 2022) | **4,000** experimentally derived organocatalysts (990 literature seed + 3,010 CSD-extracted); **8,622** combinatorial NHCs; **1,573,015** dual-H-bond donors at xTB level (1,000 per motif at DFT); IP/EA for a **2,060**-catalyst subset at DLPNO-CCSD/cc-pVTZ | computed stereoelectronic descriptors incl. conceptual-DFT reactivity indices, %V_bur, LUMO energy, nucleophilicity N-index. **No experimental catalytic outcomes.** | yes — CSV descriptors + XYZ + Chemiscope JSON | open (Materials Cloud, open access article) | https://pubs.rsc.org/en/content/articlehtml/2022/sc/d2sc04251g ; https://archive.materialscloud.org/records/x12qf-mr670 (doi:10.24435/materialscloud:gy-3h) |
| **iSynth / OrgAIcat** (Yang, Liu, Zhang, Zhang, Chen, Luo, *CCS Chem.*, 2026) | **>22,000** aminocatalytic reactions: **13,360 aldol** from 478 publications + **9,194 Michael** from 366 publications, 2000–2020, 844 papers total | ee → ΔΔG‡, plus catalyst, substrates, solvent, additive, temperature | yes — .xlsx inputs + CSV descriptors in repo | **MIT** (code repo) | https://github.com/deepsynthesis/orgaicat ; summary https://sciencesources.eurekalert.org/news-releases/1143161 |
| **Asymmetric organocatalytic Mannich dataset** (Vasechkin, Nikitin, Fedorov, 2026) | **4,015 reactions**, 153 publications (2000–2025), **468 unique organocatalysts** | ee and ΔΔG‡, plus conditions | yes — explicitly "machine-readable form suitable for training AI-models" | per publication (data paper) | https://exa.ai/library/publication/qcvdhhc8slq |
| **Pictet–Spengler organocatalysis set** (NaviCatGA + OSCAR generality study) | **820 reactions** curated | selectivity and activity | dataset deposited | open (B2FIND record) | https://b2find.eudat.eu/dataset/66b9055f-437e-5bad-b0c0-814f2d743a72 |
| **ΔΔG‡-vs-ee benchmark collection** (Giessen, 2023) | **10 datasets, 2,761 data points** across hydrogenation, Suzuki, Heck and organocatalytic reactions | ee and ΔΔG‡ | yes — CSVs (582 KB zip) + notebooks | institutional repository, open | https://jlupub.ub.uni-giessen.de/items/84928dcc-079e-43ba-a200-d9dedecbe179/full |
| **Denmark chiral phosphoric acid (CPA) set** (Zahrt et al., *Science* 363, eaau5631, 2019) | **1,075 reactions** = 43 CPA catalysts × 25 imine/thiol products, each run in duplicate (**2,150 experiments**); in-silico catalyst library >800 | ee → ΔΔG‡ | SI data files S1–S3 (tabular) | journal SI; also redistributed in ML papers | https://www.science.org/doi/10.1126/science.aau5631 ; open text https://pmc.ncbi.nlm.nih.gov/articles/PMC6417887/ |
| **Doyle/Dreher Buchwald–Hartwig HTE** (Ahneman et al., *Science* 360, 186–190, 2018) | **4,608 reactions** run; **3,960** full-factorial combinations modelled; the ORD copy holds **4,312** | yield | yes | open via ORD | ORD dataset "Ahneman", https://open-reaction-database.org ; paper https://pubmed.ncbi.nlm.nih.gov/29449509/ |
| **Perera Suzuki flow HTE** (*Science* 359, 429–434, 2018) | **5,760 reactions** (7,392 possible combinations; two-stage design → 4,608 in some analyses) | yield | yes | open via ORD | ORD "Rapid material-sparing screening of 5760 Suzuki-Miyaura coupling reactions" |
| **Curated BH HTE superset** (Neves et al., *Nat. Comput. Sci.*, 2026) | **~27,500** Pd C–N couplings (11,300 new J&J + 16,200 curated open) | yield, binary success | yes — CSV/pickle | open (Zenodo) | https://zenodo.org/records/19636649 ; https://github.com/schwallergroup/bh-hte-ood |
| **Asymmetric hydrogenation of olefins (AHO) database** (Xu, Zhang, Li, Tang, Xie, Hong, *Angew. Chem.* 2021) | **>12,000** literature transformations | ee | yes — zip in repo + web server | repo (see repo terms) | https://github.com/licheng-xu-echo/AHO/ ; http://asymcatml.net/ |
| **N,N′-dioxide/metal asymmetric Michael platform** (2025) | **>2,000** reactions curated | ee | yes — "chemically annotated, machine-readable dataset" | per publication | https://exa.ai/library/publication/9g19bf98kp8 |
| **Superbase basicity (experimental)** — phosphazenes, guanidines, guanidinophosphazenes, proton sponges | self-consistent MeCN basicity scale spanning **28 pK_a units** (Kaljurand et al., *J. Org. Chem.* 70, 1019, 2005); gas-phase superbasicity scale spanning **20 orders of magnitude**, most basic GB = 273.9 kcal/mol (Kaljurand et al., 2016) | pK_a (MeCN, THF), gas-phase GB/PA | SI tables + a maintained PDF compilation | academic, free to view | https://analytical.chem.ut.ee/HA_UT/Basicities_GB_pKa_of_Superbases.pdf ; https://exa.ai/library/publication/ysdfy4zmk8r (2016) |
| **Superbase basicity (computed)** | e.g. DFT model sets of phosphorus/nitrogen/carbon bases with R² > 0.99 vs experiment, MUE ≤ 1.0 pK_a unit (Glasovac et al., *IJMS* 23, 10576, 2022); basicity-limit predictions for whole families (Leito, Koppel et al., *Angew. Chem.* 2015) | computed pK_a(MeCN), PA, GB | tables in papers/SI | journal | https://www.mdpi.com/1422-0067/23/18/10576 ; https://onlinelibrary.wiley.com/doi/10.1002/ange.201503345 |

### 4.2 The gap that matters for this project

- **There is no public dataset of phosphazene- or guanidine-superbase *catalysis*.** Searches surface two disjoint literatures: (i) quantitative **basicity** scales (Leito/Koppel/Kaljurand, Glasovac/Maksić), which are thorough, multi-decade and partly tabulated, and (ii) scattered synthetic papers using P1/P2/P4 phosphazenes or TMG/BTMG/TBD as catalysts, with yields and ee in SI tables only. Nothing was found that pairs superbase structure with measured rate or ee in machine-readable form.
- **Peptide-catalysis data is small and SI-bound.** The largest single consistent peptide-catalyst dataset located is **200 measurements** (Wennemers/Denmark 2024); the Miller libraries are **24–54 peptides** with models fitted on as few as **9**. For a benchmark oracle this is training-set territory, not evaluation territory.
- **The large organocatalysis datasets that do exist are aminocatalysis (enamine/iminium) literature compilations** — iSynth >22,000, Mannich 4,015, Pictet–Spengler 820 — all carrying the publication bias their own authors document ("a strong skew toward high-enantioselectivity outcomes and a significant underrepresentation of negative results", OrgAIcat/iSynth, 2026).

### 4.3 How much organocatalysis is actually in the ORD — checked directly [measured here]

Queried `https://open-reaction-database.org/api/datasets` on 2026-10-02 and counted:

- **53 datasets, 2,428,291 reactions total.**
- The distribution is dominated by one text-mined set: `uspto-grants`, **1,771,032** reactions; then 409,035 + 40,000 + 30,000 train/test/validation sets from doi:10.1039/C8SC04228D; then 50,688 Pd/Ni/Cu C–N couplings; 47,015 AIChemEco amide couplings; 39,347 Pfizer HTE; 9,632 Chan–Lam; 5,760 Suzuki; 4,312 Ahneman C–N.
- **Exactly one dataset mentions organocatalysis or asymmetric catalysis**: `ord_dataset-c5b00523487a4211a194160edf45e9ab`, "Dataset from `Ultra-high-throughput mapping of the chemical space of asymmetric catalysis enables accelerated reaction discovery`" — "Full factorial of the alpha-asymmetric alkylation of aldehydes with photoredox and organocatalysis [10 bromoacetophenones x 13 photocatalysts x 11 organocatalysts, 1430 total reactions]. Enantioselectivity determined by ion mobility mass spectrometry." **1,430 reactions.**
- One further enantioselective entry exists with **3** reactions (Cu-catalysed enantioselective hydroamination of alkenes).
- Keyword counts across dataset names and descriptions: `proline` 0, `phosphoric` 0, `thiourea` 0, `chiral` 0, `guanidine` 0, `phosphazene` 0, `peptide` 0.
- So: **organocatalysis in the ORD is ~1,430 of 2.43 M reactions, i.e. 0.06%, in a single dataset.** The ORD schema does support ee (the JACS 2021 paper lists "enantiomeric excess by chiral SFC" among analytical outcomes, https://pubs.acs.org/doi/10.1021/jacs.1c09820 ); the data simply is not there. Licence/access: the full database is CC-BY-SA-4.0 on Hugging Face (`open-reaction-database/ord-data`, parquet, 1M<n<10M, last modified 2026-08-30, https://huggingface.co/datasets/open-reaction-database/ord-data ) and on GitHub via Git LFS (https://docs.open-reaction-database.org/en/stable/overview.html ).

---

## SECTION 5 — Low-data and transfer learning in catalysis ML

| approach | dataset size | reported performance | source |
|---|---|---|---|
| **Descriptor model + DFT features, single reaction class** (RF on DFT descriptors of 5 binaphthyl catalyst families, asymmetric hydrogenation) | **368** substrate–catalyst combinations | RMSE **8.4 ± 1.8 %ee**; RF beat CNN and other methods | Singh & Sunoj et al., *PNAS* (2020), https://www.pnas.org/doi/abs/10.1073/pnas.1916392117 |
| **Grid-descriptor + SVM/DNN, in-domain** (Denmark CPA) | 1,075 reactions; 600 train / 475 test, and a 384-train / 691-test split | MAD **0.161 / 0.211 / 0.238 kcal/mol** on substrate / catalyst / substrate-catalyst test sets, q² = 0.748; extrapolating from <80% ee training data to >80% ee: MAD **0.33 kcal/mol** | *Science* 363, eaau5631 (2019) |
| **Pretrain a SMILES language model on a large corpus, then fine-tune (NLP-style transfer learning)** | three reaction classes, >5,000 reactions total: BH coupling >4,000; N,S-acetal **1,027**; asymmetric hydrogenation **368** | yield RMSE **4.89 ± 0.33**; ee RMSE **8.65 ± 0.80** and **8.38 ± 1.40** %ee. Baseline without transfer (TL-m0) **11.83 ± 1.75** on one task; fine-tuning strategy (gradual unfreezing on/off) changed RMSE from 6.02 → 4.89 | Singh & Sunoj, *Digital Discovery* 1, 2022, https://pubs.rsc.org/en/content/articlehtml/2022/dd/d1dd00052g |
| **Meta-learning / few-shot (prototypical networks)** on asymmetric hydrogenation of olefins | pool split 9,032 train / 2,400 test, query sets of 128 | AUPRC **0.9117 ± 0.0026 with only 64 training reactions**, beating RF (0.8369) and GNN (0.8259) trained on the *full* data | *Nat. Commun.* 16, 2025, https://www.nature.com/articles/s41467-025-58854-8 |
| **Hierarchical learning from a large literature database down to a new substrate** | **>12,000** literature asymmetric hydrogenations | "predictive machine learning model using only dozens of enantioselectivity data with the target olefin" | Xu, …, Hong, *Angew. Chem.* 60, 2021, https://onlinelibrary.wiley.com/doi/10.1002/anie.202106880 |
| **Cross-metal transfer + delta-learning** (Pd literature → Ni/Sadphos atroposelective Suzuki) | **21** Ni/Sadphos data points for the target task; base model on Pd data; separately 10–15 "delta" points | Pearson **R = 0.811**, MAE **0.107 kcal/mol**, beating both the small-sample Ni model and direct transfer of the Pd model; 10 vs 15 delta points gave correlated predictions (R = 0.961) | 2023, https://exa.ai/library/publication/zcw0xj0tsq5 |
| **Physical-organic descriptor space + linear SVR** (cobalt/chiral carboxylic acid, indoles) | 108-dimensional descriptor space | Pearson **R = 0.859**, MAE **0.179 kcal/mol** (10-fold CV) | same source |
| **Transition-state-informed descriptors for transferability across ligand/substrate classes** | sparse Ni C(sp³) coupling data; validated on BINOL-CPA additions to imines | models transfer to unseen ligands and reaction partners; explicitly motivated by failure of "simple stereoelectronic parameters" | *Nature*, 2026, https://www.nature.com/articles/s41586-026-10239-7 ; precursor: Reid & Sigman, "Holistic prediction of enantioselectivity in asymmetric catalysis", *Nature* 571, 343 (2019), https://www.nature.com/articles/s41586-019-1384-z |
| **HTE-scale in-house library, honest negative result** | **192** chiral Rh catalysts, up to **3,552** data points | "In-domain modeling for enantioselectivity was unsuccessful, as evidenced by R2 scores not exceeding 0.2 on the test set"; usable models (R² up to 0.859) only emerged on small catalyst subsets (25 train / 4 test) | *Chem. Sci.* 2024, https://pubs.rsc.org/en/content/articlehtml/2024/sc/d4sc03647f |
| **Δ-learning, semiempirical → DFT (Δ² model)** | trained on **~167,000** reactions; GFN2-xTB geometry+energy in, B3LYP-D3/TZVP out | raw GFN2-xTB barriers: MAE **11.9** kcal/mol, RMSE 15.3, R² 0.74. Δ² model: MAE **1.32** kcal/mol, RMSE 2.26, R² 0.99; 99.3% of 36,358 test reactions within 10 kcal/mol. Fine-tuned with **5,142 G4** points: MAE 3.06 vs 4.68 for DFT itself, at xTB cost | *Chem. Sci.* 14, 2023, https://pubs.rsc.org/en/content/articlehtml/2023/sc/d3sc02408c ; code https://github.com/zhaoqy1996/Delta2ML |
| **Horizontal / diagonal transfer learning for barriers from SQM features** | as few as **33 and 39** new data points for the target task | MAE **below 1 kcal/mol**, versus >5 kcal/mol before transfer; near-chemical-accuracy retained down to 122 reactions in the intramolecular Diels–Alder set | *Digital Discovery* 2, 2023, https://pubs.rsc.org/en/content/articlehtml/2023/dd/d3dd00085k |
| **Gaussian processes / Bayesian optimisation over descriptor spaces** | used in closed-loop optimisation rather than static prediction; OrgAIcat + BO took an aldol from **41% to 96% ee in 12 experiments** | — | OrgAIcat, *CCS Chem.* 2026, https://github.com/deepsynthesis/orgaicat ; general method: Shields et al., "Bayesian reaction optimization as a tool for chemical synthesis", *Nature* 590, 89 (2021) |
| **Benchmarking study of representations under realistic small-data conditions** | curated **103** Pd DAAA reactions + **19** external validation reactions; also tested on AH and C–H activation sets | proposes a systematic protocol for small training sets and chemical-space extrapolation | 2026, https://exa.ai/library/publication/9g76hdjwqbh |

What this says about using a transfer-learned catalysis model as a benchmark oracle:

- In the **tens-to-hundreds** regime, the published successes are *within a reaction class* and *with a large same-class donor dataset* (Pd→Ni, 12,000-reaction AHO database → one new olefin, Δ²-ML's 167k reactions). For bifunctional superbase peptidomimetics there is no donor dataset of that kind (§4.2), so the published recipes do not transfer directly.
- The error floor on ee prediction is **~0.1–0.3 kcal/mol MAE in-domain and ~8 %ee RMSE out-of-domain**. An oracle with 8 %ee noise cannot score candidates that differ by a few %ee, which is precisely the regime generative optimisers drive into.
- The *Chem. Sci.* 2024 negative result (R² ≤ 0.2 on 3,552 in-house points) is the most important datum here: large, clean, single-class data is not sufficient for a stable ee model. A benchmark whose oracle is a learned ee model inherits that instability, and the optimiser will find its extrapolation artefacts. The 2026 polyimide study names this failure mode directly — "surrogate labeling failure", where the search "converges toward fictitious high-scoring candidates" (https://exa.ai/library/publication/4vwtsbm586x ) — and its answer is to use a physics-based in-loop score (GFN2-xTB) rather than a learned surrogate.

---

## SECTION 6 — Existing catalyst-design benchmarks

**No generative-design benchmark for catalysts exists that is analogous to GuacaMol or PMO.** Stated plainly, because that is the finding.

What does exist, in descending order of closeness:

- **TARTARUS** (Nigam, Pollice, Aspuru-Guzik et al., NeurIPS 2023 Datasets & Benchmarks, https://proceedings.neurips.cc/paper_files/paper/2023/file/09f8b2469a3d1089a7c60d9ef1983271-Paper-Datasets_and_Benchmarks.pdf , preprint https://arxiv.org/html/2209.12487v4 ) is the only simulation-grounded inverse-design benchmark suite with a reactivity task. Its four families are organic photovoltaics, organic emitters, protein ligands (docking), and **design of chemical reaction substrates**. The reactivity task is "the intramolecular concerted double hydrogen transfer reaction of syn-sesquinorbornenes", with ground states at **GFN0-xTB** and transition states via the **SEAM** force-field method, optimising activation energy and/or reaction energy. It designs **substrates**, not catalysts, and it does not report per-evaluation wall time. It benchmarks REINVENT, SMILES/SELFIES-VAE, MoFlow, SMILES/SELFIES-LSTM-HC, GB-GA and JANUS.
- **GuacaMol** (https://arxiv.org/abs/1811.09621 , 2018; code https://github.com/benevolentAI/guacamol ) and MOSES are drug-oriented; PMO (Gao et al., NeurIPS 2022) adds a call-budget axis (its URL was not re-verified in this search). TARTARUS's own motivation is that these are saturated/oversimplified.
- **One-off catalyst inverse-design studies, not benchmarks**: NaviCatGA + OSCAR genetic optimisation targeting *generality* in asymmetric organocatalysis, trained on a curated **820-reaction** Pictet–Spengler set over a combinatorial space of millions of catalysts (https://b2find.eudat.eu/dataset/66b9055f-437e-5bad-b0c0-814f2d743a72 ); CatDRX, a reaction-conditioned generative model for catalyst design (*Commun. Chem.*, 2025, https://www.nature.com/articles/s42004-025-01732-7 ); "Improving generative inverse design of molecular catalysts in small data regime" (*Mach. Learn.: Sci. Technol.* 6, 025057, 2025, https://backend.orbit.dtu.dk/ws/files/404484481/Cornet_2025_Mach._Learn._Sci._Technol._6_025057.pdf ); and, in materials rather than molecules, generative inverse design of high-entropy catalysts (*Nat. Synth.* 5, 730, 2026, https://www.nature.com/articles/s44160-025-00983-5 ).
- **Descriptor libraries that could supply a benchmark's chemical space but contain no performance labels**: OSCAR (4,000 organocatalysts + 1.57 M combinatorial DHBDs, §4.1) and kraken (1,558 ligands computed, ML predictions for >300,000; https://pubs.acs.org/doi/pdf/10.1021/jacs.1c09718 , 2022).

So the slot the researcher is contemplating — a GuacaMol-style suite whose oracles are catalysis-relevant and whose call budget is thousands per run — is, as far as this search reaches, empty.

---

## GAPS

- **No public dataset pairs superbase (phosphazene/TMG) catalyst structure with measured rate, yield or ee.** Basicity scales exist and are good; catalysis data exists only as per-paper SI tables. Any superbase oracle must therefore be grounded in pK_a/PA (well validated, R² > 0.99 vs experiment) rather than in catalytic performance.
- **No published correlation between global conceptual-DFT indices (ω, N, η, gap) and organocatalyst ee/rate.** Their validation is against small-molecule reactivity (Mayr scales, cycloaddition barriers). Using them as catalysis oracles would be a novel, unvalidated assumption, and they are level-of-theory-dependent (scales defined at B3LYP/6-31G(d); xTB offsets are eV-scale).
- **The dominant cost and the dominant irreproducibility are the same step: conformer handling.** No study found quantifies ETKDG *seed* sensitivity of a *computed QM descriptor* end-to-end; the closest evidence is indirect (kraken's ~1%–>75% descriptor drift from conformer-search randomness; TS-ensemble selection worth >1 kcal/mol).
- **CREST per-molecule costs are published only as log-scale figures or informal HPC guidance** (15–20 min GFN-FF / 5–8 h GFN2 for a 65-atom molecule on 24 cores). A clean, tabulated CREST cost-vs-size benchmark for 30–120 heavy atoms does not appear to exist.
- **GFN-FF gives no orbital energies**, so the cheapest tier (~0.1 s) cannot deliver HOMO/LUMO/ω/N on its own; the practical floor for an electronic descriptor is a GFN2 single point at ~0.3–1.5 s (25–65 heavy atoms, 1 core, measured here).
- **Learned ee oracles are not yet stable enough to be benchmark targets**: 3,552 clean in-house points gave R² ≤ 0.2 in-domain (*Chem. Sci.* 2024), and out-of-domain ee RMSE sits near 8 %ee even with transfer learning.
