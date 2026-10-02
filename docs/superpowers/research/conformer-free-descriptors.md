# Conformer-free electronic descriptors: what is measured, what exists, what is untested

Research note, compiled 2026-10-02. Facts and measured numbers only. Every claim carries a URL and a year.
Where a number is my own arithmetic or inference rather than a published measurement, it is marked **[inference]**.

---

## SECTION 1 — How conformer-sensitive are these descriptors?

### 1.1 The kraken claim, located exactly

The "~1% to >75%" figure is **not in the JACS paper**. It is in the README of the updated Kraken workflow repository, `SigmanGroup/kraken`:

- Run-to-run reproducibility: *"The conformer search produces slightly different conformers each run, so results vary slightly (around 1%) between runs"*
  (https://github.com/SigmanGroup/kraken, accessed 2026)
- Worst case: *"For a small number of properties (typically octant volumes), the difference between the old and new values for a given monophosphine exceeds 75% of the original value."* Cause: *"These differences likely arise from randomness in the conformer search."* This was measured by recomputing DFT properties for **more than 60 monophosphines** with the original vs. updated workflow (xTB 6.2.2/CREST 2.8 vs. xTB 6.4.0/CREST 2.12).
  (https://github.com/SigmanGroup/kraken, accessed 2026)

**Which descriptors are worst: steric, specifically octant volumes.** Octant volumes (`ovbur_max`, `ovtot_min`, `qvbur_max`, …) partition the space around phosphorus into octants; a single torsional flip moves a substituent between octants and changes the value discontinuously. This is the >75% class. The ~1% figure is the typical run-to-run drift across all properties.

A separate, **non-conformational** warning in the same README: with **xTB 6.7.0 or greater**, *"EA/IP descriptors, nucleophilicity"* show substantial variation *"because IPEA-xTB is not used for vertical IP/EA calculations"*. This matters directly: if you build an xTB-based IP/EA/nucleophilicity oracle, the xTB version you pin changes the answer independently of any conformer effect. (https://github.com/SigmanGroup/kraken)

### 1.2 The kraken paper's own statement on electronic vs. steric

Gensch, dos Passos Gomes, Friederich et al., *J. Am. Chem. Soc.* **2022**, 144, 1205–1217 (https://pubs.acs.org/doi/10.1021/jacs.1c09718; open text https://par.nsf.gov/servlets/purl/10352530):

- *"This is particularly relevant for steric properties that vary significantly with conformation, whereas electronic properties are generally less sensitive."*
- Operational consequence, and this is the strongest single piece of evidence in the literature for the robustness ranking: *"All five descriptor variants are then used for properties that were found to vary strongly across conformers (mostly steric properties) whereas **only the Boltzmann-weighted average is used for those properties that are less sensitive to conformation (mostly electronic properties)** to avoid overly redundant descriptors, resulting in a total of 190 descriptors at the DFT level."* (78 per-conformer properties → 190 ligand descriptors.)
- Scale: 1558 ligands, **21,437 DFT conformers, average 13.8 conformers per ligand**. Kraken v2 reports **33,489 unique conformers, average 20.8 per ligand** (https://descriptor-libraries.molssi.org/kraken/).
- kraken's electronic descriptor set already includes exactly the quantities in question — frontier-orbital-derived properties, IP/EA, nucleophilicity, Vmin, Fukui, NBO charges, dipole (SI section 2.2, table of contents reproduced at https://pubs.acs.org/doi/10.1021/jacs.1c09718).

### 1.3 Measured eV spread for a flexible catalyst: Rh–bisphosphine complexes

Baidun, Kalikadien, Lefort, Pidko, *J. Phys. Chem. C* **2024**, 128, 7987–7998 (doi 10.1021/acs.jpcc.4c01631; full text https://repository.tudelft.nl/file/File_4814a578-1c9d-47ce-b169-322ba8d9287e).

System: 44 L–Rh–substrate and 11 L–Rh–NBD complexes, bidentate phosphine ligands; CREST 2.12 / xtb 6.6.1 at GFN2-xTB//GFN-FF, ensembles taken within 6 kcal/mol; descriptors at PBE0-D3(BJ)/def2-SVPP, Boltzmann-weighted with standard deviations reported.

Measured **maximum standard deviation within a conformer ensemble**:

| Descriptor | Max std within ensemble | Unit |
|---|---|---|
| HOMO–LUMO gap | **0.2** | eV |
| NBO charge on Rh | 0.03 | au |
| Buried volume at Rh | 0.008 | (fraction) |
| Buried volume at donors | 0.003 | (fraction) |
| (third steric descriptor) | 0.001 | — |

Quoted: *"With a maximum standard deviation for the HOMO−LUMO gap of 0.2 eV within a conformer ensemble, a difference exceeding 0.5 eV across different substrate coordination modes is significant"*. Concrete example: L1-Rh-S major gap **3.31 eV** vs. L1-Rh-S minor gap **3.86 eV** — a 0.55 eV swing from binding mode, against 0.2 eV of conformer noise. The paper concludes electronic descriptors *"exhibited more substantial variations"* than the buried volumes it tested — note this is the opposite ordering to kraken's, because the steric descriptors here (buried volume) are the conformationally *stable* ones, unlike octant volumes/Sterimol.

**Reading of the two together:** the robust/fragile split is not "electronic vs. steric" but "scalar integral over the whole molecule vs. direction-resolved local quantity". Buried volume and HOMO–LUMO gap are both reasonably robust; octant volumes, Sterimol B5/L and the dipole vector are not. **[inference]**

### 1.4 Large-scale measurement on flexible drug-sized organics

Hamakawa & Miyao, *J. Chem. Inf. Model.* **2025**, 65, 3388–3404 (doi 10.1021/acs.jcim.5c00018; https://pmc.ncbi.nlm.nih.gov/articles/PMC12004525/).

Dataset built specifically to isolate conformational effects: 97,695 CHNO molecules (MW < 500) selected as the **top 100,000 by rotatable-bond count** from PubChemQC PM6; **644,642 conformers** from OMEGA (max 40/molecule, 2.0 Å RMSD threshold), each PM6-optimised in MOPAC. Mean of the per-compound **maximum pairwise RMSD = 3.82 Å** (mean pairwise 2.98 Å) — genuinely diverse ensembles.

Six properties profiled per compound by the standard deviation across its conformers, scaled by the interquartile range of that property (Figure 2):

- *"The dipole moment was highly conformation-dependent with significant variations due to conformational coordinates. On the other hand, the energy and enthalpy showed low conformational dependency."*
- The per-property numeric std values exist only inside Figure 2 (violin plots + histograms with "mean-range" annotations); they are not tabulated in the text, so I cannot give the eV/Debye numbers from this paper.

What *is* tabulated is the downstream consequence, which is the decision-relevant number. A purely 2D model (ECFP4 count vector) on this deliberately-flexible set:

| Target | R² from 2D ECFP4 alone |
|---|---|
| HOMO–LUMO gap | best of the set (text: *"particularly HOMO, HOMO–LUMO gap, and LUMO"*) |
| LUMO | high |
| HOMO | **0.867** |
| Dipole moment | **0.516** |

Quoted: *"The structural-formula representation of ECFP4 sufficiently predicted the QC properties of compounds with high conformational diversity, particularly HOMO, HOMO–LUMO gap, and LUMO. … The R² value for the dipole moment was 0.516, and that for HOMO was 0.867. This may be due to the conformational dependence of the dipole moment and the HOMO energy."*

This is the cleanest answer to the question asked: on ~98k highly flexible molecules, **a model with no 3D information at all reproduces the gap and LUMO well and the dipole badly**.

### 1.5 Independent corroboration that conformer geometry adds nothing for electronic targets

Cheng, Jin, Zhang, *"When Does Conformer Geometry Help?"*, arXiv:2606.08825 (**2026 preprint, not peer-reviewed — treat as weak evidence**; https://arxiv.org/html/2606.08825). ~1,000 experiments across MoleculeNet, QM9 and MARCEL, 50 ETKDG conformers per molecule. Finding: conformer-ensemble statistics give significant RMSE reductions on solvation-dependent properties (ESOL −11.0%, p<10⁻⁹; FreeSolv −13.5%) *"while providing no benefit for electronic or steric tasks"*; their property taxonomy states *"electronic properties show no improvement"* from conformer geometry.

### 1.6 Ranking, most to least conformer-robust

| Rank | Descriptor | Evidence |
|---|---|---|
| 1 | **HOMO–LUMO gap** | kraken: electronic → Boltzmann-only (JACS 2022); Baidun: max std **0.2 eV** within ensemble vs 0.5 eV chemistry signal (JPCC 2024); Hamakawa: best 2D-predictable property (JCIM 2025) |
| 2 | **LUMO** | Hamakawa: among best 2D-predictable (JCIM 2025) |
| 3 | **HOMO** | Hamakawa: R²=0.867 from 2D, explicitly attributed to *"conformational dependence of … the HOMO energy"* (JCIM 2025) |
| 4 | **IP, EA** | No direct ensemble-spread measurement found. Proxy: AIMNet-NSE reaches 0.10 eV RMSE on optimised geometries but **0.15 eV on off-equilibrium geometries** — a 50% error increase purely from geometry perturbation (Nat. Commun. 2021, https://www.nature.com/articles/s41467-021-24904-0) |
| 5 | **ω (electrophilicity), N (nucleophilicity)** | No measurement found. **[inference]**: ω = χ²/2η with χ=(IP+EA)/2 and η=IP−EA; η is a difference of two similar large numbers, so the relative conformer noise in η is larger than in IP or EA, and ω inherits it amplified. Treat ω and N as *less* robust than IP/EA until measured |
| 6 | **Dipole moment** | Worst in every source: Hamakawa *"highly conformation-dependent with significant variations"*, R²=0.516 from 2D (JCIM 2025). It is a vector sum over bond moments, so torsions reorient contributions that cancel or add |

**Bottom line for the question asked: yes — the gap is the robust one and the dipole is not.** The gap's within-ensemble σ on a real flexible catalyst is 0.2 eV, which is smaller than the 0.5 eV differences that carry chemical meaning in the same system. ω and N are the unmeasured risk in the user's descriptor list.

---

## SECTION 2 — Boltzmann averaging and cheap approximations to it

### 2.1 The four conventions in use

1. **Single lowest-energy conformer.** Still the default in catalysis feature libraries. The Doyle-group bisphosphine library's earlier protocol: force-field ensemble → RMSD-prune to five → *"From these five conformers, the lowest energy conformer defined through DFT was used for downstream featurization and modelling efforts"* (*Chem. Sci.* **2025**, DOI 10.1039/D5SC04691B, https://pubs.rsc.org/en/content/articlehtml/2025/sc/d5sc04691b).
2. **Boltzmann-weighted ensemble average** at 298 K. kraken's choice for all electronic properties (JACS 2022).
3. **min / max / delta as separate features** — the kraken and Sterimol convention. kraken keeps five variants for conformationally-variable properties: Boltzmann, min, max, delta (max−min), and the property of the conformer with smallest V_bur; at the xTB level it also stores `std` (https://descriptor-libraries.molssi.org/kraken/library_details).
4. **Range-as-uncertainty.** wSterimol (Brethomé, Fletcher, Paton, *ACS Catal.* **2019**, 9, 2313; https://doi.org/10.1021/acscatal.8b04043): Boltzmann-weighted Sterimol parameters, plus *"the range of values from conformers within 3 kcal/mol of the most stable structure provides a visual way to capture a possible source of uncertainty arising in the resulting models"*.

### 2.2 Concrete conformer counts

| Source | Count / window | Context |
|---|---|---|
| kraken, JACS 2022 | up to **20** conformers within **3 kcal/mol** of the minimum, plus every conformer that minimises or maximises any of 19 steric descriptors; RMSD pruning if >20 in window. Realised: **13.8 DFT conformers/ligand** (v1), **20.8** (v2) | P(III) ligands, 23.7 heavy atoms avg |
| Bisphosphine library, Chem. Sci. 2025 | **10** conformers, chosen at equidistant GFN2-xTB energies (min, max, evenly spaced between). Explicitly: *"The use of five conformers was also investigated but gave inadequate representation – especially with larger conformer ensembles"* | Pd(II) bisphosphine complexes |
| Baidun, JPCC 2024 | all CREST conformers within **6 kcal/mol** | Rh complexes |
| Laplaza et al., *Chem. Sci.* **2022**, 13, 6858 (https://pubs.rsc.org/en/content/articlepdf/2022/sc/d2sc01714h) | *"considering about 20 transition state conformations per catalyst"*; *"Given the availability of these data, Boltzmann weighting of all conformers is recommended"* | flexible chiral Cp Rh(III) catalysts |
| MARCEL benchmark, arXiv:2310.00115 (2023/24) | set encoders capped at **20** conformers/molecule: *"we cap the number of encoded conformers per molecule to a maximum of 20, which empirically improves training stability"* | organocatalysts + drugs |
| Conformer generation for SBDD, *JCIM* **2024**, 64, 1 (https://pubs.acs.org/doi/full/10.1021/acs.jcim.3c01245) | *"a reduced conformational ensemble of less than 25 conformers … is likely sufficient for most structure-based screening tasks"* | docking/pharmacophore, not descriptors |

### 2.3 Does more averaging actually help? One measured A/B test

*Chem. Sci.* **2025**, 10.1039/D5SC04691B built both feature sets for the same 12 ligands and the same Hayashi–Heck regioselectivity model (MLR, 22 train / 30 test ligands):

| Feature set | Train R² | Train MAE | Test R² | Test MAE |
|---|---|---|---|---|
| Conformer-weighted (lowest-E + min + max + Boltzmann + mean, 2088 features, pruned to <1300) | **0.76** | **0.30** kcal/mol | **0.57** | **0.51** kcal/mol |
| Lowest-energy conformer only | 0.46 | 0.48 kcal/mol | 0.29 | 0.69 kcal/mol |

Their pruning rule is directly reusable: *"This set was reduced to features that were conformationally dependent … For features that showed low variance across an ensemble, only the Boltzmann-weighted average value was retained"* — i.e. measure the variance first, then decide per-descriptor whether min/max/delta earn their place.

### 2.4 A warning about the word "Boltzmann"

Ryde, *"How Many Conformations Need To Be Sampled To Obtain Converged QM/MM Energies? The Curse of Exponential Averaging"* (2017, https://exa.ai/library/publication/2453jq5h8yl; published *J. Chem. Theory Comput.*) shows that **exponential** averaging of *energies* is badly conditioned: with σ = 10–17 kJ/mol between snapshots, 20 samples give errors of 8–33 kJ/mol, and ~4,200 to ~10⁹ samples are needed for 4 kJ/mol at 95% confidence. The cumulant (Gaussian) approximation cuts this to 230–1700.

This does **not** transfer directly to descriptor averaging: a Boltzmann-weighted *descriptor* average is a weighted arithmetic mean of a bounded quantity, whose conditioning is far better than an exponential average of energies. **[inference]** But the related failure does transfer: the *weights* come from an exponential of energies, so a conformer-energy error of 1 kcal/mol redistributes weight by a factor e^1.7 ≈ 5.5 at 298 K, and the weighted average of a 0.2 eV-spread descriptor can shift by a meaningful fraction of that spread. For enzyme barriers the analogous count is ~20 conformations (Li et al., *Int. J. Mol. Sci.* **2016**, https://pmc.ncbi.nlm.nih.gov/articles/PMC5000767/ — Boltzmann-weighted barrier converged after 13 and 9 conformers, disproportionate-effect criterion after 18 and 17).

**No study was found that measures the conformer count at which a Boltzmann-averaged HOMO, LUMO, gap, IP, EA, ω or N converges.** The counts above (10, 13.8, 20, 20.8) are conventions chosen for cost and for *steric* descriptor coverage, not convergence measurements on electronic descriptors.

---

## SECTION 3 — Neural surrogates, with element coverage checked

**Headline finding before the table: of the models asked about, only two output frontier-orbital or Parr/Domingo quantities at all.** AIMNet2, ANI-2x, MACE-OFF23, SpookyNet, SO3LR, AIQM1 and the OMol25/UMA family are *interatomic potentials*: energy, forces, sometimes charges and dipoles. They do not emit HOMO, LUMO, gap, IP, EA, ω or N. The two that do are **AIMNet-NSE** (IP, EA, χ, η, ω, Fukui — directly) and **OrbNet-Equi** (HOMO, LUMO, gap).

| Model | Predicts | Training data | **Elements — P?** | Max size validated | Accuracy | Speed | 3D in? | Licence / weights |
|---|---|---|---|---|---|---|---|---|
| **AIMNet-NSE** (Zubatyuk et al., *Nat. Commun.* 2021, 12, 4870) | energy, spin-polarised charges, and derived **`ip, ea, f_el, f_nuc, f_rad, chi, eta, omega, omega_el, omega_nuc, omega_rad`** in eV | ~200k UNICHEM molecules ≤16 heavy atoms; **Ions-12 = 6.44M structures ≤12 heavy atoms** used for train/val; PBE0/ma-def2-SVP | **H, C, N, O, F, Si, P, S, Cl — P YES. No Br, no I.** | **Tested on Ions-16 (13–16 heavy atoms) and ChEMBL-20 (≤20 heavy atoms). Nothing above ~20.** | IP/EA RMSE **~0.10 eV** optimised geometries, **~0.15 eV** off-equilibrium; energies 3–4 kcal/mol on 13–20 heavy atoms | not benchmarked in paper; NN inference | **Yes**, coords + charge | Public weights, 5-model ensemble, github.com/isayevlab/aimnetnse |
| **AIMNet2** (Anstine, Zubatyuk, Isayev, *Chem. Sci.* 2025, 16, 10228) | energy, forces, Hirshfeld charges. **No orbital energies.** | 2×10⁷ ωB97M-D3/def2-TZVPP, distilled from 1.2×10⁸ B97-3c; source molecules **<20 heavy atoms** from PubChem/ChEMBL + ANI-1x/2x + OrbNet sets | **H, B, C, N, O, F, Si, P, S, Cl, As, Se, Br, I (14) — P YES** | 2025 extension dataset up to **193 atoms/system** (colabfit/AIMNet2 on HF); charged + zwitterionic supported | *"outperforms semi-empirical GFN2-xTB and is on par with reference DFT"* for interaction energies, conformer search, torsions | not stated numerically | **Yes** | Public (isayevlab) |
| **AIMNet2-NSE** (Angew. Chem. 2025, doi 10.1002/anie.202516763) | open-shell energies/forces, spin states | unrestricted ωB97M-D3/def2-TZVPP | **14 elements incl. P** | gas-phase molecular only; no transition metals | — | — | Yes | isayevlab.github.io/aimnetcentral |
| **ANI-2x** (Devereux et al., *JCTC* 2020, 16, 4192) | energy, forces | 9,651,712 conformers at ωB97X/6-31G(d) | **H, C, N, O, F, Cl, S — NO PHOSPHORUS** ⚠️ | small organics (COMP6) | sub-chemical accuracy on most of COMP6 | *"~10⁶ factor speedup"* vs DFT | Yes | Open (TorchANI) |
| **ANI-1x / ANI-1ccx** | energy, forces | 5M / 0.5M | **H, C, N, O only — NO P, NO S, NO HALOGENS** ⚠️ | ≤8 heavy atoms training | ANI-1ccx at CCSD(T)*/CBS | — | Yes | Open |
| **MACE-OFF23** (Kovács et al., arXiv:2312.15211, 2023/25) | energy, forces | SPICE subset (~647k PubChem + 263k DES dimers + dipeptides) + QMugs 50–90-atom molecules; ωB97M-D3(BJ)/def2-TZVPPD | **H, C, N, O, F, P, S, Cl, Br, I (10) — P YES** | QMugs augmentation covers **50–90 atoms**; solvated protein demonstrated | — | — | Yes | ⚠️ **Academic Software License — academic use only, not commercial** (github.com/ACEsuit/mace-off). Note: 808 configurations removed for >2 eV/Å force error, *"many of which contained heavy elements, in particular phosphorus and iodine"* |
| **OrbNet / OrbNet-Equi** (Qiao et al., *PNAS* 2022, 119, e2205221119) | energies, forces, densities, **HOMO, LUMO, HOMO–LUMO gap** (needs the energy-weighted density-matrix features D^β_h, D^β_p), dipoles; zero-shot IPs | SDC21 = **235,834 geometries**, ωB97X-D3/def2-TZVP, delta-learned on GFN-xTB | **C, O, N, F, S, Cl, Br, I, P, Si, B, Na, K, Li, Ca, Mg — P YES**. Entos Sierra product page lists C,H,B,O,N,F,P,S,Cl,Si,Br,I, neutral closed-shell only | drug-like + biofragments; conformer ranking on ~700-molecule Hutchison set | On QM9 FMO energies, beats SphereNet trained on 110k samples; torsion barrier MAE 0.173 kcal/mol; S66x10 binding 0.35 kcal/mol | *"neural network inference time … on par with the GFN-xTB QM featurizer"*, overall **100–1000× faster than composite DFT** | **Yes — plus a full GFN-xTB SCF** | ⚠️ *"The software used for computing input features and gradients is proprietary to Entos, Inc."* Training set + NN code on Zenodo 6568518, featurizer is not |
| **AIQM1** (Zheng et al., *Nat. Commun.* 2021, 12, 7022) | energies, geometries, heats of formation | ANI-1x (5M) + ANI-1ccx (0.5M) | **H, C, N, O ONLY — NO P, NO S, NO HALOGENS** ⚠️ MLatom docs: *"This method is currently limited to compounds only containing H, C, N, and O elements."* | C60 geometry demonstrated | CHNO ΔHf MAD 0.9 kcal/mol; GMTKN55 WTMAD-2 8.94 (benchmarks.rowansci.com/methods/aiqm1) | semiempirical speed | Yes | MLatom; local install |
| **SpookyNet** (Unke et al., *Nat. Commun.* 2021, 12, 7273) | energy, forces; takes **total charge and spin as explicit inputs** | demonstrated on QM7-X and others | **QM7-X = C, N, O, S, Cl, ≤7 heavy atoms — NO PHOSPHORUS** ⚠️ (https://www.nature.com/articles/s41597-021-00812-2). No general P-capable pretrained model released | QM7-X ≤7 heavy atoms; transfer to larger shown by RMSD of optimised geometries | — | — | Yes | Reference implementation github.com/OUnke/SpookyNet; **architecture, not a ready oracle** |
| **SO3LR** (Kabylda et al., *JACS* 2025, doi 10.1021/jacs.5c09558) | energy, forces, **dipole moments**, Hirshfeld ratios. **No orbital energies.** | **4M** structures at PBE0+MBD: 2.7M GEMS protein fragments + 1M QM7-X + 60k AQM + 33k SPICE dipeptides + 15k DES dimers | **H, C, N, O, F, P, S, Cl (8) — P YES** | **~200,000 atoms** on a single GPU; crambin, glycoprotein, lipid bilayer in explicit solvent | TorsionNet500 MAE **1.03 kcal/mol** at PBE0+MBD | **~3 µs/atom/step on one H100** | Yes | Open, github.com/general-molecular-simulations/so3lr |
| **OMol25 + UMA / eSEN** (Levine et al., arXiv:2505.08762, 2025) | energy, forces. **No released model predicts orbital energies**, although the dataset's raw outputs contain them | **140,641,161** DFT calcs at ωB97M-V/def2-TZVPD; ~83M unique systems | **83 elements — P YES**, plus transition metals | **systems up to 350 atoms**, charges −10..+10, spin 1..11 | UMA-M-1.1 energy MAE 1.38 kcal/mol averaged over splits. OMol25 evaluations include an **ionization-energy task** | — | Yes | Dataset **CC-BY-4.0**; model checkpoints under gated **FAIR Chemistry License** |

### 3.1 Models that exclude phosphorus — flagged

- **ANI-2x**: H, C, N, O, F, Cl, S. **No P.**
- **ANI-1x, ANI-1ccx, ANI-1xnr**: H, C, N, O only. **No P, no S, no halogens.**
- **AIQM1**: H, C, N, O only. **No P.**
- **SpookyNet as trained on QM7-X**: C, N, O, S, Cl. **No P.**

For peptidomimetics carrying a phosphazene P=N centre, these four are out on element coverage alone, before any other consideration.

### 3.2 The relevant dataset gap — xTB is itself a poor proxy for DFT frontier orbitals

QMugs (Isert et al., *Sci. Data* **2022**, 9, 273, https://www.nature.com/articles/s41597-022-01390-7) is the only large dataset with **both** GFN2-xTB and DFT frontier orbitals on drug-sized molecules: 665,911 molecules, ~2.0M conformers (3 per molecule, metadynamics-generated then k-means clustered), **mean 30.6 / max 100 heavy atoms**, 10 atom types, ωB97X-D/def2-SVP.

Measured Pearson correlation between the two levels of theory over the whole set:

| Property | GFN2-xTB vs ωB97X-D/def2-SVP, PCC |
|---|---|
| Formation energy | 0.998 |
| Rotational constants | 0.999 |
| **Dipole moment** | 0.969 |
| **LUMO** | 0.924 |
| **Gap** | 0.830 |
| **HOMO** | **0.769** |

This is an important constraint on the whole plan: a GFN2-xTB single point reproduces the DFT **HOMO** with only r = 0.77 across drug-sized chemistry. If the generative loop's objective is defined at xTB level, that is self-consistent; if it is a proxy for DFT, the method error already exceeds the 0.2 eV conformer noise measured in §1.3. **[inference on the comparison]**

---

## SECTION 4 — Predicting electronic properties from 2D structure alone

### 4.1 Yes, such models exist, and they are competitive

**MARCEL / Drugs-75K** (Zhu, Hwang, Adams, Liu et al., arXiv:2310.00115, 2023; NeurIPS Datasets & Benchmarks). 75,099 molecules with **≥5 rotatable bonds**, 558,002 conformers, **mean 30.56 heavy atoms, mean 7.53 rotatable bonds**, elements **H, C, N, O, F, Si, P, S, Cl** (P included). Targets: **Boltzmann-averaged ionisation potential, electron affinity and electronegativity**, labels computed with Auto3D + AIMNet-NSE over the ensembles.

MAE (eV) for predicting the Boltzmann-averaged target:

| Representation | Model | IP | EA | χ |
|---|---|---|---|---|
| 1D fingerprint | Random forest | 0.4987 | 0.4747 | 0.2732 |
| 1D SMILES | LSTM | 0.4788 | 0.4648 | 0.2505 |
| **2D graph** | **GIN** | **0.4354** | **0.4169** | **0.2260** |
| **2D graph** | **GraphGPS** | **0.4351** | **0.4085** | **0.2212** |
| 3D, lowest-energy conformer | GemNet | 0.4069 | 0.3922 | 0.1970 |
| 3D, lowest-energy conformer | LEFTNet | 0.4174 | 0.3964 | 0.2083 |
| 3D + explicit ensemble encoder | GemNet + DeepSets | **0.4066** | **0.3910** | 0.2027 |

**The 2D graph model is within 7% of the best 3D ensemble model on IP and EA, and the ensemble encoder buys essentially nothing over a single lowest-energy conformer on this dataset** (0.4069 → 0.4066 eV for IP). The paper notes that for Drugs-75K specifically, *"this does not always extend to larger datasets"* — the ensemble gains concentrated on the small Kraken steric targets.

**Hamakawa & Miyao (JCIM 2025)**: plain ECFP4 count vectors on 97,695 highly flexible CHNO molecules reach **R² = 0.867 for HOMO**, higher still for gap and LUMO, and only 0.516 for dipole (§1.4).

**kraken's own ML layer** (JACS 2022) predicts Boltzmann-averaged DFT descriptors — including the frontier-orbital-derived ones, IP, EA and nucleophilicity — from 2D structure alone (Bag-of-Substituents / fingerprints / Gaussian processes / GCN, stacked into per-descriptor metamodels), then applies them to 331,776 ligands (VL2) and on-demand to ~191M (VL3). Reported spread: **58 properties with R²_test ≥ 0.80, 45 properties with R² < 0.50**.

### 4.2 Size extrapolation: stated plainly

**The extrapolation from QM9-scale (~9 heavy atoms) to 60–80 heavy atoms has not been tested for frontier-orbital prediction.** Nothing found tests it. What exists:

- **QMALL** (NeurIPS 2022 AI4Science, https://neurips.cc/virtual/2022/57010): train on QM9, test on Alchemy, which has *larger* molecules — but Alchemy tops out around 14 heavy atoms, not 60–80. Reports *"overall performance drop, model ranking inconsistency"* and explicitly frames OOD size extrapolation as unsolved.
- **QM9 → PC9** (Glavatskikh et al., *J. Cheminform.* **2019**, 11, 69, https://link.springer.com/article/10.1186/s13321-019-0391-2): same size and element limits, different chemistry. SchNet in-domain HOMO/LUMO MAE **0.04 / 0.03 eV**; on PC9 molecules absent from QM9 (subset B), **0.33 / 0.27 eV** — an **8–9× error increase from chemical-space shift alone, at constant molecule size**. 1.88% of compounds had HOMO error >1.5 eV; outliers were non-singlet multiplicities and specific functional groups. Models trained on the more chemically diverse PC9 generalised far better.
- **QM7 → QM9** (Collins et al., cited in the same paper): U₀ MAE rose 3.4 → 106 kcal/mol, while **HOMO/LUMO accuracy held at ~0.15 eV**. Frontier orbitals are intensive and transfer better than extensive energies.
- **The 2026 review of ML for frontier-orbital energetics** (Mantilla Santa Cruz, Faina, de Souza Pereira, https://exa.ai/library/publication/3qlnr4xy2gf; pre-proof) surveyed 59 studies: pooled median MAE 25.8 meV for equivariant architectures on QM9; **only 47.5% (28/59) report any out-of-distribution evaluation**, with no shared OOD benchmark; only **13.6% (8/59)** report any uncertainty quantification. Its own stated limitation: QM9 is C/N/O/F, ≤9 heavy atoms, neutral closed-shell, *"roughly an order of magnitude smaller than typical drug-like molecules … so we cannot assume the accuracy figures in this review transfer to larger molecules without independent validation"*, and *"Several included studies that do report OOD evaluation … show measurable accuracy loss when extrapolating to larger or differently distributed molecules."* The review also notes that the handful of studies that trained on conformer ensembles rather than the single QM9 geometry found *"conformer handling is itself a source of performance variance independent of architecture, but the corpus is too small to quantify this systematically."*

**Plain statement:** QM9-trained HOMO/LUMO models have never been validated at 60–80 heavy atoms. The one dataset that *could* support such a model at the right size, with phosphorus, is **QMugs** (max 100 heavy atoms, mean 30.6, HOMO/LUMO/gap at both GFN2-xTB and ωB97X-D/def2-SVP) — but QMugs has only 3 conformers per molecule and no ensemble-averaged labels, and no published 2D model trained on it was found in this search for the full frontier-orbital set.

---

## SECTION 5 — Surrogates trained on ensemble-averaged targets

Yes, these exist, and this is the only class that removes the conformer problem rather than moving it.

1. **MARCEL / Drugs-75K and Kraken tasks** (arXiv:2310.00115, 2023/24). The label is defined as ⟨y⟩_kB = Σ p_i y_i over the ensemble, and the model may see only the 2D graph. Quoted: *"The tasks are to predict the Boltzmann-averaged value of each property across the conformer ensemble … the goal is to predict ⟨y⟩_kB from the molecular graph G, a single conformer C_i ∈ C, or the set C."* Targets: Boltzmann-averaged **IP, EA, electronegativity** (Drugs-75K, 75,099 molecules, P included) and Boltzmann-averaged **Sterimol B5, L, buried B5, buried L** (Kraken, 1,552 ligands). Best 2D numbers in §4.1. The labels themselves come from AIMNet-NSE over Auto3D ensembles, so their absolute accuracy is bounded by AIMNet-NSE's ~0.1–0.15 eV.

2. **kraken's descriptor metamodels** (JACS 2022). Trained directly on Boltzmann-averaged DFT descriptors as a function of 2D structure; 58 of the properties reach R²_test ≥ 0.80. This is the closest existing artifact to what is wanted — including Boltzmann-averaged IP, EA, nucleophilicity and frontier-orbital descriptors — but the chemistry is monodentate P(III) ligands only (phosphines, phosphites, phosphoramidites, phosphinamines), not amide/urea/sulfonamide peptidomimetics, and the Bag-of-Substituents component *"is inherently incapable of extrapolating to unseen substituents"*.

3. **Boltzmann GNN** (arXiv:2312.13110, 2023/24). Pre-trains a graph transformer by maximising the conditional marginal likelihood of a conformer-generating diffusion model, so the graph latent encodes the Boltzmann distribution. It is a *pre-training* scheme that improves downstream property prediction, not a direct ensemble-property regressor.

4. **Ensemble-derived descriptors (not learned)**: wSterimol's Boltzmann-weighted Sterimol (*ACS Catal.* 2019) and Denmark's Average Steric Occupancy, which *"simplif[ies] the conformer population information into a location-specific numerical form"* (*Science* **2019**, 363, eaau5631, https://pmc.ncbi.nlm.nih.gov/articles/PMC6417887/). Both are steric; neither has an electronic counterpart in common use.

**Nothing was found that is trained to predict a Boltzmann-averaged ω or N.** Drugs-75K's χ is the closest conceptual neighbour (χ = (IP+EA)/2, and ω = χ²/2η follows from the same two numbers), so a Boltzmann-averaged ω label could be derived from the same AIMNet-NSE pipeline without new QM. **[inference]**

---

## SECTION 6 — What catalysis ML actually does

| Paper | System | Conformer protocol, as disclosed | Disclosed? |
|---|---|---|---|
| **kraken**, *JACS* 2022, 144, 1205 | 1,558 monodentate P(III) ligands | CREST at GFN2-xTB (toluene GBSA) on free ligand **and** on the Ni(CO)₃ complex; select every conformer that min/maximises any of 19 steric descriptors, plus ≤20 within 3 kcal/mol with RMSD pruning; all re-optimised at DFT as free ligands; **Boltzmann average only for electronic properties**, five variants for steric | **Fully**, SI §1.2 and the public workflow repo |
| **Kraken v2 / SigmanGroup**, repo README | same | xTB 6.4.0 / CREST 2.12; validation folder compares ~28–30 monophosphines against the original workflow; discloses the ~1% run-to-run drift and the >75% octant-volume outliers | **Fully, including the failure modes** — unusually honest |
| **Doyle-group bisphosphine library**, *Chem. Sci.* 2025, 10.1039/D5SC04691B | 12 bisphosphine PdCl₂ complexes, extended library | **Old practice**: force field → prune to 5 by RMSD → DFT-lowest-energy conformer only. **New practice**: 10 conformers at equidistant GFN2-xTB energies → 2088 conformer-weighted features (lowest-E, min, max, Boltzmann@298 K, arithmetic mean) → prune to <1300 by keeping Boltzmann-only for low-variance features. Measured gain in §2.3 | **Fully**, and it quantifies the cost of the old protocol |
| **Baidun et al.**, *J. Phys. Chem. C* 2024, 128, 7987 | Rh–bisphosphine hydrogenation | CREST 2.12 / xtb 6.6.1, GFN2-xTB//GFN-FF, `-noreftopo`, aromatic rings on the chiral axis fixed; ensemble within 6 kcal/mol; PBE0-D3(BJ)/def2-SVPP; **Boltzmann-weighted with standard deviations reported alongside every descriptor** (paper states T = 289 K) | **Fully**, including σ error bars — the main point of the paper |
| **Laplaza et al.**, *Chem. Sci.* 2022, 13, 6858 | chiral Cp Rh(III) C–H activation, highly flexible | ~20 TS conformers per catalyst from graph-based Molassembler, four stereochemical pathways (DR/DS/UR/US); shows that restricting to the chemically "obvious" DR/DS pathways gives significant errors because UR/US conformers are often lower in energy. Recommends Boltzmann weighting of *all* conformers | **Fully**, and it is a counterexample to lowest-energy-conformer practice |
| **Denmark et al.**, *Science* 2019, 363, eaau5631 | chiral phosphoric acid catalysts | Conformer population folded into Average Steric Occupancy grid descriptors; electronic descriptors taken from a single structure | Steric protocol disclosed; electronic-conformer handling not emphasised |
| 2026 review of ML for frontier-orbital energetics (pre-proof) | 59 ML studies | *"A small number of included studies address this directly by training on conformer ensembles … rather than the single QM9-DFT geometry … the corpus is too small to quantify this systematically."* Only 47.5% report any OOD evaluation | **Mostly not disclosed** across the surveyed field |

**Pattern.** In organometallic/organocatalysis descriptor work from the Sigman, Paton, Doyle, Pidko and Corminboeuf orbits, the conformer protocol *is* disclosed, usually in detail, and the field has converged on: CREST/GFN2-xTB ensemble → 3–6 kcal/mol window → 10–20 conformers → Boltzmann average for electronic descriptors, Boltzmann + min + max + delta for steric. In the broader ML-for-QM-properties literature the protocol is frequently absent.

---

## GAPS

- **No published measurement exists of the conformer-ensemble spread of ω (electrophilicity) or N (nucleophilicity)** for any molecule class. They are the two descriptors in the user's list with zero evidence, and the η = IP − EA denominator makes them the most likely to be noisy.
- **No study measures how many conformers are needed before a Boltzmann-averaged *electronic* descriptor converges.** The 10/14/20/21 conformer counts in the literature are cost conventions and steric-coverage heuristics, not convergence measurements on HOMO/LUMO/IP/EA.
- **The QM9 → 60–80-heavy-atom extrapolation is untested.** The largest tested jumps are QM9→Alchemy (~14 heavy atoms) and QM9→PC9 (same size, different chemistry, 8–9× HOMO/LUMO error increase). QMugs (max 100 heavy atoms, with P) could close this but no frontier-orbital model trained on it was found.
- **AIMNet-NSE — the only open model that emits IP, EA, χ, η and ω directly and supports phosphorus — was trained on ≤12 heavy atoms and validated only to ~20.** Using it at 60–80 heavy atoms with phosphazene P=N and sulfonamide motifs is a 3–4× size extrapolation into chemistry not represented in Ions-12. No published validation covers this.
- **Phosphazene/superbase chemistry is absent from every training set checked.** Kraken is P(III); AIMNet-NSE/AIMNet2/MACE-OFF23 draw from PubChem/ChEMBL/SPICE; SPICE and MACE-OFF23 specifically dropped high-error phosphorus configurations. Pentavalent P=N is at best thinly sampled anywhere.
- **The GFN2-xTB reference itself is weak for the target quantity**: r = 0.769 against ωB97X-D/def2-SVP HOMO across 665k drug-sized molecules (QMugs 2022), which is a larger discrepancy than the 0.2 eV conformer spread being worried about — and xTB ≥6.7.0 changes vertical IP/EA and nucleophilicity outright by dropping IPEA-xTB.
