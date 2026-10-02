# Peptide similarity measures and similarity-based holdouts: published facts

Scope: facts, measured numbers, and stated limitations, with URL + year per claim. Quotes kept to short phrases.
Compiled 2026-10-02.

---

## SECTION 1 — Is ECFP/Morgan-Tanimoto documented as poor for peptides?

### 1.1 The MAP4 paper (Capecchi, Probst & Reymond, J Cheminform 12:43, 2020)

Source: https://jcheminf.biomedcentral.com/articles/10.1186/s13321-020-00445-4 ; full text https://pmc.ncbi.nlm.nih.gov/articles/PMC7291580/ (2020). Code https://github.com/reymond-group/map4 (MIT licence).

Explicit statements in the paper:

- Framing claim (abstract): substructure fingerprints perform best for small molecules, "atom-pair fingerprints are preferable for large molecules such as peptides", and no then-available fingerprint was good at both. (2020)
- On ECFP4/MHFP6 specifically: both "have a poor perception of the global features of molecules such as size and shape" and "fail at perceiving structural differences that may be important in larger molecules", with three named failures: regioisomers in extended ring systems (2,7- vs 2,8-dichlorodioxin), linkers of different lengths, and **"scrambled peptide sequences of identical composition and length"**. (2020)
- Mechanistic reason given for the peptide failure: substructure fingerprints span a bounded radius (MHFP6 six bonds; ECFP4 and TT four bonds), so differentiation that requires atom pairs at longer topological distance is impossible. Worked example: heptapeptides KLLKKLL vs KLKKLLL are distinguished only by MAP4 and AP, not by the substructure fingerprints. (2020)

Measured numbers in that paper:

- **Peptide benchmark design**: 60 datasets built from 30 random linear sequences (ten 10-mers, ten 20-mers, ten 30-mers, all-proteogenic); for each, 10,000 unique scrambled and 10,000 unique point-mutated sequences; BLASTp (BLOSUM62, gap open 11, extend 1, E<10) analogues labelled active. Mean actives: 500.2 ± 0.7 (5.3%) for mutated sets, 56.0 ± 27.4 (0.6%) for scrambled sets. Metrics: AUC, EF1, EF5, BEDROC20/100, RIE20/100; 5 queries per set; ranks compared by Friedman test with Wilcoxon-Nemenyi-McDonald-Thompson post hoc. (2020)
- **Benchmark outcome**: in the peptide benchmark "atom-pair fingerprints significantly outperform substructure fingerprints"; in the small-molecule benchmark the reverse (AP and MXFP significantly worse). MAP4 is the only fingerprint good in both; MAP4 vs other atom-pair fingerprints on peptides is *not* statistically significant. Note: the paper reports **average ranks**, not per-fingerprint AUC numbers, in the main text (Figs. 2-4). (2020)
- **Nearest-neighbour resolution on HMDB 4.0** (96,456 stereochemistry-stripped metabolites; exhaustive NN search; molecules whose nearest neighbour is at Jaccard distance exactly 0, i.e. indistinguishable):

  | fingerprint | indistinguishable NN (of 96,456) | % |
  |---|---|---|
  | MAP4-1024 | 0 | 0% |
  | AP (unhashed) | 1,677 | 1.7% |
  | TT (unhashed) | 68,623 | 71.1% |
  | MHFP6-1024 | 69,972 | 72.5% |
  | **ECFP4-1024** | **70,329** | **72.9%** |

  (Table 4, 2020.) This is the paper's headline "over 70% ... indistinguishable from their nearest neighbor using substructure fingerprints". It is metabolites, not peptides, but it is the only hard collision-rate number the paper gives.
- **Bin saturation / degeneracy**: for ECFP4 and MHFP6 the ten most populated fingerprint-value bins each contain thousands of HMDB molecules (hundreds for TT); for MAP4 it is one molecule per bin and for AP at most two or three. (2020)
- **Database scale context**: the SwissProt peptide set used had HAC = 237.4 ± 104.7, versus ChEMBL 30.0 ± 17.5 — i.e. the peptide regime is ~8x the heavy-atom count fingerprints were tuned on. In SwissProt TMAPs, MAP4 separates by size and groups BLAST analogues better than MHFP6. (2020)
- **Shingle frequency**: in ChEMBL, MAP4 (r=2) yields 46,430,912 unique atom-pair shingles; half occur once, but the most common shingle is present in 85% of ChEMBL structures. (2020)
- **Stated limitation**: "the current version of the MAP fingerprint is implemented in Python and therefore it is relatively slow". (2020)

### 1.2 The follow-up (MAP4C paper) restates the ECFP peptide failure with numbers

Orsi & Reymond, "One chiral fingerprint to find them all", J Cheminform 16:53 (2024), https://pmc.ncbi.nlm.nih.gov/articles/PMC11090803/ :

- States plainly that MAP4 distinguishes structures across classes "for which other fingerprints such as the classical Morgan (ECFP4) and Atom Pair (AP) fingerprints fall short". (2024)
- Distinct-fingerprint counts, 2048 bits (Table 1, 2024) — these are direct measurements of ECFP collapse on peptides:
  - ln65 (antimicrobial undecapeptide, Lys/Leu only), 2,048 possible stereoisomers: MAP2C/4C/6C all 2,048; ECFP6C 1,140; **ECFP4C 36**; APC 196.
  - ln65 scrambled sequences, 330 possible: MAPCs 330; APC 330; **ECFP6C 8; ECFP4C 4**.
  - ln65 diastereomers x scrambled, 675,840: MAPCs 675,840; APC 90,217; ECFP6C 38,500; **ECFP4C 144**.
  - nona-arginine (R9), 512 stereoisomers: MAPCs 512; APC 146; ECFP6C 88; **ECFP4C 12**.
  - Polymyxin B2 scrambled, 1,512: MAPCs 1,512; APC 1,512; ECFP6C 861; **ECFP4C 75**.
- Reasons given: ECFP collapse on ln65 is attributed to the peptide being built from only two residue types, "which reduces the number of possible substructures"; failure on scrambled sets is attributed to "the absence of long-range substructures in ECFP fingerprints". (2024)
- Chirality is size-correlated: percentage of chiral molecules rises steadily with heavy-atom count across ChEMBL, COCONUT and ZINC (Fig. 1a, 2024). Large molecules "are almost all chiral".

### 1.3 Shared amide backbone and bit-vector saturation

No paper found that measures "amide-backbone bit saturation" under that name. What is published:

- The mechanism is documented indirectly as above: peptides built from few residue types yield few distinct circular substructures, so ECFP bit sets collapse (MAP4C, 2024), and substructure fingerprints cannot see permutations of identical composition (MAP4, 2020).
- The complementary measured statement is the HMDB table above (72.9% of molecules have an ECFP4-identical nearest neighbour) and the bin-occupancy result (thousands of molecules per bin for ECFP4/MHFP6). (MAP4, 2020)
- Background-distribution evidence: pairwise ECFP4 Tanimoto for ChEMBL actives vs random ZINC compounds is centred at **Tc = 0.11**, and active-vs-active at 0.15 (0.28 for "easy" classes) — i.e. ECFP4 values live in a narrow low band, while MACCS values are centred at 0.40/0.47. Jasial, Hu, Vogt & Bajorath, F1000Research 5:591 (2016), https://pmc.ncbi.nlm.nih.gov/articles/PMC4830209/ . These are drug-like molecules, not peptides, but they establish that the *scale* of a Tanimoto value is fingerprint-specific.
- Practical corroboration: GuacaMol's own ChEMBL preprocessing pre-filters SMILES to 5-200 characters "because processing larger molecules (e.g. peptides) takes very long" — peptides are explicitly out of scope of that pipeline. https://github.com/BenevolentAI/guacamol/blob/master/guacamol/data/get_data.py (2019).

### 1.4 Size dependence of Tanimoto

- Flower, "On the properties of bit string-based measures of chemical similarity", J Chem Inf Comput Sci 38:379-386 (1998), doi:10.1021/ci970437z — empirical result that bit strings give "a nonintuitive encoding of molecular size, shape, and global similarity" and that observed bit-string search behaviour has a large non-specific component.
- Holliday, Salim, Whittle & Willett, "Analysis and display of the size dependence of chemical similarity coefficients", J Chem Inf Comput Sci 43:819-828 (2003), doi:10.1021/ci034001x — derives the upper bounds of 14 coefficients as a function of relative bit density; identifies "an additional numerical contribution to the known size bias in the Tanimoto coefficient"; shows most coefficients are biased by the relative number of bits set in reference vs database molecule. https://eprints.whiterose.ac.uk/id/eprint/9229/
- Dixon & Koehler, "The hidden component of size in two-dimensional fragment descriptors", J Med Chem 42:2887-2900 (1999), doi:10.1021/jm980708c — 1-Tanimoto as a distance makes collections of *small* compounds look spuriously more diverse; XOR/squared-Euclidean biases the other way.
- Bajorath group summary (Bender/Todeschini-adjacent review, J Cheminform 7:20, 2015, https://pmc.ncbi.nlm.nih.gov/articles/PMC4456712/): Tanimoto tends to pick small compounds in dissimilarity selection, and Godden et al. (JCICS 40:163, 2000) showed statistically preferred Tc values (~1/3) even for structurally distant molecules.
- Consequence for a 90-390 heavy-atom corpus: the reference molecule and database molecules differ enormously in bit count, which is precisely the regime Holliday et al. (2003) show the coefficient is biased in; no paper found that quantifies this specifically at 90-390 heavy atoms.

### 1.5 Does RDKit Morgan encode stereochemistry by default? No.

- Parameter name: `useChirality` in the legacy API (`AllChem.GetMorganFingerprint`, `GetMorganFingerprintAsBitVect`), `includeChirality` in the FingerprintGenerator API (`rdFingerprintGenerator.GetMorganGenerator`).
- **Default is `False`** in both. C++ signature: `getMorganGenerator(unsigned int radius, bool countSimulation=false, bool includeChirality=false, bool useBondTypes=true, ...)`; `MorganArguments(..., bool includeChirality=false, ...)` with `df_includeChirality = false`; and the bond-invariant generator `MorganBondInvGenerator(const bool useBondTypes=true, const bool useChirality=false)`. https://www.rdkit.org/docs/cppapi/MorganGenerator_8h.html , https://rdkit.org/docs/cppapi/classRDKit_1_1MorganFingerprint_1_1MorganArguments.html (RDKit docs, accessed 2026).
- Doc wording: `includeChirality` — "include chirality in atom invariants (not for all fingerprints)". https://www.rdkit.org/docs/source/rdkit.Chem.rdFingerprintGenerator.html
- So a default Morgan/ECFP pipeline gives **identical fingerprints to a peptide and its all-D enantiomer / any diastereomer**. Even with `includeChirality=True`, ECFP4C distinguished only 36 of 2,048 ln65 stereoisomers (Orsi & Reymond 2024, above).
- Known RDKit behaviour caveat for the chiral variant: issue rdkit/rdkit#7986 (2024) — with `includeChirality=True`, Morgan distinguishes chiral atoms at radii where the enclosed environments are identical, which the maintainer flags as wrong. https://github.com/rdkit/rdkit/issues/7986
- Note also that the MAP4 reference similarity-search web tool strips chirality before searching ("chirality information is removed with RDKit") — the achiral MAP4 is chirality-blind too. (MAP4 paper, 2020)

---

## SECTION 2 — Peptide-appropriate similarity measures

### 2.1 MAP4 (MinHashed Atom-Pair fingerprint, diameter 4)

- Definition: for every atom pair (j,k), write the canonical, non-isomeric, rooted SMILES of the circular substructure around each atom at radii r=1 and r=2, order the two alphabetically, and join them with the shortest topological distance TP_jk: `CS_rj | TP_jk | CS_rk`. Shingles are SHA-1 hashed and MinHashed into a fixed-length vector (1024 default; 2048 variant). https://github.com/reymond-group/map4 (2020)
- Captures: local substructure (ECFP-like) **and** unbounded topological distance (AP-like), hence sequence order, length, scrambling, and chain topology.
- Similarity: Jaccard estimated from MinHash. **Critical implementation note from the repo**: MinHashed fingerprints cannot be compared with "standard" Jaccard/Manhattan/cosine functions computed feature-wise — the order of features matters; use the distance function provided (test.py). (2020)
- Implementation/licence: Python + RDKit, `reymond-group/map4`, **MIT**, last pushed 2023-05-01 (GitHub API, 2026). Also LSH-forest search portals at http://map-search.gdb.tools/ .
- Non-canonical residues: no monomer dictionary needed — it is an atom-graph fingerprint, so NCAAs, N-methylation and macrocycles are handled natively.
- Stereochemistry: **no** (achiral by construction; SMILES written non-isomeric).
- Speed: authors state the Python implementation "is relatively slow" and suggest a C/C++ rewrite. (2020)

### 2.2 MAP4C / MAPC (chiral MAP)

- Paper: Orsi & Reymond, "One chiral fingerprint to find them all", J Cheminform 16:53 (2024), doi:10.1186/s13321-024-00849-6, https://pmc.ncbi.nlm.nih.gov/articles/PMC11090803/ . Data/code: https://zenodo.org/records/10389905 ; **repo https://github.com/reymond-group/mapchiral** (mirror https://github.com/markusorsi/mapchiral), pip `mapchiral`. GitHub API reports **no recognised SPDX licence** (NOASSERTION) on reymond-group/mapchiral, 2 stars, last pushed 2025-01-16 — licence must be confirmed with the authors before redistribution.
- Definition delta vs MAP4: at the **largest** radius only, when the central atom of a circular substructure is a stereocentre, its first atom symbol is replaced by its CIP descriptor (R, S, r, s) wrapped in `$...$`, or `?` if undefined; `@`/`@@` are stripped from the extracted SMILES while `/` and `\` (E/Z) are kept; radius 0 skipped. Allene and atropisomeric (biaryl, helicene) chirality are **not** handled. (2024)
- Design consequence (stated): chiral information appears in only a small fraction of shingles, roughly proportional to the fraction of chiral atoms, so a defined stereoisomer stays relatively similar to the undefined-stereochemistry parent, and stereochemistry's weight scales down with molecule size. (2024)
- API: `from mapchiral.mapchiral import encode, jaccard_similarity`; `encode(mol, max_radius=2, n_permutations=2048, mapping=False)`.
- Measured capability on exactly the researcher's molecule class (Table 1, 2048 bits, 2024): distinguishes all 2,048 ln65 stereoisomers, all 330 scrambled, all 675,840 combined; all 512 R9 stereoisomers; all 4,096 polymyxin B2 stereoisomers, 1,512 scrambled, 774,144 combined, 531,441 R/S/undefined assignments. For macrocycles: MAP4C resolves all 136 quinaldopeptin and all 2,080 onchidin stereoisomers (C2-symmetric), but **fails on internal rotational symmetry**: gramicidin S (C2) 504/528 for MAP4C (MAP6C 528/528); valinomycin (C3) 714/1,376; nonactin (C4) 16,176/16,456; NP213 cyclic hepta-arginine (C7) 13/20. Larger bit sizes or no MinHashing "did not increase performance significantly".
- Ordering of distances (2024): for ln65 and polymyxin B2, Jaccard distance grows with Levenshtein distance, and stereoisomer pairs are closer than sequence-isomer pairs **only** for chiral MAP fingerprints and APC; chiral ECFPs get this backwards.
- Virtual screening: chiral vs achiral differences not significant (Friedman-Nemenyi); MAP(C) significantly beats ECFP(C) and AP(C) across metrics except AP(C) on AUC. MAP4C best rank on small molecules, MAP6C best on peptides. (2024)
- Authors' recommendation: MAP4C, because it computes faster than MAP6C (fewer atom pairs).

### 2.3 MHFP6

- Probst & Reymond, J Cheminform 10:66 (2018), doi:10.1186/s13321-018-0321-8, https://link.springer.com/article/10.1186/s13321-018-0321-8 . Code https://github.com/reymond-group/mhfp .
- Definition: enumerate circular substructures around every atom up to diameter 6 as SMILES ("molecular shingling"), **plus the SMILES of each ring in the SSSR** (added precisely because "for either small radii r or macrocycles the ring information of a molecule is lost"), deduplicate, hash to 32-bit, MinHash. Jaccard/Tanimoto is the right metric; MinHash estimation error is O(1/log n) in the number of permutations.
- Captures: substructures only — **no inter-substructure distance**, so it inherits ECFP's peptide blindness. Measured: 72.5% of HMDB molecules have an MHFP6-identical nearest neighbour (vs 0% for MAP4). (MAP4 paper, 2020)
- Advantage over ECFP4 for large corpora: LSH forest ANN search is directly applicable (ECFP's folding scheme prevents it); reported ~2 orders of magnitude faster ANN with lower error. (2018)
- Non-canonical residues: yes (atom graph). Stereochemistry: no.

### 2.4 Atom-pair fingerprints (AP, Carhart)

- Carhart, Smith & Venkataraghavan, J Chem Inf Comput Sci 25:64-73 (1985), doi:10.1021/ci00046a002. RDKit `rdMolDescriptors.GetAtomPairFingerprint` / `GetHashedAtomPairFingerprint`.
- Definition: every pair of atoms, each typed by atomic number + number of heavy neighbours + number of pi electrons, together with their shortest-path bond distance. Unbounded distance.
- Captures: shape, size, and long-range arrangement; this is exactly why it wins the MAP4 peptide benchmark. It does **not** perceive atom environments, only atomic properties — measured consequence: AP fails to separate 4-phenanthrol from 9-phenanthrol, and 1,677 HMDB molecules (96.1% of them phospholipid-like) have an AP-identical nearest neighbour. (MAP4, 2020)
- Chiral variant APC exists (used in Orsi & Reymond 2024) and is much weaker than MAPC on peptide stereoisomers (196/2,048 for ln65) but did resolve all 330 and all 1,512 scrambled sets.
- Speed/licence: RDKit C++, BSD-3 — the fastest option here by a wide margin. Handles NCAAs natively.
- RDKit's AP has `includeChirality` as well; default `False`.

### 2.5 MXFP (macromolecule extended atom-pair fingerprint)

- Capecchi, Awale, Probst & Reymond, Mol Inf 38:1900016 (2019), doi:10.1002/minf.201900016, https://onlinelibrary.wiley.com/doi/10.1002/minf.201900016 . Code https://github.com/reymond-group/mxfp_python .
- Definition: 217-D "fuzzy" vector. Atoms are assigned to 7 pharmacophore categories (heavy atoms, hydrophobic, aromatic, HBA, HBD, positive, negative); for each category every atom pair contributes a Gaussian of 18% width centred at its topological (or 3D Euclidean) distance, sampled at 31 exponentially spaced bins from 0 to 317.8 bonds; the per-bin sum is normalised by N_c^1.5 **specifically to reduce size sensitivity**, x100, rounded. 7 x 31 = 217.
- Distance: **Manhattan / city-block**, not Tanimoto (as used in the MAP4 benchmark).
- Built for exactly this size regime (non-Lipinski PubChem/ChEMBL, peptides, PDB). Handles NCAAs natively; no stereochemistry.
- Measured position: significantly worse than substructure fingerprints on small molecules; in the peptide benchmark it is in the atom-pair group that significantly beats substructure fingerprints, but MAP4 ranks above it. (MAP4, 2020)
- Used as the search space of PDGA, the peptide genetic algorithm that generates linear, cyclic/polycyclic and dendritic topologies (Capecchi, Zhang & Reymond, JCIM 60:121-132, 2020, doi:10.1021/acs.jcim.9b01014, https://github.com/reymond-group/PeptideDesignGA).

### 2.6 Monomer-level fingerprints (order-aware monomer graphs)

- **PepFuNN** (Ochoa & Deibler, Novo Nordisk, J Pept Sci 31:e3666, 2025, doi:10.1002/psc.3666; https://github.com/novonordisk-research/pepfunn, **MIT**, actively maintained — last push 2026-03). Preprint PDF: https://chemrxiv.org/engage/api-gateway/chemrxiv/assets/orp/resource/item/66c5ba0020ac769e5f50bbfa/original/...pdf
  - Monomer-based fingerprint: the peptide is parsed from **BILN** into an igraph graph whose nodes are monomers (so multiple chains, cyclisation and non-canonical monomers are representable), the graph is decomposed into fragments of up to 3 consecutive monomers in all directions — "the Morgan approach ... but with the monomers playing the role of the atoms" — each fragment becomes a token built from monomer physicochemical properties, hashed into a fixed-size (e.g. 2048-bit) fingerprint compared by Tanimoto.
  - Also ships: weighted alignments using a precomputed monomer-vs-monomer similarity matrix over the open **HELM monomer dictionary** (322 monomers embedded; extensible by editing the monomer SDF), Biopython `pairwise2` for unequal lengths, Hamming distance, and RDKit Morgan (radius 4) on the full peptide SMILES.
  - Stated limitation: non-canonical residues are handled "in some cases, based on the availability of a public monomer dictionary"; in the Sequence module NCAAs are **replaced by their closest natural analogue** (chosen by RDKit fingerprint similarity) before property calculation.
  - Stereochemistry: inherited from how monomers are defined in the dictionary (D- and L- forms are different monomer symbols); not an explicit CIP encoding.
- **pyPept** (Boehringer Ingelheim, https://github.com/Boehringer-Ingelheim/pyPept, **MIT**): BILN/HELM/FASTA to RDKit mol and 3D; the standard route from a monomer notation to an atom graph you can then fingerprint with MAP4C.
- **CycPeptMPDB** monomer definition (Li et al., JCIM 63:2240-2250, 2023, doi:10.1021/acs.jcim.2c01573, http://cycpeptmpdb.com/): 7,991 cyclic peptides from 56 publications, decomposed into **385 monomer types** (312 in v1.0) classified into 21 natural-analogue categories, with capping rules for cleaved attachment points. This is the de-facto public monomer vocabulary for macrocyclic peptidomimetics. Context number: **>99.6% of CycPeptMPDB cyclic peptides contain non-natural amino acids** (CycPeptMP, Brief Bioinform 25:bbae417, 2024, https://github.com/akiyamalab/cycpeptmp).
- Fully order-independent monomer **sets** (composition-only bags) were not found as a published, named fingerprint; the published monomer fingerprints above are fragment-based and therefore partly order-aware. Treat a pure composition vector as a trivially implementable but unvalidated baseline.

### 2.7 Sequence alignment (Needleman-Wunsch, BLOSUM/PAM, BLAST)

- Plainly: **standard substitution matrices do not cover non-canonical residues at all.** BLOSUM62/PAM250 are 20x20 (plus ambiguity codes) matrices estimated from natural protein alignments; there is no entry for Aib, D-Leu, N-Me-Phe, etc. The MAP4 benchmark itself only used BLAST because its peptides were "generated with each of all 20 proteogenic amino acids" (2020).
- What the field does instead:
  - Replace each NCAA by its closest natural analogue and align normally (PepFuNN Sequence module, 2025; CycPeptMPDB's 21 natural-analogue categories, 2023) — lossy by construction.
  - Build a **custom monomer-vs-monomer scoring matrix** over a monomer dictionary and run NW/SW with it (PepFuNN Similarity module over 322 HELM monomers, 2025; the matrix is regenerable when monomers are added).
  - Treat NCAAs as **opaque tokens with no substitution matrix** (match/mismatch only) — e.g. PepSeA (cited by PepFuNN as the open HELM-aware MSA tool) and the Datagrok `@datagrok/sequenceutils` HELM MSA package (NW with affine gaps over integer-encoded monomers, cyclic rotation normalisation; https://registry.npmjs.org/@datagrok%2Fsequenceutils).
- Alignment is also blind to D/L inversion unless D-residues are distinct monomer symbols, and blind to cyclisation topology unless handled specially (Section 6).
- Standard redundancy-removal practice in peptide ML is CD-HIT identity clustering (e.g. 70% identity cluster-level train/test split, or 90% for redundancy reduction) — see ProDCARL (2026 preprint, https://github.com/HIVE-UofT/ProDCARL) — which again requires canonical alphabets.

### 2.8 Levenshtein / edit distance on HELM or BILN strings

- Used as a *reference* axis in the MAP4C paper: Levenshtein distance is defined there as "the minimum number of mutations necessary to transform one sequence into another one, considering residue type changes, **stereochemical inversions**, insertions and deletions"; Jaccard distances of all fingerprints rose monotonically with it. (2024)
- Used as an *operational novelty filter* in a real peptide design campaign: Reymond group ML-PDGA anticancer peptides (SI, https://github.com/reymond-group/ML-PDGA-anticancer-peptides): generated sequences were kept only with **LD > 5 from the classifiers' training set and LD > 4 from the test set**; the applicability domain was set from the 90% quantile of test-to-train minimum LD, excluding anything with LD >= 8 from the training set; final selection clustered with the RDKit Butina module using Levenshtein as the distance function (threshold 10).
- Properties: trivially handles NCAAs and D-residues if they are distinct tokens; cheap (O(mn), `python-Levenshtein`); **but** it is defined on a linearisation, so it is sensitive to the arbitrary start point of a macrocycle and cannot see branches/staples unless the notation is normalised first.
- Normalised variants in use: normalised LD = LD / max(len) as diversity/novelty metric in AMP generative work (PLOS Comput Biol, 2026, doi:10.1371/journal.pcbi.1014771).

### 2.9 Pharmacophore / CATS descriptors

- CATS: Schneider, Neidhart, Giller & Schmid, Angew Chem Int Ed 38:2894-2896 (1999), doi:10.1002/(SICI)1521-3773(19991004)38:19<2894::AID-ANIE2894>3.0.CO;2-F — topological pharmacophore correlation vector; atoms typed as lipophilic / H-bond donor / H-bond acceptor (+ positive/negative), all type-pair occurrences counted per topological distance, then scaled ("types scaling") to stop frequent feature types dominating; max correlation distance typically 10 bonds. CATS2 adds an aromatic type. 3D and surface variants: CATS3D, SURFCATS.
- Measured behaviour (Reutlinger et al., Mol Inf 32:133-138, 2013, https://pmc.ncbi.nlm.nih.gov/articles/PMC3743170/): radial (Morgan) and Carhart atom-pair fingerprints **retrieve more actives** (higher BEDROC), but CATS1 gives the **highest ratio of distinct scaffolds** among retrieved actives. CATS is deliberately "fuzzy".
- Implication here: fuzziness is the wrong property for a leakage filter — a descriptor built to call different scaffolds similar will over-exclude, and one built to scaffold-hop has by design low resolution between close analogues. Useful as a *second, orthogonal* filter, not as the primary one. No peptide-specific CATS benchmark found; MXFP is effectively a pharmacophore atom-pair descriptor designed for the large-molecule regime and is better supported for peptides.

### 2.10 Ranking for this use case (planting a target, excluding its neighbourhood from a peptide/peptidomimetic pretraining corpus)

Criteria: resolves scrambled sequences; resolves stereochemistry (D-AA, which the molecules have); handles NCAAs/N-methylation/macrocycles without a dictionary; published evidence in this size regime; usable implementation.

1. **MAP4C / MAPC** (2048 permutations, max_radius=2). Only measure with published evidence of resolving thousands of peptide stereoisomers *and* sequence isomers, with distances ordered sensibly against Levenshtein. Caveats: licence unclear on the repo; Python speed; fails on high internal rotational symmetry (C3+); MinHash means you must use the package's own Jaccard function.
2. **MAP4** (achiral) — same topology perception, MIT-licensed, more widely used; use when a chirality-blind filter is acceptable or as a conservative superset filter (it will call a peptide and its diastereomer identical, so it over-excludes, which is the safe direction for a holdout).
3. **RDKit Atom Pair (+ `includeChirality=True`)** — the fast screening pass. Significantly beats substructure fingerprints on peptides (MAP4 benchmark, 2020), C++ speed, BSD. Weakness: no atom environments; APC resolved only 196/2,048 ln65 stereoisomers.
4. **MXFP (Manhattan)** — purpose-built for 90-390 heavy atoms and explicitly size-normalised; good as an orthogonal second axis. Not Tanimoto-comparable; no stereochemistry.
5. **Monomer-level fingerprint (PepFuNN over BILN) + monomer-weighted alignment** — the only measure that speaks the medicinal chemist's unit of change (one residue swapped). Requires a monomer dictionary covering your NCAAs; MIT; actively maintained.
6. **Levenshtein on a normalised BILN/HELM string** — cheap, interpretable, has precedent as a published novelty criterion (LD > 5), but linearisation-dependent; needs cyclic rotation normalisation.
7. **MHFP6** — strictly dominated by MAP4 here (same substructure blindness, 72.5% HMDB NN collisions) although its SSSR ring shingles help with macrocycles relative to plain ECFP.
8. **ECFP4/Morgan-Tanimoto** — documented failure mode for exactly these molecules; keep only as the legacy comparator you report against.
9. **CATS / pharmacophore** — deliberately fuzzy, scaffold-hopping-oriented; wrong resolution for a leakage filter.
10. **BLOSUM/NW on raw sequence** — cannot represent the alphabet; excluded except via a custom monomer matrix.

---

## SECTION 3 — Threshold conventions, and whether they transfer

### 3.1 Where 0.85 came from, and what it actually buys

- Origin of the "reject anything >= 0.85 similar" rule: compound-selection / neighbourhood-behaviour work of the mid-1990s. Patterson, Cramer, Ferguson, Clark & Weinberger, "Neighborhood behavior: a useful concept for validation of molecular diversity descriptors", J Med Chem 39:3049-3059 (1996), doi:10.1021/jm960290n — validates descriptors by the trapezoidal plot of descriptor distance vs activity difference, and ranks 11 descriptors on 20 datasets.
- The claim attached to 0.85: early clustering/diversity work using **MACCS keys and UNITY fingerprints** indicated that on average "85% of compounds that yielded a Tc value of 0.85 compared to a known active molecule were also active" (as summarised by Jasial et al., 2016, https://pmc.ncbi.nlm.nih.gov/articles/PMC4830209/).
- The measurement that overturned it: Martin, Kofron & Traphagen, "Do structurally similar molecules have similar biological activity?", J Med Chem 45:4350-4358 (2002), doi:10.1021/jm020155c — using **Daylight** fingerprints and IC50 follow-ups to 115 HTS assays, there is **only a 30% chance** that a compound >= 0.85 Tanimoto to an active is itself active. The authors attribute this partly to fingerprint/Tanimoto deficiencies and partly to the fact that similar compounds need not bind similarly; they stress library design is probabilistic.
- 0.4 and 0.7 are not single-origin conventions. Documented anchors:
  - ECFP4 Tc ~0.3-0.4 is the *activity-relevant* band measured by Jasial et al. (2016): the random-vs-active ECFP4 distribution reaches baseline at ~**0.3**, which they put on a par with MACCS ~**0.8**. For "easy" activity classes, 38.2% of active-active ECFP4 pairs are >= 0.3 vs 0.03% of random-active pairs.
  - Muchmore, Debe, Metz et al. (belief theory, J Chem Inf Model 48:941-948, 2008, doi:10.1021/ci7004498), as reported by Jasial et al.: at Tc 0.85 with atom-pathway fingerprints ~30% of pairs share activity; **ECFP6 reaches the same point at Tc = 0.42**.
  - 0.7 circulates as a generic clustering/novelty cutoff and is used that way in current practice — e.g. a 2025 molecular generative supplementary defines novelty as max Tanimoto to training < 0.7 and similarity-uniqueness clusters at >= 0.7 to centroid, calling 0.7 "a common benchmark for chemical similarity searches" (https://www.repository.cam.ac.uk/bitstreams/2df67203-b01c-4579-a81e-fd19cda9e65f/download). Butina's 1999 clustering (JCICS 39:747) used Daylight fingerprints at 0.8.

### 3.2 Published guidance on choosing a threshold for a DIFFERENT fingerprint

Yes — and it is the method GuacaMol itself used.

- **Landrum, RDKit blog, "Fingerprint thresholds" (October 2013)**, http://rdkit.blogspot.com/2013/10/fingerprint-thresholds.html . Method: take 25,000 random drug-like ChEMBL pairs (MW < 600), compute the similarity distribution of *random* pairs per fingerprint, and read off the 90th/95th/99th percentile. A threshold is then "the value a random pair almost never exceeds" for that fingerprint. Published table:

  | fingerprint | 90% | 95% | 99% |
  |---|---|---|---|
  | MACCS | 0.528 | 0.573 | 0.652 |
  | Morgan0 (counts) | 0.525 | 0.568 | 0.649 |
  | Morgan1 (counts) | 0.333 | 0.365 | 0.428 |
  | **Morgan2 (counts) = ECFP4** | **0.230** | **0.255** | **0.306** |
  | RDKit4 | 0.283 | 0.325 | 0.425 |
  | Avalon (512 bits) | 0.462 | 0.504 | 0.579 |
  | Atom pairs (counts) | 0.238 | 0.266 | 0.326 |
  | Torsions (counts) | 0.165 | 0.198 | 0.264 |

  The spread across fingerprints at a fixed percentile (0.65 for MACCS vs 0.31 for Morgan2 vs 0.26 for torsions at 99%) is the direct, quantitative answer to "do thresholds transfer": **they do not transfer as numbers; they transfer as percentiles of the fingerprint's own random-pair background.**
- GuacaMol cites exactly this blog post (reference 84) for its 0.323 ECFP4 cutoff. https://ar5iv.labs.arxiv.org/html/1811.09621 (2019).
- Independent statement that no universal threshold exists: Jasial, Hu, Vogt & Bajorath (2016) — "generally applicable Tc threshold values as an indicator of activity similarity do not exist", because similarity-value distributions are compound-class- and fingerprint-dependent; and even the activity-relevant values they *do* identify (MACCS >= 0.8, ECFP4 >= 0.3) cannot be used as search thresholds because of the ratio of active-active to active-random comparisons in a real database. Their worked example: ECFP4 at 0.3 over 50 actives + 500,000 inactives yields ~19 true and ~150 false positives.

### 3.3 Measured similarity distributions for peptide datasets

Thin, but these exist:

- **MAP4C supplementary Figures S1 and S2** (2024) report, for every dataset in the extended Riniker-Landrum benchmark **including the 60 peptide sets**, the mean and standard deviation of pairwise **ECFP4C Tanimoto** (S1) and pairwise **MAP4C Jaccard** (S2) similarities of the five selected actives. This is the closest published side-by-side peptide similarity distribution under two fingerprints; the values themselves are only in the figures of Additional file 1 (doi:10.1186/s13321-024-00849-6).
- **MAP4C Figure 4b/4c** (2024): Jaccard distance vs Levenshtein distance for 330 ln65 sequence isomers vs 2,048 diastereomers, and for 1,512 polymyxin B2 sequence isomers vs 512 diastereomers, for MAP2C/4C/6C/APC/ECFP4C/ECFP6C (and Figs S17-S18). Reported qualitative result: distance rises with LD for all fingerprints; only chiral MAPs and APC place stereoisomers closer than sequence isomers; TMAPs show complete separation of the stereoisomer and sequence-isomer clouds.
- **HMDB nearest-neighbour table** (Section 1.1) is a distribution statement at distance 0 for 96,456 large/medium molecules.
- Jasial et al. (2016) distributions are drug-like, MW < 550, explicitly excluding anything above that to balance size effects — i.e. **the published Tc background distributions were measured on a size regime your molecules sit entirely outside of.**

Practical consequence stated by the sources, not opinion: to pick a cutoff for MAP4C on a 90-390-heavy-atom peptide corpus, the only published method that applies is Landrum's — compute the random-pair distribution **in your own corpus with your own fingerprint** and take a high percentile. No pre-computed table exists for MAP4/MAP4C or for peptides.

---

## SECTION 4 — How benchmarks actually define holdouts

### 4.1 GuacaMol (Brown, Fiscato, Segler & Vaucher, J Chem Inf Model 59:1096-1108, 2019)

Paper: doi:10.1021/acs.jcim.8b00839 ; arXiv:1811.09621 ; code https://github.com/BenevolentAI/guacamol (MIT).

Exact procedure (Appendix 8.1 and `guacamol/data/get_data.py`):
1. Start from **ChEMBL 24.1** chemreps.
2. Pre-filter raw SMILES to length 5-200 chars — the code comments that this excludes large molecules "e.g. peptides" because they are slow.
3. Remove salts; neutralise charges; drop SMILES > 100 chars after canonicalisation; keep only H, B, C, N, O, F, Si, P, S, Cl, Se, Br, I.
4. **Holdout exclusion**: compute ECFP4 bit vectors — concretely `AllChem.GetMorganFingerprintAsBitVect(m, radius=2, nBits=4096)` — for the 10 holdout drugs (celecoxib, aripiprazole, cobimetinib, osimertinib, troglitazone, ranolazine, thiothixene, albuterol, fexofenadine, mestranol, file `holdout_set_gcm_v1.smiles`), and **remove every ChEMBL molecule whose maximum Tanimoto to any of them exceeds `TANIMOTO_CUTOFF = 0.323`**. Stereocentres are dropped during canonicalisation (`include_stereocenters=False`).
5. Shuffle with `np.random.seed(1337)`, split 5% valid / 15% test / 80% train, verify md5 hashes.

Key facts about it: the holdout is **10 molecules**, the criterion is **one fingerprint at one number**, 0.323 is taken from Landrum's random-pair percentile table (Section 3.2), and the rediscovery benchmarks (Celecoxib / Troglitazone / Thiothixene) then score `top-1 sim(target, ECFC4)` with threshold 1.0 — i.e. **the same fingerprint family is used to build the holdout and to score the task**. Similarity benchmarks score against FCFP4 (albuterol), AP (mestranol), ECFP4/FCFP4 (aripiprazole, version-dependent) at Thresholded(0.75).

Measured consequence: Graph GA and SMILES LSTM both score **1.000** on all three rediscovery benchmarks, while "Best of Dataset" scores 0.505 / 0.419 / 0.456. The authors themselves note of the trivial distribution benchmarks that "the ChEMBL database alone achieves excellent scores already".

### 4.2 MOSES (Polykovskiy et al., Front Pharmacol 11:565644, 2020)

https://pmc.ncbi.nlm.nih.gov/articles/PMC7775580/ ; https://github.com/molecularsets/moses (MIT).

- Dataset: 1,936,962 molecules from ZINC Clean Leads; 448,854 unique Bemis-Murcko scaffolds; IntDiv1 = 0.857.
- Split: train 1,584,664 / test 176,075 / **scaffold test (TestSF) 176,226**. Procedure in `scripts/prepare_dataset.py`: count scaffold frequencies, sort by (-count, scaffold), take **every 10th scaffold (`scaffolds[9::10]`)** as the scaffold-test set; all molecules carrying those scaffolds go to test_scaffolds; the test set is then a random 10% of what remains.
- Criterion is therefore **exact Bemis-Murcko scaffold identity** (RDKit implementation, which also counts carbonyls attached to rings), not similarity. "Scaffold similarity from the training set to the scaffold test set (TestSF) is zero by design."
- Metrics include Novelty (fraction of valid unique generated molecules not in train — explicitly "measures memorization") and SNN (mean Tanimoto of each generated molecule to its nearest neighbour in the reference set). Reported: CharRNN rediscovered 11% novel scaffolds; VAE has the highest SNN (0.626) and Scaff (0.939) but low novelty (0.695), which the authors flag as possible overfitting to training prototypes.

### 4.3 TDC / PMO

- PMO (Gao, Fu, Sun & Coley, NeurIPS 2022 Datasets & Benchmarks, https://proceedings.neurips.cc/paper_files/paper/2022/file/8644353f7d307baaf29bc1e56fe8e0ec-Paper-Datasets_and_Benchmarks.pdf ; https://github.com/wenhao-gao/mol_opt, MIT): 25 algorithms x 23 oracles (QED, DRD2, GSK3b, JNK3 + **19 GuacaMol oracles**), 10,000-oracle-call budget, AUC top-10. Data rule: "**We restrict all our methods to using the ZINC 250K dataset only**" for pretraining, screening, and fragment extraction.
- **No similarity-based exclusion of the target molecules from ZINC 250k is described** in PMO. The holdout logic that GuacaMol applied to ChEMBL is not reproduced when the GuacaMol oracles are re-run against ZINC 250k. This is a documented gap, not an inference about intent.
- TDC oracles (https://tdcommons.ai/functions/oracles/, MIT) are scoring functions; they do not ship a paired, similarity-filtered pretraining corpus.

### 4.4 Peptide-specific splits

- No peptide analogue of GuacaMol/MOSES with a published similarity holdout was found. What exists is bioinformatics practice carried over: **CD-HIT identity clustering followed by cluster-level splitting** — e.g. clustering AMP sequences at **70% identity** and assigning whole clusters to train or test "to reduce homologous leakage", and 90% identity for redundancy reduction in a toxicity set (ProDCARL, 2026 preprint, https://github.com/HIVE-UofT/ProDCARL). This only works on canonical alphabets.
- The closest peptide-design precedent for an explicit novelty gate is the Reymond ML-PDGA campaign's **Levenshtein distance > 5 to the training set, > 4 to the test set, and < 8 for applicability-domain validity** (Section 2.8).
- The review "Deep generative models for peptide design" (Digital Discovery, 2022) calls for exactly this and states it does not yet exist: benchmarking "may use standardized datasets and metrics", with the dataset "split into a training set (on which all generative models will be trained) and a test set for evaluation", and lists the required heuristics for de novo peptides as (1) non-redundancy, (2) novelty with respect to the training set, (3) diversity with respect to general sequence space. (Retrieved via PMC, 2026.)

### 4.5 Scaffold splits vs similarity splits, and the published criticism

- **Guo, Hernández-Hernández & Ballester, "Scaffold Splits Overestimate Virtual Screening Performance", arXiv:2406.00873 / ICANN 2024** (https://arxiv.org/abs/2406.00873, https://link.springer.com/chapter/10.1007/978-3-031-72359-9_5): 3 model types x 60 NCI-60 datasets (~30k-50k molecules each) x 3 split types (scaffold, Butina, UMAP), 2,100 models per algorithm. Finding: "molecules with different chemical scaffolds are often similar", so scaffold splits leave unrealistically high train-test similarity; performance is much worse under UMAP splits. Recommendation: avoid scaffold splits. Follow-up: J Cheminform 17:94 (2025), doi:10.1186/s13321-025-01039-8.
- **Fooladi, Vu, Mathea & Kirchmair, J Chem Inf Model (2025)**, https://pmc.ncbi.nlm.nih.gov/articles/PMC12529777/ : 14 models, 8 datasets, 10 splitting strategies. Measured: models perform on Bemis-Murcko scaffold splits "not substantially different from random splitting"; UMAP clustering on ECFP4 is the hardest split; ID-OOD performance correlation is Pearson r ~0.9 for scaffold splits but drops to ~0.4 for cluster-based splits.
- **Reimer et al., "Scaffold and UMAP splits are closer to random than we think"** (https://people.cs.umass.edu/~annagreen/publication/reimer_chemspace/reimer_chemspace.pdf): across 8 MoleculeNet datasets and 7 models, after Bonferroni correction only 5/34 scaffold-split and 5/34 UMAP-split model-dataset combinations differed significantly from random; in 3/8 datasets scaffold and UMAP splits were not harder than random by cross-split overlap (CSO). Proposes reporting CSO — an explicit train-test similarity measurement — for any split.
- **DataSAIL** (Joeres, Blumenthal & Kalinina, Nat Commun 16:3337, 2025, doi:10.1038/s41467-025-58606-8, https://www.nature.com/articles/s41467-025-58606-8): formalises leakage-minimising splitting as a combinatorial optimisation (proved NP-hard), with a clustering + ILP heuristic and an explicit leakage score L(pi); shows lower L(pi) correlates with larger performance drops. Also reports that LoHi and DeepChem fingerprint splitting only partly reduce leakage. Important caveat stated by the authors: L(pi) depends entirely on the chosen similarity function, and the user must ensure "the selected similarity function sim indeed captures the intended generalization task" — i.e. a leakage-free split under ECFP can still be leaky under a peptide-appropriate measure.
- Coverage-bias angle: Nat Commun 16 (2025), doi:10.1038/s41467-024-55462-w — even under scaffold splits, "almost all molecular structures in the test split have a somewhat close molecular structure in the train split".
- Protein-side analogue worth noting for method: Bushuiev et al. / "Revealing data leakage in protein interaction benchmarks" (arXiv:2404.10457, 2024) — strict 30%-sequence-identity MMseqs2 splitting still left a **30% structural leakage rate**; family-level splits of DIPS left ~53%. The general lesson, measured: splitting on a *proxy* similarity leaves large residual leakage under the similarity that actually matters.

---

## SECTION 5 — Leakage detection after the fact

Concrete, citable procedures for verifying a generated winner was not memorised:

1. **Exact-match novelty.** Fraction of valid unique generated molecules whose canonical SMILES is absent from the training set (MOSES "Novelty", which the paper describes as measuring memorisation; GuacaMol's novelty benchmark penalises generating training molecules). MOSES: https://pmc.ncbi.nlm.nih.gov/articles/PMC7775580/ (2020); GuacaMol: doi:10.1021/acs.jcim.8b00839 (2019).
   **Documented failure of this test**: Renz, Van Rompaey, Wegner, Hochreiter & Klambauer, "On failure modes in molecule generation and optimization", Drug Discov Today Technol 32-33:55-63 (2019), doi:10.1016/j.ddtec.2020.09.003, https://epub.jku.at/obvulioa/download/pdf/5687408 — a trivial **"AddCarbon"** model that adds a single carbon to a random training molecule achieves near-perfect scores on most GuacaMol distribution-learning benchmarks and beats all baselines but the LSTM on FCD. The authors name this the **copy problem**: current metrics check only exact matches, so minimal edits evade novelty detection. Recommendation: report **test-set likelihood** for likelihood-based models, as in NLP.
2. **Nearest-neighbour similarity to the training set (SNN / aSNN).** MOSES defines SNN(G,R) = mean over generated molecules of max Tanimoto to the reference set; high SNN with low novelty is flagged as probable memorisation (VAE: SNN 0.626, novelty 0.695). (2020)
   - Operationalised for goal-directed campaigns as **aSNN** — average single-nearest-neighbour Tanimoto (Morgan/RDKit) of generated compounds to a named reference subset, implemented in **MolScore**. Used in "On the difficulty of validating molecular generative models realistically: a case study on public and proprietary data", J Cheminform (2023), doi:10.1186/s13321-023-00781-1, https://pmc.ncbi.nlm.nih.gov/articles/PMC10664602/ : measured aSNN of REINVENT output to low/middle/high/ultra-high activity compounds = 0.304/0.367/0.420/0.408 for public projects, 0.431/0.425/0.427/0.348 for in-house projects; rediscovery of middle/late-stage compounds was 1.60%/0.64%/0.21% of the top 100/500/5000 for public projects vs 0.00%/0.03%/0.04% in-house.
   - That paper also states the structural criticism of GuacaMol-style holdouts directly: GuacaMol "just removes the target compound from the training dataset ... However, analogues may still remain in the training dataset", because ChEMBL is assembled from papers that publish whole SAR series. (2023)
3. **Threshold-based novelty on nearest-neighbour similarity.** Define a design novel iff max Tanimoto to any training molecule < T. Current practice uses T = 0.7 on RDKit/Morgan fingerprints, and reports the inverse-novelty curve (fraction of generated molecules above several thresholds) rather than a single number. https://www.repository.cam.ac.uk/bitstreams/2df67203-b01c-4579-a81e-fd19cda9e65f/download (2025).
4. **Control scores (the strongest published audit for goal-directed generation).** Renz et al. (2019): train the scoring model on data split 1 (**optimization score, OS**); train a second model on the same data with a different seed (**model control score, MCS**); train a third on an independent split 2 (**data control score, DCS**). Divergence OS-MCS quantifies model-specific bias; divergence OS-DCS quantifies data-specific bias. Measured on EGFR, DRD2, JAK2 with GA, particle swarm and LSTM: OS rises monotonically while DCS stagnates or falls; optimised molecules drift toward split-1 actives in nearest-neighbour Tanimoto and in t-SNE of ECFP4 distance. Recommendation: stop optimising when control scores plateau.
5. **Likelihood-decile audit.** Hoffmann/Grisoni-style analysis in "How evaluation choices distort the outcome of generative drug discovery" (https://pmc.ncbi.nlm.nih.gov/articles/PMC12613558/): generate 1,000,000 designs, bin by model likelihood decile, and per decile report validity, novelty (not in pre-training or fine-tuning sets), max Tanimoto on ECFP (radius 2, 2048 bits) to the fine-tuning set, and number of distinct substructures. Measured: high-likelihood deciles have high validity and high training-set similarity but low novelty — the exploitation end is where memorised output lives. The same paper identifies a **"size trap"**: FCD and Fréchet descriptor distance both fall with library size and only plateau above ~10,000 (recommended >= 100,000) designs, so similarity/novelty comparisons must be made at equal library size.
6. **Distributional distances as supporting evidence**: FCD (Fréchet ChemNet Distance) and FDD; both in MOSES/GuacaMol. Treat as corroboration only — the AddCarbon result shows FCD alone does not detect copying.

For the researcher's specific question — "was the generated winner memorised?" — the published recipe that actually answers it is: (a) exact-match check against the full pretraining corpus; (b) nearest-neighbour **under the same peptide-appropriate measure used to build the holdout**, reported as the full NN-similarity distribution plus the single winner's value; (c) the Renz control-score experiment if a learned oracle is in the loop; (d) report at a fixed, large library size.

---

## SECTION 6 — Similarity for macrocycles specifically

- **The problem is explicitly stated in the literature.** CyclicPepedia / cyclicpeptide (Liu, Cao, Liu, Zhu & Wu, Brief Bioinform, 2024/2025; https://github.com/dfwlab/cyclicpeptide, https://www.biosino.org/iMAC/cyclicpepedia/): their GraphAlignment tool "overcomes the limitations of traditional linear sequence alignment algorithms when applied to cyclic structures". It converts a cyclic peptide sequence into a NetworkX graph (nodes = residues, edges = bonds) and scores either **graph edit distance** or **maximum common subgraph**; a GCN trained on ~80,000 aligned cyclic-peptide pairs approximates the GED score for speed. Their web tool offers Smith-Waterman local alignment only for one-letter-code peptides and graph alignment (graph isomorphism) otherwise.
- **Rotation invariance is the core failure mode.** Cyclome930 (2026 preprint, https://exa.ai/library/publication/k4788s26c4n ; cyclome930.studio): "Two identical (or highly similar) cyclic peptides may appear dissimilar if one is 'rotated' relative to the other in sequence form", citing Grossi et al. on cyclic permutation invariance. Their fix: for head-to-tail (e2e) cyclised peptides, duplicate and concatenate the template N->C to make a double-length "T-template", then slide the query across all rotational windows and keep the maximum local-alignment score; for side-to-end and side-to-side cyclisation the linear order is well defined and no duplication is used. They also introduce cyclicity-aware ESM-C embeddings with a cyclic offset vector so that rotations of the same peptide get the same representation. Stated motivation: conventional linear calculators "systematically underestimate" cyclic-peptide similarity.
- **Production implementation of the same idea**: `@datagrok/sequenceutils` HELM MSA (https://registry.npmjs.org/@datagrok%2Fsequenceutils) — Needleman-Wunsch with affine gaps over integer-encoded monomers, with explicit macrocycle detection and rotation normalisation. Published rules: head-to-tail (`R2->R1` self-connection) always rotatable; lariat (`R1->R3`/`R2->R3` self) rotatable if the span is >= 70% of the chain; backbone CHEM-bridge rotatable if ring span >= 50%; **side-chain staple (`R3->R3` through CHEM) never rotatable**. Handles non-canonical monomers as opaque tokens with no substitution matrix. Reported throughput: 2,800 sequences in under 1 second.
- **Fingerprint side.** MHFP6 adds the SMILES of each SSSR ring to the shingling precisely because "for either small radii r or macrocycles the ring information of a molecule is lost" (Probst & Reymond, 2018) — the only explicit macrocycle accommodation found inside a mainstream fingerprint.
- **Measured macrocycle results** (Orsi & Reymond, 2024, Table 1 — see Section 2.2): MAP4C resolves all stereoisomers of C2-symmetric quinaldopeptin (136) and onchidin (2,080); MAP6C is needed for gramicidin S (528); **all tested chiral fingerprints fail** on C3 valinomycin (1,250/1,376 best), C4 nonactin (16,425/16,456 best) and C7 NP213 (13/20 best, MAP4C). The stated cause is internal rotational symmetry, not bit size — larger bit vectors did not help. This is a hard, published limitation for symmetric cyclic peptides.
- **CycPeptMPDB-adjacent representation work**: CycPeptMP (Li, Yanagisawa & Akiyama, Brief Bioinform 25:bbae417, 2024; https://github.com/akiyamalab/cycpeptmp) builds features at **atom, monomer and peptide levels** and fuses them, explicitly because whole-molecule small-molecule descriptors "ignore the unique structural properties of cyclic peptides"; monomer-level features are motivated by the fact that permeability SAR is explored one residue at a time, and that "even a minor alteration in a single residue can lead to substantial changes". Reported performance: MAE 0.355, r = 0.883 for log permeability. Monomer capping rules (methylate cleaved amide/ester, aldehyde for carboxyl) are specified. Context: >99.6% of CycPeptMPDB peptides contain NCAAs.
- **Mass-spec dereplication literature** is the oldest statement of the problem: a cyclic peptide of length k gives k ion series rather than 2, so linear-peptide methods do not transfer, and database search only works "if an identical or very close variant is present" (Mohimani et al., https://pmc.ncbi.nlm.nih.gov/articles/PMC3398611/). Relevant as precedent, not as a similarity measure for this task.
- No published study was found that measures *distributions* of macrocycle-vs-macrocycle similarity under competing fingerprints, nor any threshold recommendation for macrocycles.

---

## GAPS

- No published random-pair similarity background distribution, and therefore no percentile-derived threshold table, exists for MAP4/MAP4C/MXFP or for any peptide/peptidomimetic corpus; Landrum's 2013 table is drug-like, MW<600, ECFP/MACCS/AP/TT only. The threshold must be measured in-corpus.
- No molecular generative benchmark was found that builds a similarity holdout with a non-ECFP fingerprint; GuacaMol's 0.323/ECFP4/10-molecule procedure is the only fully specified one, PMO inherits GuacaMol's oracles without reproducing any holdout filter, and MOSES uses exact scaffold identity rather than similarity.
- No peptide-specific generative benchmark with a published similarity-based holdout exists; peptide practice is CD-HIT identity clustering (canonical alphabets only) or ad-hoc Levenshtein gates (LD>5 in the one documented campaign).
- "Amide-backbone bit saturation" is not measured under that name anywhere found; the evidence is indirect (ECFP4 identical-NN rates on HMDB, ECFP4C resolving 36/2,048 stereoisomers of a two-residue-type peptide, bin-occupancy counts).
- All chiral fingerprints tested, including MAP4C/MAP6C, fail on macrocycles with C3 or higher internal rotational symmetry, and larger bit sizes do not fix it — unresolved for symmetric cyclic peptides.
- The mapchiral repository carries no machine-readable licence (GitHub reports NOASSERTION); the achiral map4 repo is MIT. Licence for MAP4C must be confirmed with the authors before use in a released benchmark.
