# Project Specification — hv/CRKP Convergence Prediction
**Status:** DRAFT — commit with a dated hash before any real data exploration begins.
**Frozen on:** 2026-08-29 (UTC)
**Version control:** [git commit hash / OSF registration link — fill in once committed]

---

## 1. Primary Research Question

Can genomic information available *before* convergence identify *K. pneumoniae* lineages
at elevated risk of subsequently acquiring the complementary trait — either
hypervirulence-associated determinants (in a resistant background) or carbapenem
resistance (in a hypervirulent background)?

## 2. Secondary Research Questions

1. Do the genomic predictors of convergence differ between the two evolutionary routes —
   (a) hypervirulent → acquires resistance vs. (b) resistant → acquires
   hypervirulence-associated determinants?
2. Does an evolutionary accumulation model (HyperTraPS-CT) outperform simple baselines
   (ST alone, ST + geography, current genomic state) at this prediction task?
3. Which genomic and mobile-element features (plasmid replicons, insertion sequences,
   capsule type, accessory-genome content) are associated with elevated predicted risk?
4. Which currently non-convergent lineages are predicted to be at elevated risk of
   future convergence?

## 3. Operational Definitions

**Resistance (target trait, primary focus: carbapenem resistance)**
A genome is called carbapenem-resistant if Kleborate reports presence of any gene or
mutation conferring resistance to the carbapenem class (`Bla_Carb` character, per
Kleborate's convention: presence of *any* qualifying gene/mutation = resistant).
Secondary/sensitivity resistance characters (for the fine-grained character set) follow
the same presence/absence convention for other drug classes reported by Kleborate.

**Hypervirulence-associated (target trait)**
- *Coarse definition:* Kleborate virulence score ≥ 3 (established convention from the
  Kleborate reference framework — indicates yersiniabactin plus aerobactin and/or
  virulence plasmid detection).
- *Fine-grained definition:* individual presence/absence of `ybt`, `iuc`, `iro`, `clb`,
  `rmpA`, `rmpA2`.
- Both definitions are carried forward as Model A (coarse) and Model B (fine-grained) —
  see Section 4. Neither is discarded until Phase 1 Step 5 / Phase 2 empirically
  justifies dropping one.

**Convergent (for descriptive/labeling purposes only, not the prediction target itself)**
Virulence score ≥ 3 **and** at least one high-priority resistance determinant present
(carbapenemase, 16S rRNA methyltransferase, plasmid-mediated quinolone resistance,
ESBL, or colistin resistance) — following the convention already established in the
Kleborate genomic surveillance literature, so results remain comparable to prior work.

**Acquisition event (working definition — final version locked in Phase 4, mentor-reviewed)**
A character transitions from absent in an ancestral state to present in a descendant
state, as inferred via ancestral state reconstruction on the phylogeny, where the
transition is not attributable to:
(a) within-clade inheritance from a common ancestor (pseudoreplication),
(b) loss followed by reappearance (handled via the reversibility sensitivity analysis
    inherited from the reference HyperTraPS methodology),
(c) an assembly or genotyping artifact (handled via genome QC exclusion, Section 6).
This definition is provisional pending Phase 4 mentor review and is not to be modified
after Phase 1 data exploration begins.

## 4. Candidate Character Sets

| | Resistance | Virulence |
|---|---|---|
| **Model A (coarse)** | Single carbapenem-resistance binary character | Single virulence-score-≥3 binary character |
| **Model B (fine-grained)** | Individual major resistance classes (carbapenemase, ESBL, aminoglycoside, fluoroquinolone, colistin, etc. — following Kleborate's ~22-class convention) | Individual loci: `ybt`, `iuc`, `iro`, `clb`, `rmpA`, `rmpA2` |

**Decision rule (to be applied in Phase 1 Step 5 / Phase 2, not before):** Model B is
adopted as primary only if (i) total character count stays within the computational
tractability bound established in Phase 1 Step 6, and (ii) no individual character in
Model B falls below a minimum prevalence threshold to be set from the Phase 1 metadata
audit (candidate default: 1% of the qualifying genome set). Otherwise Model A is primary
and Model B is retained only as a sensitivity analysis.

## 5. Intended Prediction Unit

**Default:** clonal complex / sequence type (ST/CC)-level lineage, to correct for the
non-independence of near-identical isolates.
**Fallback:** individual isolate, phylogenetically weighted, if Phase 1 Step 3 shows
insufficient lineage-level longitudinal density to support ST/CC-level modeling.
**Decision rule:** the fallback is triggered only if fewer than [N — to be set from
Phase 1 lineage-density plots] major lineages show ≥ [Y — to be set] years of
longitudinal coverage with observable acquisition events.

## 6. Candidate Temporal Validation Design

Three-window split by collection year:

```
|---------------|------------|------------------|
     TRAIN          VALIDATE           TEST
```

**Decision rule for setting T1 (train/validate boundary) and T2 (validate/test
boundary):** chosen as the earliest boundaries such that each window contains at least
[threshold — to be set from Phase 1 Step 3's acquisition-event-over-time plot] candidate
acquisition events across the primary character set. Exact years are NOT to be chosen
before this plot exists (i.e., not arbitrarily 2022/2023).

**Hard rule (non-negotiable, applies throughout):** no information from after a given
timepoint is used to predict transitions occurring before or at that timepoint.

**Additional holdouts (robustness suite, Phase 9):** lineage holdout, geographic
holdout, repeated subsampling — run in addition to, not instead of, the temporal split.

## 7. Primary Evaluation Metrics

- AUROC and AUPRC for the primary "acquires complementary trait within test window"
  classification task
- Sensitivity at a pre-defined false-positive rate (threshold to be fixed *before*
  looking at test-window results)
- Calibration (calibration curve / Brier score)
- Rank-based accuracy for "which character is acquired next," following the reference
  HyperTraPS validation convention (predicted rank of the true next acquisition vs. a
  prevalence-only null model)
- All of the above reported broken down by major ST and by GBD-style geographic
  superregion, not only in aggregate

## 8. Explicit Exclusion Criteria

- Genomes with collection-date quality classified as "missing" (Phase 1 Step 2 tiers)
- Genomes with collection-date quality classified as "inferred/uncertain," **unless**
  Phase 1 Step 3 shows this exclusion drops usable data below a viable threshold — in
  which case inclusion is permitted only as a documented, explicit sensitivity arm, never
  silently merged into the primary analysis
- Genomes failing genome-quality QC (candidate threshold, following reference
  methodology: >20% unaligned bases against the *K. pneumoniae* reference assembly —
  to be confirmed during Phase 1 Step 7 reproduction)
- Genomes not confirmed *K. pneumoniae* sensu stricto (species complex members excluded)
- Genomes lacking either AMR **or** virulence calls (both are required for a genome to
  be usable in the joint model)
- Countries/regions with fewer than [threshold — to be set from Phase 1 Step 3] total
  qualifying genomes (mirrors the reference paper's exclusion of countries with <2
  genomes, adjusted upward once real density is known)

---

## Sign-off

- [ ] Reviewed against Phase 1 Step 1 checklist
- [ ] Committed to version control with timestamp: ________________
- [ ] Mentor review scheduled for acquisition-event definition (Section 3) and temporal
      design (Section 6): ________________
