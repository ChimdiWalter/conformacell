# Molecular evidence bridge contract

Status: original synthetic engineering implementation inside InvariantCell. The actual ConformaCell PoC has not been located in the reviewed scope. This contract does not claim to reproduce that source or its reported 1CRN example.

`src/invariantcell/molecular.py` adds `MolecularConstruct`, `ResidueMapping`, `MolecularIntervention`, `ConformerEvidence`, `ExperimentalConstraint`, and `MolecularStateEnsemble`. It uses existing `DependencyRef`, `EvidenceRecord`, `Snapshot`, `PhenotypeMeasurement`, and temporal certificate/replay interfaces. It imports no original research project and performs no network access, training, or structural simulation.

## Identity and sequence

A construct stores immutable identity, protein, organism, context, explicit canonical amino-acid sequence, hashed sequence source, versioned residue mapping, and synthetic status. Mapping records distinguish a one-based sequence position from chain, author residue number, and insertion code. Missing positions may remain unmapped; a partial map is never represented as a complete structure. Duplicate structural residues or repeated sequence positions within a chain are rejected. Multiple chains may map to the same sequence position. Noncanonical residues, insertions/deletions, multi-substitution interventions, post-translational modifications, ligands and heteromers are outside this minimal contract.

A supported intervention is one explicit substitution between distinct WT/variant construct IDs with the same protein/context/organism/synthetic status. The reference residue must match the WT sequence, and the complete variant sequence must equal the declared substitution. Every ensemble conformer must bind to the exact construct and sequence digest.

`coordinate_hash` hashes ordered numeric XYZ values independent of residue labels. The same coordinate hash for WT and variant yields `SAME_COORDINATES_NOT_MUTATION_AWARE`, even if the variant label and sequence have changed. Different coordinates and a declared sequence-sensitive method only produce `MUTATION_SPECIFIC_INPUTS_DECLARED_NOT_VALIDATED`. Noise, an ANM coordinate perturbation, or a self-declared method flag does not validate mutation sensitivity. Coordinate units, atom ordering, mapping completeness, and physical method behavior require independent source review. Exact hashes do not detect structures related by rotation, translation, or atom reordering.

## Record semantics and validation

Each conformer and constraint retains `OBSERVED`, `PREDICTED`, `INFERRED`, or `SYNTHETIC` declarations. A computed ensemble cannot be `OBSERVED`; synthetic input ancestry must remain synthetic. Weights are finite, nonnegative and normalized, and uncertainty is either nonnegative finite or explicitly unknown. These are data-shape checks, not calibrated probability or uncertainty claims.

A constraint identifies its original acquisition and endpoint, construct/context and hashed source. `validate_independent_endpoints` rejects fitting/validation overlap by acquisition-plus-endpoint, source identity, or source hash; changing a filename or transformed-data hash alone does not hide a retained original acquisition ID. Reuse of a common source file is conservatively rejected even if a caller claims different endpoints. Truly independent acquisitions require reviewed provenance. The validator does not establish statistical independence, solve alias/entity resolution, or detect deliberately falsified declarations. Nonobserved validation remains engineering-only. Unknown original acquisition identity must not be replaced with invented real metadata.

Real evaluation is NOT RUN. Before any real promotion, freeze a matched-budget comparison including sequence only, sequence plus static structure, sequence plus WT flexibility, a mutation-sensitive ensemble, and direct observable features. NMA/BioEmu/pooled or reweighted variants are candidates only if actually available and permitted. The direct-observable comparator cannot be omitted from an ensemble-value claim. Select a genuinely independent molecular endpoint; a fitted SAXS curve cannot also be the independent validation measurement.

## Shared evidence and phenotype gate

`molecular_snapshot` creates only `has_molecular_ensemble_record`, with the full ensemble digest and source dependencies. It always uses predicted evidence semantics, preserving synthetic status. It does not assert an observed molecular effect, a cell response, a mechanism, or causality. Event time and biological validity remain unknown unless supplied separately; fictional demo availability dates do not become experiment dates.

`phenotype_bridge_status` checks declared intervention, construct, context and organism matching to the existing phenotype record. Unmatched records abstain. Matched synthetic records return `SYNTHETIC_BRIDGE_ENGINEERING_ONLY`; predicted/inferred phenotypes cannot stand in for observations. Even declared matched observations return `MATCHED_DECLARATIONS_CAUSALITY_NOT_ESTABLISHED`. A shared protein/gene identifier, a maize phenotype, and an unrelated human context cannot establish a molecular-to-cell bridge. The measurement's intervention/construct link is an explicit caller declaration requiring upstream provenance review, not an independently verified join.

A mechanistic mediator model, real cellular prediction bridge and reveal-only assay-selection benchmark are NOT IMPLEMENTED. No causal predicate or automatic causal promotion is exposed. Future work needs matched intervention/context/assay records, controls, confounding assessment, independent validation and a frozen decision protocol.

## Reproduction and limits

Run from this repository with the already installed Python 3.12.3 environment:

```bash
source ~/.venvs/lesegenv/bin/activate
PYTHONPATH=src python -B -m unittest discover -s tests -p 'test_molecular.py' -v
PYTHONPATH=src python -B scripts/molecular_demo.py
```

The demo emits `reports/conformacell/synthetic_demo.json` and `receipt.json`. Eighteen engineering tests cover sequence/residue mapping, malformed weights, provenance, fitting/validation reuse, relabeling, synthetic/observed separation, context matching, snapshot roundtrip and independent certificate replay/tampering. These invented records are not biological data. Some negative tests deliberately exercise declared non-synthetic/observed metadata; they do not establish that any real observation was loaded.

The trusted-input boundary is explicit: typed declarations and hashes must be verified against source assets before real ingestion. These checks are software safeguards, not an OS security boundary or scientific proof. No source rights, commercial rights, NIH resource eligibility, or public disclosure clearance follows from passing tests.
