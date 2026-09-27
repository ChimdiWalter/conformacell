"""Engineering assertions over invented records; no biological validation data."""
from dataclasses import replace
import unittest

from invariantcell.core import DependencyRef
from invariantcell.io import digest
from invariantcell.molecular import (
    ConformerEvidence, ExperimentalConstraint, MolecularConstruct,
    MolecularIntervention, MolecularStateEnsemble, ResidueMapping,
    coordinate_hash, molecular_snapshot, mutation_sensitivity,
    phenotype_bridge_status, validate_independent_endpoints,
)
from invariantcell.phenotype import PhenotypeFeature, PhenotypeMeasurement
from invariantcell.replay import verify_certificate
from invariantcell.temporal import Snapshot, TemporalQuery, issue_certificate


def reference(name, value):
    return DependencyRef(name, digest(value), "2026-01-01T00:00:00Z")


def fixture(synthetic=True):
    def construct(name, sequence):
        mapping = tuple(ResidueMapping("A", 100 + i, "", i, aa)
                        for i, aa in enumerate(sequence, 1))
        return MolecularConstruct(name, "SYNTHETIC:protein", "SYNTHETIC:context", "synthetic",
            sequence, reference("sequence:" + name, sequence), "fixture-v1", mapping, synthetic)
    wt, variant = construct("SYNTHETIC:WT", "ACDE"), construct("SYNTHETIC:C2V", "AVDE")
    intervention = MolecularIntervention("SYNTHETIC:mutation", wt, variant, 2, "C", "V")
    coordinates = [[0, 0, 0], [1, 1, 1], [2, 2, 2], [3, 3, 3]]
    kind = "SYNTHETIC" if synthetic else "PREDICTED"
    conformer = ConformerEvidence("SYNTHETIC:conformer", variant.construct_id,
        variant.sequence_sha256, coordinate_hash(coordinates), reference("coordinates:variant", coordinates),
        "fixture-no-physical-model", "v1", kind, True)
    ensemble = MolecularStateEnsemble("SYNTHETIC:ensemble", variant, (conformer,), (1.0,), kind,
        "fixture-weighting", "v1", public_available_at="2026-01-02T00:00:00Z",
        ingested_at="2026-01-02T00:00:00Z")
    return intervention, conformer, ensemble


def measurement(kind="SYNTHETIC"):
    return PhenotypeMeasurement("SYNTHETIC:measurement", "SYNTHETIC:subject", "SYNTHETIC:context",
        "synthetic", "fixture", "SYNTHETIC:assay", digest("fixture-assay"), None, None, None,
        (PhenotypeFeature("fixture", 1.0, "arbitrary", "morphology"),), (), "fixture", "v1", kind,
        None, (), None, None, None, ("Invented engineering fixture; not biological evidence.",))


class MolecularTests(unittest.TestCase):
    def test_wt_variant_sequence_identity(self):
        intervention, _, _ = fixture()
        self.assertEqual(intervention.variant.sequence, "AVDE")
        with self.assertRaisesRegex(ValueError, "reference residue"):
            replace(intervention, reference_amino_acid="A")
        with self.assertRaisesRegex(ValueError, "declared substitution"):
            replace(intervention, alternate_amino_acid="W")
        with self.assertRaises(ValueError):
            replace(intervention, position=True)

    def test_mapping_author_number_is_not_sequence_position(self):
        intervention, _, _ = fixture()
        mapped = replace(intervention.reference, residue_mapping=(
            ResidueMapping("A", 10, "", 1, "A"), ResidueMapping("A", 10, "A", 2, "C")))
        self.assertEqual(mapped.residue_mapping[1].sequence_position, 2)
        with self.assertRaisesRegex(ValueError, "disagrees"):
            replace(mapped, residue_mapping=(ResidueMapping("A", 10, "", 2, "A"),))

    def test_duplicate_mapping_and_out_of_range_rejected(self):
        intervention, _, _ = fixture()
        row = ResidueMapping("A", 101, "", 1, "A")
        for mapping in ((row, row), (ResidueMapping("A", 999, "", 5, "A"),)):
            with self.assertRaises(ValueError):
                replace(intervention.reference, residue_mapping=mapping)

    def test_mapping_same_sequence_position_in_distinct_chains_allowed(self):
        intervention, _, _ = fixture()
        c = replace(intervention.reference, residue_mapping=(ResidueMapping("A", 1, "", 1, "A"),
            ResidueMapping("B", 1, "", 1, "A")))
        self.assertEqual(len(c.residue_mapping), 2)

    def test_same_structure_relabel_not_mutation_aware(self):
        intervention, variant, _ = fixture()
        wt = replace(variant, conformer_id="SYNTHETIC:WT-conformer",
            construct_id=intervention.reference.construct_id,
            sequence_sha256=intervention.reference.sequence_sha256)
        self.assertEqual(mutation_sensitivity(intervention, wt, variant),
                         "SAME_COORDINATES_NOT_MUTATION_AWARE")

    def test_changed_coordinates_alone_do_not_prove_mutation_awareness(self):
        intervention, variant, _ = fixture()
        wt = replace(variant, construct_id=intervention.reference.construct_id,
            sequence_sha256=intervention.reference.sequence_sha256, coordinate_sha256=digest("other-xyz"))
        self.assertEqual(mutation_sensitivity(intervention, wt, replace(variant, uses_variant_sequence=False)),
                         "MUTATION_SENSITIVITY_UNSUPPORTED")
        self.assertEqual(mutation_sensitivity(intervention, wt, variant),
                         "MUTATION_SPECIFIC_INPUTS_DECLARED_NOT_VALIDATED")

    def test_coordinate_hash_has_numeric_normalization_and_rejects_nonfinite(self):
        self.assertEqual(coordinate_hash([[1, 2, 3]]), coordinate_hash([[1.0, 2.0, 3.0]]))
        for coords in ([], [[1, 2]], [[1, float("nan"), 3]], [[True, 2, 3]]):
            with self.assertRaises(ValueError): coordinate_hash(coords)

    def test_prediction_uses_correct_construct(self):
        intervention, conformer, ensemble = fixture()
        with self.assertRaisesRegex(ValueError, "construct sequence"):
            replace(ensemble, conformers=(replace(conformer, construct_id=intervention.reference.construct_id),))
        with self.assertRaises(ValueError):
            replace(ensemble, conformers=(replace(conformer, sequence_sha256=digest("different-sequence")),))

    def test_predicted_inferred_synthetic_never_promoted_to_observed(self):
        for synthetic in (True, False):
            _, _, ensemble = fixture(synthetic)
            self.assertEqual(molecular_snapshot(ensemble).evidence[0].state_kind, "predicted")
            self.assertEqual(molecular_snapshot(ensemble).evidence[0].synthetic, synthetic)
            with self.assertRaises(ValueError): replace(ensemble, record_kind="OBSERVED")
        _, _, ensemble = fixture(False)
        self.assertEqual(molecular_snapshot(replace(ensemble, record_kind="INFERRED")).evidence[0].state_kind,
                         "predicted")

    def test_synthetic_provenance_cannot_be_removed(self):
        _, _, ensemble = fixture()
        with self.assertRaisesRegex(ValueError, "synthetic"):
            replace(ensemble, record_kind="PREDICTED")

    def test_normalized_weights_and_uncertainty(self):
        _, _, ensemble = fixture()
        for weights in ((-1.0,), (0.5,), (float("nan"),), (True,)):
            with self.assertRaises(ValueError): replace(ensemble, weights=weights)
        self.assertIsNone(ensemble.uncertainty)
        with self.assertRaises(ValueError): replace(ensemble, uncertainty=-1)

    def test_saxs_fitting_not_reused_as_independent_validation(self):
        _, _, ensemble = fixture()
        fit = ExperimentalConstraint("SYNTHETIC:acquisition", "SAXS:I(q)", ensemble.construct.construct_id,
            ensemble.construct.context_id, reference("SYNTHETIC:SAXS", "synthetic-curve"), "SYNTHETIC")
        for validation in (fit, replace(fit, source=reference("aliased-source", "transformed-curve")),
                           replace(fit, acquisition_id="alias", endpoint="renamed-endpoint")):
            with self.assertRaisesRegex(ValueError, "independent validation"):
                validate_independent_endpoints((fit,), (validation,))

    def test_disjoint_synthetic_endpoint_is_engineering_only(self):
        _, _, ensemble = fixture()
        validation = ExperimentalConstraint("SYNTHETIC:acquisition", "fixture-endpoint",
            ensemble.construct.construct_id, ensemble.construct.context_id,
            reference("SYNTHETIC:endpoint", "fixture"), "SYNTHETIC")
        self.assertEqual(validate_independent_endpoints((), (validation,)),
                         "NONOBSERVED_VALIDATION_ENGINEERING_ONLY")
        self.assertEqual(validate_independent_endpoints((), ()), "NO_INDEPENDENT_VALIDATION")
        with self.assertRaises(ValueError):
            replace(ensemble, independent_validation=(replace(validation, construct_id="different"),))

    def test_unknown_event_time_not_invented(self):
        _, _, ensemble = fixture()
        evidence = molecular_snapshot(ensemble).evidence[0]
        self.assertIsNone(evidence.event_time)
        self.assertIsNone(evidence.valid_from)
        self.assertEqual(evidence.validity_kind, "unknown")

    def test_dependency_and_snapshot_roundtrip(self):
        _, _, ensemble = fixture()
        snapshot = molecular_snapshot(ensemble)
        self.assertEqual(Snapshot.from_dict(snapshot.to_dict()), snapshot)
        self.assertEqual(len(snapshot.dependencies), 3)
        self.assertEqual(snapshot.evidence[0].predicate, "has_molecular_ensemble_record")

    def test_temporal_certificate_replays_and_tampering_fails(self):
        _, _, ensemble = fixture()
        snapshot = molecular_snapshot(ensemble)
        query = TemporalQuery(ensemble.construct.construct_id, "has_molecular_ensemble_record",
            ensemble.construct.context_id, "2026-01-03T00:00:00Z", "2026-01-03T00:00:00Z",
            "2026-01-03T00:00:00Z", allow_predictions=True)
        cert = issue_certificate(snapshot, query)
        self.assertTrue(verify_certificate(cert, snapshot, digest(snapshot.to_dict()), query).valid)
        altered = replace(snapshot, evidence=(replace(snapshot.evidence[0], state_kind="observed"),))
        self.assertFalse(verify_certificate(cert, altered, digest(snapshot.to_dict()), query).valid)

    def test_unrelated_phenotypes_do_not_form_bridge(self):
        intervention, _, ensemble = fixture()
        for change in ({"context_id": "unrelated-cell"}, {"organism": "Zea mays"}):
            status = phenotype_bridge_status(ensemble, intervention, replace(measurement(), **change),
                measurement_intervention_id=intervention.intervention_id,
                measurement_construct_id=ensemble.construct.construct_id)
            self.assertEqual(status, "ABSTAIN_UNMATCHED_INTERVENTION_CONSTRUCT_OR_CONTEXT")
        self.assertEqual(phenotype_bridge_status(ensemble, intervention, measurement(),
            measurement_intervention_id="same-protein-but-unmatched-intervention",
            measurement_construct_id=ensemble.construct.construct_id),
            "ABSTAIN_UNMATCHED_INTERVENTION_CONSTRUCT_OR_CONTEXT")

    def test_matched_ids_do_not_certify_causality(self):
        intervention, _, ensemble = fixture(False)
        status = phenotype_bridge_status(ensemble, intervention, measurement("OBSERVED"),
            measurement_intervention_id=intervention.intervention_id,
            measurement_construct_id=ensemble.construct.construct_id)
        self.assertEqual(status, "MATCHED_DECLARATIONS_CAUSALITY_NOT_ESTABLISHED")
        self.assertEqual(phenotype_bridge_status(ensemble, intervention, measurement("PREDICTED"),
            measurement_intervention_id=intervention.intervention_id,
            measurement_construct_id=ensemble.construct.construct_id), "ABSTAIN_NO_OBSERVED_PHENOTYPE")


if __name__ == "__main__": unittest.main()
