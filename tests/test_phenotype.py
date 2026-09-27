from dataclasses import replace
import unittest

from invariantcell.io import digest
from invariantcell.phenotype import (PathologyImageAdapter, PhenotypeRegion, as_of,
                                    measurement_decision, measurement_snapshot)


def fixture():
    image = [[1, 2], [3, 4]]
    mask = [[1, 1], [0, 2]]
    manifest = dict(measurement_id='SYNTHETIC:mask:v1', subject_id='SYNTHETIC:subject',
        context_id='SYNTHETIC:cell_image', organism='synthetic', image_ref='SYNTHETIC:image',
        mask_ref='SYNTHETIC:mask', image_content_sha256=digest(image), mask_content_sha256=digest(mask),
        source_kind='SYNTHETIC', annotation_kind='synthetic', annotation_version='v1',
        measurement_version='v1', area_support='supported', count_support='unsupported',
        instance_labels=True, uncertainty=None, quality_flags=['possible_merge'],
        observed_at='2026-01-01T00:00:00Z', public_available_at='2026-01-02T00:00:00Z',
        ingested_at='2026-01-02T00:00:00Z')
    return image, mask, manifest


class PhenotypeTests(unittest.TestCase):
    def test_prediction_and_synthetic_never_promoted(self):
        image, mask, m = fixture()
        for kind in ('PREDICTED', 'SYNTHETIC', 'INFERRED'):
            r = PathologyImageAdapter().measure(image, mask, dict(m, source_kind=kind))
            s = measurement_snapshot(r)
            self.assertEqual(s.evidence[0].state_kind, 'predicted')
            self.assertEqual(s.evidence[0].synthetic, kind == 'SYNTHETIC')
        r = PathologyImageAdapter().measure(image, mask, dict(m, source_kind='OBSERVED', annotation_kind='model'))
        self.assertEqual(r.record_kind, 'INFERRED')

    def test_count_and_area_support_separate(self):
        r = PathologyImageAdapter().measure(*fixture())
        self.assertEqual(r.phenotype_features[0].value, .75)
        self.assertEqual(r.phenotype_features[1].value, 2)
        self.assertEqual(measurement_decision(r, 'count'), 'ABSTAIN_UNSUPPORTED_QUANTITY')
        self.assertEqual(measurement_decision(r, 'area'), 'DEFER_DERIVED_OR_SYNTHETIC_MEASUREMENT')

    def test_semantic_labels_not_claimed_instances(self):
        image, mask, m = fixture(); m.update(count_support='supported', instance_labels=False)
        r = PathologyImageAdapter().measure(image, mask, m)
        self.assertEqual(measurement_decision(r, 'count'), 'ABSTAIN_UNSUPPORTED_QUANTITY')

    def test_region_hierarchy_rejects_dangling_and_cycles(self):
        r = PathologyImageAdapter().measure(*fixture())
        for regions in [(PhenotypeRegion('a','b','cell'),),
                        (PhenotypeRegion('a','b','cell'),PhenotypeRegion('b','a','cell'))]:
            with self.assertRaises(ValueError): replace(r, regions=regions)

    def test_annotation_update_preserves_event_and_history(self):
        image, mask, m = fixture(); old = PathologyImageAdapter().measure(image, mask, m)
        newmask = [[1,1],[2,2]]
        new = PathologyImageAdapter().measure(image, newmask, dict(m,
            measurement_id='SYNTHETIC:mask:v2', annotation_version='v2', measurement_version='v2',
            mask_content_sha256=digest(newmask), public_available_at='2026-02-01T00:00:00Z',
            ingested_at='2026-02-01T00:00:00Z'))
        self.assertEqual(old.observed_at, new.observed_at)
        self.assertNotEqual(old.annotation_sha256, new.annotation_sha256)
        self.assertEqual(as_of([old,new], '2026-01-20T00:00:00Z', '2026-01-20T00:00:00Z'), (old,))
        self.assertEqual(old.annotation_version, 'v1')
        with self.assertRaises(ValueError): as_of([old,old], m['public_available_at'], m['ingested_at'])

    def test_provenance_replay_rejects_tampering(self):
        image, mask, m = fixture(); adapter = PathologyImageAdapter(); r = adapter.measure(image,mask,m)
        cert = adapter.certificate(r,m)
        self.assertTrue(adapter.replay(cert,r,image,mask,m,digest(m)))
        self.assertFalse(adapter.replay(cert,replace(r,record_kind='OBSERVED'),image,mask,m,digest(m)))
        forged = dict(m, source_kind='OBSERVED')
        self.assertFalse(adapter.replay(cert,r,image,mask,forged,digest(m)))
        with self.assertRaises(ValueError): adapter.measure(image,[[0,0],[0,0]],m)

    def test_unknown_time_not_invented(self):
        r = PathologyImageAdapter().measure(*fixture())
        unknown = replace(r, observed_at=None, public_available_at=None)
        self.assertEqual(as_of([unknown], '2027-01-01T00:00:00Z', '2027-01-01T00:00:00Z'), ())
        self.assertIsNone(measurement_snapshot(unknown).evidence[0].event_time)

    def test_missing_annotation_version_rejected(self):
        r = PathologyImageAdapter().measure(*fixture())
        with self.assertRaises(ValueError): replace(r, annotation_version=None)


if __name__ == '__main__': unittest.main()
