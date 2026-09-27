import inspect
import unittest
import numpy as np
from invariantcell.identifiability_v10 import (MixtureSet, ModeledAssay, InfeasibleEvidence,
    outcome_bounds, worst_case_remaining_width, reveal, choose_assay)


class IdentifiabilityTests(unittest.TestCase):
    def simplex(self, n=2):
        return MixtureSet(np.empty((0,n)), np.empty(0), np.empty(0), np.ones(n,dtype=bool))

    def test_simplex_analytic_extrema_and_duality(self):
        r=outcome_bounds(self.simplex(3), [-2,4,1])
        self.assertEqual((r['lower'],r['upper']),(-2,4))
        for side in ('lower','upper'):
            self.assertLess(r[side+'_diagnostics']['absolute_primal_dual_gap'],1e-8)

    def test_interior_optima_not_vertex_enumeration(self):
        p=MixtureSet([[0,1]],[.2],[.3],[True,True])
        r=outcome_bounds(p,[0,10])
        self.assertAlmostEqual(r['lower'],2);self.assertAlmostEqual(r['upper'],3)

    def test_outcome_model_uncertainty_is_propagated(self):
        r=outcome_bounds(self.simplex(),[0,.8],[.2,1.2])
        self.assertAlmostEqual(r['upper'],1.2);self.assertEqual(r['decision'],'UNRESOLVED')
        self.assertFalse(r['confidence_claim'])

    def test_empty_evidence_cannot_certify(self):
        with self.assertRaises(InfeasibleEvidence):
            outcome_bounds(MixtureSet([[0,1]],[2],[3],[True,True]),[0,1])

    def test_qc_exclusions(self):
        r=outcome_bounds(MixtureSet(np.empty((0,2)),[],[],[True,False]),[.2,1])
        self.assertAlmostEqual(r['upper'],.2);self.assertEqual(r['decision'],'NEGATIVE')

    def test_threshold_contact_abstains(self):
        self.assertEqual(outcome_bounds(self.simplex(),[.5,.5])['decision'],'UNRESOLVED')

    def test_exact_two_copy_noise_oracle(self):
        for eps in (0,.05,.2,1):
            a=ModeledAssay('assay',(0.,1.),eps)
            r=worst_case_remaining_width(self.simplex(),a,[0,1])
            self.assertAlmostEqual(r['width'],min(1,2*eps))

    def test_uninformative_assay_does_not_shrink(self):
        a=ModeledAssay('constant',(1.,1.),0)
        self.assertAlmostEqual(worst_case_remaining_width(self.simplex(),a,[0,1])['width'],1)

    def test_reveal_narrows_and_noise_contains_truth(self):
        p=self.simplex();a=ModeledAssay('modeled',(0.,1.),.1)
        r=outcome_bounds(reveal(p,a,.7),[0,1])
        self.assertAlmostEqual(r['lower'],.6);self.assertAlmostEqual(r['upper'],.8)
        self.assertEqual(r['decision'],'POSITIVE')

    def test_final_endpoint_and_observed_assay_rejected(self):
        for a in (ModeledAssay('endpoint',(0,1),0,role='final_endpoint'),
                  ModeledAssay('observed',(0,1),0,evidence_kind='OBSERVED')):
            with self.assertRaises(ValueError): choose_assay(self.simplex(),[a],[0,1],policy='fixed',budget=1)
        self.assertNotIn('truth',inspect.signature(choose_assay).parameters)

    def test_equal_access_budget_and_decision_criterion(self):
        p=self.simplex(4);f=[0,0,1,1]
        assays=[ModeledAssay('nuisance',(0,10,0,10),0),ModeledAssay('decision',(0,.1,.9,1),.05)]
        self.assertEqual(choose_assay(p,assays,f,policy='uncertainty',budget=1),'nuisance')
        self.assertEqual(choose_assay(p,assays,f,policy='worst_case_decision',budget=1),'decision')
        self.assertIsNone(choose_assay(p,assays,f,policy='fixed',budget=.9))
        self.assertEqual(choose_assay(p,assays,f,policy='fixed',budget=1,already=('nuisance',)),'decision')

    def test_validation_and_no_silent_axis_broadcast(self):
        with self.assertRaises(ValueError):MixtureSet([[0,1]],[0],[1],[1,1])
        with self.assertRaises(ValueError):outcome_bounds(self.simplex(),[0])
        with self.assertRaises(ValueError):outcome_bounds(self.simplex(),[1,1],[0,0])
        with self.assertRaises(ValueError):ModeledAssay('bad',(0,1),-1).validate(2)


if __name__=='__main__':unittest.main()
