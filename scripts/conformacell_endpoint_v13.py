#!/usr/bin/env python3
"""Metadata-only endpoint admission receipt; never fits or reads endpoint tables.

Extends existing molecular endpoint-use records with conservative declared-use
checks. Passing declarations alone cannot establish scientific independence.
"""
import argparse
import hashlib
import json
from pathlib import Path

PROJECT = Path(__file__).resolve().parents[1]
MANIFEST = 'data/post_confirmation_v13/conformacell/endpoint_manifest.json'
MATCH_FIELDS = ('sequence_definition_matched', 'chemical_construct_matched',
                'condition_compatibility_resolved', 'sample_lineage_matched',
                'quantitative_endpoint_and_uncertainty_available',
                'measurement_acquisition_lineage_resolved', 'evaluation_protocol_frozen')
DISALLOWED_USES = {'generation', 'fitting', 'selection', 'interpretation'}


def assess_endpoint(record):
    reasons = ['NOT_CONFIRMED_' + name.upper() for name in MATCH_FIELDS if record.get(name) is not True]
    if record.get('prior_use_review_complete') is not True:
        reasons.append('PRIOR_USE_UNKNOWN')
    if record.get('outcome_values_exposed') is not False:
        reasons.append('OUTCOME_VALUES_EXPOSED_OR_UNKNOWN')
    uses = record.get('dependency_roles')
    if not isinstance(uses, list):
        reasons.append('DEPENDENCY_ROLES_MISSING')
        uses = []
    for use in uses:
        if (not isinstance(use, dict) or not use.get('acquisition_id')
                or not use.get('endpoint_family') or not use.get('source_id')
                or not use.get('role') or not use.get('model_id')):
            reasons.append('INCOMPLETE_USE_EDGE')
            continue
        same_measurement = (use['acquisition_id'] == record.get('acquisition_id')
                            or use['endpoint_family'] == record.get('endpoint_family')
                            or use['source_id'] == record.get('endpoint_source_id'))
        if same_measurement and use['role'] in DISALLOWED_USES:
            reasons.append('ENDPOINT_USED_FOR_' + use['role'].upper())
        if use['role'] not in DISALLOWED_USES | {'heldout_validation', 'metadata_only', 'unknown'}:
            reasons.append('UNRECOGNIZED_USE_ROLE')
        if same_measurement and use['role'] == 'unknown':
            reasons.append('PRIOR_USE_UNKNOWN')
    if not all(record.get(k) for k in ('candidate_id','acquisition_id','endpoint_family','endpoint_source_id')):
        reasons.append('ENDPOINT_IDENTITY_MISSING')
    return {'candidate_id': record.get('candidate_id'),
            'status': 'NOT_ADMITTED' if reasons else 'DECLARATIONS_PASS_REQUIRES_INDEPENDENT_REVIEW',
            'reasons': sorted(set(reasons)),
            'scientific_independence_certified': False,
            'endpoint_access_authorized': False}


def check(project=PROJECT):
    project=Path(project).resolve()
    manifest=json.loads((project/MANIFEST).read_text())
    if manifest.get('schema_version') != 1 or manifest.get('boundary') != 'post_confirmation_v13':
        raise ValueError('Wrong endpoint metadata boundary')
    if manifest.get('contains_endpoint_values') is not False:
        raise ValueError('Metadata-only manifest required')
    bindings=manifest['source_sha256']
    if not bindings:
        raise ValueError('Empty source bindings')
    for name, expected in bindings.items():
        path=(project/name).resolve()
        if Path(name).is_absolute() or '..' in Path(name).parts or not path.is_relative_to(project):
            raise ValueError('Source path outside project')
        if path.suffix.lower() in {'.parquet','.h5ad','.npy','.npz','.h5'}:
            raise ValueError('Outcome container cannot be opened by metadata checker')
        if hashlib.sha256(path.read_bytes()).hexdigest() != expected:
            raise ValueError('Changed metadata dependency: '+name)
    pten=manifest['PTEN_schema_only_intake']
    receipt=json.loads((project/pten['receipt_path']).read_text())
    schema=json.loads((project/pten['schema_path']).read_text())
    if pten['receipt_path'] not in bindings or pten['schema_path'] not in bindings:
        raise ValueError('PTEN receipt/schema must be hash-bound')
    if receipt.get('rows_read') is not False or receipt.get('identity_and_split_ready') is not False:
        raise ValueError('PTEN boundary changed; independent review required')
    if len(schema)!=3628 or schema.get('aaChanges')!='String' or schema.get('Metadata_Object_Count')!='Int64':
        raise ValueError('PTEN schema no longer matches reviewed metadata')
    if any(name in schema for name in pten['missing_identity_columns']):
        raise ValueError('PTEN identity schema changed; review before opening outcomes')
    decisions=[assess_endpoint(row) for row in manifest['endpoint_candidates']]
    return {'status':'METADATA_CHECK_COMPLETE_NO_REAL_BENCHMARK_ADMITTED',
            'candidate_assessments':decisions,'metadata_dependencies_verified':len(bindings),
            'PTEN_columns':len(schema),'PTEN_rows_read_by_checker':False,
            'endpoint_values_read_by_checker':False,'fits':0,
            'source_consistency_only_not_scientific_certification':True,
            'next_command':'PYTHONPATH=src ~/.venvs/lesegenv/bin/python -B scripts/conformacell_endpoint_v13.py --check'}


def main():
    p=argparse.ArgumentParser();p.add_argument('--check',action='store_true',required=True)
    p.parse_args();print(json.dumps(check(),indent=2))


if __name__=='__main__':main()
