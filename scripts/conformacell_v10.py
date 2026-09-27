#!/usr/bin/env python3
"""New v10 metadata audit and explicitly synthetic finite-library pilot.

--check does no model fit or LP solve. --run requires its pinned task plan and
an unused output directory. It never opens a cell/clinical response store.
"""
import argparse
import csv
import hashlib
import json
import os
from pathlib import Path
import shlex
import time

for variable in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS'):
    os.environ[variable]='1'

P=Path(__file__).resolve().parents[1]
D=P/'data/post_confirmation_v10/conformacell'
R=P/'reports/conformacell'
POLICIES=['random','cheapest','fixed','uncertainty','worst_case_decision']
CODE=['scripts/conformacell_v10.py','src/invariantcell/identifiability_v10.py']


def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def read(path):return json.loads(Path(path).read_text())
def write(path,value):Path(path).write_text(json.dumps(value,indent=2,sort_keys=True)+'\n')


def scalar_metadata(path):
    result={}
    # Read deposition metadata only, stop before intensity/model blocks.
    with Path(path).open() as f:
        for line in f:
            if line.startswith('_sas_scan_intensity.'):break
            if line.startswith('_'):
                values=shlex.split(line)
                if len(values)==2:result[values[0]]=values[1]
    return result


def verify_sources():
    sources={}
    manifest=P/'data/public/molecular/download_manifest.json'
    for row in read(manifest):
        if Path(row['path']).name in {'SASDAU4.dat','SASDAT4.dat','SASDAU4.html','SASDAT4.html','SASBDB_rights.html'}:
            path=P/row['path']
            if sha(path)!=row['sha256']:raise ValueError('Changed cached public input: '+str(path))
            sources[str(path.relative_to(P))]=row['sha256']
    previous=P/'reports/conformacell/construct_review.json';r=read(previous)
    sources[str(previous.relative_to(P))]=sha(previous)
    for name,value in r['source_hashes'].items():
        path=P/'data/public/conformacell_construct_review'/name
        if sha(path)!=value:raise ValueError('Changed resolved construct source: '+name)
        sources[str(path.relative_to(P))]=value
    acquisition=read(D/'acquisition.json')
    if acquisition['status']!='DOWNLOADED' or acquisition['bytes']>20*1024*1024:
        raise ValueError('Verified small primary deposition metadata unavailable')
    sascif=P/acquisition['path']
    if sascif.resolve()!=D/'SASDAU4.sascif' or sha(sascif)!=acquisition['sha256']:
        raise ValueError('Unexpected new asset binding')
    sources[str(sascif.relative_to(P))]=acquisition['sha256']
    meta=scalar_metadata(sascif)
    expected={'_sas_result.SASBDB_code':'SASDAU4','_sas_scan.title':'uPAR H47C/N259C',
              '_struct_ref.db_code':'Q03405','_sas_scan.cell_temperature':'10.0',
              '_sas_buffer.concentration':'25.00'}
    if any(meta.get(k)!=v for k,v in expected.items()):raise ValueError('Deposition metadata changed')
    # Reuse the resolved sequence result; do not infer a mutant from archive FASTA.
    resolved=r['sequence_resolution']
    if resolved['mature_sequence_length']!=283 or [(x['author_position'],x['reference'],x['alternate']) for x in resolved['mutations']]!=[(47,'H','C'),(259,'N','C')]:
        raise ValueError('Resolved mapping no longer matches prior review')
    return r,meta,sources,acquisition


def audit():
    D.mkdir(parents=True,exist_ok=True)
    r,meta,sources,acquisition=verify_sources()
    paper='data/public/conformacell_construct_review/PMC3464537.html'
    dep='data/post_confirmation_v10/conformacell/SASDAU4.sascif'
    archive='data/public/molecular/SASDAU4.html'
    rows=[]
    def add(subject,field,value,status,source,locator,use,limit):
        rows.append(dict(subject=subject,field=field,value=value,status=status,source_path=source,
                         source_sha256=sources[source],source_locator=locator,endpoint_use=use,limitation=limit))
    add('SASDAU4','construct','uPAR H47C/N259C; Q03405; mature 283 aa mapping reused','RESOLVED_POLYPEPTIDE',dep,'_sas_scan.title; _struct_ref.db_code','identity','No exact sample chemistry/batch claim; precursor archive sequence not mutant construct')
    for field,tag in [('temperature_C','_sas_scan.cell_temperature'),('phosphate_mM','_sas_buffer.concentration'),('buffer','_sas_buffer.composition'),('pH','_sas_buffer.pH'),('concentration_mg_ml','_sas_sample.specimen_concentration'),('concentration_method','_sas_sample.concentration_method'),('instrument','_sas_beam.instrument_name'),('detector','_sas_detc.name'),('measurement_date','_sas_scan.measurement_date')]:
        add('SASDAU4',field,meta[tag],'ARCHIVE_DECLARATION',dep,tag,'SAXS conditioning metadata','Archive declaration, not a corrected paper value or newly measured result')
    add('uPAR WT/mutant SAXS','temperature_C','25','CONFLICT',paper,'SAXS Data Collection','conditioning','Archive declares 10 C; raw acquisition log absent')
    add('uPAR WT/mutant SAXS','phosphate_mM','20','CONFLICT',paper,'SAXS Data Collection','conditioning','Archive declares 25 mM; no primary adjudication acquired')
    add('uPAR WT/mutant SAXS','buffer','pH 7.2; 50 mM NaSO4; 5% w/v glycerol','AGREEMENT_WITH_UNIT_QUALIFIER',paper,'SAXS Data Collection','conditioning','Archive gives glycerol percent without w/v basis')
    add('uPAR WT/mutant','preparation','S2 secreted 1–283; affinity purified; SEC; SAXS within 48 h','DECLARED',paper,'Experimental Procedures / Materials and SAXS Data Collection','sampling','No cross-assay batch identifiers')
    add('uPAR WT/mutant','glycosylation','Glycosylated preparations; reported masses 35.1–35.4 kDa','PARTIAL',paper,'Table 1 legend','interpretation','No site-resolved glycoform matching; no chemical identity inference')
    add('SASDAU4','ligand','Entry names unliganded double mutant; distinct liganded records exist','DECLARED_NOT_CHEMICALLY_PROVEN',archive,'Title and study collection','conditioning','No assay-specific ligand-occupancy assay acquired')
    add('uPAR HDX','conditions','33 uM starting PBS; 20-fold D2O-PBS dilution; 25 C; 10/100/1000 s; duplicates','DECLARED',paper,'Hydrogen Deuterium Measurements','structure selection / interpretation','17 peptides/44% coverage; no machine-readable uptake replicates/raw spectra acquired')
    add('uPAR HDX','endpoint_independence','HDX discriminated among SAXS-compatible domain models','NOT_INDEPENDENT_FOR_SELECTED_MODEL',paper,'HDX and Discussion','selection','Orthogonal assay name does not remove selection dependence')
    add('uPAR SPR','conditions','20 C; HEPES/NaCl/EDTA/P20; DI1–92 + DIIDIII88–283 and ATF comparison','UNMATCHED_ENDPOINT',paper,'Surface Plasmon Resonance; Figure 5','interpretation','Not exact WT/H47C-N259C paired endpoint')
    add('SASDAU4','curve','Processed merged q/I/sigma with 2339 rows; 84 negative intensities retained','MEASURED_PROCESSED',dep,'_sas_result.type_of_curve; _sas_scan_intensity','fitting','Not raw detector frames; sigma does not specify covariance or confidence threshold')
    fields=list(rows[0])
    with (R/'condition_and_independence_ledger.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=fields);w.writeheader();w.writerows(rows)
    graph={'schema_version':1,'extension_of':'molecular.ExperimentalConstraint + construct_review orthogonal_endpoints',
           'nodes':{'SAXS_processed':'observed processed curve','3U74':'inferred static crystal model',
                    'CORAL_models':'SAXS-fitted structural models','HDX_published':'exposed summarized measurements; raw missing',
                    'HDX_selected_interpretation':'SAXS+HDX selected interpretation','SPR_published':'different construct comparison'},
           'edges':[{'from':'3U74','to':'CORAL_models','use':'sampling/structural template'},
                    {'from':'SAXS_processed','to':'CORAL_models','use':'fitting/selection'},
                    {'from':'CORAL_models','to':'HDX_selected_interpretation','use':'candidate interpretation'},
                    {'from':'HDX_published','to':'HDX_selected_interpretation','use':'selection/interpretation'}],
           'hyperparameter_independence':'unknown outside described publication procedure',
           'excluded_as_independent_validation':['SAXS fit residuals','SAXS-derived Rg/p(r)','HDX used for selection','B factors'],
           'prior_eligibility':r['orthogonal_endpoints'],'source_sha256':sources,
           'no_new_independence_claim':True}
    write(D/'evidence_use_graph.json',graph)
    plan={'schema_version':1,'boundary':'post_confirmation_v10','status':'READY_SYNTHETIC_ONLY',
          'source_sha256':sources,'code_sha256':{x:sha(P/x) for x in CODE},
          'task':'finite-library modeled-assay identifiability engineering pilot',
          'resource_caps':{'cpu_threads':1,'ram_bytes':2*1024**3,'wall_seconds':180,'gpu':False},
          'estimated_wall_seconds':15,'policies':POLICIES,'budgets':[0,1,2],'random_seed':0,
          'case_families':['decision-relevant and nuisance observables','uninformative observations'],
          'actual_assay_cost':None,'modeled_assay_cost':1.0,'final_endpoint_in_query_pool':False,
          'noise':'deterministic absolute-error intervals; no distribution or confidence level',
          'nuisance':'Synthetic units calibrated; fixed scale=1/background=0. Real nuisance scale/background/covariance unresolved.',
          'real_uPAR_fit':False,'new_public_download_bytes':acquisition['bytes'],
          'claim_limit':'Synthetic correctness and conditional numerical bounds; no real assay savings, global molecular certification, novelty or causal claim'}
    write(D/'task_plan.json',plan)
    readiness={'status':'CHECK_PASS_SYNTHETIC_READY_REAL_INDEPENDENT_ENDPOINT_PENDING',
               'mature_sequence_mapping':'REUSED_RESOLVED_283_H47C_N259C','SASDAU4_lead_verified':True,
               'new_metadata':meta,'source_sha256':sources,'condition_conflict_resolved':False,
               'independent_endpoint_available':False,'forward_conformer_library_available':False,
               'real_uPAR_decision_or_assay_policy_result':False,'LP_solved_by_check':False,
               'new_download_bytes':acquisition['bytes'],'task_plan_sha256':sha(D/'task_plan.json')}
    write(D/'readiness.json',readiness)
    (R/'observed_bridge_status.md').write_text('''# ConformaCell v10 observed bridge status

The SASDAU4 lead is verified against cached primary inputs and a newly acquired official SASCIF (648,583 bytes). It identifies the H47C/N259C uPAR sample and Q03405. The resolved 283-residue mapping is reused, not reconstructed from the shared precursor FASTA.

The structured deposition repeats **10 C / 25 mM phosphate**, reports measurement date 2010-06-10 and A280 concentration estimation. It does not resolve the paper's **25 C / 20 mM** declaration. Buffer additives, concentration range and instrument metadata are recorded alongside the disagreement in `condition_and_independence_ledger.csv`. Batch, exact glycans/oxidation and cross-assay chemical matching remain absent. Date and conditions are depositor declarations, not independently verified acquisition logs.

The evidence-use graph preserves SAXS fitting and HDX-guided interpretation. HDX is an orthogonal modality but not independent validation for the published HDX-selected models. No numerical peptide/time/replicate table, raw HDX spectra, exact mutant SPR endpoint or matched cell phenotype was acquired. Prior supplement 403/metadata failures remain preserved; no restriction was bypassed. The scoped new asset is deposition metadata with already exposed SAXS/model content, not new independent evidence.

Real molecular fitting/ambiguity: **not executed**. Existing static/mixed-representation fitted models do not constitute an independently specified compatible forward-curve library with justified scale/background/error covariance. The new LP is exercised only on explicitly synthetic calibrated observables. There is no outcome function validated for uPAR, so it supplies no real assay recommendation or saved-experiment claim.

Milestones: public sequence/static-coordinate intake and processed SAXS QC are available; independent observable prediction, intervention/context transfer, incremental cell prediction, mediator evidence and rescue validation remain unfulfilled. PTEN VIS-seq is a separate lead, not adopted or inspected for endpoints in this bounded task; it first needs rights/schema and residue-group partition/protocol review. Original ConformaCell PoC remains unreproduced.

Primary sources: [SASDAU4](https://www.sasbdb.org/data/SASDAU4/), [official structured deposition](https://www.sasbdb.org/media/sascif/sascif_files/SASDAU4.sascif), [2012 paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC3464537/). SASBDB's previously cached deposited-data use declaration permits noncommercial/commercial use with author attribution; PMC publisher rights are not generalized. NIH/institutional eligibility remains unresolved. Exact source and acquisition hashes are in `data/post_confirmation_v10/conformacell/{readiness,acquisition,evidence_use_graph}.json`.
''')
    (R/'decision_identifiability_pilot.md').write_text(method_report()+'\nStatus: implementation and engineering tests; queued synthetic pilot not yet run by this command.\n')
    return readiness


def method_report():
    return '''# Decision-relevant identifiability v10

This is an original small finite-library LP adapter, not an original ConformaCell reproduction, a causal model or a novelty claim. For admitted library states, W={w>=0, 1^T w=1, lower<=A w<=upper}; inadmissible states have weight zero. Error intervals are deterministic assumptions, not confidence intervals. Linear outcome coefficients may themselves range between f_lower and f_upper, giving min_W f_lower.w and max_W f_upper.w. The solver records primal residual, primal-dual objective gap, stationarity residual, status and witnesses; floating-point diagnostics are not a formal global certificate.

Decision is positive only if the lower bound exceeds threshold 0.5 by 1e-7; negative only if the upper bound lies below it by 1e-7. Boundary contact and disagreement return unresolved. An empty feasible set raises an inconsistency rather than certifying a decision.

For a modeled scalar assay a with |z-a.w|<=epsilon, the largest possible post-reveal interval width is a two-copy LP: maximize f_upper.w2-f_lower.w1 over w1,w2 in W with |a.w1-a.w2|<=2epsilon. The condition is necessary and sufficient for some shared z. This optimizes worst-case residual outcome width, not expected utility or a posterior probability. Ranking uses width reduction per modeled cost; it does not prove optimal multi-step planning or minimize threshold decision error.

Random (seed 0), cheapest-first, fixed order, uncertainty-range-per-cost and worst-case-decision-width policies see identical evidence and assay models. Query assays exclude the final endpoint by typed role; policy APIs have no hidden truth argument. Two explicitly designed synthetic families include useful/nuisance assays and an uninformative null family, with seven fixed state/mixture truths each. Budgets 0/1/2 use unit modeled cost; actual laboratory costs are unknown. Every policy spends the same available budget, even after resolving a decision. Outcome coefficients have +/-0.02 uncertainty. Synthetic reveals use noiseless values within declared +/-0.05 intervals. This establishes finite-fixture algorithm behavior, not biological utility or calibrated generalization.

All intervals depend on the finite proposal library and linear forward model. Missing states, nonlinear outcomes, nuisance scale/background and correlated errors are outside this implementation. Real SAXS sigma columns do not justify these assumptions automatically. No real uPAR forward library/outcome is substituted with synthetic numbers. Conditions/endpoint-use audit is separate from the synthetic policy experiment.

Reproduce metadata readiness with `PYTHONPATH=src ~/.venvs/lesegenv/bin/python -B scripts/conformacell_v10.py --check`; then root queues `PYTHONPATH=src ~/.venvs/lesegenv/bin/python -B scripts/conformacell_v10.py --run --plan data/post_confirmation_v10/conformacell/task_plan.json --output data/post_confirmation_v10/conformacell/pilot-001`. The output directory must not exist. CPU=1, RAM<=2 GiB, wall<=180 s, no GPU/network/training.
'''


def run(plan_path, output):
    import resource
    import signal
    started=time.monotonic();plan=read(plan_path)
    if Path(plan_path).resolve()!=D/'task_plan.json' or plan['status']!='READY_SYNTHETIC_ONLY':raise ValueError('Pinned checked plan required')
    for name,value in {**plan['source_sha256'],**plan['code_sha256']}.items():
        if sha(P/name)!=value:raise ValueError('Changed pinned dependency: '+name)
    output=Path(output).resolve()
    if output.parent!=D or output.exists():raise ValueError('Unique direct child of v10 molecular output directory required')
    os.sched_setaffinity(0,{min(os.sched_getaffinity(0))})
    resource.setrlimit(resource.RLIMIT_AS,(2*1024**3,2*1024**3))
    resource.setrlimit(resource.RLIMIT_CPU,(180,180));signal.alarm(180)
    import numpy as np
    from invariantcell.identifiability_v10 import MixtureSet, ModeledAssay, choose_assay, reveal, outcome_bounds
    output.mkdir()
    write(output/'task_plan.json',plan)
    f=np.array([.1,.2,.8,.9]);lo=f-.02;hi=f+.02
    truths=[*np.eye(4),np.full(4,.25),np.array([.6,0,.4,0]),np.array([.3,.2,.4,.1])]
    families={'informative':[('A_nuisance',(0,10,0,10)),('B_decision_proxy',(0,.1,.9,1)),('C_mixed',(.1,.6,.5,.9))],
              'null':[('A_constant',(1,1,1,1)),('B_constant',(2,2,2,2)),('C_constant',(3,3,3,3))]}
    rows=[];maxgap=0.0
    for family,items in families.items():
        assays=[ModeledAssay(name,tuple(v),.05) for name,v in items];lookup={x.assay_id:x for x in assays}
        for case,truth in enumerate(truths):
            for budget in plan['budgets']:
                for policy in plan['policies']:
                    problem=MixtureSet(np.empty((0,4)),[],[],np.ones(4,dtype=bool));chosen=[];spent=0.
                    while spent<budget:
                        name=choose_assay(problem,assays,lo,hi,policy=policy,budget=budget-spent,already=chosen,seed=plan['random_seed']+len(chosen))
                        if name is None:break
                        assay=lookup[name];observed=float(np.asarray(assay.forward)@truth)
                        problem=reveal(problem,assay,observed);chosen.append(name);spent+=assay.cost
                    bounds=outcome_bounds(problem,lo,hi)
                    for side in ('lower','upper'):maxgap=max(maxgap,bounds[side+'_diagnostics']['absolute_primal_dual_gap'])
                    actual=float(f@truth);reference='POSITIVE' if actual>.5 else 'NEGATIVE' if actual<.5 else 'BOUNDARY'
                    rows.append({'family':family,'case':case,'policy':policy,'budget':budget,'spent':spent,'assays':chosen,
                                 'bounds':bounds,'synthetic_reference':actual,'reference_decision':reference,
                                 'resolved':bounds['decision']!='UNRESOLVED',
                                 'correct_if_resolved':bounds['decision']==reference if bounds['decision']!='UNRESOLVED' else None})
    aggregate=[]
    for family in families:
        for budget in plan['budgets']:
            for policy in plan['policies']:
                rr=[r for r in rows if (r['family'],r['budget'],r['policy'])==(family,budget,policy)]
                aggregate.append({'family':family,'budget':budget,'policy':policy,'cases':len(rr),
                                  'mean_width':float(np.mean([r['bounds']['width'] for r in rr])),
                                  'resolved':sum(r['resolved'] for r in rr),'wrong_resolved':sum(r['correct_if_resolved'] is False for r in rr),
                                  'mean_cost':float(np.mean([r['spent'] for r in rr]))})
    result={'status':'COMPLETE_SYNTHETIC_ENGINEERING_PILOT','real_molecular_result':False,'plan_sha256':sha(plan_path),
            'aggregate':aggregate,'rows':rows,'largest_absolute_primal_dual_gap':maxgap,
            'wall_seconds':time.monotonic()-started,'max_rss_KiB':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
            'source_sha256':plan['source_sha256'],'code_sha256':plan['code_sha256']}
    write(output/'result.json',result)
    write(output/'receipt.json',{'status':'COMPLETE','result_sha256':sha(output/'result.json'),'wall_seconds':result['wall_seconds'],
                               'max_rss_KiB':result['max_rss_KiB'],'new_download_bytes':0,'paid_cost_usd':0})
    table='\n| Family | Budget | Policy | Mean width | Resolved / 7 | Wrong resolved |\n|---|---:|---|---:|---:|---:|\n'
    table+=''.join(f"| {r['family']} | {r['budget']} | {r['policy']} | {r['mean_width']:.6f} | {r['resolved']} | {r['wrong_resolved']} |\n" for r in aggregate)
    (R/'decision_identifiability_pilot.md').write_text(method_report()+'\nExecuted synthetic pilot: `'+str(output.relative_to(P))+'`. Actual artifacts and solver diagnostics are pinned by its receipt.\n'+table)
    return {k:v for k,v in result.items() if k not in {'rows','source_sha256','code_sha256','aggregate'}}


def main():
    p=argparse.ArgumentParser();mode=p.add_mutually_exclusive_group(required=True);mode.add_argument('--check',action='store_true');mode.add_argument('--run',action='store_true')
    p.add_argument('--plan');p.add_argument('--output');args=p.parse_args()
    if args.check:result=audit();print(json.dumps({k:v for k,v in result.items() if k not in {'new_metadata','source_sha256'}},indent=2))
    else:
        if not args.plan or not args.output:p.error('--run requires --plan and --output')
        print(json.dumps(run(args.plan,args.output),indent=2))


if __name__=='__main__':main()
