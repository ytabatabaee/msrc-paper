import csv
import importlib.util
import json
from pathlib import Path

BASE = Path(__file__).parents[1]
SPEC = importlib.util.spec_from_file_location('stage4d', BASE/'scripts/stage4d_core.py')
s = importlib.util.module_from_spec(SPEC); SPEC.loader.exec_module(s)

def test_published_outlier_exact_window_start_rule():
    assert s.PNAS_INTERVALS == ((25030000,32670000),(33510000,34470000),(44130000,56810000))
    assert s.published_outlier('chr4',25030000)
    assert s.published_outlier('chr4',56810000)
    assert not s.published_outlier('chr4_random',25030000)
    assert not s.published_outlier('chr4',56810001)

def test_chr4_rule_matches_stage4c_random_scaffolds():
    assert s.is_chr4('chr4') and s.is_chr4('chr4_AADN_random')
    assert not s.is_chr4('chr14')

def test_zero_anchored_representatives_are_deterministic_and_scaffold_local():
    rows=[{'locus_id':'z','chromosome':'chr4','midpoint':250000},
          {'locus_id':'a','chromosome':'chr4','midpoint':250000},
          {'locus_id':'b','chromosome':'chr4','midpoint':750000},
          {'locus_id':'r','chromosome':'chr4_random','midpoint':250000},
          {'locus_id':'x','chromosome':'chr5','midpoint':250000}]
    expected={('chr4',0):'a',('chr4',1):'b',('chr4_random',0):'r'}
    assert s.representatives(rows)==expected
    assert s.representatives(list(reversed(rows)))==expected

def test_block_selector_has_no_focal_topology_inputs():
    import inspect
    signature=str(inspect.signature(s.representatives))
    source=inspect.getsource(s.representatives)
    assert signature == '(rows, width=500000)'
    for forbidden in ('q1','q2','q3','N61','N62','Columbea','S2024','J2014'):
        assert forbidden not in source

def test_treatments_and_random_control_membership():
    rows=[{'locus_id':'n','chromosome':'chr5','midpoint':1,'is_chr4':False,'published_outlier_region':False,'frozen_structural_primary':False},
          {'locus_id':'a','chromosome':'chr4','midpoint':250000,'is_chr4':True,'published_outlier_region':True,'frozen_structural_primary':False},
          {'locus_id':'b','chromosome':'chr4','midpoint':260000,'is_chr4':True,'published_outlier_region':False,'frozen_structural_primary':True},
          {'locus_id':'c','chromosome':'chr4','midpoint':750000,'is_chr4':True,'published_outlier_region':False,'frozen_structural_primary':False}]
    assert s.treatment_ids(rows,'NONCHR4')=={'n'}
    assert s.treatment_ids(rows,'T0')=={'n','a','b','c'}
    assert s.treatment_ids(rows,'T_PNAS')=={'n','b','c'}
    assert s.treatment_ids(rows,'T_STRUCT')=={'n','a','c'}
    block=s.treatment_ids(rows,'T_BLOCK'); assert block=={'n','a','c'}
    random=s.random_matched_ids(rows,2,431729)
    assert len(random-{'n'})==2 and 'n' in random
    assert random==s.random_matched_ids(list(reversed(rows)),2,431729)

def test_collapse_threshold_and_path_lengths():
    bits={'A':1,'B':2,'C':4,'D':8}
    n=s.parse_newick('((A:1,B:1)/0.94:2,(C:1,D:1)/0.95:3);')
    c=s.collapse_named(n)
    allbits,splits,leaves,depths=s.tree_signature(c,bits)
    assert len(splits)==1
    assert depths['A']==3 and depths['C']==4

def test_taxon_pruning_suppresses_unary_paths():
    n=s.parse_newick('((A:1,B:1)1:2,(C:1,D:1)1:3);')
    p=s.prune_tree(n,{'A','C','D'})
    _,_,leaves,depths=s.tree_signature(p,{'A':1,'C':2,'D':4})
    assert set(leaves)=={'A','C','D'} and depths['A']==3

def test_rf_identity_and_single_nni():
    assert s.rf_comparison('((A,B),(C,D));','((A,B),(C,D));')['rf']==0
    result=s.rf_comparison('((A,B),(C,D));','((A,C),(B,D));')
    assert result['rf']==2 and result['normalized_rf']==1

def test_frozen_prior_inventory_unchanged():
    inventory=BASE/'results/stage4d/frozen_prior_inventory.json'
    if not inventory.exists(): return
    for row in json.loads(inventory.read_text()):
        p=s.ROOT/row['path']
        assert p.exists() and s.sha256(p)==row['sha256'], row['path']

def test_manifest_treatment_membership_if_audit_complete():
    p=BASE/'results/stage4d_locus_manifest.tsv'
    if not p.exists(): return
    with p.open() as f: rows=list(csv.DictReader(f,delimiter='\t'))
    assert len(rows)==63430 and len({r['locus_id'] for r in rows})==63430
    assert sum(r['is_chr4']=='True' for r in rows if not s.is_chr4(r['chromosome']))==0

def checkpoint_runner():
    import sys
    sys.modules['stage4d_core'] = s
    spec = importlib.util.spec_from_file_location('checkpoint_runner', BASE/'scripts/04d_run_checkpoint.py')
    runner = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(runner)
    return runner

def test_checkpoint_completion_requires_success_log_and_tree(tmp_path):
    runner = checkpoint_runner()
    tree, log = tmp_path/'tree.nwk', tmp_path/'stderr.log'
    tree.write_text('')
    log.write_text('ASTRAL finished in 1 secs\n')
    assert not runner.completed(0, tree, log)
    tree.write_text('((A,B),(C,D));\n')
    assert not runner.completed(1, tree, log)
    assert runner.completed(0, tree, log)
    log.write_text('Building set of clusters (X) from gene trees\n')
    assert not runner.completed(0, tree, log)

def test_checkpoint_inference_settings_are_fixed(tmp_path):
    runner = checkpoint_runner()
    commands = [runner.command('/java', tmp_path, tmp_path/name, tmp_path/'out', '8g', 4)
                for name in runner.INPUTS.values()]
    for cmd in commands:
        assert cmd[cmd.index('-s')+1] == '692'
        assert cmd[cmd.index('-t')+1] == '3'
        assert '-C' in cmd
        assert not any(x in cmd for x in ('-x','-p','-q'))
    assert all(cmd[:cmd.index('-i')] == commands[0][:commands[0].index('-i')] for cmd in commands)

def test_focal_split_extraction_handles_alternatives_polytomy_and_missing_role():
    groups={'example':{'C1':{'A'},'C2':{'B'},'S':{'C'}}}
    for tree, state in [('((A,B),(C,D));','q1'),('((A,C),(B,D));','q2'),
                        ('((B,C),(A,D));','q3'),('(A,B,C,D);','mixed_or_unresolved'),
                        ('(A,B,C);','missing_role')]:
        assert s.focal_split_states(tree,groups)[0]['state']==state

def test_focal_split_extraction_does_not_choose_a_favorable_representative():
    groups={'example':{'C1':{'A','E'},'C2':{'B'},'S':{'C'}}}
    assert s.focal_split_states('((A,B),(E,(C,D)));',groups)[0]['state']=='q1'
    assert s.focal_split_states('((A,B),(C,(E,D)));',groups)[0]['state']=='mixed_or_unresolved'
