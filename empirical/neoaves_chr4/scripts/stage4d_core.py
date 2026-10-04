"""Stage 4D input audits and topology-neutral rules. No inference on import."""
import csv
import hashlib
import importlib.util
import re
import random
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
BASE = ROOT / 'empirical/neoaves_chr4'
EXTERNAL = BASE / 'external/stiller2024'
DATA = BASE / 'data/stage4d'
RESULTS = BASE / 'results'
PNAS_INTERVALS = ((25030000, 32670000), (33510000, 34470000), (44130000, 56810000))

def sha256(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for b in iter(lambda: f.read(1048576), b''): h.update(b)
    return h.hexdigest()

def load_stage(name):
    spec = importlib.util.spec_from_file_location(name, BASE / 'scripts' / name)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod

def read_tsv(path):
    with Path(path).open() as f: return list(csv.DictReader(f, delimiter='\t'))

def write_tsv(path, rows, fields=None):
    rows = list(rows)
    with Path(path).open('w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=fields or list(rows[0]), delimiter='\t', lineterminator='\n')
        w.writeheader(); w.writerows(rows)

def is_chr4(chromosome):
    stem = chromosome.split('_', 1)[0]
    return (stem[3:] if stem.startswith('chr') else stem) == '4'

def published_outlier(chromosome, window_start):
    return chromosome == 'chr4' and any(a <= window_start <= b for a, b in PNAS_INTERVALS)

def representatives(rows, width=500000):
    """Zero-anchored, half-open bins on each named sequence, including random scaffolds.

    Uses only chromosome, midpoint and locus_id. Sequence-local coordinates must
    never be pooled across scaffolds. Ties are broken by lexical locus identifier.
    """
    if width <= 0: raise ValueError('positive bin width required')
    chosen = {}
    for row in rows:
        chrom, midpoint = row['chromosome'], float(row['midpoint'])
        if midpoint < 0: raise ValueError('negative coordinate')
        if not is_chr4(chrom): continue
        key = (chrom, int(midpoint // width))
        rank = (abs(midpoint - (key[1] + .5) * width), row['locus_id'])
        if key not in chosen or rank < chosen[key][0]: chosen[key] = (rank, row['locus_id'])
    return {key: value[1] for key, value in sorted(chosen.items())}

def treatment_ids(rows, treatment, width=500000):
    """Return locus IDs for a predeclared treatment from manifest fields only."""
    if treatment == 'T0': return {r['locus_id'] for r in rows}
    if treatment == 'NONCHR4': return {r['locus_id'] for r in rows if not _bool(r['is_chr4'])}
    if treatment == 'T_PNAS': return {r['locus_id'] for r in rows if not _bool(r['published_outlier_region'])}
    if treatment == 'T_STRUCT': return {r['locus_id'] for r in rows if not _bool(r['frozen_structural_primary'])}
    if treatment == 'T_BLOCK':
        reps=set(representatives(rows,width).values())
        return {r['locus_id'] for r in rows if not _bool(r['is_chr4']) or r['locus_id'] in reps}
    raise ValueError(treatment)

def random_matched_ids(rows, n_chr4, seed):
    non={r['locus_id'] for r in rows if not _bool(r['is_chr4'])}
    chrom=sorted(r['locus_id'] for r in rows if _bool(r['is_chr4']))
    if not 0 <= n_chr4 <= len(chrom): raise ValueError('invalid chr4 count')
    return non | set(random.Random(seed).sample(chrom,n_chr4))

def _bool(value):
    return value is True or str(value).lower() in {'true','1','yes'}

# The archive has plain unquoted Newick labels and no comments. Reject other
# dialects instead of silently misparsing them. Species-tree comparisons below
# use DendroPy, independently of this fast streaming gene-tree parser.
TOKEN = re.compile(r'[(),:;]|[^\s(),:;]+')

def parse_newick(text):
    if any(x in text for x in "[]'\""): raise ValueError('unsupported Newick dialect')
    tokens = TOKEN.findall(text)
    i = 0
    def node():
        nonlocal i
        children = []
        if tokens[i] == '(':
            i += 1
            children.append(node())
            while tokens[i] == ',':
                i += 1; children.append(node())
            if tokens[i] != ')': raise ValueError('unbalanced Newick')
            i += 1
        label = ''
        if tokens[i] not in ',):;': label = tokens[i]; i += 1
        length = 0.0
        if tokens[i] == ':': i += 1; length = float(tokens[i]); i += 1
        if not children and not label: raise ValueError('unnamed leaf')
        return [children, label, length]
    root = node()
    if tokens[i:] != [';']: raise ValueError('trailing Newick content')
    return root

def tree_signature(root, bits, threshold=None):
    edges, leaves, depths = [], [], {}
    def walk(n, depth):
        depth += n[2]
        if not n[0]:
            if n[1] in depths: raise ValueError('duplicate taxon')
            leaves.append(n[1]); depths[n[1]] = depth
            return bits[n[1]]
        mask = 0
        for c in n[0]: mask |= walk(c, depth)
        support = float(n[1].lstrip('/')) if n[1] else None
        if n is not root and (threshold is None or (support is not None and support >= threshold)):
            edges.append((mask, support))
        return mask
    allbits = walk(root, 0)
    splits = {}
    for b, support in edges:
        other = allbits ^ b
        if b.bit_count() > 1 and other.bit_count() > 1:
            splits[min(b, other)] = support
    return allbits, splits, leaves, depths

def collapse_named(root, threshold=.95):
    """Contract support<threshold; retain IDs externally. Preserve path lengths."""
    def visit(n, is_root=False):
        children=[]
        for child in n[0]: children.extend(visit(child))
        n=[children, n[1].lstrip('/') if children else n[1], n[2]]
        if children and not is_root:
            # An unlabeled multifurcation is already unresolved and is retained.
            if n[1] and float(n[1]) < threshold:
                for c in children: c[2] += n[2]
                return children
        return [n]
    return visit(root, True)[0]

def emit_newick(n):
    def emit(x):
        return ('('+','.join(emit(c) for c in x[0])+')' if x[0] else '')+x[1]+':'+format(x[2],'.12g')
    return emit(n)+';'

def prune_tree(root, keep):
    """Prune leaves and suppress unary nodes while preserving path lengths."""
    def visit(n, is_root=False):
        if not n[0]: return n if n[1] in keep else None
        children=[x for c in n[0] if (x:=visit(c)) is not None]
        if not children: return None
        if len(children)==1:
            child=children[0]; child[2]+=n[2]
            return child
        return [children,n[1],n[2]]
    return visit(root,True)

def rf_comparison(a, b):
    """Unrooted split symmetric difference on identical taxon sets."""
    import dendropy
    from dendropy.calculate import treecompare
    ns = dendropy.TaxonNamespace()
    ta = dendropy.Tree.get(data=a, schema='newick', taxon_namespace=ns, rooting='force-unrooted', preserve_underscores=True)
    tb = dendropy.Tree.get(data=b, schema='newick', taxon_namespace=ns, rooting='force-unrooted', preserve_underscores=True)
    if {n.taxon.label for n in ta.leaf_node_iter()} != {n.taxon.label for n in tb.leaf_node_iter()}:
        raise ValueError('RF requires identical taxa; pruning must be explicit')
    ta.encode_bipartitions(); tb.encode_bipartitions()
    distance = treecompare.symmetric_difference(ta, tb)
    sa = {e.bipartition.split_bitmask for e in ta.postorder_edge_iter() if not e.bipartition.is_trivial()}
    sb = {e.bipartition.split_bitmask for e in tb.postorder_edge_iter() if not e.bipartition.is_trivial()}
    assert distance == len(sa ^ sb)
    return {'rf': distance, 'normalized_rf': distance / (len(sa) + len(sb)) if sa or sb else 0,
            'branches_only_a': len(sa - sb), 'branches_only_b': len(sb - sa)}

def named_splits(newick):
    """Canonical unrooted nontrivial splits, retaining supplied edge labels."""
    import dendropy
    tree = dendropy.Tree.get(data=newick, schema='newick', rooting='force-unrooted',
                             preserve_underscores=True)
    taxa = frozenset(n.taxon.label for n in tree.leaf_node_iter())
    if len(taxa) != len(list(tree.leaf_node_iter())):
        raise ValueError('duplicate species-tree tips')
    splits = {}
    for node in tree.postorder_node_iter():
        side = frozenset(n.taxon.label for n in node.leaf_iter())
        other = taxa - side
        if min(len(side), len(other)) <= 1:
            continue
        key = min(tuple(sorted(side)), tuple(sorted(other)))
        splits[key] = node.label
    return taxa, splits


def frozen_focal_groups():
    """Read the already-frozen Stage-3B/4C role definitions without edits."""
    import zipfile
    roles = read_tsv(RESULTS / 'stage3b_v2_quartet_role_mapping.tsv')
    archive = ROOT / 'data/neoaves_chr4/raw/genetreesupport/clade-analysis/examined-clades.zip'
    groups = {}
    with zipfile.ZipFile(archive) as z:
        for row in roles:
            if row['role'] == 'O':
                continue
            groups.setdefault(row['clade'], {})[row['role']] = set(
                z.read(row['source_definition']).decode().split())
    return groups


def focal_split_states(newick, groups):
    """Conservatively classify exact frozen quadripartition splits.

    No best representative taxon or majority quartet is chosen. If no exact
    focal split exists, report mixed_or_unresolved, not a preferred topology.
    """
    taxa, splits = named_splits(newick)
    out = []
    for clade, definition in groups.items():
        parts = {r: set(definition[r]) & taxa for r in ('C1', 'C2', 'S')}
        if any(parts[a] & parts[b] for a, b in (('C1','C2'), ('C1','S'), ('C2','S'))):
            raise ValueError('overlapping frozen focal roles')
        parts['O'] = set(taxa) - set.union(*parts.values())
        matches = []
        for state, pair in (('q1', ('C1','C2')), ('q2', ('C1','S')), ('q3', ('C2','S'))):
            side = parts[pair[0]] | parts[pair[1]]
            key = min(tuple(sorted(side)), tuple(sorted(taxa-side)))
            if all(parts.values()) and key in splits:
                matches.append((state, splits[key]))
        if len(matches) > 1:
            raise ValueError('incompatible focal splits in species tree')
        out.append(dict(clade=clade,
                        state=matches[0][0] if matches else ('mixed_or_unresolved' if all(parts.values()) else 'missing_role'),
                        branch_support=matches[0][1] if matches else None,
                        **{f'n_{r}': len(v) for r, v in parts.items()}))
    return out
