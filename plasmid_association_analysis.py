import os, math, collections
import numpy as np, pandas as pd
from scipy.stats import chi2

PHASE_DIR = 'phase10'
OUT_DIR = os.path.join(PHASE_DIR, 'assoc_v2')
os.makedirs(OUT_DIR, exist_ok=True)

master = pd.read_csv(os.path.join(PHASE_DIR, 'phase10_master_table_v2.csv'), low_memory=False, dtype=str)
N = len(master)
RES_COLS = ['AGly_acquired','Bla_acquired','Bla_ESBL_acquired','Bla_Carb_acquired','Col_acquired','Fcyn_acquired','Flq_acquired','MLS_acquired','Phe_acquired','Rif_acquired','Sul_acquired','Tet_acquired','Tgc_acquired','Tmt_acquired']
VIR_COLS = ['Yersiniabactin','Aerobactin','Salmochelin','Colibactin','rmpA','rmpA2']
missing = [c for c in RES_COLS + VIR_COLS + ['base_ST','genome_id'] if c not in master.columns]
if missing:
    raise SystemExit('missing columns: %s' % missing)


def is_present(col):
    return (~master[col].fillna('.').isin(['.','','0','0.0','nan'])).to_numpy()


skip_cols = ('NUM_FOUND', '#FILE')
char_pres = {c: is_present(c) for c in RES_COLS + VIR_COLS}                                    # AMR classes + virulence
card_pres = {c[5:]: is_present(c) for c in master.columns if c.startswith('CARD_') and c[5:] not in skip_cols}  # CARD genes
rep_pres = {c[3:]: is_present(c) for c in master.columns if c.startswith('PF_') and c[3:] not in skip_cols}    # plasmid replicons
print('rows', N, '| CARD genes', len(card_pres), '| replicons', len(rep_pres))

# strata = base_ST
st_code, st_levels = pd.factorize(master['base_ST'].fillna('NA'))
n_st = len(st_levels)
st_n = np.bincount(st_code, minlength=n_st).astype(float)


def mh_test(x, y):
    # Mantel-Haenszel OR with 95% CI and chi-square test, stratified by base_ST
    x = x.astype(float); y = y.astype(float)
    a = np.bincount(st_code, weights=x*y, minlength=n_st)
    nx = np.bincount(st_code, weights=x, minlength=n_st)
    ny = np.bincount(st_code, weights=y, minlength=n_st)
    ok = st_n > 1
    a, nx, ny, n = a[ok], nx[ok], ny[ok], st_n[ok]
    b = nx - a; c = ny - a; d = n - nx - ny + a
    Ri = a*d/n; Si = b*c/n
    Rs = Ri.sum(); Ss = Si.sum()
    if Rs == 0 or Ss == 0:
        mh_or = np.inf if (Ss == 0 and Rs > 0) else np.nan
        lo = hi = np.nan
    else:
        P = (a+d)/n; Q = (b+c)/n
        var = (P*Ri).sum()/(2*Rs**2) + (P*Si + Q*Ri).sum()/(2*Rs*Ss) + (Q*Si).sum()/(2*Ss**2)
        mh_or = Rs/Ss; se = math.sqrt(var)
        lo, hi = mh_or*math.exp(-1.96*se), mh_or*math.exp(1.96*se)
    E = nx*ny/n
    V = nx*(n-nx)*ny*(n-ny)/(n**2*(n-1))
    p = chi2.sf(max(0.0, abs((a-E).sum())-0.5)**2/V.sum(), 1) if V.sum() > 0 else np.nan
    return mh_or, lo, hi, p


def bh_adjust(pvals):
    # Benjamini-Hochberg q-values, NaNs left as NaN
    pvals = np.array(pvals, float)
    qvals = np.full(len(pvals), np.nan)
    valid = ~np.isnan(pvals)
    pv = pvals[valid]
    if len(pv) == 0:
        return qvals
    order = np.argsort(pv)
    ranked = pv[order]*len(pv)/(np.arange(len(pv))+1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    adj = np.empty(len(pv))
    adj[order] = np.minimum(ranked, 1)
    qvals[valid] = adj
    return qvals


def family(x_pres, x_names, y_pres, y_names, name):
    rows = []
    for xk in x_names:
        x = x_pres[xk]
        if x.sum() < 50:
            continue
        for yk in y_names:
            y = y_pres[yk]
            nb = int((x & y).sum())
            if nb < 10:
                continue
            a = nb; b = int(x.sum()) - nb; c = int(y.sum()) - nb; d = N - a - b - c
            crude = ((a+.5)*(d+.5))/((b+.5)*(c+.5))
            mh_or, lo, hi, p = mh_test(x, y)
            rows.append(dict(x=xk, y=yk, n_both=nb, n_x=int(x.sum()), n_y=int(y.sum()), crude_OR=crude, MH_OR=mh_or, CI_lo=lo, CI_hi=hi, p=p))
    df = pd.DataFrame(rows)
    df['q_BH'] = bh_adjust(df['p'])
    df.to_csv(os.path.join(OUT_DIR, name + '.csv'), index=False)
    print(name, 'eligible tests', len(df), '| q<0.05:', int((df['q_BH'] < 0.05).sum()))
    return df


genes_ok = [g for g in card_pres if 50 <= card_pres[g].sum() <= 0.8*N]
reps_ok = [r for r in rep_pres if rep_pres[r].sum() >= 50]
F1 = family(char_pres, RES_COLS, card_pres, genes_ok, 'F1_class_x_CARDgene')
F2 = family(char_pres, RES_COLS, rep_pres, reps_ok, 'F2_class_x_replicon')
F3 = family(char_pres, VIR_COLS, rep_pres, reps_ok, 'F3_virulence_x_replicon')

print('\nTop 15 F1 by MH OR (q<0.001, n_both>=50):')
top = F1[(F1.q_BH < 0.001) & (F1.n_both >= 50)].sort_values('MH_OR', ascending=False).head(15)
print(top[['x','y','n_both','MH_OR','CI_lo','CI_hi','q_BH']].to_string(index=False))


def lookup(df, x_name, y_name):
    r = df[(df.x == x_name) & (df.y == y_name)]
    return (round(float(r.MH_OR.iloc[0]), 2), int(r.n_both.iloc[0])) if len(r) else ('NA', 'NA')


# (class, CARD gene, MH OR reported in the paper)
paper_card = [
    ('Tet_acquired','tet(A)',757.26),
    ('Sul_acquired','sul1',360.84),
    ('Bla_ESBL_acquired','CTX-M-15',346.30),
    ('Sul_acquired','APH(6)-Id',34.84),
    ('Flq_acquired',"AAC(6')-Ib10",27.68),
    ('Flq_acquired','OXA-1',25.26),
    ('Sul_acquired','qacEdelta1',24.52),
    ('AGly_acquired','sul1',23.48),
    ('Rif_acquired','floR',15.90),
    ('MLS_acquired','dfrA12',15.44),
    ('MLS_acquired','sul1',11.85),
    ('Tmt_acquired','sul1',9.79),
    ('Bla_acquired',"AAC(6')-Ib10",9.44),
    ('Bla_ESBL_acquired','TEM-1',7.21),
]
print('\nPaper Table 3 vs this run (MH OR, n_both):')
for x_name, y_name, paper_or in paper_card:
    print(x_name, y_name, '| paper', paper_or, '| this run', lookup(F1, x_name, y_name))

# (class or virulence locus, replicon, MH OR reported in the paper, table to look in)
paper_pf = [
    ('Bla_Carb_acquired','IncX3_1',17.68,F2),
    ('Bla_Carb_acquired','IncL/M(pOXA-48)_1_pOXA-48',27.85,F2),
    ('AGly_acquired','IncA/C2_1',13.13,F2),
    ('rmpA','repB_KLEB_VIR',162.97,F3),
    ('Salmochelin','repB_KLEB_VIR',136.10,F3),
    ('rmpA2','repB_KLEB_VIR',86.04,F3),
    ('Aerobactin','repB_KLEB_VIR',58.40,F3),
    ('rmpA','IncHI1B_1_pNDM-MAR',35.43,F3),
]
print('\nPaper Table 4 ORs vs this run:')
for x_name, y_name, paper_or, tab in paper_pf:
    print(x_name, y_name, '| paper', paper_or, '| this run', lookup(tab, x_name, y_name))

print('\nGene-level context counts (reading raw abricate and MOB files)...')
mob = pd.read_csv(os.path.join(PHASE_DIR, 'mob_consolidated', 'contig_to_cluster.tsv'), sep='\t', dtype=str, usecols=['genome_id','molecule_type','primary_cluster_id','contig_id'])
mob['contig'] = mob['contig_id'].str.split().str[0]
# (genome, contig) -> (molecule_type, primary_cluster_id)
mob_map = {(gid, contig): (mol_type, cluster)
           for gid, contig, mol_type, cluster in zip(mob.genome_id, mob.contig, mob.molecule_type, mob.primary_cluster_id)}
focal_reps = ['IncX3_1','IncL/M(pOXA-48)_1_pOXA-48','IncFII_1_pKP91','IncA/C2_1','repB_KLEB_VIR','IncHI1B_1_pNDM-MAR']
focal_reps = [r for r in focal_reps if r in rep_pres]
gene_set = set(genes_ok); rep_set = set(focal_reps)


def read_hits(path, wanted):
    # abricate .tab -> {gene/replicon name: set of contigs it was found on}
    hits = collections.defaultdict(set)
    try:
        with open(path) as f:
            for line in f:
                if line.startswith('#'):
                    continue
                fields = line.rstrip('\n').split('\t')
                if len(fields) > 5 and fields[5] in wanted:
                    hits[fields[5]].add(fields[1])
    except FileNotFoundError:
        pass
    return hits


# counts[(gene, replicon)] = [n_both, n_mob_plasmid, n_same_contig, n_same_unit]
counts = collections.defaultdict(lambda: [0, 0, 0, 0])
for gid in master['genome_id']:
    gene_hits = read_hits(os.path.join(PHASE_DIR, 'abricate_out', 'card', gid + '.tab'), gene_set)
    rep_hits = read_hits(os.path.join(PHASE_DIR, 'abricate_out', 'plasmidfinder', gid + '.tab'), rep_set)
    if not gene_hits or not rep_hits:
        continue
    for rep_name, rep_contigs in rep_hits.items():
        rep_mob = [mob_map.get((gid, c)) for c in rep_contigs]
        rep_on_plasmid = any(x and x[0] == 'plasmid' for x in rep_mob)
        rep_units = {x[1] for x in rep_mob if x and x[0] == 'plasmid' and x[1] not in ('-', '')}
        for gene_name, gene_contigs in gene_hits.items():
            cell = counts[(gene_name, rep_name)]
            cell[0] += 1
            if rep_on_plasmid: cell[1] += 1
            if gene_contigs & rep_contigs: cell[2] += 1
            gene_units = set()
            for c in gene_contigs:
                x = mob_map.get((gid, c))
                if x and x[0] == 'plasmid' and x[1] not in ('-', ''):
                    gene_units.add(x[1])
            if rep_units & gene_units: cell[3] += 1

rows = [dict(gene=gene_name, replicon=rep_name, n_both=v[0], n_mob_plasmid=v[1], n_same_contig=v[2], n_same_unit=v[3],
             pct_mob=100*v[1]/v[0], pct_same_contig=100*v[2]/v[0], pct_same_unit=100*v[3]/v[0])
        for (gene_name, rep_name), v in counts.items() if v[0] >= 30]
ctx = pd.DataFrame(rows).sort_values('n_both', ascending=False)
ctx.to_csv(os.path.join(OUT_DIR, 'context_gene_x_replicon.csv'), index=False)
print(ctx.head(40).round(1).to_string(index=False))
print('\nDone. Files in', OUT_DIR)