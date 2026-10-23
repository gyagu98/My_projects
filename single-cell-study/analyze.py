"""Reanalyse the checksum-verified processed Kang PBMC export; no cell-level inference."""
import argparse
import hashlib
import json
import re
import time
import warnings
from pathlib import Path

import anndata as ad
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats
from sklearn.decomposition import PCA

SHA = 'e6a5adac64dcdeb36eaba27db49b63e0c64bb0ed4a64c6705971506b41c39830'
SOURCE = 'https://api.figshare.com/v2/file/download/34464122'


def aggregate(x, obs, genes, cell_type, min_cells=20):
    """Return integer sample-by-gene counts, aligned metadata and a stratum audit."""
    needed = ['donor', 'condition', 'cell_type']
    if obs[needed].isna().any().any() or any(obs[c].astype(str).str.strip().eq('').any() for c in needed):
        raise ValueError('Missing sample metadata')
    if not obs.index.is_unique or not pd.Index(genes).is_unique or x.shape != (len(obs), len(genes)):
        raise ValueError('Counts/metadata dimensions or identifiers are invalid')
    values = x.data if hasattr(x, 'tocsr') else np.asarray(x).ravel()
    if not np.isfinite(values).all() or (values < 0).any() or (values != np.round(values)).any():
        raise ValueError('Raw nonnegative integer counts required')
    if set(obs.condition) != {'ctrl', 'stim'}:
        raise ValueError('Expected ctrl/stim conditions')
    mask = obs.cell_type.eq(cell_type)
    audit = (obs.loc[mask].groupby(['donor','condition'], observed=True).size()
             .unstack(fill_value=0).reindex(columns=['ctrl','stim'], fill_value=0))
    donors = audit.index[(audit >= min_cells).all(axis=1)].tolist()
    audit['retained'] = audit.index.isin(donors)
    if len(donors) < 3:
        return None, None, audit
    rows, metadata = [], []
    for donor in donors:
        for condition in ['ctrl','stim']:
            selected = mask & obs.donor.eq(donor) & obs.condition.eq(condition)
            rows.append(np.asarray(x[np.asarray(selected)].sum(axis=0)).ravel())
            metadata.append({'donor':donor,'condition':condition,'n_cells':int(selected.sum())})
    meta = pd.DataFrame(metadata, index=[f'pb{i}' for i in range(len(rows))])
    counts = pd.DataFrame(np.asarray(rows, dtype=np.int64), index=meta.index, columns=genes)
    return counts, meta, audit


def fit(counts, meta):
    keep = (counts >= 10).sum(axis=0) >= 3
    if keep.sum() < 2:
        raise ValueError('Insufficient genes after filtering')
    dds = DeseqDataSet(counts=counts.loc[:,keep], metadata=meta,
                      design='~ donor + condition', n_cpus=2, quiet=True)
    with warnings.catch_warnings(record=True) as captured:
        warnings.simplefilter('always')
        dds.deseq2()
        stats = DeseqStats(dds, contrast=['condition','stim','ctrl'], n_cpus=2, quiet=True, alpha=.05)
        stats.summary()
    result=stats.results_df
    result.attrs['fit_warnings']=sorted(set(str(w.message) for w in captured))
    return result, pd.DataFrame(dds.layers['normed_counts'], index=meta.index, columns=dds.var_names)


def run(input_path, outdir):
    started = time.time()
    out = Path(outdir); out.mkdir(parents=True, exist_ok=True)
    checksum = hashlib.sha256(Path(input_path).read_bytes()).hexdigest()
    if checksum != SHA:
        raise ValueError('Input checksum differs from the reviewed processed export')
    with warnings.catch_warnings():
        warnings.simplefilter('ignore', FutureWarning)
        a = ad.read_h5ad(input_path)
    if not a.obs_names.is_unique or not a.var_names.is_unique:
        raise ValueError('Nonunique identifiers')
    obs = a.obs[['replicate','label','cell_type']].rename(columns={'replicate':'donor','label':'condition'}).copy()
    if obs.isna().any().any():
        raise ValueError('Source contains missing donor, condition or annotation metadata')
    for col in obs:
        obs[col] = obs[col].astype(str)
    total = np.asarray(a.X.sum(axis=1)).ravel()
    detected = np.asarray((a.X > 0).sum(axis=1)).ravel()
    if not np.allclose(total, a.obs['nCount_RNA']):
        raise ValueError('Matrix does not match source raw-count library totals')
    valid = detected >= 200
    qc = obs.assign(total_counts=total, detected_genes=detected, retained=valid)
    qc.groupby(['donor','condition'], observed=True).agg(cells=('retained','size'),
        retained=('retained','sum'), median_counts=('total_counts','median'),
        median_genes=('detected_genes','median')).to_csv(out/'sample_qc.csv')
    obs = obs.loc[valid]; x = a.X[valid]
    provenance = {'source':SOURCE, 'sha256':checksum, 'source_dimensions':list(a.shape),
        'retained_cells':len(obs), 'donors':sorted(obs.donor.unique()),
        'count_matrix_matches_source_totals':True, 'detected_gene_minimum':200,
        'design':'~ donor + condition', 'contrast':['condition','stim','ctrl'],
        'gene_filter':'count >=10 in at least 3 pseudobulk samples',
        'mitochondrial_qc':'Unavailable: processed export has no MT- genes',
        'annotations':'Inherited cell_type labels; not independently validated',
        'multiplets':'Flags unavailable in this processed export; no new doublet calling',
        'p_adjustment':'BH within cell type; independent filtering and Cooks filtering enabled',
        'cross_cell_type_claims':'Descriptive; no treatment-by-cell-type interaction test'}
    (out/'provenance.json').write_text(json.dumps(provenance, indent=2))
    summary, sensitivity, fit_warnings = [], [], {}
    for cell_type in sorted(obs.cell_type.unique()):
        slug = re.sub(r'[^a-z0-9]+','_',cell_type.lower()).strip('_')
        counts, meta, audit = aggregate(x, obs, a.var_names, cell_type)
        audit.to_csv(out/f'{slug}_donor_audit.csv')
        if counts is None:
            summary.append({'cell_type':cell_type,'status':'insufficient complete pairs'})
            continue
        print(f'Fitting {cell_type}: {len(meta)//2} pairs', flush=True)
        result, norm = fit(counts, meta)
        fit_warnings[cell_type+' / primary']=result.attrs['fit_warnings']
        result.index.name = 'gene'
        result.sort_values('padj').to_csv(out/f'{slug}_differential_expression.csv')
        meta.to_csv(out/f'{slug}_samples.csv')
        counts.to_csv(out/f'{slug}_pseudobulk_counts.csv.gz', compression='gzip')
        norm.to_csv(out/f'{slug}_normalized_counts.csv.gz', compression='gzip')
        sig = result.padj.lt(.05)
        top = result.sort_values('padj').head(10)
        top.to_csv(out/f'{slug}_top_genes.csv')
        summary.append({'cell_type':cell_type,'status':'analyzed','donor_pairs':len(meta)//2,
            'cells':int(meta.n_cells.sum()),'tested_genes':len(result),
            'finite_padj':int(result.padj.notna().sum()), 'padj_below_005':int(sig.sum()),
            'up':int((sig & result.log2FoldChange.gt(0)).sum()),
            'down':int((sig & result.log2FoldChange.lt(0)).sum())})
        for threshold in [10,20,50]:
            c,m,au = aggregate(x,obs,a.var_names,cell_type,threshold)
            same = m is not None and m.donor.tolist()==meta.donor.tolist()
            row = {'cell_type':cell_type,'min_cells':threshold,
                   'pairs':0 if m is None else len(m)//2,'same_samples_as_primary':same}
            if m is not None and not same:
                alt,_ = fit(c,m)
                fit_warnings[cell_type+f' / threshold {threshold}']=alt.attrs['fit_warnings']
                both = result.add_suffix('_primary').join(alt.add_suffix('_sensitivity'))
                row['log2fc_correlation'] = both.log2FoldChange_primary.corr(both.log2FoldChange_sensitivity)
                alt.to_csv(out/f'{slug}_threshold_{threshold}_de.csv')
            sensitivity.append(row)
        # Descriptive donor-level and leave-one-donor-out mean differences;
        # deliberately NOT represented as refitted statistical significance.
        lognorm = np.log2(norm.loc[:,top.index]+1)
        deltas=[]
        for donor in sorted(meta.donor.unique()):
            stim=lognorm.loc[(meta.donor==donor)&(meta.condition=='stim')].iloc[0]
            ctrl=lognorm.loc[(meta.donor==donor)&(meta.condition=='ctrl')].iloc[0]
            deltas.append(stim-ctrl)
        deltas=pd.DataFrame(deltas,index=sorted(meta.donor.unique()));deltas.index.name='donor'
        deltas.to_csv(out/f'{slug}_top_gene_donor_log2_differences.csv')
        pd.DataFrame({d:deltas.drop(index=d).mean() for d in deltas.index}).to_csv(out/f'{slug}_leave_one_out_mean_differences.csv')
        z=np.log2(norm+1)
        variable=z.var().nlargest(min(500,z.shape[1])).index
        pc=PCA(n_components=2).fit_transform(z[variable])
        fig,axs=plt.subplots(1,2,figsize=(11,4),layout='constrained')
        for condition,color in [('ctrl','#397ca7'),('stim','#d16532')]:
            sel=meta.condition.eq(condition)
            axs[0].scatter(pc[sel,0],pc[sel,1],label=condition,color=color)
        for donor in meta.donor.unique():
            sel=meta.donor.eq(donor);axs[0].plot(pc[sel,0],pc[sel,1],color='gray',alpha=.35)
        axs[0].set(xlabel='PC1',ylabel='PC2',title='Pseudobulk log-count PCA (descriptive)');axs[0].legend()
        im=axs[1].imshow(deltas.to_numpy(),aspect='auto',cmap='RdBu_r',vmin=-5,vmax=5)
        axs[1].set_xticks(range(len(top)),top.index,rotation=80,fontsize=8)
        axs[1].set_yticks(range(len(deltas)),deltas.index,fontsize=7)
        axs[1].set_title('Top genes: donor-level log2(count+1) differences')
        fig.colorbar(im,ax=axs[1],label='stim minus ctrl, normalized counts')
        fig.suptitle(cell_type);fig.savefig(out/f'{slug}_diagnostics.png',dpi=150)
        fig.savefig(out/f'{slug}_diagnostics.svg');plt.close(fig)
    pd.DataFrame(summary).to_csv(out/'cell_type_summary.csv',index=False)
    pd.DataFrame(sensitivity).to_csv(out/'cell_threshold_sensitivity.csv',index=False)
    (out/'fit_warnings.json').write_text(json.dumps(fit_warnings,indent=2))
    provenance['runtime_seconds']=round(time.time()-started,2)
    (out/'provenance.json').write_text(json.dumps(provenance,indent=2))
    print(pd.DataFrame(summary).to_string(index=False),flush=True)


if __name__ == '__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--input',required=True);parser.add_argument('--outdir',required=True)
    args=parser.parse_args();run(args.input,args.outdir)
