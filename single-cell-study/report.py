"""Create Markdown and standalone HTML reports from the actual analysis outputs."""
import argparse
import base64
import html
import json
import re
from pathlib import Path
import pandas as pd


def markdown_table(df):
    rows=[[str(x) for x in df.columns]]+[[str(x) for x in row] for row in df.itertuples(index=False,name=None)]
    return '\n'.join(['| '+' | '.join(rows[0])+' |', '| '+' | '.join(['---']*len(df.columns))+' |']+
                     ['| '+' | '.join(r)+' |' for r in rows[1:]])


def render(outdir):
    out=Path(outdir)
    prov=json.loads((out/'provenance.json').read_text())
    summary=pd.read_csv(out/'cell_type_summary.csv').astype(object).fillna('—')
    compact=summary[['cell_type','donor_pairs','tested_genes','padj_below_005']].copy()
    compact.columns=['Cell type','Complete donor pairs','Genes tested','Within-type adjusted p < 0.05']
    for col in compact.columns[1:]:
        compact[col]=compact[col].map(lambda x:str(int(x)) if isinstance(x,(float,int)) else x)
    sensitivity=pd.read_csv(out/'cell_threshold_sensitivity.csv')
    examples=[]
    for ct in ['CD14+ Monocytes','CD4 T cells']:
        slug=re.sub(r'[^a-z0-9]+','_',ct.lower()).strip('_')
        r=pd.read_csv(out/f'{slug}_differential_expression.csv',index_col=0)
        for gene in ['ISG15','IFIT1','MX1']:
            examples.append({'Cell type':ct,'Gene':gene,'Log2 fold change (stim/ctrl)':f'{r.loc[gene,"log2FoldChange"]:.2f}'})
    example=pd.DataFrame(examples)
    sections=[
      ('Question and scope', 'How does the transcriptional response to ex vivo interferon-beta stimulation vary across PBMC cell types, and is the response consistent across donors? This is a reanalysis of a published processed dataset, not a new clinical discovery or an independent replication cohort.'),
      ('Data and QC', f'The checksum-verified processed Kang export contains {prov["source_dimensions"][0]:,} cells and {prov["source_dimensions"][1]:,} genes. A minimum of 200 detected genes retained {prov["retained_cells"]:,} cells (111 excluded). All eight donors have ctrl and stim observations. The nonnegative integer matrix matches recorded RNA library totals. Cell-type labels are inherited from the processed export. Mitochondrial genes and multiplet flags are absent, so mitochondrial QC, original doublet filtering and annotation accuracy cannot be independently checked here. No integrated or SCT matrix is used for count modeling.'),
      ('Methods', 'Sum raw counts within each donor × condition × cell type. Retain only complete donor pairs with at least 20 cells in both conditions; require three or more pairs. Fit PyDESeq2 0.5.4 with design ~ donor + condition, contrasting stim versus ctrl. Keep genes with count >=10 in at least three pseudobulk samples. Use BH adjusted p-values within each cell type at 0.05, independent filtering and Cooks filtering. Effect sizes are unshrunken. Donors, not cells, supply biological replication. Python and R implementations are separate; numerical identity with R DESeq2 is not asserted.'),
      ('Observed results', 'Seven cell types met the paired-sample rule; megakaryocytes were excluded for insufficient complete pairs. The table describes model outputs, not a ranking of biological response strength: donor counts, cell abundance and statistical power differ across cell types. Multiplicity across cell types is not controlled by these within-type adjusted p-values.'),
      ('Illustrative expression changes', 'ISG15, IFIT1 and MX1 show positive fitted stimulated/control effects in CD14+ monocytes and CD4 T cells. These examples were selected for illustration after the analysis; no gene-set enrichment or pathway significance test was performed. Among each of these two cell types’ ten lowest-adjusted-p genes, all ten have positive descriptive log-normalized count differences in all eight donors. Selection from the same dataset means this is a consistency diagnostic, not independent validation.'),
      ('Sensitivity and model diagnostics', 'Repeat donor inclusion at 10, 20 and 50 cells per stratum. CD14+ monocyte and CD4 T-cell inclusion is unchanged at all three thresholds. Rarer populations lose pairs as the threshold rises. Where samples changed and fitting remained possible, log2-effect correlations versus primary results ranged from 0.969 to 0.994 on intersecting tested genes. Correlation does not establish identical significant-gene lists. Dendritic cells have fewer than three complete pairs at threshold 50. The threshold-10 dendritic fit falls back to a mean dispersion trend; threshold-50 CD8 and FCGR3A monocyte fits warn about low residual degrees of freedom. Primary fits emitted dependency deprecation warnings but no recorded statistical fitting warning. See fit_warnings.json for the complete record. Leave-one-donor-out tables summarize mean descriptive log-count differences; they do not refit differential-expression models or demonstrate significance robustness.'),
      ('Interpretation and limits', 'The observed response is consistent across donors for the selected prominent genes in the two abundant populations shown. It does not establish clinical efficacy, disease specificity, a biomarker, or that one lineage responds more strongly than another. The experiment compares treated and untreated cells from this lupus cohort; it is not lupus versus healthy. Condition-associated processing effects cannot be separated without additional experimental metadata. Large global responses can affect relative normalization: apparent decreases are relative, not absolute RNA changes. Processed gene filtering, inherited labels, missing upstream QC, small donor numbers, untested normalization sensitivity and lack of an external cohort limit inference. A formal cell-type interaction model and pathway enrichment remain future extensions.'),
      ('Reproducibility', 'Source: https://api.figshare.com/v2/file/download/34464122. SHA-256: '+prov['sha256']+'. The checksum matches the upstream LIANA dataset registry. The input is downloaded separately and is not redistributed with this report. The bundle includes scripts, a pinned package snapshot, complete gene tables, pseudobulk counts, donor audits, diagnostics and sensitivity outputs. Synthetic-effect, aggregation and invalid-input tests pass locally. The full observed study ran on macOS ARM64 with Python 3.12; other platforms are not yet verified. A reproduction run must be interpreted with its warnings and environment record.'),
    ]
    refs=[('Original study','https://doi.org/10.1038/nbt.4042'),
          ('Upstream processed-data registry','https://github.com/saezlab/liana-py/blob/main/src/liana/datasets/registry.yaml'),
          ('Pertpy source documenting the Figshare export','https://github.com/scverse/pertpy/blob/0.9.5/pertpy/data/_datasets.py'),
          ('PyDESeq2 methods','https://pydeseq2.readthedocs.io/en/stable/auto_examples/plot_minimal_pydeseq2_pipeline.html')]
    md=['# Donor-paired PBMC interferon-response reanalysis','Run date: 2026-10-23. Published-data portfolio case study.']
    body=['<header><p class="eyebrow">COMPUTATIONAL BIOLOGY · PORTFOLIO CASE STUDY</p><h1>Donor-paired PBMC interferon response</h1><p>Published-data reanalysis · 8 donors · 7 analyzed cell types · 23 October 2026</p></header>']
    for title,text in sections:
        md += ['','## '+title,'',text]
        body += ['<section><h2>'+title+'</h2><p>'+html.escape(text)+'</p>']
        if title=='Observed results':
            md += ['',markdown_table(compact)];body += [compact.to_html(index=False,escape=True)]
        if title=='Illustrative expression changes':
            md += ['',markdown_table(example)];body += [example.to_html(index=False,escape=True)]
        if title=='Sensitivity and model diagnostics':
            for ct in ['cd14_monocytes','cd4_t_cells']:
                img=out/f'{ct}_diagnostics.png'
                md += ['',f'![{ct} diagnostics]({ct}_diagnostics.svg)']
                data=base64.b64encode(img.read_bytes()).decode()
                body += [f'<img alt="{ct} donor-level diagnostics" src="data:image/png;base64,{data}">']
        body += ['</section>']
    md += ['','## References','']+[f'- [{name}]({url})' for name,url in refs]
    body += ['<section><h2>References</h2><ul>']+[f'<li><a href="{url}">{name}</a></li>' for name,url in refs]+['</ul></section>']
    (out/'REPORT.md').write_text('\n'.join(md)+'\n')
    css='body{font:17px/1.65 system-ui,sans-serif;color:#17313b;background:#f5f7f8;margin:0}main{max-width:1050px;margin:auto;padding:35px}header,section{background:white;padding:28px 34px;margin:18px 0;border-radius:12px}header{border-top:7px solid #247b83}h1{font-size:37px;line-height:1.15}h2{font-size:24px;color:#17666e}.eyebrow{font-size:12px;letter-spacing:.12em}table{border-collapse:collapse;width:100%;font-size:14px}th,td{text-align:left;padding:10px;border-bottom:1px solid #dce5e7}img{width:100%;height:auto;margin:20px 0}a{color:#17666e}p{overflow-wrap:anywhere}@media(max-width:650px){main{padding:8px}header,section{padding:18px}table{font-size:11px}h1{font-size:29px}}'
    (out/'REPORT.html').write_text('<!doctype html><html lang="en"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>Donor-paired PBMC interferon response</title><style>'+css+'</style><main>'+''.join(body)+'</main></html>')


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--outdir',required=True);a=p.parse_args();render(a.outdir)
