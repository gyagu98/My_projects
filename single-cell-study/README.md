# Executed donor-paired PBMC case study

This directory extends the original PBMC3k demonstration with an executed
multi-donor reanalysis of a processed Kang et al. interferon-beta dataset.
Start with [the observed-data report](REPORT.md).

## Reproduce (Python 3.12)

From this directory, create an isolated environment and install the exact
tested package snapshot. The snapshot was executed on macOS ARM64;
other operating systems have not yet been verified.

```bash
python3.12 -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements-lock.txt
python -m unittest -v test_analyze.py
python download.py --output data/kang_2018.h5ad
python analyze.py --input data/kang_2018.h5ad --outdir results
python report.py --outdir results
```

The downloader uses Figshare's public API and checks SHA-256 against the
upstream registry. The original dataset is not committed. Re-running creates
full differential-expression and count tables locally in results, along with
diagnostic images. GitHub includes compact summaries, audits and figures;
the full local bundle also includes complete gene tables and compressed
pseudobulk/normalized counts.

## Scientific boundaries

Donor pairing is preserved in `~ donor + condition`; cells are not independent
replicates. Conditions are `ctrl` and `stim`, and positive effects mean higher
relative expression after stimulation. Raw-count integers and source library
totals are checked. Cell labels are inherited, not newly validated.
Mitochondrial genes and multiplet flags are missing from the processed export.
This is a published-data reanalysis, not a new discovery or a clinical model.

Threshold sensitivity and descriptive donor/leave-one-donor-out summaries
are included. The latter do not refit DE models. No pathway enrichment,
external cohort validation, treatment-by-cell-type interaction testing or
normalization sensitivity analysis is claimed. Consult the report and the
recorded fit warnings before interpreting differences.
