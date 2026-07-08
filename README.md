# DecoyPYratPlus

DecoyPYratPlus generates decoy protein FASTA databases for proteomics workflows. It keeps the original DecoyPYrat pseudo-reverse approach and adds overlap-aware peptide shuffling, gzip and multi-file input, existing-decoy reuse, optional target normalization, deduplication, concatenation for multi-stage searches, and an optional Cython fast-digest path.

## Highlights

- Pseudo-reverse decoy generation that respects enzyme cleavage orientation.
- Optional `--checkSimilar` mode to reduce target/decoy peptide collisions under rules such as `I=L`, `N=D`, `Q=E`, and `GG=N`.
- Supports `.fasta`, `.fa`, `.txt`, and their `.gz` variants, with multiple input files in one run.
- Can reuse an existing decoy database with `--existing-decoy`.
- Optional target export with `--target`, `--target_I2L`, `--concat`, and `--dedup`.
- Optional Cython acceleration with `--fast_digest`.

## Installation

Clone the repository and install the Python package:

```bash
git clone https://github.com/ATPs/DecoyPYratPlus.git
cd DecoyPYratPlus
pip install -e .
```

Required dependency:

```bash
pip install numpy
```

Recommended dependency for progress reporting:

```bash
pip install tqdm
```

Conda or mamba also works:

```bash
mamba install -c conda-forge numpy tqdm
```

## Usage

Installed entry point:

```bash
DecoyPYratPlus -h
```

Module execution without installation:

```bash
python -m DecoyPYratPlus.decoyPYratPlus -h
```

Direct script execution:

```bash
python DecoyPYratPlus/decoyPYratPlus.py -h
```

## Quick Start

Generate a decoy FASTA from one or more target FASTA files:

```bash
DecoyPYratPlus proteins.fa contaminants.fa.gz -o decoy.fa
```

Enable similarity-aware overlap checking:

```bash
DecoyPYratPlus proteins.fa \
  --checkSimilar "GG=N,N=D,Q=E" \
  --min_peptide_length 7 \
  --max_peptide_length 40 \
  -o decoy.fa
```

Reuse an existing decoy database and keep target names:

```bash
DecoyPYratPlus proteins.fa \
  --existing-decoy old_decoy.fa \
  --keep_names \
  -o decoy.fa
```

Write a processed target FASTA and concatenate target plus decoy sequences:

```bash
DecoyPYratPlus proteins.fa \
  --target target.fa \
  --target_I2L \
  --concat "*" \
  -o decoy_concat.fa
```

## Key Options

- `--checkSimilar "GG=N,N=D,Q=E"` enables similarity-aware overlap checks. Replacements are applied from left to right in the order provided.
- `--all_shuffle_mimic` randomizes all decoy peptides after pseudo-reversal instead of only repairing collisions.
- `--existing-decoy FILE` reuses matching decoy entries from a prior database and reports unmatched records.
- `--target FILE` writes the target database used for comparison; combine with `--target_I2L` when target output should normalize `I` to `L`.
- `--concat SEP` appends target sequence to decoy sequence using a separator such as `*` or the cleavage residue expected by the downstream search engine.
- `--dedup` removes duplicated target sequences before decoy generation.
- `--fast_digest` enables the Cython digestion backend when it is available.

Run `DecoyPYratPlus -h` for the full argument list.

## Optional Fast Digest

`--fast_digest` uses a Cython implementation when it can be imported. Build it in place with:

```bash
python setup_fast_digest.py build_ext --inplace
```

If the extension is not available, DecoyPYratPlus falls back to the pure Python implementation automatically.

## Digestion CLI

The digestion helpers are also available as a standalone CLI:

```bash
python -m DecoyPYratPlus.digestion --sequence AKRPQK -l 2
```

In this standalone CLI, `I` and `L` are kept distinct by default. Add `--isobaric` if you want `I -> L` normalization before digestion.

Digest one or more FASTA files with missed-cleavage expansion:

```bash
python -m DecoyPYratPlus.digestion proteins.fa --method trypsin -L 2 -l 6 -M 40
```

Use Comet-style enzyme presets:

```bash
python -m DecoyPYratPlus.digestion proteins.fa --enzyme Trypsin
python -m DecoyPYratPlus.digestion proteins.fa --enzyme No_cut
python -m DecoyPYratPlus.digestion proteins.fa --enzyme Cut_everywhere
```

Manual cleavage flags still work and override the preset:

```bash
python -m DecoyPYratPlus.digestion proteins.fa --enzyme Trypsin -a ""
```

Output formats:

- `tsv` writes `header<TAB>peptide` rows such as `>sp|P12345<TAB>AK`.
- `peptide` writes one peptide per line such as `AK`.
- `fasta` writes each peptide as FASTA using `--header-template`, default `{protein_id}_{index}`.

Examples:

```bash
python -m DecoyPYratPlus.digestion proteins.fa --output-format tsv
python -m DecoyPYratPlus.digestion proteins.fa --output-format peptide
python -m DecoyPYratPlus.digestion proteins.fa --output-format fasta
python -m DecoyPYratPlus.digestion proteins.fa --output-format fasta --header-template "{protein_id}|pep{index}"
```

Supported `--header-template` placeholders are `{protein_id}`, `{index}`, and `{header}`.

## Behavior Notes

- By default the output contains decoy proteins only, written to `decoy.fa`.
- In the main `DecoyPYratPlus` decoy-generation CLI, `I` is treated as `L` unless `--no_isobaric` is set.
- With `--keep_names`, output headers keep the original target identifier as `DECOY_<target_header>`.
- With `--existing-decoy`, matching decoys are preserved, missing targets get new decoys, and decoys without targets are dropped.
- The legacy single-file DecoyPYrat script is still present as `DecoyPYratPlus/decoyPYrat.py`.

## Citation

If you use the original method, cite:

[DecoyPyrat: Fast Non-redundant Hybrid Decoy Sequence Generation for Large Scale Proteomics. J Proteomics Bioinform. 2016 Jun 27;9(6):176-180.](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4941923/)
