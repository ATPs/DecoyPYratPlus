# Repository Guidance

## Scope

This repository now uses the `DecoyPYratPlus/` package directory. Keep packaging, README examples, and command references aligned with that module path.

## Environment

- Use `zsh` for commands.
- When validating changes, activate the base conda environment first:
  `source /data/p/anaconda3/bin/activate base`
- Add the local tools path before running checks:
  `export PATH=/data/p/bin:$PATH`

## Code Changes

- Prefer small, targeted edits over cleanup passes.
- Keep both execution modes working:
  `python -m DecoyPYratPlus.decoyPYratPlus`
  `python DecoyPYratPlus/decoyPYratPlus.py`
- Treat `numpy` as required. `tqdm` and the Cython fast-digest extension are optional.
- If packaging changes touch the module name, update `setup.py`, `setup_fast_digest.py`, `meta.yaml`, and `README.md` together.

## Verification

- For CLI changes, run `-h` and at least one minimal FASTA smoke test.
- For import-path changes, verify both module execution and direct script execution.
- If Cython or packaging changes are involved, check whether `python setup_fast_digest.py build_ext --inplace` still makes sense before editing the build instructions.

## Repo Hygiene

- Do not commit caches, `*.so`, `build/`, `dist/`, or `*.egg-info/`.
- Keep release notes under `release/YYYYMMDD.md`.
