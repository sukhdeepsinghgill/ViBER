# ViBEr: release-readiness review

Source snapshot reviewed 2026-11-02. Scientific C source is retained unchanged. This audit distinguishes release defects from the published method's validation.

## Blocking findings

1. **Allocation uses an uninitialized size.** `4_fbe/viber.c:171–175` allocates `num`, `binu`, `binps`, `logU`, `logPU` with `NB`, but `NB=15` is not assigned until line 266. This is undefined behavior and may allocate the wrong size before later loops write to the arrays. Move `NB=15` ahead of all five allocations; also initialize the associated range before use. Rebuild and run with address/undefined-behavior sanitizers on a small compatible fixture.
2. **No active argument-count check for viber.** The old usage block is commented out, but `argv[1]` and `argv[2]` are read at lines 87–88. Add `argc==3`, validate the realization and use bounded string copying. A no-argument launch can crash.
3. **Input data missing.** `2b_visfitsgrid/make_copies.sh` expects `single_chan.fits`, absent from the archive. `2a_MWA_baselines/` is empty. The code cannot generate telescope sampling from nothing; provide a redistributable template or an exact data-acquisition/preparation recipe.
4. **Broken top-level script and inconsistent ranges.** `pipeline.sh` enters `1_ngrf` and then attempts `cd 2b_visfitsgrid` without returning to the root. Generator produces realizations 1–50; injection runs 1–100; template copying runs 201–500; gridding and estimator scripts request up to 500. The pipeline also does not run the final bispectrum stage. Replace with a path-anchored, fail-fast single-realization driver after input validation.
5. **Hard-coded binary geometry.** `viber.c` fixes `N1=N2=1176`, `Du=1.0254384902642049`; the active binary reader has no header. Any changed gridding parameter needs a matching estimator geometry. Add a versioned header or at least validate file sizes against grid dimensions before reading.

## Numerical interpretation and robustness

- `Bispec /= Num` is performed without a zero-support guard around line 414. Record empty bins explicitly instead of dividing by zero.
- The output is the FFT-ratio visibility statistic; the brightness-temperature angular bispectrum additionally needs the appropriate beam/Q factors. Do not label raw output as mK³ without checking the selected analysis path.
- The analytic comparison `P_T` is a fixed test model. It is not an automatic prediction for observed foregrounds or an arbitrary input field.
- `Area` and some diagnostic metadata are printed without initialization in the active path; misleading logs are possible even if these variables do not feed the bispectrum result.
- Check all fopen/fread/allocation returns. Current shared intermediate names are unsuitable for concurrent runs in one directory.
- Source comments describe several gridding/tapering paths. Confirm which routine is active before claiming a particular taper kernel is applied.

## Packaging defects handled by the root build

Legacy default targets include absent sources and cleanup operations; `3_vcg`'s first target is merely one object, not the executable. Compiled ELF binaries and `.o` objects were bundled from another machine. The release copy removes these and adds a root build with explicit source lists. The stage scripts remain historical examples.

Notebook sources were rescued from checkpoint directories into `legacy_notebooks/`. At least one imports `cldnufuncs`, which is not supplied. These are not self-contained, tested examples.

## Acceptance test before release

Supply a tiny UVFITS fixture with known sampling, verify dimensions and visibility/sampling binary sizes, then compare a few ring-triplet sums against direct enumeration. Test filled and empty cells, zero signal, and a missing input. Only after that compare the recovered angular bispectrum and unit conversion to the known input model. Verify the general and equilateral implementations agree on identical bins.
