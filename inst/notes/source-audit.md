# Source-code audit

## Finding

No publicly indexed copy of the original FORTRAN source has yet been located.
The recoverable 2005 manual describes a compiled standalone distribution named
`PERMANOVA.exe` for Windows and `PERMANOVA.mac` for Macintosh.

The manual explicitly states that:

- intellectual property and copyright in the source code and program remained
  with Marti J. Anderson;
- the program depended on Numerical Recipes routines `tred2`, `tqli`,
  `piksr2`, and `ran2`;
- it also used Applied Statistics algorithms AS91, AS111, AS239, AS66, and
  AS245.

This matters because a redistributable reconstruction should be clean-room and
should not reproduce proprietary Numerical Recipes source.

## Confirmed legacy scope

The standalone 2005 version documented in the manual is PERMANOVA v1.6. It:

- accepts raw response data or an N x N distance matrix;
- handles balanced designs only: equal observations per cell and no missing
  cells;
- permits fixed/random and crossed/nested factor declarations;
- constructs term-specific pseudo-F ratios using numerator and denominator
  terms;
- permutes units identified by the denominator term;
- offers raw-data, reduced-model residual, and full-model residual
  permutations;
- reports permutation and Monte Carlo P-values.

This legacy scope must be kept separate from the later PERMANOVA+ for PRIMER
implementation, which includes unbalanced-design machinery and additional
sums-of-squares choices.

## Search trail

1. University of Auckland 2005 annual report: identifies the program and the
   original `Programs.htm` distribution page.
2. Recoverable 24-page 2005 user manual: documents executable names, v1.6,
   algorithms, limitations, and dependencies.
3. Searches of indexed GitHub, Zenodo, Figshare, Archive.org results, and the
   open web did not return a source archive.

## Next source-recovery actions

- Query Internet Archive CDX records for every historical URL beneath
  `stat.auckland.ac.nz/~mja/` and `~mja3/`.
- Recover the old `Programs.htm` link targets and filenames.
- Search institutional mirrors and old course/workshop bundles for the Windows
  and Macintosh distributions.
- Contact the copyright holder before requesting or redistributing source.
- If only an executable is recovered, use it solely as a black-box validation
  oracle unless explicit permission allows deeper inspection.
