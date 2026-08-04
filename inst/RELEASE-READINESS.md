# Release readiness

## Technical candidate

The source tree is on the 0.9.2 pre-release development line; see
`RELEASE-HANDOFF.md`. The intended first stable release is 1.0.0. No release tag
or publication is created from development. Its local technical gate is:

```sh
Rscript tools/release-audit.R
Rscript -e 'testthat::test_local(reporter = "progress")'
R CMD build .
R CMD check --no-manual distmer_0.9.2.tar.gz
R CMD check --as-cran --no-manual -o /tmp/distmer-cran distmer_0.9.2.tar.gz
```

The first four commands must have no failures or warnings. An offline
`--as-cran` run may report only environment notes that it could not contact
CRAN or verify the system clock; those are not package findings. Online CI is
the authoritative incoming/dependency check.

The GitHub Actions workflow runs the structural audit and `--as-cran` check on
release R for Linux, macOS, and Windows. A local green run does not substitute
for all three hosted jobs.

## Publication metadata

The owner-supplied repository metadata is present:

- `DESCRIPTION` contains the canonical project `URL`;
- `DESCRIPTION` contains the public `BugReports` URL.

The strict metadata audit must remain green:

```sh
STRICT_RELEASE_METADATA=true Rscript tools/release-audit.R
```

## Scope declaration

The development tree implements only the support envelopes documented in the
README and vignettes. Random-slope and simultaneous-inference guarantees apply
only within their explicitly documented exchangeability and geometry
conditions; unsupported extensions must continue to fail explicitly.

The package is an independent clean-room implementation. It does not contain
or redistribute proprietary PERMANOVA or Numerical Recipes source.
