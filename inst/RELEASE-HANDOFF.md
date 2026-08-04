# Development release handoff

The active source is the 0.9.2 development tree, currently in pre-release.
The intended first stable release is 1.0.0. It is not yet a tagged or submitted
release candidate. Historical sealed artifacts retain their original internal
checkpoint labels; their checksums do not validate this development version.
The release numbering rules are recorded in `VERSIONING.md`.

## Owner-supplied values

Two owner identity fields are present:

1. `Authors@R` names Roff George and a working maintainer email;
2. `LICENSE` names Roff George as copyright holder.

The repository owner has supplied:

1. the canonical project `URL` in `DESCRIPTION`;
2. the public `BugReports` URL in `DESCRIPTION`.

These values are present and validated by the strict audit, which must pass:

```sh
STRICT_RELEASE_METADATA=true Rscript tools/release-audit.R
```

## Candidate reconstruction

To reconstruct the candidate, run:

```sh
Rscript tools/release-audit.R
Rscript -e 'testthat::test_local(reporter = "progress")'
R CMD build .
R CMD check --no-manual distmer_0.9.2.tar.gz
R CMD check --as-cran --no-manual distmer_0.9.2.tar.gz
```

Dispatch `.github/workflows/R-CMD-check.yaml` and require green release-R jobs
on Linux, macOS, and Windows. Record their URLs in `cran-comments.md`.

## Tag and publication boundary

Only the repository owner should create and push the 1.0.0 tag, publish a
hosted release, or submit to CRAN. Before doing so, set `Version: 1.0.0`, verify the source artifact
checksum, the four metadata fields, all hosted checks, and the final release
notes. Development work can continue on the 0.9.2 tree without modifying
historical sealed artifacts.
