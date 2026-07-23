# Version policy

`distmer` follows semantic versioning for public package releases and the
four-component `.9000` convention for development snapshots.

## Current path to the first stable release

- `0.9.0.9000` is the active, pre-release development tree.
- `1.0.0` will be the first stable public release, after metadata, hosted
  platform checks, adversarial review, and the release gate are complete.
- A `.9000` build is never tagged or submitted as a stable release.

The source was renumbered from the former internal `2.9.1.9000` checkpoint to
`0.9.0.9000`. No published versions were withdrawn or superseded because the
earlier numbers were implementation checkpoints rather than package releases.
They remain labelled as such in `NEWS.md` and `VALIDATION.md` where needed for
artifact and checksum traceability.

## After 1.0.0

- increment **PATCH** for compatible bug fixes (`1.0.1`);
- increment **MINOR** for backward-compatible capabilities (`1.1.0`);
- increment **MAJOR** for deliberate breaking changes (`2.0.0`).

Development after a release uses a `.9000` suffix until the next release
number is chosen. Internal phases, validation batches, and statistical feature
milestones belong in issue, NEWS, or validation labels and must not drive the
package `Version` field.

## Release mechanics

For a stable release, set the three-component version in `DESCRIPTION` and the
first `NEWS.md` heading, rebuild documentation and vignettes, run the complete
release gate, and tag the exact checked source as `v<version>` (for example,
`v1.0.0`). Immediately after release, development may resume on a `.9000`
version.
