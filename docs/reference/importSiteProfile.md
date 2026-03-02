# Import a Site Profile from CSV

Reads a site profile from CSV files previously written by
[`exportSiteProfile`](exportSiteProfile.md). Reconstructs the profile
list object that can be passed to
[`poolLikelihoodProfiles`](poolLikelihoodProfiles.md).

## Usage

``` r
importSiteProfile(profilePath, metadataPath = NULL)
```

## Arguments

- profilePath:

  Character. Path to the profile CSV file (the file containing beta and
  log_likelihood columns).

- metadataPath:

  Character. Path to the metadata CSV file. If NULL (default), the
  function infers the path by replacing "\_profile\_" with
  "\_metadata\_" in `profilePath`.

## Value

A list with the same structure as
[`fitOutcomeModel`](fitOutcomeModel.md) output, suitable for passing to
[`poolLikelihoodProfiles`](poolLikelihoodProfiles.md).

## See also

[`exportSiteProfile`](exportSiteProfile.md),
[`poolLikelihoodProfiles`](poolLikelihoodProfiles.md)

## Examples

``` r
if (FALSE) { # \dontrun{
profile <- importSiteProfile("medusa_profile_site_A.csv")
} # }
```
