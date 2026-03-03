<div id="main" class="col-md-9" role="main">

# Import a Site Profile from CSV

<div class="ref-description section level2">

Reads a site profile from CSV files previously written by
`exportSiteProfile`. Reconstructs the profile list object that can be
passed to `poolLikelihoodProfiles`.

</div>

<div class="section level2">

## Usage

<div class="sourceCode">

``` r
importSiteProfile(profilePath, metadataPath = NULL)
```

</div>

</div>

<div class="section level2">

## Arguments

-   profilePath:

    Character. Path to the profile CSV file (the file containing beta
    and log\_likelihood columns).

-   metadataPath:

    Character. Path to the metadata CSV file. If NULL (default), the
    function infers the path by replacing "\_profile\_" with
    "\_metadata\_" in `profilePath`.

</div>

<div class="section level2">

## Value

A list with the same structure as `fitOutcomeModel` output, suitable for
passing to `poolLikelihoodProfiles`.

</div>

<div class="section level2">

## See also

<div class="dont-index">

`exportSiteProfile`, `poolLikelihoodProfiles`

</div>

</div>

<div class="section level2">

## Examples

<div class="sourceCode">

``` r
if (FALSE) { # \dontrun{
profile <- importSiteProfile("medusa_profile_site_A.csv")
} # }
```

</div>

</div>

</div>
