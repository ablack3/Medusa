# Default FeatureExtraction settings for Medusa

Creates a FeatureExtraction covariate settings object tailored for
Mendelian Randomization analysis. Includes demographics, conditions,
drug exposures, and measurements in standard lookback windows.

## Usage

``` r
createDefaultMRCovariateSettings()
```

## Value

A FeatureExtraction `covariateSettings` object.

## Details

Create Default Covariate Settings for MR Analysis

## See also

[`buildMRCovariates`](buildMRCovariates.md)

## Examples

``` r
settings <- createDefaultMRCovariateSettings()
```
