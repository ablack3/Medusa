<div id="main" class="col-md-9" role="main">

# Consistent plotting theme for Medusa package

<div class="ref-description section level2">

A clean ggplot2 theme applied to all plots produced by the Medusa
package. Based on `theme_minimal` with consistent font sizes, axis
formatting, and color scheme.

</div>

<div class="section level2">

## Usage

<div class="sourceCode">

``` r
mrTheme(baseSize = 12)
```

</div>

</div>

<div class="section level2">

## Arguments

-   baseSize:

    Base font size in points. Default is 12.

</div>

<div class="section level2">

## Value

A ggplot2 theme object.

</div>

<div class="section level2">

## Details

Medusa ggplot2 Theme

</div>

<div class="section level2">

## Examples

<div class="sourceCode">

``` r
library(ggplot2)
#> Warning: package ‘ggplot2’ was built under R version 4.5.2
ggplot(data.frame(x = 1:10, y = rnorm(10)), aes(x, y)) +
  geom_point() +
  mrTheme()

```

</div>

</div>

</div>
