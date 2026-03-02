# Consistent plotting theme for Medusa package

A clean ggplot2 theme applied to all plots produced by the Medusa
package. Based on `theme_minimal` with consistent font sizes, axis
formatting, and color scheme.

## Usage

``` r
mrTheme(baseSize = 12)
```

## Arguments

- baseSize:

  Base font size in points. Default is 12.

## Value

A ggplot2 theme object.

## Details

Medusa ggplot2 Theme

## Examples

``` r
library(ggplot2)
#> Warning: package ‘ggplot2’ was built under R version 4.5.2
ggplot(data.frame(x = 1:10, y = rnorm(10)), aes(x, y)) +
  geom_point() +
  mrTheme()

```
