<div id="main" class="col-md-9" role="main">

# Load and Render SQL from Package

<div class="ref-description section level2">

Convenience wrapper around SqlRender functions that loads SQL from the
package inst/sql directory, renders parameters, and translates to the
target dialect.

</div>

<div class="section level2">

## Usage

<div class="sourceCode">

``` r
loadRenderTranslateSql(sqlFileName, dbms, ...)
```

</div>

</div>

<div class="section level2">

## Arguments

-   sqlFileName:

    Name of the SQL file in inst/sql/sql\_server/.

-   dbms:

    Target database management system (e.g., "postgresql", "sql
    server").

-   ...:

    Named parameters to substitute into the SQL template.

</div>

<div class="section level2">

## Value

Translated SQL string.

</div>

</div>
