# Load and Render SQL from Package

Convenience wrapper around SqlRender functions that loads SQL from the
package inst/sql directory, renders parameters, and translates to the
target dialect.

## Usage

``` r
loadRenderTranslateSql(sqlFileName, dbms, ...)
```

## Arguments

- sqlFileName:

  Name of the SQL file in inst/sql/sql_server/.

- dbms:

  Target database management system (e.g., "postgresql", "sql server").

- ...:

  Named parameters to substitute into the SQL template.

## Value

Translated SQL string.
