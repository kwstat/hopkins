# Calculate the p-value for Hopkins statistic

Calculate the p-value for Hopkins statistic

## Usage

``` r
hopkins.pval(x, m)
```

## Arguments

- x:

  Observed value of Hopkins statistic

- m:

  Number of events/points sampled.

## Value

A p-value between 0 and 1.

## Details

Under null hypothesis of spatial randomness, Hopkins statistic has a
Beta(m,m) distribution, where 'm' is the number of events/points
sampled. This function calculates the p-value for the statistic.

## References

Michael T. Gastner (2005). Spatial distributions: Density-equalizing map
projections, facility location, and two-dimensional networks. Ph.D.
dissertation, Univ. Michigan (Ann Arbor, 2005).
http://hdl.handle.net/2027.42/125368

## Author

Kevin Wright

## Examples

``` r
hopkins.pval(0.21, 10) # .00466205
#> [1] 0.00466205
```
