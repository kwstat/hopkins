# A protocol for Hopkins statistic

Hopkins statistic is used to test a Null Hypothesis of spatial
randomness. Under the null distribution of spatial randomness, the value
of the statistic should be 0.5.

## A protocol for using Hopkins statistic

The important point of this protocol is to raise awareness of potential
problems. We leave it to the practitioner to decide what do with the
answers to these questions.

1.  \*\*Does the data need to be centered and scaled to unit variance?
2.  **Is the number of events $`n > 100`$ and the number of
    randomly-sampled events at most 10% of $`n`$?** This is recommended
    by Cross and Jain (1982).
3.  **Is spatial randomness of the events even possible?** If the events
    are known or suspected to be correlated, this violates the null
    hypothesis of spatial uniformity, and may also mean that the
    sampling frame is not shaped like a hypercube.
4.  **Could nearest-neighbor events have occurred outside the boundary
    of the sampling frame?** If yes, it may be appropriate to calculate
    nearest-neighbor distances using a torus geometry.
5.  **Is the sampling frame non-rectangular?** If yes, then be extremely
    careful with the use of Hopkins statistic in how points are samples
    from $`U`$.
6.  **Is the dimension of the data much greater than 2?** Edge effects
    are more common in higher dimensions.
7.  **Is the dimension of the data large?** It is difficult to detect
    spatial randomness in a large number of dimensions.

#### Example 1

We can simulate 1000 points uniformly in a unit square and then
calculate Hopkins statistic.

``` r

library(hopkins)
set.seed(42)
dat1 <- matrix(runif(2000), ncol=2)
plot(dat1)
```

![](hopkins_protocol_files/figure-html/unnamed-chunk-1-1.png)

``` r

hopkins(dat1) # .52
```

    ## [1] 0.5217802

``` r

hopkins.pval(.52, 200)
```

    ## [1] 0.4238149

The value of the statistic is 0.52, which is not far from 0.5. The
p-value of the test is non-significant. Spatial uniformity of the points
is possible.

#### Example 2

Simulate 1000 points from a bivariate normal distribution (with 0
covariance). The sampling frame for generating new points $`u`$ is from
the minimum value to maximum value of the events in each axis. Roughly
-3 to 3 for Normal data. The points form a circular “cluster” within
this bounding box.

``` r

set.seed(42)
dat2 <- matrix(rnorm(1000), ncol=2)
plot(dat2)
```

![](hopkins_protocol_files/figure-html/unnamed-chunk-2-1.png)

``` r

hopkins(dat2) # .89
```

    ## [1] 0.8927865

``` r

hopkins.pval(.89, 100)
```

    ## [1] 0

The value of Hopkins statistic is 0.89. This value is far from 0.5. The
p-value of the test is significant. The distribution of points is
non-uniform.

Cross, GR, and AK Jain. 1982. “Measurement of Clustering Tendency.” In
*Theory and Application of Digital Control*. Elsevier.
<https://doi.org/10.1016/B978-0-08-027618-2.50054-1>.
