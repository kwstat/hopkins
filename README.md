# hopkins <img src="man/figures/logo.png" align="right" />

[![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/hopkins)](https://cran.r-project.org/package=hopkins)
[![CRAN_Downloads](https://cranlogs.r-pkg.org/badges/hopkins)](https://cranlogs.r-pkg.org/badges/hopkins)

Homepage: https://kwstat.github.io/hopkins

Repository: https://github.com/kwstat/hopkins

## Key features

The `hopkins()` function:

+ Can assess the spatial distribution of points for randomness.
+ Correctly calculates Hopkins statistic using exponents d=D.
+ Is very fast for `method="simple"`.
+ Can use `method="torus"` to correct for edge effects in higher dimensions.

We believe that most software implementations of Hopkins statistics are incorrect because they use "d=1" as the exponent of the distances in the calculation instead of "d=2" for 2-dimensional data or the more general "d=D" for D-dimensional data.

For full details, see Wright (2023), "Will the Real Hopkins Statistic Please Stand Up", https://journal.r-project.org/articles/RJ-2022-055/

## Installation

```R
# Install the released version from CRAN:
install.packages("hopkins")

# Install the development version from GitHub:
install.packages("devtools")
devtools::install_github("kwstat/hopkins")
```

## Usage

```R
library(hopkins)
hopkins(iris[, -5], m=15) # .9952293
hopkins.pval(0.21, 10) # .00466205
```

## About the logo

The Hopkins statistic attempts to assess the spatial distribution of points. Prairie dogs cluster together in burrows, while the burrows are spaced somewhat evenly apart from each other. Also, prairie dogs (sort of) "hop", and "kins" is both "endearing" and "a group with similar characteristics". So "hopkins" could mean "a cute group of hoppers".

Thanks to TW for the artwork.
