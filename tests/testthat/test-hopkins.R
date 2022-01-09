# test-hopkins.R

require(hopkins)

X1 <- matrix(c(7,5,
              3,3,
              8,4,
              1,2,
              4,3,
              7,4), ncol=2, byrow=TRUE)
X2 <- matrix(c(3,3,
               3,5,
               5,2,
               5,4,
               5,6,
               7,4), ncol=2, byrow=TRUE)
U <- matrix(c(4,2,
              3,5,
              7,2), ncol=2, byrow=TRUE)

if(FALSE){
  # Hand-calculate Hopkins statistic
  plot(NA, xlim=c(0,10), ylim=c(0,10))
  text(X1, labels=1:nrow(X1))
  text(U, col="red")
  debugonce(hopkins(X1, m=3, U=U))
  ## U               nearest X1 dist^2
  ## [1,]   4  2     4  3       sqrt((4-4)^2+(2-3)^2)^2 = 1
  ## [2,]   3  5     3  3       sqrt((3-3)^2+(5-3)^2)^2 = 4
  ## [3,]   7  2     7  4       sqrt((7-7)^2+(2-4)^2)^2 = 4
  ## W=X1[c(1,5,6),] nearest X1  dist^2
  ## [1,]   7  5     7  6       sqrt((7-7)^2+(5-6)^2)^2 = 1
  ## [2,]   4  3     3  3       sqrt((4-3)^2+(3-3)^2)^2 = 1
  ## [3,]   7  4     7  5       sqrt((7-7)^2+(4-5)^2)^2 = 1
  # sum(dux^d) / sum( dux^d + dwx^d ) = (1+4+4)/(1+4+4 + 1+1+1) = 0.75
}

test_that("Hopkins statistic", {
  set.seed(42) # this seed samples rows 1,5,6
  expect_equal(hopkins(X1, m=3, U=U),
               0.75)
  set.seed(1)
  expect_equal(hopkins(iris[, -5], m=15),
               0.9952293289 )

})

test_that("Hopkins p-value", {
  # These two tests adapted from Gastner 2005, page 121
  expect_equal(hopkins.pval(0.21, 10),
               0.00466205,  tolerance=1e-6)
  expect_equal(hopkins.pval(.84, 10),
               0.0004981939, tolerance=1e-6)
})
