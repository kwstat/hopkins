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
  text(X1[c(1,5,6),], labels=c("W      ", "W       ", "W       "), col="blue")
  text(U, labels=c("U1","U2","U3"), col="red")
  debugonce(hopkins)
  hopkins(X1, m=3, U=U)
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
  # Would be nice to have a hand-calculated torus check
  #expect_equal(hopkins(X1, m=3, U=U,method="torus"),
  #             0.75)
  #set.seed(2)
  #expect_equal(hopkins(iris[, -5], m=15, method="simple"),
  #             0.9952293289 )
})

test_that("japanesepines data", {
  # This example is used in the paper.  We use it to compare simple/torus.
  
  jpines <- data.frame(
    x = c(
    0.09, 0.29, 0.38, 0.39, 0.48, 0.59, 0.65, 
    0.67, 0.73, 0.79, 0.86, 0.89, 0.98, 0.02, 0.11, 0.42, 0.48, 0.62, 
    0.73, 0.89, 0.02, 0.03, 0.07, 0.52, 0.64, 0.08, 0.08, 0.12, 0.12, 
    0.17, 0.31, 0.32, 0.42, 0.52, 0.91, 0.94, 0.34, 0.37, 0.47, 0.52, 
    0.59, 0.66, 0.76, 0.73, 0.89, 0.94, 0.98, 0.97, 0.12, 0.11, 0.17, 
    0.21, 0.29, 0.32, 0.35, 0.39, 0.52, 0.58, 0.69, 0.77, 0.36, 0.36, 
    0.39, 0.43, 0.62),
    y = c(
      0.09, 0.02, 0.03, 0.18, 0.03, 0.02, 
      0.16, 0.13, 0.13, 0.03, 0.13, 0.08, 0.02, 0.18, 0.31, 0.22, 0.13, 
      0.21, 0.23, 0.23, 0.41, 0.44, 0.42, 0.42, 0.43, 0.59, 0.63, 0.63, 
      0.66, 0.58, 0.53, 0.52, 0.49, 0.52, 0.52, 0.58, 0.68, 0.68, 0.67, 
      0.67, 0.67, 0.68, 0.66, 0.73, 0.74, 0.78, 0.79, 0.86, 0.84, 0.94, 
      0.95, 0.79, 0.84, 0.83, 0.86, 0.79, 0.93, 0.83, 0.93, 0.93, 0.97, 
      0.96, 0.96, 0.96, 0.97))
  
  set.seed(28)
  # (.023^2+.076^2+.07^2) /
  #   ((.023^2+.076^2+.07^2) + (.104^2+.1^2+.058^2)) # .32
  expect_equal( round(hopkins(jpines, m=3, method="simple"),6),
               0.313107)
  set.seed(28)
  # (.0374^2+.0225^2+.0697^2) /
  #   ((.0374^2+.0225^2+.0697^2)+(.0583^2+.1044^2+.1000^2)) = 0.218
  expect_equal(round(hopkins(jpines, m=3, method="torus"), 6),
               0.217831)
})

test_that("Hopkins p-value", {
  # These two tests adapted from Gastner 2005, page 121
  expect_equal(hopkins.pval(0.21, 10),
               0.00466205,  tolerance=1e-6)
  expect_equal(hopkins.pval(.84, 10),
               0.0004981939, tolerance=1e-6)
})
