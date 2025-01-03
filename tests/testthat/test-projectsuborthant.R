# the exact method described in Amaral et al 2007 "Pivotal Bootstrap Methods for k-Sample Problems in Directional Statistics and Shape Analysis", Section 3.2.1
rotationmat_amaral  <- function(a, b){ #assumes a and b are unit vectors
  ab <- (a %*% b)[[1]]
  alpha <- acos(ab)
  c <- b - a*ab
  c <- c/sqrt(c %*% c)[[1]]
  A <- a%o%c - c%o%a
  Q = diag(length(a)) + sin(alpha)*A + (cos(alpha) - 1)*(a%o%a + c%o%c)
  return(Q)
}

test_that("projectsuborthant() behaves with simulated data", {
  # Plan: start with some data, increase the D by using a given rotation matrix, apply projection to get back to original data with scores
  ###
  # starting data ON SPHERE
  set.seed(34)
  X_5 <- matrix(runif(5*20), nrow = 20)
  X_5 <- X_5 / sqrt(rowSums(X_5^2))

  # rotating into a higher dimension, need to set up the basis for higher dimensions
  # for simplicity lets suppose 6
  V6 <- diag(6)
  V5 <- rbind((0.7 * V6[1, ] + 0.3 * V6[2, ]) / sqrt(0.3^2 + 0.7^2),
              V6[3:6, ])
  X <- X_5 %*% V5 #X_5 expressed in the ambient coordinates given by V6

  #I'll rotate X a smidgeon off its plane
  torotvec <- (V5[1, ] + 1E-2 * V6[2, ])
  torotvec <- torotvec/ sqrt(sum(torotvec^2))

  Q <- rotationmat_amaral(V5[1, ], torotvec)
  #check Q
  expect_equal(t(Q %*% t(V5))[2:5, ], V5[2:5, ])
  expect_equal(torotvec, drop(t(Q) %*% V5[1, ]))

  rotX <- t(t(Q) %*% t(X))
  angles <- acos(rowSums(rotX * X))
  stopifnot(angles < drop(acos(V5[1, ] %*% torotvec)))

  proj <- projectsuborthant(X = divL1(rotX), V = V6, v1 = 1, v2 = 2, w = 0.7)
  expect_equal(proj$X, divL1(X_5), tolerance = 1E-4, ignore_attr = TRUE) #expect the projection of a small rotation to be very close to the original
  expect_equal(abs(proj$scores), angles, tolerance = 1E-4, ignore_attr = TRUE) #expect the projection of a small rotation to be very close to the original
  expect_equal(proj$V, V5)
  expect_equal(proj$mergedirection, V6[2, ] - V6[1, ])

  #expect when passed X, the result to be X, with scores of zero
  proj <- projectsuborthant(X = divL1(X), V = V6, v1 = 1, v2 = 2, w = 0.7)
  expect_equal(proj$X, divL1(X_5), ignore_attr = TRUE)
  expect_equal(abs(proj$scores), rep(0, nrow(X)), tolerance = 1E-7, ignore_attr = TRUE)
})

test_that("project.subsphere works on fake data", {
  # starting data ON SPHERE
  lowX <- matrix(runif(4*3), nrow = 3)
  lowX <- lowX / sqrt(rowSums(lowX^2))
  #treat these points as if they have basis vectors of
  basis <- rbind(c(0.2, 0.8, 0, 0,0) /(sqrt(0.2^2 + 0.8^2)),
                 diag(5)[-(1:2), ])
  # so in the natural basis of 1,0,0,0, etc lowX is
  lowX_5 <- lowX %*% basis
  w <- 0.2
  sizew1mw <- sqrt((w-1)^2 + w^2)
  nvec <- ((w-1)/sizew1mw) * diag(5)[1, ] + (w/sizew1mw) * diag(5)[2, ]
  stopifnot(sum(nvec*basis[1, ]) < 1E-5)

  expect_equal(project.subsphere(lowX_5, nvec), lowX_5)

  #if I rotate everything to be parallel to nvec, then the projection should be zero in the basis 1 direction
  expect_equal(project.subsphere(t(rotationmat_amaral(basis[1, ], nvec) %*% t(lowX_5)),
                                 nvec)[, 1:2], 0*lowX_5[, 1:2])

  #if I rotate to  eps * nvec + basis[1,], I expect the projection to be close to the original
  newdirection <- basis[1, ] + 1E-2 * nvec
  newdirection <- newdirection / sqrt(sum(newdirection^2))
  expect_equal(project.subsphere(
    t(rotationmat_amaral(basis[1, ], newdirection) %*% t(lowX_5)), nvec),
    lowX_5, tolerance = 1E-3)
})


test_that("divL2 reversed by div by rowSums",{
  X <- matrix(runif(5*20), nrow = 20)
  X <- X / rowSums(X)

  Xsph <- divL2(X)
  expect_equal(Xsph / rowSums(Xsph), X)
})
