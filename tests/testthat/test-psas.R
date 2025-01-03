test_that("psas() on 5-component vertices", {
  # the two components with the smallest number of 1s should merge first and so forth
  # 1,2,3,4,5
  X <- matrix(0, ncol = 5, nrow = 15)
  X[1, 1] <- 1
  X[2:3, 2] <- 1
  X[4:6, 3] <- 1
  X[7:10, 4] <- 1
  X[11:15, 5] <- 1
  colSums(X)
  # since 1s become 1s after a merge. I expect psas to get pairs
  # 1,2 (weight of 1/3 ish)
  # 1,2 (weight of 1/2 ish)
  # 2,3 (weight of 4/9 ish)
  # 1,2 (weight of 9/15 ish)
  out <- suppressMessages(psas(X, testweights = seq(0, 1, length.out = 100)))
  expect_true(out$merges[2, ]["w"] %in% seq(0, 1, length.out = 100))
  expect_equal(out$merges[2, c("v1", "v2")] |> unlist(), c(2, 1), ignore_attr = TRUE)
  expect_equal(out$merges[2, "w"], 2/3, tolerance = 1E-1, ignore_attr = TRUE) #I thought new vertex would be closer to v2 because 2x more points at v2 - AND IT IS! weight on v1 is half that of v2
  expect_equal(out$merges[3, c("v1", "v2")] |> unlist(), c(2, 1), ignore_attr = TRUE)
  expect_equal(out$merges[3, "w"], 1/2, tolerance = 1E-1, ignore_attr = TRUE)
  expect_equal(out$merges[4, c("v1", "v2")] |> unlist(), c(3, 2), ignore_attr = TRUE)
  expect_equal(out$merges[4, "w"], 1 - 4/9, tolerance = 1E-1, ignore_attr = TRUE)
  expect_equal(out$merges[5, c("v1", "v2")] |> unlist(), c(1, 2), ignore_attr = TRUE)
  expect_equal(out$merges[5, "w"], 9/15, tolerance = 1E-1, ignore_attr = TRUE)
})

test_that("psas() can be applied using provided merges", {
  # 1,2,3,4,5
  X <- matrix(0, ncol = 5, nrow = 15)
  X[1, 1] <- 1
  X[2:3, 2] <- 1
  X[4:6, 3] <- 1
  X[7:10, 4] <- 1
  X[11:15, 5] <- 1
  merges <- data.frame(v1 = c(NA_integer_, rep(1, 4)),
                       v2 = c(NA_integer_, rep(2, 4)),
                       w = c(NA_real_, rep(0.5, 4)))
  out <- suppressMessages(psas(X, merges = merges))
  expect_equal(out$merges, merges, ignore_attr = TRUE)

  expect_true(all(round(abs(do.call(cbind, out$scores[-1])), 10) %in% round(c(sqrt(2)/2, 0), 10))) #data points are vertices and the merges occur at halfway - along the sqrt(2) distance between vertices
})
