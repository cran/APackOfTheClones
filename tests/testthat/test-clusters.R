source("testdata/cluster_lists.R")
source("testdata/SCIP.R")

test_that("find_centroids() works", {
  input <- data.frame(
    "UMAP_1" = c(-1, 1, 1, 2, 3, 4, 4, 5, 9, 9),
    "UMAP_2" = c(-1, 3, 5, 4, 2, 7, 9, 0, 1, 10),
    "clstrs" = c(0, 9, 5, 5, 4, 0, 0, 9, 4, 4)
    )

  expected_output <- list(
    `0` = c(2.333333, 5),
    `4` = c(7, 4.333333),
    `5` = c(1.5, 4.5),
    `9` = c(3, 1.5)
  )

  expect_equal(find_centroids(input),
               expected_output,
               tolerance=3e-7)
})

test_that("get_cluster_centroids() works", {
  expect_equal(
    get_cluster_centroids(test_pbmc),
    list(`1` = c(5.72926706448197, 8.48345941543579),
         `2` = c(-9.54877844080329, -14.1390990257263)
         ))
})

test_that("trans_coord() works withOUT new_coord", {
  c1_new <- c1
  c1_new$centroid <- c(4, 5)
  expect_equal(trans_coord(c1_new), c1_shifted_by_4_5)
})

test_that("trans_coord() works WITH new_coord",{
  expect_equal(trans_coord(c1, c(9, 0)), c1_shifted_to_9_0)
})

# need to test if it doesn't work when the input is null/na!

test_that("move_cluster() works", {
  expect_equal(move_cluster(c1,c(4,5)), c1_shifted_by_4_5)
  expect_equal(move_cluster(c1,c(9,0)), c1_shifted_to_9_0)
})