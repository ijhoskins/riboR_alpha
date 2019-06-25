context("coverage functions")
library(ribor)

ribo.object <- ribo("sample.ribo")

cov_1   <- get_coverage(ribo.object,
                        name = "MYC",
                        range.lower = 2,
                        range.upper = 5,
                        experiments = "Hela_1")

cov_2   <- get_coverage(ribo.object,
                        name = "VEGFA",
                        range.lower = 2,
                        range.upper = 5,
                        experiments = "Hela_1")

cov_3   <- get_coverage(ribo.object,
                        name = "GAPDH",
                        range.lower = 2,
                        range.upper = 5,
                        experiments = c("Hela_1"))

actual <- sum(colSums(cov_1[, -1])) + 
          sum(colSums(cov_2[, -1])) + 
          sum(colSums(cov_3[, -1]))

expected <- 118

test_that("get_coverage- total reads of an experiment",
          expect_equal(actual, expected))

cov_4   <- get_coverage(ribo.object,
                        name = "GAPDH",
                        range.lower = 2,
                        range.upper = 2,
                        experiments = c("Hela_1"))

actual   <- rowSums(cov_4[, -1])
expected <- 14 

test_that("get_coverage- test individual read length",
          expect_equal(actual, expected))

cov_5   <- get_coverage(ribo.object,
                        name = "MYC",
                        range.lower = 3,
                        range.upper = 3,
                        experiments = c("Hela_2"))

actual <- rowSums(cov_5[, -1])
expected <- 3 

test_that("get_coverage- test individual read length",
          expect_equal(actual, expected))
