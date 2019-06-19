context("rnaseq functions")
library(ribor)

experiments <- c("Hela_1", "Hela_2", "WT_1")
ribo.object <- ribo("sample.ribo")

rnaseq_tidy <- get_rnaseq(ribo.object,
                          experiments = experiments)

rnaseq <- get_rnaseq(ribo.object,
                     tidy = FALSE,
                     experiments = experiments)

actual <- c(nrow(rnaseq_tidy), ncol(rnaseq_tidy))
expected <- c(45, 4)

test_that("get_rnaseq- tidy version size",
          expect_equal(actual, expected))

hela_1_rnaseq <- get_rnaseq(ribo.object,
                            tidy = FALSE,
                            experiments = "Hela_1")

actual   <- rnaseq[rnaseq$experiment == "Hela_1", ]
expected <- hela_1_rnaseq

test_that("get_rnaseq- subsetting experiments",
          expect_equal(actual, expected))

actual   <- rnaseq_tidy[rnaseq_tidy$experiment == "Hela_1" &
                        rnaseq_tidy$region == "UTR5", 4]
actual   <- unname(actual)
actual   <- unlist(actual)

expected <- c(2.04, 3.2, 0)
error <- 1e-2

test_that("get_rnaseq- correct values",
          expect_equal(all(abs(expected - actual) <= error), TRUE))

actual   <- rnaseq[rnaseq$experiment == "WT_1" &
                   rnaseq$transcript == "GAPDH", c(-1,-2)]
actual   <- unname(actual)
actual   <- unlist(actual)

expected <- c(2.04, 5.3, 2.89, 1.78, 3.0)

result <- all(abs(expected - actual) <= error)

test_that("get_rnaseq- correct values",
          expect_equal(all(abs(expected - actual) <= error), expected = TRUE))

actual <- c(nrow(rnaseq), ncol(rnaseq))
expected <- c(9, 7)

test_that("get_rnaseq- non-tidy version size",
          expect_equal(actual, expected))

