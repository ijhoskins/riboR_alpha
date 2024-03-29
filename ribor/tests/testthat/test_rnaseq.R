context("rnaseq functions")
library(ribor)

experiments <- c("Hela_1", "Hela_2", "WT_1")
ribo.object <- ribo("sample.ribo")

rnaseq_tidy <- get_rnaseq(ribo.object,
                          experiments = experiments)

rnaseq <- get_rnaseq(ribo.object,
                     tidy = FALSE,
                     experiments = experiments)

actual <- get_rnaseq(ribo.object,
                     regions = c("UTR3j", "cds"),
                     tidy = FALSE,
                     experiments = experiments)

expected <-rnaseq[, -c("UTR5", "UTR5J", "UTR3")]

test_that("get_rnaseq = non-tidy CDS and UTR3J",
          expect_true(all(actual == expected)))

actual <- get_rnaseq(ribo.object,
                     regions = c("UTR3", "cds", "utR5"),
                     tidy = FALSE,
                     experiments = experiments)

expected <-rnaseq[, -c("UTR3J", "UTR5J")]

test_that("get_rnaseq = non-tidy CDS and UTR3J",
          expect_true(all(actual == expected)))

actual   <- get_rnaseq(ribo.object,
                         tidy = TRUE,
                         regions = c("UTR5"),
                         experiments = experiments)

expected <- rnaseq_tidy[rnaseq_tidy$region == "UTR5", ]


test_that("get_rnaseq- tidy UTR5 only",
          expect_true(all(actual == expected)))

actual   <- get_rnaseq(ribo.object,
                       tidy = TRUE,
                       regions = c("CDS", "utR5j"),
                       experiments = experiments)

expected <- rnaseq_tidy[rnaseq_tidy$region == "UTR5J" | 
                        rnaseq_tidy$region == "CDS"]

test_that("get_rnaseq- tidy UTR5 and CDS",
          expect_true(all(actual == expected)))

actual   <- c(nrow(rnaseq_tidy), ncol(rnaseq_tidy))
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

