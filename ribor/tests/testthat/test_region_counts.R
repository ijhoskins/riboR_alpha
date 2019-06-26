context("region count functions")
library(ribor)

ribo.object<- ribo("sample.ribo")

rc_1 <- get_region_counts(ribo.object,
                          regions = c("UtR5", "cDs", "utr3"),
                          2,
                          5,
                          experiments = "Hela_1")

actual <- c(nrow(rc_1), ncol(rc_1))
expected <- c(3, 3)

test_that("get_region_counts- size",
          expect_equal(actual, expected))


rc_all <- get_region_counts(ribo.object,
                          regions = c("UtR5", "UTR5j", "cDs", "utr3j", "utr3"),
                          2,
                          5,
                          experiments = "Hela_1")

actual <- c(nrow(rc_all), ncol(rc_all))
expected <- c(5, 3)

test_that("get_region_counts- size",
          expect_equal(actual, expected))

actual   <- sum(rc_all[, 3])
expected <- 118 

test_that("get_region_counts- total count",
           expect_equal(actual, expected))

rc_3 <- get_region_counts(ribo.object,
                          regions = c("CDS"),
                          2,
                          5,
                          length = TRUE,
                          transcript = FALSE,
                          experiments = c("Hela_1"))

actual <- c(nrow(rc_3), ncol(rc_3))
expected <- c(3, 4)

test_that("get_region_counts- size",
          expect_equal(actual, expected))

actual <- rc_3[1, ][["count"]]
expected <- 13

test_that("get_region-count- preserving transcripts",
          expect_equal(actual, expected))

rc_4 <- get_region_counts(ribo.object,
                          regions = c("UtR5j"),
                          2,
                          5,
                          length = FALSE,
                          transcript = TRUE)

actual <- c(nrow(rc_4), ncol(rc_4))
expected <- c(20, 4)

test_that("get_region-count- size",
          expect_equal(actual, expected))

actual <- sum(rc_4$count)
expected <- 102
test_that("get_region_count- preserving lengths",
          expect_equal(actual, expected))
