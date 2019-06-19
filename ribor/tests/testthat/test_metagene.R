context("metagene functions")
library(ribor)

ribo.object   <- ribo("sample.ribo")

meta_1 <- get_metagene(ribo.object,
                       "start",
                       range.lower = 2,
                       range.upper = 5,
                       length = TRUE,
                       transcript = TRUE,
                       experiments = c("Hela_1", "Hela_2"))

actual   <- c(nrow(meta_1), ncol(meta_1)) 
expected <- c(2, 6)

test_that("get_metagene size",
          expect_equal(actual, expected))



expected <- 17
test_that("get_metagene, checks the default case",
          expect_equal(meta_1[meta_1$experiment == "Hela_1", ][["-1"]], 
                       expected))


meta_2 <- get_metagene(ribo.object,
                       "start",
                       range.lower = 2,
                       range.upper = 3,
                       length = TRUE,
                       transcript = FALSE,
                       experiments = c("Hela_1"))

actual   <- c(nrow(meta_2), ncol(meta_2)) 
expected <- c(3, 7)

test_that("get_metagene- size",
          expect_equal(actual, expected))

expected <- 4
actual   <- meta_2[meta_2$transcript == "GAPDH", ][["-1"]]
test_that("get_metagene, only condense lengths",
          expect_equal(actual, 
                       expected))



meta_3 <- get_metagene(ribo.object,
                       "stop",
                       range.lower = 3,
                       range.upper = 4,
                       length = FALSE,
                       transcript = TRUE)

expected <- 4
test_that("get_metagene, only condense transcripts",
          expect_equal(meta_3[meta_3$experiment == "Hela_1" &
                              meta_3$length == 3, ][["-2"]], 
                       expected))



meta_4 <- get_metagene(ribo.object,
                       "stop",
                       range.lower = 4,
                       range.upper = 4,
                       length = FALSE,
                       transcript = FALSE)

actual <- meta_4[meta_4$experiment == "Hela_1" &
                   meta_4$transcript == "GAPDH", ][["1"]]
expected <- 1
test_that("get_metagene, preserve everything",
          expect_equal(actual, 
                       expected))


actual   <- c(nrow(meta_4), ncol(meta_4)) 
expected <- c(15, 8)

test_that("get_metagene- size",
          expect_equal(actual, expected))

tidy_meta_1 <- get_tidy_metagene(ribo.object,
                                 "start",
                                 2,
                                 3,
                                 length = TRUE)
actual   <- c(nrow(tidy_meta_1), ncol(tidy_meta_1)) 
expected <- c(25, 3)

test_that("get_metagene- size",
          expect_equal(actual, expected))


actual <- tidy_meta_1[tidy_meta_1$experiment == "Hela_1" &
                        tidy_meta_1$position == -1, ][["count"]]
expected <- 9

test_that("get_metagene, preserve everything",
          expect_equal(actual, 
                       expected))



tidy_meta_2 <- get_tidy_metagene(ribo.object,
                                 "start",
                                 2,
                                 3,
                                 length = FALSE)
actual <- c(nrow(tidy_meta_2), ncol(tidy_meta_2))
expected <- c(50, 4)

test_that("get_metagene- size",
          expect_equal(actual, 
                       expected))




  
actual <- tidy_meta_2[tidy_meta_2$experiment == "Hela_2" &
                        tidy_meta_2$length == 3 &
                        tidy_meta_2$position == -1, ][["count"]]

expected <- 5
test_that("get_metagene, preserve everything",
          expect_equal(actual, 
                       expected))

