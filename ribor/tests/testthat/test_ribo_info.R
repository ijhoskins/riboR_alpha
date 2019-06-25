context("file info functions")
library(ribor)

ribo.object   <- ribo("sample.ribo")
expected.exp <- c("Hela_1", "Hela_2", "WT_1", "WT_2", "WT_3")

#test get_experiments
test_that("get_experiments returns the correct output",
          expect_equal(get_experiments(ribo.object = ribo.object), expected.exp))

actual <- get_info(ribo.object)

names.list     <- expected.exp
reads.list     <- list(118, 118, 12, 12, 12)
coverage.list  <- list(TRUE, TRUE, TRUE, TRUE, FALSE)
rna.seq.list   <- list(TRUE, TRUE, TRUE, FALSE, FALSE)
metadata.list  <- list(TRUE, FALSE, TRUE, TRUE, TRUE)
expected.table <- data.table(names    = names.list,
                             reads    = reads.list,
                             coverage = coverage.list,
                             rna.seq  = rna.seq.list,
                             metadata = metadata.list)

properties     <- list(1, 3, 5, 2, 2, "appris_human", 2)
unnamed.actual <- actual[1:7]
unname(unnamed.actual)

#test get_info
test_that("get_info- returns a list of the correct length",
           expect_equal(length(actual), 8))

test_that("get_info- checks first seven elements",
           expect_equal(as.character(unnamed.actual),
                        as.character(properties)))


#test get_rnaseq
abundance.list  <- c(2.89, 15.46, 8.0, 6.42, 1.37, 0.06, 2.89, 15.46, 8.0)

experiment.list <- c("Hela_1", "Hela_1", "Hela_1",
                     "Hela_2", "Hela_2", "Hela_2",
                     "WT_1"  , "WT_1"  , "WT_1")

transcript.list <- c("GAPDH", "VEGFA", "MYC",
                     "GAPDH", "VEGFA", "MYC",
                     "GAPDH", "VEGFA", "MYC")

expected.rnaseq <- data.table(transcript = transcript.list,
                              experiment = experiment.list,
                              abundance  = abundance.list)

hela.only <- expected.rnaseq[expected.rnaseq$experiment == "Hela_1" |
                             expected.rnaseq$experiment == "Hela_2",]

hela.list     <- c("Hela_1", "Hela_2")
hela.no.match <- c("Hela_3", "Hela_4")


#test get_metadata
hela.meta <- c("HeLa", "5 min", "HindIII", "https://www.encodeproject.org/")

wt.meta <- c("WT", "5 min", "HindIII", "https://www.encodeproject.org/")

actual.hela.meta <- get_metadata(ribo.object, "Hela_1")
actual.hela.meta <- unname(actual.hela.meta)

test_that("get_metadata- checks correct input for Hela_1",
          expect_equal(as.character(actual.hela.meta), as.character(hela.meta)))

actual.wt.meta <- get_metadata(ribo.object, "WT_1")
actual.wt.meta <- unname(actual.wt.meta)

test_that("get_metadata- checks correct input for WT_1",
          expect_equal(as.character(actual.wt.meta), as.character(wt.meta)))

test_that("get_metadata- checks output for file with no metadata",
          expect_warning(get_metadata(ribo.object, "Hela_2")))