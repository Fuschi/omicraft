test_that("omics constructor creates empty object", {
  obj <- omics()
  expect_s4_class(obj, "omics")
  expect_length(obj@omics, 0)
})

test_that("omics constructor accepts named list of omic objects", {
  dummy_abun <- matrix(1:4, nrow = 2, dimnames = list(c("s1", "s2"),
                                                      c("t1", "t2")))
  dummy_obj <- omic(abun = dummy_abun)
  
  obj <- omics(A = dummy_obj, B = dummy_obj)
  expect_s4_class(obj, "omics")
  expect_named(obj@omics, c("A", "B"))
  expect_true(all(sapply(obj@omics, function(x) inherits(x, "omic"))))
})

test_that("omics constructor accepts list of omic objects if named", {
  dummy_abun <- matrix(1:4, nrow = 2, dimnames = list(c("s1", "s2"),
                                                      c("t1", "t2")))
  dummy_obj <- omic(abun = dummy_abun)
  
  named_list <- list(A = dummy_obj, B = dummy_obj)
  obj <- omics(named_list)
  expect_s4_class(obj, "omics")
  expect_equal(length(obj@omics), 2)
  expect_named(obj@omics)
})

test_that("omics constructor fails on unnamed input", {
  dummy_abun <- matrix(1:4, nrow = 2, dimnames = list(c("s1", "s2"),
                                                      c("t1", "t2")))
  dummy_obj <- omic(abun = dummy_abun)
  
  expect_error(
    omics(dummy_obj, dummy_obj),
    regexp = "must be named"
  )
})

test_that("omics constructor fails on invalid content", {
  expect_error(
    omics(A = "not an omic"),
    regexp = "must be instances of"
  )
})

test_that("omics constructor fails on duplicated names", {
  dummy_abun <- matrix(1:4, nrow = 2, dimnames = list(c("s1", "s2"),
                                                      c("t1", "t2")))
  dummy_obj <- omic(abun = dummy_abun)
  
  duplicated_list <- list(A = dummy_obj, A = dummy_obj)
  expect_error(
    omics(duplicated_list),
    regexp = "must be unique"
  )
})

test_that("omics `$` accessor works", {
  dummy_abun <- matrix(1:4, nrow = 2, dimnames = list(c("s1", "s2"),
                                                      c("t1", "t2")))
  dummy_omic <- omic(dummy_abun)
  dummy_omics <- omics(test = dummy_omic)
  expect_equal(dummy_omics$test, dummy_omic)
})

test_that("omics `$<-` setter works and validates", {
  dummy_abun <- matrix(1:4, nrow = 2, dimnames = list(c("s1", "s2"),
                                                     c("t1", "t2")))
  dummy_omic <- omic(abun = dummy_abun)
  obj <- omics()
  obj$test <- dummy_omic
  expect_s4_class(obj$test, "omic")
  expect_identical(obj$test, dummy_omic)
})

test_that("omics `$<-` fails on invalid object", {
  obj <- omics()
  expect_error(
    obj$invalid <- "not an omic",
    regexp = "Assigned value must be"
  )
})
