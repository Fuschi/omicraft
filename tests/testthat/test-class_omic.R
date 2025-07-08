# test/test-class_omic.R

test_that("valid omic object passes validation", {
  abun <- matrix(1:9, nrow = 3, dimnames = list(paste0("S", 1:3), paste0("T", 1:3)))
  rela <- abun / rowSums(abun)
  norm <- log2(abun + 1)
  meta <- data.frame(group = c("A", "B", "A"), row.names = paste0("S", 1:3))
  taxa <- data.frame(class = c("c1", "c2", "c3"), row.names = paste0("T", 1:3))
  netw <- igraph::make_empty_graph(n = 3)
  igraph::V(netw)$name <- paste0("T", 1:3)
  comm <- igraph::cluster_optimal(netw)
  
  omic <- new("omic", abun = abun, rela = rela, norm = norm,
              meta = meta, taxa = taxa, netw = netw, comm = comm)
  expect_true(validObject(omic))
})

test_that("duplicate row names in abun triggers validation error", {
  abun <- matrix(1:4, nrow = 2, dimnames = list(c("S1", "S1"), c("T1", "T2")))
  expect_error(new("omic", abun = abun), "Duplicate row names")
})

test_that("abun must be a numeric matrix", {
  abun <- matrix("a", nrow = 2, ncol = 2,
                 dimnames = list(c("s1", "s2"), c("f1", "f2")))
  meta <- data.frame(group = c("A", "B"),
                     row.names = c("s1", "s2"))
  
  expect_error(
    new("omic", abun = abun, meta = meta),
    regexp = "must be a numeric matrix"
  )
})

test_that("netw and comm with mismatched vertices triggers validation error", {
  abun <- matrix(1:4, nrow = 2, dimnames = list(c("S1", "S2"), c("T1", "T2")))
  netw <- igraph::make_empty_graph(n = 2)
  igraph::V(netw)$name <- c("T1", "T2")
  comm <- structure(list(membership = 1:3), class = "communities")
  
  expect_error(
    new("omic", abun = abun, netw = netw, comm = comm),
    "number of community memberships.*vertices"
  )
})

test_that("abun cannot contain negative values", {
  abun <- matrix(c(1, -2, 3, 4), nrow = 2,
                 dimnames = list(c("s1", "s2"), c("f1", "f2")))
  rela <- matrix(c(.1, -.2, .3, .4), nrow = 2,
                 dimnames = list(c("s1", "s2"), c("f1", "f2")))
  meta <- data.frame(group = c("A", "B"),
                     row.names = c("s1", "s2"))
  
  expect_error(
    new("omic", abun = abun, meta = meta),
    regexp = "must be >= 0"
  )
  expect_error(
    new("omic", rela = rela, meta = meta),
    regexp = "must be >= 0"
  )
})

test_that("meta must be a data.frame", {
  abun <- matrix(1:4, nrow = 2,
                 dimnames = list(c("s1", "s2"), c("f1", "f2")))
  meta <- c("A", "B")
  
  expect_error(
    new("omic", abun = abun, meta = meta),
    regexp = "meta slot must be a <data.frame>"
  )
})

test_that("meta must not contain reserved keywords", {
  abun <- matrix(1:4, nrow = 2,
                 dimnames = list(c("s1", "s2"), c("f1", "f2")))
  meta <- data.frame(sample_id = c("A", "B"),
                     row.names = c("s1", "s2"))
  
  expect_error(
    new("omic", abun = abun, meta = meta),
    regexp = "must not match reserved keyword"
  )
})
