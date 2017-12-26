library(megadaph.mtdna)
context("Merge test table indel rows")

test.tab <- readRDS("../dat/test_table.Rds")

find.test.tab <- data.frame(
  species = rep("pulex", 3),
  genotype = rep("L", 3),
  isolate = rep("L", 3),
  sample = rep(1, 3),
  pos = c(1:3),
  coverage = rep(2000, 3),
  ref = rep("A", 3),
  alt = c("A", "-", "-"),
  class = c("snv", "deletion", "deletion"),
  af = rep(1, 3),
  af.diff = c(0, 1, 1),
  strand.bias = rep(0, 3),
  unique = rep(TRUE, 3),
  coverage.prop = rep(1, 3),
  p = c(1, 0, 0)
)

# Test merge found indel

# Test merges indels simple case

# Test ignores non-significant

# 
