
#' @export
af_by_sample <- function(var_table, line_info) {
  sapply(line_info$sample, function(smp) {
    var_table[which(var_table$sample == smp),"af"]})
}

#' @export
het_vs_hom <- function(var_table, line_info) {
  var_by_spp <- split(var_table, var_table$species)
  li_by_spp <- split(line_info, line_info$species)

  by_spp <- lapply(1:(length(li_by_spp)), function(i) {
    hom <- dplyr::filter(var_by_spp[[i]], af > 0.99)
    het <- dplyr::filter(var_by_spp[[i]], af <= 0.99)
    prop_hom <- nrow(hom)/nrow(het)
    hom_af <- lapply(li_by_spp[[i]]$sample, function(x) {
      hom[which(hom$sample == x), "af"]
})
    het_af <- lapply(li_by_spp[[i]]$sample, function(x) {
      het[which(het$sample == x), "af"]
})
    hom_mutation_rate <- quantile_mutation_rate(hom_af, li_by_spp[[i]]$generations,
      sum(li_by_spp[[i]]$bp))
    het_mutation_rate <- quantile_mutation_rate(het_af, li_by_spp[[i]]$generations,
      sum(li_by_spp[[i]]$bp))
    list(prop_hom=prop_hom, hom_u = hom_mutation_rate, het_u = het_mutation_rate)
    })
  fill <- c("Homoplasmic", "Heteroplasmic", "Homoplasmic", "Heteroplasmic")
  dat <- as.data.frame(rbind(by_spp[[1]]$hom_u, by_spp[[1]]$het_u,
      by_spp[[2]]$hom_u, by_spp[[2]]$het_u))
  group <- c("D. magna", "D. magna", "D. pulex", "D. pulex")
  dat <- cbind(group, dat)
  colnames(dat) <- c("group", "V1", "ci1", "q25", "q50", "q75", "ci2")
  custom_boxplot(dat, "Species", "Mutation Rate", expression(italic("D. magna"), italic("D. pulex")),
    fill=fill)
}
#' @export
#' @importFrom Hmisc capitalize
mutation_rate_by <- function(var_table, line_info, by, quant=TRUE) {
  groups <- unique(line_info[, by])
  mu_table <- t(sapply(groups, function(g) {
      var_table_by <- var_table[which(var_table[, by]==g),]
      line_info_by <- line_info[which(line_info[, by]==g),]
      af <- af_by_sample(var_table_by, line_info_by)
      gen <- line_info_by$generations
      nsamples <- length(gen)
      nuc <- sum(line_info_by$bp)
      mu <- mutation_rate(unlist(af), mean(gen), nuc)
      if (quant) {
        quantile_mutation_rate(af, gen, nuc)
      } else {
        mu
      }
      }))
  if (quant) {
    mu_table <- cbind(as.data.frame(capitalize(groups),
        stringsAsFactors=FALSE),
      mu_table)
    colnames(mu_table) <- c("group", "mutation_rate", "ci1", "q25",
      "q50", "q75", "ci2")
    mu_table
  } else {
    mu_table
  }
}

# CANT DO PERMUTATION TESTS WITH UNEQUAL SAMPLE SIZES

#' @export
sample_mu <- function(line_info, var_table) {
  sample_mutation_rates <- sapply(line_info$sample, function(sample) {
    sample_afs <- var_table[which(var_table$sample == sample), "af"]
    sample_gen <- line_info[which(line_info$sample == sample),
      "generations"]
    sample_nuc <- line_info[which(line_info$sample == sample),
      "bp"]
    mutation_rate(sample_afs, sample_gen, sample_nuc)
      })
  barplot(sample_mutation_rates, col=as.factor(line_info$genotype))
  sample_mutation_rates
}

#' @export
#' @importFrom dplyr filter
confounding_factors <- function(line_info, sample_mutation_rates, merged_table) {
  leveneTest(sample_mutation_rates~as.factor(line_info$genotype))
  leveneTest(sample_mutation_rates~as.factor(line_info$isolate))

  pulex_table <- filter(merged_table, species=="pulex")
  magna_table <- filter(merged_table, species=="magna")
  plot(pulex_table$pos~log(pulex_table$af), xlim=c(-8, -2))
  plot(magna_table$pos~log(magna_table$af), xlim=c(-8, -2))
  plot(density(pulex_table$af, n=10000, from=0, to=1), xlim=c(0, 0.02), ylim=c(0, 50))
  lines(density(magna_table$af, n=10000, from=0, to=1), col="red")
  plot(density(pulex_table$af, n=10000, from=0, to=1), xlim=c(0.01, 0.05), ylim=c(0, 2))
  lines(density(magna_table$af, n=10000, from=0, to=1), col="red")
  nmat <- sapply(sort(unique(merged_table$af)), function(n) {
    c(length(which(pulex_table$af==n)), length(which(magna_table$af==n)))
      })
  nmat <- apply(nmat, 1, function(r) r/sum(r))
  colnames(nmat) <- ""
  plot(nmat[,1]-nmat[,2], type="l", ylim=c(-0.2, 0.2))
  plot(log(nmat[1:16,1])-log(nmat[1:16,2]), type="l", xlim=c(0, 16))

  arg_mat <- readRDS("../consensus_calling/arg_mat.rds")
  nsamples <- lapply(arg_mat[1,], length)
  ntables <- c(rep(14948, 9), rep(15333, 2))


  ctabs <- lapply(1:ncol(arg_mat), function(i) {
    simulate_contingency_tables(ntables[i], nsamples[[i]], unlist(arg_mat[3,i]),
      c(1000, 1000))
      })
  pulex_ctabs <- unlist(ctabs[10:11], recursive=FALSE)
  magna_ctabs <- unlist(ctabs[1:9], recursive=FALSE)
  mut_wt <- readRDS("../tables/mut_wt_matrices.sub.Rds")
  pulex_mut_wt <- mut_wt[(length(mut_wt)-(15333*2)+1):length(mut_wt)]
  magna_mut_wt <- mut_wt[1:(length(mut_wt)-(15333*2))]
  equals1000 <- function(ctab) {
    all(rowSums(ctab)==1000)
  }
  pulex_inc_mat <- which(sapply(pulex_mut_wt, equals1000))
  magna_inc_mat <- which(sapply(magna_mut_wt, equals1000))
  pulex_mut_wt <- pulex_mut_wt[pulex_inc_mat]
  pulex_ctabs <- pulex_ctabs[pulex_inc_mat]
  magna_mut_wt <- magna_mut_wt[magna_inc_mat]
  magna_ctabs <- magna_ctabs[magna_inc_mat]
  pulex_afs <- unlist(lapply(pulex_mut_wt, function(x) x[,1]/(x[,1]+ x[,2])))
  magna_afs <- unlist(lapply(magna_mut_wt, function(x) x[,1]/(x[,1]+ x[,2])))
  sim_pulex_afs <- unlist(lapply(pulex_ctabs, function(x) x[,1]/(x[,1]+ x[,2])))
  sim_magna_afs <- unlist(lapply(magna_ctabs, function(x) x[,1]/(x[,1]+ x[,2])))
  plot(density(pulex_afs-mean(sim_pulex_afs)), xlim=c(0,1), ylim=c(0, 0.05))
  lines(density(magna_afs-mean(sim_magna_afs)), col="red")
  nmat <- sapply(sort(unique(c(magna_afs, pulex_afs))), function(n) {
    c(length(which(pulex_afs==n)), length(which(magna_afs==n)))
      })
  nmat <- t(apply(nmat, 1, function(r) r/sum(r)))
  plot(nmat[1,], xlim=c(0, 40), ylim=c(0, 0.001))
  lines(nmat[1,])
  points(nmat[2,], col="red")
  lines(nmat[2,], col="red")

  line_info <- cbind(line_info, sample_mutation_rates)
  plot(density(line_info$sample_mutation_rates))
  shapiro.test(sample_mutation_rates)
  cor.test(line_info$generations,
    line_info$sample_mutation_rates,
    method="spearman")
  ## Non-significant after removal of pulex
  cor.test(line_info$generations[which(line_info$species=="magna")],
    line_info$sample_mutation_rates[which(line_info$species=="magna")],
    method="spearman")
  non_zero <- line_info[which(line_info$sample_mutation_rates != 0),]
  plot(log(non_zero$generations), log(non_zero$sample_mutation_rates),
    xlim=c(-1, 6), ylim=c(-20, -10))
  abline(lm(log(non_zero$sample_mutation_rates)~log(non_zero$generations)), col="red")
  tab <- cbind(as.data.frame(log(sample_mutation_rates)),
    log(line_info$generations),
    line_info$genotype)
  tab <- tab[which(!is.infinite(tab[,1])),]
  require(car)
  scatterplot(tab[,1]~tab[,2] | tab[,3], smoother=FALSE,
    xlab="log(generations)", ylab="log(Mutation Rate)")
  cor.test(line_info$sequencing_error,
    line_info$sample_mutation_rates,
    method="spearman")
  cor.test(line_info$coverage,
    line_info$sample_mutation_rates,
    method="spearman")
}

#' @export
by_isolate <- function(var_table, line_info) {
  # BY POPULATION
  isolate_mutation_rates <- mutation_rate_by(var_table, line_info, "isolate")
  custom_boxplot(isolate_mutation_rates, "Isolate", "Mutation Rate",
    unique(line_info$isolate))

  isolate_diffs <- boot_mutation_rate_diff_within(var_table, line_info, "isolate", "genotype")
}

#' @export
by_population <- function(line_info, var_table) {
  population_mutation_rates <- mutation_rate_by(var_table, line_info, "genotype")
  custom_boxplot(population_mutation_rates, "Genotype", "Mutation Rate",
    unique(line_info$genotype))
  population_diffs <- boot_mutation_rate_diff_within(var_table, line_info, "genotype",
    "species")
}

#' @export
by_species <- function(var_table, line_info) {
  species_mutation_rates <- mutation_rate_by(var_table, line_info, "species")
  species_diffs <- boot_mutation_rate_diff_within(var_table, line_info, "species")
  custom_boxplot(species_mutation_rates, "Species", "Mutation Rate",
    c("D. magna", "D. pulex"))
}

#' @export
by_reproduction <- function(var_table, line_info) {
  repro_mutation_rates <- mutation_rate_by(var_table, line_info, "reproduction")
  custom_boxplot(repro_mutation_rates, "Reproduction", "Mutation Rate",
    c("Asexual", "Cyclical Parthenogenetic"))
  repro_diffs <- boot_mutation_rate_diff_within(var_table, line_info, "reproduction")
  repro_diffs
}

#' @export
diff_p_value <- function(boot_diff) {
  boot_diff <- unlist(boot_diff)
  avg <- mean(boot_diff)
  if (avg > 0) {
    p <- (length(which(unlist(boot_diff) <= 0)) + 1) / (length(boot_diff)+ 1)
  } else {
    p <- (length(which(unlist(boot_diff) >= 0)) + 1) / (length(boot_diff)+ 1)
  }
  p
}

ci <- function(x) {
  error <- qnorm(0.975)*sd(x)/sqrt(length(x))
  ci1 <- mean(x)-error
  ci2 <- mean(x)+error
  c(ci1, ci2)
}

#' @export
by_indel_snv <- function(line_info, var_table) {
  indel_var_table <- var_table[which(var_table$class %in%
    c("insertion", "deletion")),]
  snv_var_table <- var_table[which(var_table$class == "snv"),]
  snv_mu_by_species <- mutation_rate_by(snv_var_table, line_info, "species")
  indel_mu_by_species <- mutation_rate_by(indel_var_table, line_info, "species")
  mut_class_df <- rbind(snv_mu_by_species, indel_mu_by_species)
  fill <- c("SNV", "SNV", "Indel", "Indel")
  snv_boot_diff <- boot_mutation_rate_diff_within(snv_var_table, line_info, "species")
  snv_p <- diff_p_value(snv_boot_diff)
  indel_boot_diff <- boot_mutation_rate_diff_within(indel_var_table, line_info, "species")
  indel_p <- diff_p_value(indel_boot_diff)
  custom_boxplot(mut_class_df, "Mutation Class", "Mutation Rate",
    c("D. magna", "D. pulex"), fill=fill)
  ninsertions <- length(which(indel_var_table$class == "insertion"))
  ndeletions <- length(which(indel_var_table$class == "deletion"))
  fisher.test(matrix(c(ndeletions, ninsertions, ninsertions, ndeletions), ncol=2))
  ins_lengths <- nchar(filter(indel_var_table, class=="insertion")$alt) - 1
  del_lengths <- nchar(filter(indel_var_table, class=="deletion")$ref) - 1
  wilcox.test(ins_lengths, del_lengths)
  indel_afs <- var_table[which(var_table$class
    %in% c("insertion", "deletion")),"af_diff"]
  snv_afs <- var_table[which(var_table$class=="snv"),"af_diff"]
  wilcox.test(indel_afs, snv_afs)
}

#' @export
by_effect <- function(var_table, by) {

  split_by <- split(var_table, var_table[, by])

  ts_tv_by <- t(sapply(split_by, function(t) {
      nts <- length(which(t$ts_tv=="transition"))
      ntv <- length(which(t$ts_tv=="transversion"))
      c(nts, ntv)
      }))

  nts <- length(which(var_table$ts_tv=="transition"))
  ntv <- length(which(var_table$ts_tv=="transversions"))
  prop.test(nts, nts+ntv, p=1/3)

  a_t <- length(which((var_table$ref %in% c("A", "T")) &
      (var_table$alt %in% c("A", "T"))))
  g_t <- length(which((var_table$ref %in% c("G", "T")) &
      (var_table$alt %in% c("G", "T"))))
  g_c <- length(which((var_table$ref %in% c("G", "C")) &
      (var_table$alt %in% c("G", "C"))))
  a_c <- length(which((var_table$ref %in% c("A", "C")) &
      (var_table$alt %in% c("A", "C"))))
  a_g <- length(which((var_table$ref %in% c("A", "G")) &
      (var_table$alt %in% c("A", "G"))))
  c_t <- length(which((var_table$ref %in% c("C", "T")) &
      (var_table$alt %in% c("C", "T"))))

  snv_df <- data.frame(id=c("A->T", "G->T", "G->C", "A->C", "A->G",
      "C->T"), n=c(a_t, g_t, g_c, a_c, a_g, c_t),
    stringsAsFactors = FALSE)
  snv_df <- snv_df[order(snv_df$n),]
}

#' @export
#' @importFrom coin independence_test
by_mutation_effect <- function(var_table) {
  independence_test(var_table$af~as.factor(var_table$effect))
  independence_test(var_table$af~as.factor(var_table$gene))
  synon <- as.factor(var_table$effect=="synonymous_variant")
  independence_test(var_table$af~synon)
  independence_test(var_table$p_value~synon)
}





#' @export
locations <- function(var_table, n) {
  mag_mut_by_loc <- sapply(1:(n-1), function(i) {
    nrow(dplyr::filter(var_table, species=="magna", pos > (i*(14948/n)),
        pos <= ((i+1)*(14948/n))))
    })
  pul_mut_by_loc <- sapply(1:(n-1), function(i) {
    nrow(dplyr::filter(var_table, species=="pulex", pos > (i*(15333/n)),
        pos <= ((i+1)*(15333/n))))
    })
}


#' @export
read_tables <- function() {
  var_table <<- read.csv("../data/tables/annot_table.csv", stringsAsFactors = F)
  line_info <<- read.csv("../data/tables/line_info.csv", stringsAsFactors = F)
  merged_table <<- read.csv("../data/tables/merged_table.csv",
    stringsAsFactors = F)
}

#' @export
species_mu <- function(line_info, var_table) {
  species_mutation_rates <- mutation_rate_by(var_table, line_info, "species")
  box <- custom_boxplot(species_mutation_rates, "Species", "Mutation Rate",
    c("D. magna", "D. pulex"))
  box

}

#' @export
species_diffs <- function(line_info, var_table) {
  species_diffs <- boot_mutation_rate_diff_within(var_table, line_info, "species")
  indel_diffs <- boot_mutation_rate_diff_within(dplyr::filter(var_table, class != "snv"),
    line_info, "species")

  c(species_diffs, indel_diffs)
}
