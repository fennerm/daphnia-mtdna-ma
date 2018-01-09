#!/usr/bin/env Rscript

##### MAIN #####

#' @export
generate_params <- function() {
    fdr_table <- readRDS("../tables/fdr_table.rds")
    line_info <- read.csv("../tables/line_info.csv")
    files <- list.files("../consensus_calling/vcf", full.names=TRUE)
    isolates <- get_isolate(files)
    outfile <- sapply(isolates, function(x) paste0(x, "_mle.rds"))
    threads <- 1
    nsamples <- sapply(unique(line_info$isolate), function(x) {
        nrow(line_info[line_info$isolate==x,])
    })
    generations <- sapply(unique(line_info$isolate), function(x) {
        mean(line_info[line_info$isolate==x, "generations"])
    })
    bases <- sapply(unique(line_info$isolate), function(x) {
        mean(line_info[line_info$isolate==x, "bases_surveyed"])
    })
    coverage <- sapply(unique(line_info$isolate), function(x) {
        mean(line_info[line_info$isolate==x, "coverage"])
    })
    arg_mat <- cbind(as.data.frame(files), nsamples,
                     outfile, generations, coverage,
                     bases)
    arg_mat
}

#' @export
produce_mle_table <- function(models) {

    rows <- sapply(models, function(mod) {
        f0 <- sapply(mod, function(m) {
            m$estimate["f0"]
        })
        ve <- sapply(mod, function(m) {
            m$estimate["VE"]
        })
        loglik <- sapply(mod, function(m) {
            m$maximum
        })
        N_hat <- which(vr==max(vr))
        VE_hat <- ve[N_hat]
        f0_hat <- f0[N_hat]
    })
}

# if (interactive()) {
#     if(exists("models")) {
#         model <- models[[j]]
#         isolates <- get_isolate(as.character(arg_mat$files))
#         vr <<- vector(mode="numeric", length=length(model))
#         f0r <<- vector(mode="numeric", length=length(model))
#         ver <<- vector(mode="numeric", length=length(model))
#         for (i in 1:length(model)) {
#             f0r[i] <- model[[i]]$estimate[1]
#             ver[i] <- model[[i]]$estimate[2]
#             vr[i] <- model[[i]]$maximum
#         }
#
#         plot(unlist(f0r), type="l")
#         plot(unlist(ver), type="l")
#         plot(unlist(vr), type="l")
#         N_hat <- which(vr==max(vr))
#         VE_hat <- ver[N_hat]
#         sd_hat <- sqrt(VE_hat)
#         f0_hat <- f0r[N_hat]
#         mu <- (1-f0_hat)/(N_hat*arg_mat$generations[j])
#         ### USE TO SHOW GOODNESS OF FIT
#         xs <- density(d)$x
#         ys <- density(d)$y
#         ys <- ys/(sum(ys))
#         #v_hat <- scale_vprime(generate_vprime(N_hat, param$generations), f0_hat,
#                           # N_hat)
#         # lik_est <- pdf(xs, gen=param$gen, v=v_hat, N=N_hat, sd=sd_hat)
#
#         # qqplot(ys, lik_est)
#         # abline(0,1)
#
#     }
# } else {
#     isolates <- sapply(as.character(arg_mat[,"files"]), get_isolate)
#     models <- lapply(1:nrow(arg_mat), function(i) {
#         vcf <- read_vcf(as.character(arg_mat[i, "files"]))
#         d <- get_vcf_afs(vcf)
#         cov <- round(arg_mat[i, "coverage"])
#         replicated_d <- rbinom(length(d), cov, d)/cov
#         ds <- cbind(d, replicated_d)
#         mut_mle(d=ds, gen=arg_mat[i, "generations"],
#                 bp=arg_mat[i, "bases"]*arg_mat[i, "nsamples"], threads=1,
#                 N_max=400)
#     })
#     names(models) <- isolates
# }

