# #!/bin/bash
# args <- commandArgs(TRUE)
# args <- as.numeric(args)
# 
# RowOfDesign <- args[1]
# Replication <- args[2]

## Note: Simulations are very computationally intensive and were run on a cluster

args <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))

Repetition <- args

library("partykit")
library("glmertree")
library("Formula")
library("LongCART")
library("semtree")
library("lavaan")

##################################################
##
## Generate training data
##

## Data-generating design
N <- c(100, 250) # number of clusters
sigma_int <- c(1, 2) # SD of random intercept
sigma_slope <- c(sqrt(.1), sqrt(.4)) # SD of random slope
p_noise <- c(5, 25) # number of noise variables
rho <- c(0, 0.3) # correlation between partitioning variables
design_matrix <- expand.grid(N = N, sigma_int = sigma_int, 
                             sigma_slope = sigma_slope, 
                             p_noise = p_noise, rho = rho)

set.seed(Repetition + 41)

## Assume same number of measurement occasions for all datasets
n_timepoints <- 5
## Assume same error variance for all datasets:
sigma_epsilon <- sqrt(5) # SD of error
## Assume same variance for all partitioning variables
sigma_x <- 5

for (RowOfDesign in 1:NROW(design_matrix)) {

  Design <- design_matrix[RowOfDesign,]
  n_rep <- Repetition

  ## Generate the actual dataset

  ## Generate partitioning variables (including noise variables)
  nvars <- design_matrix[RowOfDesign, "p_noise"] + 3
  ## Create correlation matrix
  Rho <- matrix(design_matrix[RowOfDesign, "rho"], nrow = nvars, ncol = nvars)
  diag(Rho) <- 1
  ## Generate nvars variables according to correlation matrix
  data <- data.frame(t(t(chol(Rho)) %*% matrix(
    rnorm(nvars * design_matrix[RowOfDesign, "N"], sd = sigma_x), 
    nrow = nvars, ncol = design_matrix[RowOfDesign, "N"])))
  names(data) <- paste0("x", 1:nvars)
  data$x1 <- ifelse(data$x1 > 0, 1, 0)

  ## Generate intercepts and slopes, with means conditional on x1, x2, and x3,
  ## and SDs the same for all
  intercepts <- ifelse(data$x1 == 0, -1, 1)
  intercepts <- rnorm(design_matrix[RowOfDesign, "N"], mean = intercepts, 
                      sd = design_matrix[RowOfDesign, "sigma_int"])
  slopes <- ifelse(data$x1 == 0 & data$x2 < 0, -1, 0)
  slopes <- slopes + ifelse(data$x1 == 1 & data$x3 >= 0, 1, 0)
  slopes <- rnorm(design_matrix[RowOfDesign, "N"], mean = slopes, 
                  sd = design_matrix[RowOfDesign, "sigma_slope"])
  data <- data.frame(data, intercepts, slopes, 
                     cluster_id = factor(1:design_matrix[RowOfDesign, "N"]))
  
  ## Repeat the data for each timepoint
  data <- data[rep(seq_len(design_matrix[RowOfDesign, "N"]), n_timepoints), ]
  ## Add column with time for each timepoint
  data$time <- rep(0:(n_timepoints - 1), each = design_matrix[RowOfDesign, "N"])
  ## then y = intercept + time*slope + error
  data$y <- with(data, intercepts + slopes * time + 
                   rnorm(design_matrix[RowOfDesign, "N"] * n_timepoints, 
                         sd = sigma_epsilon)) 

  data <- data[ , -which(names(data) %in% c("intercepts", "slopes"))]
  data$x1 <- factor(data$x1)
  datasets <- data
  
  setwd("/exports/fsw/claramuntj/Marjolein")
  #save(datasets, file = paste0("datasets",RowOfDesign, "Repetition", Repetition, ".Rdata"))
  
  
  ################################
  ##
  ## Fit LM(M) trees
  ##
  
  inds <- c("lmer", "lmer_c", "lmer_r", "lmer_cr",
            "lmer_s", "lmer_s_c", "lmer_s_r", "lmer_s_cr",
            "lm", "lm_c")
  
  get_1stsplitvar <- function(object) {
    if (length(object) > 1L) {
      return(names(object$data)[object[[1]]$node$split$varid])
    } else {
      return(NA)
    }
  }
  
  random_effects <- list()
  splitvars <- times <- tree_sizes <- as.data.frame(
    matrix(NA, nrow = 1, ncol = length(inds)))
  names(tree_sizes) <- names(times) <- inds
  
  trees <- list()
  p <- if (ncol(datasets) > 12L) 28L else 8L
  
  ## Fit LMM trees
  form <- Formula(formula(paste0("y ~ time | cluster_id | ",
                                 paste0("x", 1L:p, collapse = " + "))))
  form.s <- Formula(formula(paste0("y ~ time | ((1+time)|cluster_id) | ",
                                   paste0("x", 1L:p, collapse = " + "))))
  times[1, 1L] <- system.time(trees[[1L]] <- try(
    lmertree(form, data = datasets)))["elapsed"]
  times[1, 2L] <- system.time(trees[[2L]] <- try(
    lmertree(form, data = datasets, cluster = cluster_id)))["elapsed"]
  times[1, 3L] <- system.time(trees[[3L]] <- try(
    lmertree(form, data = datasets, ranefstart = TRUE)))["elapsed"]
  times[1, 4L] <- system.time(trees[[4L]] <- try(
    lmertree(form, data = datasets, cluster = cluster_id,
             ranefstart = TRUE)))["elapsed"]
  times[1, 5L] <- system.time(trees[[5L]] <- try(
    lmertree(form.s, data = datasets)))["elapsed"]
  times[1, 6L] <- system.time(trees[[6L]] <- try(
    lmertree(form.s, data = datasets, cluster = cluster_id)))["elapsed"]
  times[1, 7L] <- system.time(trees[[7L]] <- try(
    lmertree(form.s, data = datasets, ranefstart = TRUE)))["elapsed"]
  times[1, 8L] <- system.time(trees[[8L]] <- try(
    lmertree(form.s, data = datasets, cluster = cluster_id,
             ranefstart = TRUE)))["elapsed"]
  for (k in 1L:8L) tree_sizes[1, k] <- try(length(trees[[k]]$tree))
  
  ## Extract first splitting variable
  for (k in 1L:8L) {
    splitvars[1, k] <- get_1stsplitvar(trees[[k]]$tree[[1]])       
  }
  
  ## Fit LM trees
  form_lmtree <- Formula(formula(paste0(
    "y ~ time | ", paste0("x", 1L:p, collapse = " + "))))
  times[1, 9L] <- system.time(trees[[9L]] <- try(
    lmtree(form_lmtree, data = datasets)))["elapsed"]
  times[1, 10L] <- system.time(trees[[10L]] <- try(
    lmtree(form_lmtree, data = datasets,
           cluster = cluster_id)))["elapsed"]
  tree_sizes[1, 9L] <- try(length(trees[[9L]]))
  tree_sizes[1, 10L] <- try(length(trees[[10L]]))
  splitvars[1, 9L] <- get_1stsplitvar(trees[[9L]]) 
  splitvars[1, 10L] <- get_1stsplitvar(trees[[10L]]) 
  
  save(tree_sizes, file = paste0("LMM_sizes", RowOfDesign, "Repetition", Repetition, ".Rdata"))
  save(times, file = paste0("LMM_times", RowOfDesign, "Repetition", Repetition, ".Rdata"))
  save(splitvars, file = paste0("LMM_splitvars", RowOfDesign, "Repetition", Repetition, ".Rdata"))
  
  
  
  
  
  ################################
  ##
  ## Fit LongCART
  ##

  inds <- c("LongCART")
  
  trees <- list()
  splitvars <- times <- tree_sizes <- as.data.frame(
    matrix(NA, nrow = 1, ncol = length(inds)))
  names(tree_sizes) <- names(times) <- inds
  p <- if (ncol(datasets) > 12L) 28L else 8L
  data <- datasets
  times[1, 1L] <- system.time(trees[[1L]] <- try(
    LongCART(data = data, patid = "cluster_id", fixed = y ~ time,
             gvars = paste0("x", 1L:p),
             tgvars = c(0, rep(1L, times = p-1L)))))["elapsed"]
  tree_sizes[1, 1L] <- try(nrow(trees[[1L]]$Treeout))
  splitvars[1, 1L] <- try(trees[[1L]]$Treeout$var[1])
  
  save(tree_sizes, file = paste0("Lcrt_sizes", RowOfDesign, "Repetition", Repetition, ".Rdata"))
  save(times, file = paste0("Lcrt_times", RowOfDesign, "Repetition", Repetition, ".Rdata"))
  save(splitvars, file = paste0("Lcrt_splitvars", RowOfDesign, "Repetition", Repetition, ".Rdata"))
  
  

  
  
  #################################
  ##
  ## Fit SEM trees
  ##
  
  #n_rep <- 100
  inds <- c("f_sem", "f_sem_i", "f_sem_s", "s_sem", "s_sem_i", "s_sem_s")
  
  get_1stsplitvar <- function(object) {
    split <- semtree:::getNodeList(object)[[1L]]$caption
    if (split == "TERMINAL") return(NA) else return(strsplit(split, " ")[[1L]][1L])
  }
  
  trees <- list()
  splitvars <- times <- tree_sizes <- as.data.frame(
    matrix(NA, nrow = 1, ncol = length(inds)))
  names(tree_sizes) <- names(times) <- inds
  p <- if (ncol(datasets) > 12L) 28L else 8L
  data <- reshape(data = datasets, v.names = "y",
                  timevar = "time", idvar = "cluster_id",
                  direction = "wide")
  
  mod <- '
    i =~ 1*y.0 + 1*y.1 + 1*y.2 + 1*y.3 + 1*y.4
    s =~ 1*y.1 + 2*y.2 + 3*y.3 + 4*y.4
    s ~~ 0*s
    i ~~ 0*i
    i ~~ 0*s
    '
  fit <- growth(mod, data = data, do.fit = FALSE)
  
  mod.i <- '
    i =~ 1*y.0 + 1*y.1 + 1*y.2 + 1*y.3 + 1*y.4
    s =~ 1*y.1 + 2*y.2 + 3*y.3 + 4*y.4
    s ~~ 0*s
    i ~~ 0*s
  '
  fit.i <- growth(mod.i, data = data, do.fit = FALSE)
  
  mod.s <- '
    i =~ 1*y.0 + 1*y.1 + 1*y.2 + 1*y.3 + 1*y.4
    s =~ 1*y.1 + 2*y.2 + 3*y.3 + 4*y.4
  '
  fit.s <- growth(mod.s, data = data, do.fit = FALSE)
  
  times[1, 1L] <- system.time(trees[[1L]] <- try(
    semtree(fit, predictors = paste0("x", 1L:p), data = data,
            control = semtree.control(method = "naive", 
                                      bonferroni = TRUE))))["elapsed"]
  times[1, 2L] <- system.time(trees[[2L]] <- try(
    semtree(fit.i, predictors = paste0("x", 1L:p), data = data,
            control = semtree.control(method = "naive", 
                                      bonferroni = TRUE))))["elapsed"]
  times[1, 3L] <- system.time(trees[[3L]] <- try(
    semtree(fit.s, predictors = paste0("x", 1L:p), data = data,
            control = semtree.control(method = "naive", 
                                      bonferroni = TRUE))))["elapsed"]
  
  times[1, 4L] <- system.time(trees[[4L]] <- try(
    semtree(fit, predictors = paste0("x", 1L:p), data = data,
            control = semtree.control(method = "score", 
                                      bonferroni = TRUE))))["elapsed"]
  times[1, 5L] <- system.time(trees[[5L]] <- try(
    semtree(fit.i, predictors = paste0("x", 1L:p), data = data,
            control = semtree.control(method = "score", 
                                      bonferroni = TRUE))))["elapsed"]
  times[1, 6L] <- system.time(trees[[6L]] <- try(
    semtree(fit.s, predictors = paste0("x", 1L:p), data = data,
            control = semtree.control(method = "score", 
                                      bonferroni = TRUE))))["elapsed"]
  
  for (k in 1L:6L) {
    tree_sizes[1, k] <- try(semtree:::getNumNodes(trees[[k]]))
    splitvars[1, k] <- try(get_1stsplitvar(trees[[k]]))
  }
  
  ## Save results
  save(tree_sizes, file = paste0("SEM_sizes", RowOfDesign, "Repetition", Repetition, ".Rdata"))
  save(times, file = paste0("SEM_times", RowOfDesign, "Repetition", Repetition, ".Rdata"))
  save(splitvars, file = paste0("SEM_splitvars", RowOfDesign, "Repetition", Repetition, ".Rdata"))
  
}