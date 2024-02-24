### R code from vignette source '_Partitioning_GCMs_with_GLMM_trees.Rnw'

library("glmertree")
load("Science ability data.Rdata")
load("sizesScience.Rda")
Ns <- c(250, 1000)
nreps <- 100L

## Get sample (design has changed over time, need to recover same sample previously used)
set.seed(42)
for (resp in c("Math", "Reading", "Science")) {
  for (N in Ns) {
      
    ## Read in data and prepare objects for saving results
    load(paste0(resp, " ability data.Rdata"))
    resp_short <- ifelse(resp == "Math", "math", ifelse(resp == "Reading", "read", "scie"))
    data <- get(paste0(resp_short, "data"))
    data$CHILDID <- factor(data$CHILDID)
    
    ## Generate samples for each repetition
    bag_ids <- sapply(1L:nreps, function(x) sample(unique(data$CHILDID), size = N))
    if (N == 250L && resp == "Science") {
      ## Get training sample
      intro_sample <- bag_ids[ , 39L]
      ## Get equally-sizes sample of test observations
      tmp <- bag_ids
    }
  }
}
sciedata$CHILDID <- factor(sciedata$CHILDID)
sciedata$GENDER <- factor(sciedata$GENDER)
sciedata$RACE <- factor(sciedata$RACE)
## simplify variable names
names(sciedata) <- gsub("P1", "", gsub("T1", "", gsub("C1", "", names(sciedata))))
names(sciedata)[names(sciedata) == "WKSESL"] <- "SES"
data <- sciedata[sciedata$CHILDID %in% intro_sample, ]
library("lmerTest")
LMM <- lmer(score ~ months + (1|CHILDID), data = data)

y <- t(sapply(unique(data$CHILDID), function(x) data$score[data$CHILDID == x]))
x <- t(sapply(unique(data$CHILDID), function(x) data$months[data$CHILDID == x]))
plot(data$months, data$score, type = "n", xlab = "Months after baseline", 
     ylab = "Science ability", xaxt = "n", cex.lab = 1, cex.axis=1,
     xlim = c(0.5, 57^(2/3)), main = "Full sample (N = 250)")
axis(1, at=c(0, 12, 36, 60)^(2/3), labels=c(0, 12, 36, 60), cex.axis = .8)
sapply(1:nrow(x), function(row) lines(x[row, ], y[row, ], col = gray(0.2, alpha = 1/6), lwd = 2))
lines(x = c(min(data$months), max(data$months)), 
     y = c(fixef(LMM)[1], fixef(LMM)[1] + max(data$months)*fixef(LMM)[2]),
     lwd = 3, col = 2)

LMM_form <- score ~ months | CHILDID | GENDER + RACE + SES + GMOTOR + FMOTOR + 
  INTERN + EXTERN + INTERP + CONTRO + FIRKDG + AGEBASELINE

## Fit models
lmm_tree <- lmertree(LMM_form, data = data)
lmm_tree_c <- lmertree(LMM_form, data = data, cluster = CHILDID)
lmm_tree_r <- lmertree(LMM_form, data = data, ranefstart = TRUE)

## Plot tree
plot(lmm_tree_r, which = "growth", ip_args = list(pval = FALSE), 
     nodesize_level = "CHILDID", tp_args = 
       list(xscale = c(0, 60^(2/3)), xaxis.at = c(0, 12, 36, 60)^(2/3),
            xaxis.labs = c(0, 12, 36, 60)), gp = gpar(cex = 1.5))


library("kableExtra")
coefs <- coef(lmm_tree_r)
coefs <- data.frame(cbind(Node = as.numeric(rownames(coefs)), coefs))
colnames(coefs) <- c("Node", "Intercept", "Slope")
coefs$Slope <- paste0("$", formatC(coefs$Slope, digits = 3, format = "f"), "$")
coefs$Intercept <- paste0("$", formatC(coefs$Intercept, digits = 3, format = "f"), "$")
coefs$Intercept[!grepl("-", coefs$Intercept)] <- paste("$\\:\\:$",
                                                      coefs$Intercept[!grepl("-", coefs$Intercept)])

## function nodestats returns estimated coefficients, standard errors and number of observations per node
## if newdata is specified, estimates are based on the new data, otherwise based on training data
nodestats <- function(object, newdata = NULL, drop = FALSE, ...) {
  
  merMod_type <- ifelse(inherits(object, "lmertree"), "lmer", "glmer")
  
  if (is.null(newdata)) warning("Standard errors computed on training data are likely overly optimistic, because they do not account for the searching of the tree structure.") 
    
  ## Preparations
  if (!is.null(newdata)) {
    newdata$.tree <- factor(predict(object, newdata = newdata, type = "node"), 
                            levels = levels(object$data$.tree))
    if (merMod_type == "lmer") {
      object$lmer <- lmer(update(formula(object$lmer), ~ .), 
                          data = newdata)
    } else {
      object$glmer <- glmer(update(formula(object$glmer), ~ .), 
                            data = newdata, family = family(object$glmer)$family)
    }
  } else {
    newdata <- object$data
  }
  
  ## Get coefficients and SEs
  coefs <- fixef(object, which = "tree", drop = drop)
  ses <- glmertree:::get_merMod_SEs(object, which = "tree", global_intercept = TRUE)
  
  ## Get node counts
  ns <- as.data.frame(table(predict(object, newdata = newdata, type = "node")))
  ns <- matrix(ns[ , 2L], dimnames = list(ns[ , 1L], "Freq"))
  
  return(list(coef = coefs, se = ses, n = ns))
}
set.seed(44)
test_sample <- sample(unique(sciedata$CHILDID)[!unique(sciedata$CHILDID) %in% unique(intro_sample)], size = 250)
test_data <- sciedata[sciedata$CHILDID %in% test_sample, ]

## Add columns: training N, test intercept, test intercept se, test slope, test slope se, test n
nodestats_train <- nodestats(lmm_tree_r)
nodestats_test <- nodestats(lmm_tree_r, newdata = test_data)
coefs <- cbind(coefs,
cbind(cbind(n = formatC(nodestats_train$n[,"Freq"], digits = 0, format = "f"), 
      Intercept = formatC(nodestats_test$coef[ , "(Intercept)"], digits = 3, format = "f"),
      s.e. = formatC(nodestats_train$se[ , "(Intercept)"], digits = 3, format = "f"),
      Slope = formatC(nodestats_test$coef[ , "months"], digits = 3, format = "f"),
      s.e. = formatC(nodestats_test$se[ , "months"], digits = 3, format = "f"),
      n = formatC(nodestats_test$n[,"Freq"], digits = 0, format = "f"))))

kable_styling(add_footnote(add_header_above(
  kable(coefs, format = "latex", booktabs = TRUE, label = "local_coefs",
                    align = "ccccccccc", escape=FALSE, linesep = "", # linesep="" suppressess addlinesep every 5 rows
                    caption = "Estimated fixed-effects coefficients for the terminal nodes of Figure 2.",
                    row.names = FALSE, digits = 3L),
header = c(" " = 1, "Training data" = 3, "Test data" = 5)),
"\\footnotesize \\\\ \\textit{Note.} Training data refers to the observations also depicted in Figures 1 and 2, test data refers to a new sample, used to (re)compute node-specific statistics and obtain valid uncertainty estimates. s.e. = standard error, n = node-specific sample size.", notation = "none", threeparttable = TRUE, escape = FALSE),
              font_size = 11, full_width=FALSE, position = "left")


lmm_tree <- lmertree(LMM_form, data = data, maxdepth = 6)
plot(lmm_tree, which = "growth", ip_args = list(pval = FALSE), 
     nodesize_level = "CHILDID", tp_args = 
       list(xscale = c(0, 60^(2/3)), xaxis.at = c(0, 12, 36, 60)^(2/3),
            xaxis.labs = c(0, 12, 36, 60)), gp = gpar(cex = 1))

vcov_d <- as.data.frame(VarCorr(lmm_tree))
vcov_r <- as.data.frame(VarCorr(lmm_tree_r))
vcov_c <- as.data.frame(VarCorr(lmm_tree_c))
vcv <- as.data.frame(VarCorr(LMM))
tab <- data.frame(Figure = 1:4,
                  Model = c("LMM", rep("LMM tree", times = 3)), 
                  Estimation = c("default", "random effects initialization", 
                                 "default", "clustered covariances"),
                  n_subs = c(1L, (length(lmm_tree_r$tree)+1)/2, 
                             (length(lmm_tree$tree)+1)/2, (length(lmm_tree_c$tree)+1)/2))
tab$sig_b <- c(vcv[vcv$grp == "CHILDID", "sdcor"], vcov_r[vcov_r$grp == "CHILDID", "sdcor"], 
               vcov_d[vcov_d$grp == "CHILDID", "sdcor"], vcov_c[vcov_c$grp == "CHILDID", "sdcor"])
tab$sig_epsb <- c(vcv[vcv$grp == "Residual", "sdcor"], vcov_r[vcov_r$grp == "Residual", "sdcor"], 
               vcov_d[vcov_d$grp == "Residual", "sdcor"], vcov_c[vcov_c$grp == "Residual", "sdcor"])


plot(lmm_tree_c, which = "growth", ip_args = list(pval = FALSE), 
     nodesize_level = 2, tp_args = list(xscale = c(0, 60^(2/3)), 
                                        xaxis.at = c(0, 12, 36, 60)^(2/3),
                                        xaxis.labs = c(0, 12, 36, 60)),
     gp = gpar(cex = 1.5))


library("kableExtra")
colnames(tab)[4:6] <- c("$J$", "$\\hat{\\sigma}_b$", "$\\hat{\\sigma}_\\epsilon$")
kable_styling(add_footnote(
  kable(tab, format = "latex", linesep = "", row.names = FALSE, 
        label = "ranef", digits = 3, align = c("cllccc"), booktabs = TRUE, escape = FALSE, 
        caption = "Estimated random-effects parameters for the mixed-effects models in Figures 1 through 4."),
  "\\footnotesize \\\\ \\textit{Note.} $J$ is the number of subgroups; $\\hat{\\sigma}_b$ is the estimated standard deviation of the random intercept; $\\hat{\\sigma}_\\epsilon$ is the estimated residual standard deviation.", notation="none", threeparttable = TRUE, escape = FALSE),
  font_size = 11, full_width=FALSE, position = "left")


library("partykit")
set.seed(42)
u1 <- sample(0:1, size = 1250, replace = TRUE) 
u2 <- round(rnorm(1250))
u3 <- round(rnorm(1250))
time <- rep(0:4, each = 250)
## left-most node: x1 == 0 & x2 <= 0
## intercept = -1, slope = -1
node3 <- as.numeric(u1 == 0 & u2 <= 0)
## left-middle node: x1 == 0 & x2 > 0
## intercept = -1, slope = 0
node4 <- as.numeric(u1 == 0 & u2 > 0)
## right-middle node: x1 == 1 & x3 =< 0
## intercept = 1, slope = 0
node6 <- as.numeric(u1 == 1 & u3 <= 0)
## right-most node: x1 == 1 & x3 > 0
## intercept = 1, slope = 1
node7 <- as.numeric(u1 == 1 & u3 > 0)
y <- -(node3 + node3*time) - node4 + node6 + (node7 + node7*time) 
data <- data.frame(u1 = as.factor(u1), u2, u3, time, y)   
tree <- lmtree(y ~ time | u1 + u2 + u3, data = data, maxdepth = 3)

fig <- party(
  partynode(1L,
            split = partysplit(1L, breaks = 1),
            kids = list(
              partynode(2L,
                        split = partysplit(2L, breaks = 0),
                        kids = list(partynode(3L, info = c(
                            expression(beta[j0] == '-1.0'),
                            expression(''),
                            expression(beta[j1] == '-1.0'))),
                          partynode(4L, info = c(
                            expression(beta[j0] == '-1.0'),
                            expression(''),
                            expression(beta[j1] == '0.0'))))),
              partynode(5L,
                        split = partysplit(3L, breaks = 0),
                        kids = list(
                          partynode(6L, info = c(
                            expression(beta[j0] == '1.0'),
                            expression(''),
                            expression(beta[j1] == '0.0'))),
                          partynode(7L, info = c(
                            expression(beta[j0] == '1.0'),
                            expression(''),
                            expression(beta[j1] == '1.0'))))))),
  data.frame(data)
)

## Panel-combining function:
combine_panels <- function(panel1, panel2, party2 = NULL) {
  function(node) {
    nid <- id_node(node)
    pushViewport(viewport(
      layout = grid.layout(nrow = 2, ncol = 1, heights = c(1, 1)),
      name = paste("node_mob_mypanel", nid, sep = "")))
    grid.rect(gp = gpar(fill = "white", col = 0))
    pushViewport(viewport(layout.pos.col = 1, layout.pos.row = 1))
    panel1(node)
    popViewport()
    pushViewport(viewport(layout.pos.col = 1, layout.pos.row = 2))
    node2 <- if(is.null(party2)) node else node_party(party2[nid])
    panel2(node2)
    popViewport(2)
  }
}


plot(tree, tnex = 3.5, terminal_panel = combine_panels(
    panel1 = glmertree:::node_growthplot(tree, xaxis.at = c(0, 4), xscal = c(0, 4), ylines = 2),
    panel2 = node_terminal(fig, FUN = identity, height = 4, width = 11, id = FALSE, fill = "white",
                           align = "right", just = "top", top = 0.85, gp = gpar(lty = 0, fill = "white")),
    party2 = fig), 
  gp = gpar(cex = .7), ip_args = list(pval = FALSE))


#library("devtools")
#install_github("brandmaier/semper")
library("semper")

## low to medium slope reliability for the model with the 
## small random slope:

## low sigma_b0 and low sigma_b1 
model1 <- lgcm(timepoints=0:4, slope.variance = 0.1,
              intercept.slope.covariance = 0, intercept.variance = 1,
              residual.variance = 5)
ecr(model1)
##[1] 0.2857143

## high sigma_b0 and low sigma_b
model2 <- lgcm(timepoints=0:4, slope.variance = 0.1,
              intercept.slope.covariance = 0, intercept.variance = 4,
              residual.variance = 5)
ecr(model2)
##[1] 0.21875

## medium to good slope reliability for the condition with large slope variance:
model3 <- lgcm(timepoints=0:4, slope.variance = 0.4, intercept.slope.covariance = 0,
              intercept.variance = 1, residual.variance = 5)
ecr(model3)
##[1] 0.6153846

model4 <- lgcm(timepoints=0:4, slope.variance = 0.4, intercept.slope.covariance = 0,
              intercept.variance = 4, residual.variance = 5)
ecr(model4)
##[1] 0.5283019


###########################
##
## Load tree size data
##

load("tree_sizes.Rda")
tree_size <- data.frame(sizes, true = 7)
N <- c(100, 250)
sigma_int <- c(1, 2)
sigma_slope <- c(sqrt(.1), sqrt(.4))
p_noise <- c(5, 25)
rho <- c(0, .3)
design_matrix <- expand.grid(N = N, sigma_int = sigma_int,
                             sigma_slope = sigma_slope,
                             p_noise = p_noise, rho = rho)
design_matrix[] <- lapply(design_matrix, factor)
tree_sizes_long <- data.frame(stack(tree_size),
                              dataset_id = factor(rep(1:nrow(tree_size), 
                                                      times = ncol(tree_size))))
names(tree_sizes_long)[1:2] <- c("tree_size", "method_id") 
tree_sizes_long <- cbind(tree_sizes_long, design_matrix)
tree_sizes_long$method_id <- relevel(tree_sizes_long$method_id, ref = "true")
LMM_ids <- which(tree_sizes_long$method_id %in% c("true", "lm", "lm_c", "lmer", "lmer_c", 
                                                  "lmer_cr", "lmer_r", "lmer_s",
                                                  "lmer_s_c", "lmer_s_cr", "lmer_s_r"))
stdCoef <- function(object) {
  sdy <- sd(getME(object,"y"))
  sdx <- apply(getME(object,"X"), 2, sd)
  sc <- fixef(object)*sdx/sdy
  se.fixef <- coef(summary(object))[,"Std. Error"]
  se <- se.fixef*sdx/sdy
  return(data.frame(stdcoef=sc, stdse=se))
}


library("lmerTest")

## Compare means between settings
tmp <-  tree_sizes_long[LMM_ids, ]
tmp$tree_size <- (tmp$tree_size-1)/2
tmp$default_ids <- tmp$method_id %in% c("lm", "lmer", "lmer_s")
tmp$cluster_ids <- tmp$method_id %in% c("lm_c", "lmer_c", "lmer_s_c")
tmp$ranefstart_ids <- tmp$method_id %in% c("lmer_r", "lmer_s_r")
tmp$both_ids <- tmp$method_id %in% c("lmer_cr", "lmer_s_cr")

M_default <- mean(tmp$tree_size[tmp$default_ids])
M_cluster <- mean(tmp$tree_size[tmp$cluster_ids])
M_ranefstart <- mean(tmp$tree_size[tmp$ranefstart_ids])
M_both <- mean(tmp$tree_size[tmp$both_ids])

means <- tapply(tmp$tree_size, tmp$method_id, mean)
sds <- tapply(tmp$tree_size, tmp$method_id, sd)
mads <- tapply(tmp$tree_size, tmp$method_id, function(x) mean(abs(x-3)))
props_more_splits <- tapply(tmp$tree_size, tmp$method_id, function(x) mean(x > 3L))
props_less_splits <- tapply(tmp$tree_size, tmp$method_id, function(x) mean(x < 3L))
props_corr_splits <- tapply(tmp$tree_size, tmp$method_id, function(x) mean(x == 3L))



## Prepare data
tmp <-  tree_sizes_long
tmp$tree_size <- (tmp$tree_size-1)/2 ## count splits instead of nodes
tmp <- tmp[-which(tmp$method_id == "true") , ]
tmp$method_id <- factor(tmp$method_id)

## Create indicator for algorithm and random effects specification
tmp$ranef <- NA
tmp$ranef[tmp$method_id %in% c("lm", "lm_c")] <- "{hat(sigma)[b[0]] == hat(sigma)[b[1]]} == 0"
tmp$ranef[tmp$method_id %in% c("lmer", "lmer_c", "lmer_r", "lmer_cr")] <- "list(hat(sigma)[b[0]] > 0, hat(sigma)[b[1]] == 0)"
tmp$ranef[tmp$method_id %in% c("lmer_s", "lmer_s_c", "lmer_s_r", "lmer_s_cr")] <- "list(hat(sigma)[b[0]] > 0, hat(sigma)[b[1]] > 0)"
tmp$ranef[tmp$method_id %in% c("LongCART")] <- "LongCART"
tmp$ranef[tmp$method_id %in% c("s_sem", "f_sem")] <- "SEM trees \n(no random effects)"
tmp$ranef[tmp$method_id %in% c("s_sem_i", "f_sem_i")] <- "SEM tree \n(random intercept)"
tmp$ranef[tmp$method_id %in% c("s_sem_s", "f_sem_s")] <- "SEM tree \n(random intercept + slope)"

## Create indicator for additional settings
tmp$setting <- NA
tmp$setting[tmp$method_id %in% c("lmer", "lmer_s", "lm")] <- "default"
tmp$setting[tmp$method_id %in% c("lmer_c", "lmer_s_c", "lm_c")] <- "cl.cov."
tmp$setting[tmp$method_id %in% c("lmer_r", "lmer_s_r")] <- "ran.eff."
tmp$setting[tmp$method_id %in% c("lmer_cr", "lmer_s_cr")] <- "cl.cov.+\nran.eff."
tmp$setting[tmp$method_id %in% c("LongCART")] <- "default"
tmp$setting[tmp$method_id %in% c("s_sem", "s_sem_i", "s_sem_s")] <- "score"
tmp$setting[tmp$method_id %in% c("f_sem", "f_sem_i", "f_sem_s")] <- "fair"
tmp$setting <- factor(tmp$setting, levels = c("cl.cov.+\nran.eff.", "cl.cov.", 
                                              "default", "fair",
                                              "ran.eff.", "score"))
levels(tmp$sigma_int) <- c("1", "4")
levels(tmp$sigma_slope) <- c("0.1", "0.4")
tmp$setting <- factor(tmp$setting, levels = c("default", "cl.cov.",  "ran.eff.", 
                       "cl.cov.+\nran.eff.", "fair", "score"))

## Create the plot
library("ggplot2")
library("gridExtra")
library("colorspace")
library("scales")

## Set up data.frame of summary stats
summ_stats <- data.frame(M = round(means[c("lm", "lm_c", "lmer", "lmer_c", "lmer_cr", "lmer_r", 
                   "lmer_s", "lmer_s_c", "lmer_s_cr", "lmer_s_r")], digits = 2), 
                   SD = format(sds[c("lm", "lm_c", "lmer", "lmer_c", "lmer_cr", "lmer_r", 
                   "lmer_s", "lmer_s_c", "lmer_s_cr", "lmer_s_r")], digits = 2),
                   MAD = format(mads[c("lm", "lm_c", "lmer", "lmer_c", "lmer_cr", "lmer_r", 
                   "lmer_s", "lmer_s_c", "lmer_s_cr", "lmer_s_r")], digits = 2))
summ_stats$method <- rownames(summ_stats)
summ_stats$y <- 50
## Note!!! Labels of ran.eff. and cl.cov.+ran.eff. have been switched to have text labels
## above the right boxplot / countplot. Setting it as a factor and releveling does not help. 
summ_stats$setting <- c("default", "cl.cov.", 
                        "default", "cl.cov.", "cl.cov.+\nran.eff.", "ran.eff.", 
                        "default", "cl.cov.", "cl.cov.+\nran.eff.", "ran.eff.")

summ_stats$ranef <- NA
summ_stats$ranef[summ_stats$method %in% c("lm", "lm_c")] <- "{hat(sigma)[b[0]] == hat(sigma)[b[1]]} == 0"
summ_stats$ranef[summ_stats$method %in% c("lmer", "lmer_c", "lmer_r", "lmer_cr")] <- "list(hat(sigma)[b[0]] > 0, hat(sigma)[b[1]] == 0)"
summ_stats$ranef[summ_stats$method %in% c("lmer_s", "lmer_s_c", "lmer_s_r", "lmer_s_cr")] <- "list(hat(sigma)[b[0]] > 0, hat(sigma)[b[1]] > 0)"
summ_stats$all <- with(summ_stats, paste(M, SD, MAD, sep = "\n"))

## LM(M)
LMM_ids <- which(tmp$method_id %in% c("lm", "lm_c", "lmer", "lmer_c", "lmer_cr", 
                                      "lmer_r", "lmer_s", "lmer_s_c", 
                                      "lmer_s_cr", "lmer_s_r"))
breaks <- c(0, 1, 3, 7, 20, 38, 50, 64)
labels <- c("0", "1", "3", "7", "20", "MAD", "SD", "M")
lim <- c(0, 60)

## Plot
theme_set(theme_bw())
plot0 <- ggplot(tmp[LMM_ids, ]) +
  geom_boxplot(aes(x=setting, y=tree_size), 
               position=position_dodge(1), show.legend = FALSE, 
               alpha = 0.5, width = .6, outlier.shape = NA, coef = NULL) + 
  geom_count(aes(x=setting, y=tree_size), col = "gray", alpha = .75) +
  scale_y_continuous(trans = pseudo_log_trans(sigma = 1, base = exp(1)), 
                     breaks = breaks, labels = labels, lim = lim) +
  geom_text(data = summ_stats, aes(x=setting, y=y, label = all),
            size = 3, col = "gray44") +
  labs(x = "", y = "no. of splits") + theme(legend.position = "none") +
  geom_hline(yintercept=3, col = "darkgray") + 
  facet_grid(~ranef, scales = "free", space = "free", 
             labeller = label_parsed)  +
  coord_cartesian(clip = 'off') 
plot0


load(file = "splitvars.Rda")
splitvar <- apply(splits, 2, unlist)
rownames(splitvar) <- NULL
Lcrt_splitvars <- splitvar[, "Lcrt"] ## save for later use
splitvar[tree_size$Lcrt == 1, "Lcrt"] <- NA

## Get splitting proportions
func <- function(x) table(x, useNA = "ifany") / sum(table(x, useNA = "ifany"))
probs <- apply(splitvar, 2, func)
probs <- sapply(probs, function(x) {
  names(x)[is.na(names(x))] <- "No split"
  return(x)
})

## Have column for u1, u2, u3, u4-up, no split
splits <- matrix(sapply(probs, function(x) x["x1"]), ncol = 1,
                 dimnames = list(names(probs), "x1"))
splits[is.na(splits)] <- 0
for (u_set in list("x2","x3", paste0("x", 4:30), "No split")) {
  splits <- cbind(splits, sapply(probs, function(x) sum(x[names(x) %in% u_set])))
}
splits <- data.frame(splits[c(9:10, 1:8, 11:17), ]) ## reorder rows so LM tree is 1st
rows <- cbind(c("LM tree", " ", 
                "LMM tree", " ", " ", " ", " ", " ", " ", " ", 
                "SEM tree", " ", " ", " ", " ", " ", 
                "LongCART"),
                c("$\\hat{\\sigma}_{b_0} = \\hat{\\sigma}_{b_1} = 0$",
                  "$\\hat{\\sigma}_{b_0} = \\hat{\\sigma}_{b_1} = 0$",
                  "$\\hat{\\sigma}_{b_0} > 0$, $\\hat{\\sigma}_{b_1} = 0$",
                  "$\\hat{\\sigma}_{b_0} > 0$, $\\hat{\\sigma}_{b_1} = 0$",
                  "$\\hat{\\sigma}_{b_0} > 0$, $\\hat{\\sigma}_{b_1} = 0$",
                  "$\\hat{\\sigma}_{b_0} > 0$, $\\hat{\\sigma}_{b_1} = 0$",
                  "$\\hat{\\sigma}_{b_0} > 0$, $\\hat{\\sigma}_{b_1} > 0$",
                  "$\\hat{\\sigma}_{b_0} > 0$, $\\hat{\\sigma}_{b_1} > 0$",
                  "$\\hat{\\sigma}_{b_0} > 0$, $\\hat{\\sigma}_{b_1} > 0$",
                  "$\\hat{\\sigma}_{b_0} > 0$, $\\hat{\\sigma}_{b_1} > 0$",
                  "$\\hat{\\sigma}_{b_0} = \\hat{\\sigma}_{b_1} = 0$",
                  "$\\hat{\\sigma}_{b_0} > 0$, $\\hat{\\sigma}_{b_1} = 0$",
                  "$\\hat{\\sigma}_{b_0} > 0$, $\\hat{\\sigma}_{b_1} > 0$",
                  "$\\hat{\\sigma}_{b_0} = \\hat{\\sigma}_{b_1} = 0$",
                  "$\\hat{\\sigma}_{b_0} > 0$, $\\hat{\\sigma}_{b_1} = 0$",
                  "$\\hat{\\sigma}_{b_0} > 0$, $\\hat{\\sigma}_{b_1} > 0$",
                  "$\\hat{\\sigma}_{b_0} > 0$, $\\hat{\\sigma}_{b_1} = 0$"),
             c("default", "cl.cov.", "default", "cl.cov.", 
               "ran.eff.", "cl.cov. + ran.eff.", "default", "cl.cov.", 
               "ran.eff.", "cl.cov. + ran.eff.",
               "LRT", "LRT", "LRT", "score-based", "score-based", "score-based", ""))
splits <- cbind(rows, splits)
colnames(splits) <- c("Algorithm", "Random effects", "Fitting approach", "$u_1$", 
                      "$u_2$", "$u_3$", "$u_4$--$u_{25}$", "No split")

## Print the table
kable_styling(add_footnote(row_spec(kable(splits[1:10, -1], format = "latex", booktabs = TRUE, 
                                 label = "first_splits", align = "llccccc", 
        escape=FALSE, linesep = "", # linesep command suppressess addlinesep every 5 rows
        caption = "Variables selected for the first split by each LM tree (top two rows) and LMM tree (bottom eight rows) estimation approach.", row.names = FALSE, digits = 3L), row = c(2, 6), hline_after =  TRUE),
  "\\footnotesize \\\\ \\textit{Note.} $u_1$ is the true first splitting variable and is binary; all other partitioning variables are continuous, with $u_2$ and $u_3$ being true splitting variables (nodes 2 and 3). $\\hat{\\sigma}_{b_0}$ and $\\hat{\\sigma}_{b_1}$ are the estimated standard deviations of the random intercept and slope, respectively.", 
  notation="none", threeparttable = TRUE, escape = FALSE),
  font_size = 11, full_width=FALSE)


tmp2 <- tmp
LMM_c_ids <- which(tmp2$method_id %in% c("lm_c", "lmer_c", "lmer_s_c"))
SEM_ids <- which(tmp2$method_id %in% c("f_sem", "f_sem_i", "f_sem_s",
                                      "s_sem", "s_sem_i", "s_sem_s"))
Lcrt_ids <- which(tmp2$method_id %in% "Lcrt")
tmp2 <-  tmp2[c(LMM_c_ids, SEM_ids, Lcrt_ids), ]
tmp2$method_id <- factor(tmp2$method_id)


## Note!!! White spaces added to LM(M) and SEM tree, b/c with use of geom_text panels get ordered
## alphabetically
tmp2$method <- ifelse(
  tmp2$method_id %in% c("lm_c", "lmer_c", "lmer_s_c"), " LM(M) tree ", 
      ifelse(tmp2$method_id %in% c("f_sem", "f_sem_i", "f_sem_s"), " SEM tree (LRT-based) ", 
             ifelse(tmp2$method_id %in% c("s_sem", "s_sem_i", "s_sem_s"),
                             " SEM tree (score-based) ", "LongCART")))
tmp2$method <- factor(tmp2$method, levels = c(" LM(M) tree ", " SEM tree (LRT-based) ", 
                                              " SEM tree (score-based) ", "LongCART"))

tmp2$random <- ifelse(tmp2$method_id %in% c("lm_c", "f_sem", "s_sem"), 
                      "atop(hat(sigma)[b[0]] == 0, hat(sigma)[b[1]] == 0)", 
                      ifelse(tmp2$method_id %in% c("lmer_c", "f_sem_i", "s_sem_i", "Lcrt"),
                             "atop(hat(sigma)[b[0]] > 0, hat(sigma)[b[1]] == 0)", 
                              "atop(hat(sigma)[b[0]] > 0, hat(sigma)[b[1]] > 0)"))

## Set up data.frame of summary stats
means <- tapply(tmp2$tree_size, tmp2$method_id, mean)
sds <- tapply(tmp2$tree_size, tmp2$method_id, sd)
mads <- tapply(tmp2$tree_size, tmp2$method_id, function(x) mean(abs(x - 3)))
summ_stats <- data.frame(M = round(means[c("lm_c", "lmer_c", "lmer_s_c", 
                                            "f_sem", "f_sem_i", "f_sem_s",
                                            "s_sem", "s_sem_i", "s_sem_s", "Lcrt")], digits = 2), 
                   SD = round(sds[c("lm_c", "lmer_c", "lmer_s_c", 
                                            "f_sem", "f_sem_i", "f_sem_s",
                                            "s_sem", "s_sem_i", "s_sem_s", "Lcrt")], digits = 2),
                   MAD = sprintf(mads[c("lm_c", "lmer_c", "lmer_s_c", 
                                            "f_sem", "f_sem_i", "f_sem_s",
                                            "s_sem", "s_sem_i", "s_sem_s", "Lcrt")],  fmt = '%#.2f'))
summ_stats$method_id <- factor(rownames(summ_stats))

summ_stats$method <- ifelse(
  summ_stats$method_id %in% c("lm_c", "lmer_c", "lmer_s_c"), " LM(M) tree ", 
      ifelse(summ_stats$method_id %in% c("f_sem", "f_sem_i", "f_sem_s"), " SEM tree (LRT-based) ", 
             ifelse(summ_stats$method_id %in% c("s_sem", "s_sem_i", "s_sem_s"),
                             " SEM tree (score-based) ", "LongCART")))
tmp2$method <- factor(tmp2$method, levels = c(" LM(M) tree ", " SEM tree (LRT-based) ", 
                                              " SEM tree (score-based) ", "LongCART"))
summ_stats$y <- 25.5
summ_stats$random <- ifelse(summ_stats$method_id %in% c("lm_c", "f_sem", "s_sem"), 
                      "atop(hat(sigma)[b[0]] == 0, hat(sigma)[b[1]] == 0)", 
                      ifelse(summ_stats$method_id %in% c("lmer_c", "f_sem_i", "s_sem_i", "Lcrt"),
                             "atop(hat(sigma)[b[0]] > 0, hat(sigma)[b[1]] == 0)", 
                              "atop(hat(sigma)[b[0]] > 0, hat(sigma)[b[1]] > 0)"))
summ_stats$all <- with(summ_stats, paste(M, SD, MAD, sep = "\n"))

breaks <- c(0, 1, 3, 7, 15, 20, 25.5, 32)
labels <- c("0", "1", "3", "7", "15", "MAD", "SD", "M")
lim <- c(0, 30)
plot1 <- ggplot(tmp2) +
  geom_boxplot(aes(x=random, y=tree_size), alpha = .5, width = .6, coef = NULL) +
  geom_count(aes(x=random, y=tree_size), col = "gray", alpha = .75) +
  geom_hline(yintercept=3, col = "darkgray") +
  geom_text(data = summ_stats, aes(x=random, y=y, label = all),
            size = 3, col = "gray44") +
  scale_y_continuous(trans = pseudo_log_trans(), 
                     breaks = breaks, lim = lim, labels = labels) +
  facet_grid(~ method, scales = "free", space = "free") +
  labs(x = "", y = "no. of splits") + theme(legend.position = "none",
                       text=element_text(size=9)) +
  scale_x_discrete(labels = scales::parse_format()) +
  coord_cartesian(clip = 'off')
plot1


library("kableExtra")
splits$Algorithm[2] <- "LM(M) tree"
splits$Algorithm[4] <- "  (cl.cov.)"
splits$Algorithm[8] <- " "
splits$Algorithm[12] <- "  (LRT-based)"
splits$Algorithm[14] <- "SEM tree"
splits$Algorithm[15] <- "  (score-based)"

## Print table
kable_styling(add_footnote(
  kable(splits[c(2, 4, 8, 11:17), -3], format = "latex", booktabs = TRUE, 
        label = "first_splits2", align = "llccccc", row.names = FALSE, digits = 3,
        escape=FALSE, linesep = "", # linesep command suppressess addlinesep every 5 rows
        caption = "Variables selected for the first splits by each of the partitioning approaches."),
  "\\footnotesize \\\\ \\textit{Note.} $u_1$ is the true first splitting variable and is a binary factor; all other partitioning variables are continuous, with $u_2$ and $u_3$ being true splitting variables (nodes 2 and 3). The first column indicates wheter the random intercept and/or slope were estimated or not.", 
  notation="none", threeparttable = TRUE, escape = FALSE),
  font_size = 11, full_width=FALSE)


load("comp_times.Rda")
colnames(times) <- c(rep("LMM tree", times = 8),
                     rep(" LM tree", times = 2), 
                     rep("SEM tree\n(LRT)", times = 3),
                     rep("SEM tree\n(score)", times = 3),
                     "LongCART")
times <- reshape2::melt(times)
mean_times <- tapply(times$value, times$Var2, mean)
summ_stats <- data.frame(M = sprintf(mean_times, fmt = '%#.2f'),
                         Var2 = levels(times$Var2),
                         y = 100000)

ggplot(times, aes(Var2, value)) + 
  geom_boxplot(alpha = .5, fill = "gray", width=.5) +
  labs(x = "", y = "Computation time (in seconds)") +
  geom_text(data = summ_stats, aes(x=Var2, y=y, label = M),
            size = 3, col = "gray44") +
  theme(text = element_text(size = 10)) +
  scale_y_continuous(trans = "log", limits = c(0.01, 100000), 
                     breaks = c(.01, .1, 1, 10, 100, 1000, 10000, 100000),
                     labels = c("1e-02", "1e-01", "1e+00", "1e+01", "1e+02", "1e+03", "1e+04", "M"))


## Load MSEs
load("MSEsMath.Rda")
math_MSEs <- MSEs
load("MSEsReading.Rda")
read_MSEs <- MSEs
load("MSEsScience.Rda")
scie_MSEs <- MSEs

## Load datasets
load("Reading ability data.Rdata")
var_read <- var(readdata$score)
load("Math ability data.Rdata")
var_math <- var(mathdata$score)
load("Science ability data.Rdata")
var_scie <- var(sciedata$score)

R2 <- round(data.frame(scie_min = 1 - sapply(scie_MSEs, min)/var_scie, 
                       read_min = 1 - sapply(read_MSEs, min)/var_read,
                       math_min = 1 - sapply(math_MSEs, min)/var_math,
                       scie_max = 1 - sapply(scie_MSEs, max)/var_scie,
                       read_max = 1 - sapply(read_MSEs, max)/var_read,
                       math_max = 1 - sapply(math_MSEs, max)/var_math), digits = 2)

load("vars_selectedMath.Rda")
vs_math <- vars_selected
load("vars_selectedReading.Rda")
vs_read <- vars_selected
load("vars_selectedScience.Rda")
vs_scie <- vars_selected


## Create table with M, SD, R2

## Math
MSE_df <- sprintf("%.4f", sapply(math_MSEs[c(1:5, 7, 11)], mean, na.rm = TRUE))
MSE_df <- cbind(MSE_df, sprintf("%.3f", sapply(math_MSEs[c(1:5, 7, 11)], sd, na.rm = TRUE)))
MSE_df <- cbind(MSE_df, sprintf("%.3f", sapply(math_MSEs[c(1:5, 7, 11)], \(x) 1 - mean(x, na.rm = TRUE)/var_math)))

## Read
MSE_df <- cbind(MSE_df, sprintf("%.4f", sapply(read_MSEs[c(1:5, 7, 11)], mean, na.rm = TRUE)))
MSE_df <- cbind(MSE_df, sprintf("%.3f", sapply(read_MSEs[c(1:5, 7, 11)], sd, na.rm = TRUE)))
MSE_df <- cbind(MSE_df, sprintf("%.3f", sapply(read_MSEs[c(1:5, 7, 11)], \(x) 1 - mean(x, na.rm = TRUE)/var_read)))

## Science
MSE_df <- cbind(MSE_df, sprintf("%.4f", sapply(scie_MSEs[c(1:5, 7, 11)], mean, na.rm = TRUE)))
MSE_df <- cbind(MSE_df, sprintf("%.3f", sapply(scie_MSEs[c(1:5, 7, 11)], sd, na.rm = TRUE)))
MSE_df <- cbind(MSE_df, sprintf("%.3f", sapply(scie_MSEs[c(1:5, 7, 11)], \(x) 1 - mean(x, na.rm = TRUE)/var_scie)))

## tmp3 is math
tmp3 <- data.frame(stack(math_MSEs[ , c(1:5, 7, 11)]),
                    dataset_id = factor(rep(1:nrow(math_MSEs), times = 7)))
names(tmp3)[1:2] <- c("MSE", "method") 
tmp3$method2 <- as.character(tmp3$method)
tmp3$method2[grepl("LMtree", tmp3$method2)] <- "LM tree\n\n(cluster)"
tmp3$method2[tmp3$method2 %in% c("LMMtree_c", "LMMtree_sc")] <- "LMM tree\n\n(cluster)"
tmp3$method2[tmp3$method2 == "LMMtree_r"] <- "LMM tree\n\n(ranef)"
tmp3$method2[tmp3$method2 %in% c("LMMtree_cr", "LMMtree_scr")] <- "LMM tree\n\n(cluster + ranef)"
tmp3$method2[tmp3$method2 == "LMtree_c"] <- "LM tree\n\n(cl.cov.)"
tmp3$method2[tmp3$method2 == "LMMtree"] <- "LMM tree"
tmp3$method2[tmp3$method2 %in% c("LMMtree_c", "LMMtree_sc")] <- "LMM tree\n\n(cl.cov.)"
tmp3$method2[tmp3$method2 == "LMMtree_r"] <- "LMM tree\n\n(ran.eff.)"
tmp3$method2[tmp3$method2 %in% c("LMM_s", "LMM_s3")] <- " LMM"
tmp3$random <- ifelse(tmp3$method == "LMtree_c", "atop(hat(sigma)[b[0]] == 0, hat(sigma)[b[1]] == 0)", 
                      ifelse(tmp3$method %in% c("LMMtree", "LMMtree_c", "LMMtree_r"), 
                             "atop(hat(sigma)[b[0]] > 0, hat(sigma)[b[1]] == 0)", 
                             ifelse(tmp3$method == "LMMtree_sc", 
                                "atop(hat(sigma)[b[0]] > 0, hat(sigma)[b[1]] > 0)",
                                ifelse(tmp3$method == "LMM_s", " time ", "all"))))
theme_set(theme_bw(base_size = 12))
ggplot(tmp3) +
  geom_boxplot(aes(x=random, y=MSE), 
               position=position_dodge(1), alpha = .5, width = .45, fill = "gray") +
  scale_y_continuous(#trans = log_trans(),
                     sec.axis = sec_axis(trans=~1-(./var_math))) +
  facet_grid(~method2, scales = "free", space = "free") +
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(),
        legend.position="none", axis.text.y = element_text(size = 8),
        axis.title.y = element_text(size = 10)) +
  labs(x = "", y = "MSE")


## tmp4 is reading
tmp4 <- data.frame(stack(read_MSEs[ , c(1:5, 7, 11)]),
                    dataset_id = factor(rep(1:nrow(read_MSEs), times = 7)))
names(tmp4)[1:2] <- c("MSE", "method") 
tmp4$method2 <- as.character(tmp4$method)
tmp4$method2[grepl("LMtree", tmp4$method2)] <- "LM tree\n\n(cluster)"
tmp4$method2[tmp4$method2 %in% c("LMMtree_c", "LMMtree_sc")] <- "LMM tree\n\n(cluster)"
tmp4$method2[tmp4$method2 == "LMMtree_r"] <- "LMM tree\n\n(ranef)"
tmp4$method2[tmp4$method2 %in% c("LMMtree_cr", "LMMtree_scr")] <- "LMM tree\n\n(cluster + ranef)"
tmp4$method2[tmp4$method2 == "LMtree_c"] <- "LM tree\n\n(cl.cov.)"
tmp4$method2[tmp4$method2 == "LMMtree"] <- "LMM tree"
tmp4$method2[tmp4$method2 %in% c("LMMtree_c", "LMMtree_sc")] <- "LMM tree\n\n(cl.cov.)"
tmp4$method2[tmp4$method2 == "LMMtree_r"] <- "LMM tree\n\n(ran.eff.)"
tmp4$method2[tmp4$method2 %in% c("LMM_s", "LMM_s3")] <- " LMM"
tmp4$random <- ifelse(tmp4$method == "LMtree_c", "atop(hat(sigma)[b[0]] == 0, hat(sigma)[b[1]] == 0)", 
                      ifelse(tmp4$method %in% c("LMMtree", "LMMtree_c", "LMMtree_r"), 
                             "atop(hat(sigma)[b[0]] > 0, hat(sigma)[b[1]] == 0)", 
                             ifelse(tmp4$method == "LMMtree_sc",
                                "atop(hat(sigma)[b[0]] > 0, hat(sigma)[b[1]] > 0)", 
                                ifelse(tmp4$method == "LMM_s", " time ", "all"))))
ggplot(tmp4) +
  geom_boxplot(aes(x=random, y=MSE), 
               position=position_dodge(1), alpha = .5, width = .45, fill = "gray") +
  scale_y_continuous(#trans = log_trans(),
                       sec.axis = sec_axis(trans=~1-(./var_read))) +
  facet_grid(~method2, scales = "free", space = "free") +
  theme(strip.background = element_blank(), strip.text.x = element_blank(),
        axis.title.x=element_blank(), axis.text.x=element_blank(),
        legend.position="none", axis.text.y = element_text(size = 8),
        axis.title.y = element_text(size = 10)) +
  labs(x = "", y = "MSE")


## tmp5 is science
tmp5 <- data.frame(stack(scie_MSEs[ , c(1:5, 7, 11)]),
                    dataset_id = factor(rep(1:nrow(scie_MSEs), times = 7)))
names(tmp5)[1:2] <- c("MSE", "method") 
tmp5$method2 <- as.character(tmp5$method)
tmp5$method2[tmp5$method2 == "LMtree_c"] <- "LM tree\n\n(cl.cov.)"
tmp5$method2[tmp5$method2 == "LMMtree"] <- "LMM tree"
tmp5$method2[tmp5$method2 %in% c("LMMtree_c", "LMMtree_sc")] <- "LMM tree\n\n(cl.cov.)"
tmp5$method2[tmp5$method2 == "LMMtree_r"] <- "LMM tree\n\n(ran.eff.)"
tmp5$method2[tmp5$method2 %in% c("LMM_s", "LMM_s3")] <- " LMM"
tmp5$random <- ifelse(tmp5$method == "LMtree_c", "atop(hat(sigma)[b[0]] == 0, hat(sigma)[b[1]] == 0)", 
                      ifelse(tmp5$method %in% c("LMMtree", "LMMtree_c", "LMMtree_r"), 
                             "atop(hat(sigma)[b[0]] > 0, hat(sigma)[b[1]] == 0)", 
                             ifelse(tmp5$method == "LMMtree_sc", 
                                "atop(hat(sigma)[b[0]] > 0, hat(sigma)[b[1]] > 0)", 
                                ifelse(tmp5$method == "LMM_s", " time ", "all"))))
ggplot(tmp5, aes(x=random, y=MSE)) +
  geom_boxplot(position=position_dodge(1), alpha = .5, width = .45, fill = "gray") +
  scale_y_continuous(sec.axis = sec_axis(trans=~1-(./var_scie))) +
  facet_grid(~method2, scales = "free", space = "free") +
  theme(strip.background = element_blank(), strip.text.x = element_blank(),
        legend.position="none", axis.text.y = element_text(size = 9),
        axis.text.x = element_text(size = 9), axis.title.y = element_text(size = 10)) +
  labs(x = "", y = "MSE") +
  scale_x_discrete(labels = scales::parse_format()) +
  coord_cartesian(clip = 'off') 


library("kableExtra")
rownames(MSE_df) <- c("LM tree$^c$", "LMM tree$^i$", "LMM tree$^{i,c}$", 
                      "LMM tree$^{i,r}$", "LMM tree$^{i,s,c}$", "LMM$^t$", "LMM$^a$")
colnames(MSE_df) <- rep(c("$M$", "$SD$", "$R^2$"), times = 3)

kable_styling(
  add_header_above(
  add_footnote(
  row_spec(
  kable(MSE_df[c(1:3, 5, 4, 6:7), ], format = "latex", booktabs = TRUE, label = "application_MSEs", align = "c", 
        escape=FALSE, linesep = "", # linesep command suppresses addlinesep every 5 rows
        caption = "Cross-validated mean squared errors for each of the response variables."),
  5, hline_after = TRUE,), 
  "\\footnotesize \\\\ \\textit{Note.} Means and standard deviations computed over 100 cross-validation repetitions. $R^2$ was computed as $1-\\frac{\\text{mean(MSE)}}{\\text{var}(y)}$. $^c$~cluster-level covariances; $^r$~estimation initialized with random effects; $^i$~random-intercept variance freely estimated; $^s$~random-slope variance freely estimated; $^t$~LMM with fixed effect of time; $^a$~LMM with fixed effects of time, all covariates and all time-by-covariate interactions.", 
  notation="none", threeparttable = TRUE, escape = FALSE), 
  c(" ", "Math" = 3, "Reading" = 3, "Science" = 3), align = "c"),
  font_size = 11, full_width=FALSE)


vs <- cbind(sprintf("%.2f", c(apply(sapply(vs_math, \(x) apply(x, 1, sum)), 2, mean), 0, 11)),
            sprintf("%.2f", c(apply(sapply(vs_math, \(x) apply(x, 1, sum)), 2, sd), 0, 0)),
            sprintf("%.2f", c(apply(sapply(vs_read, \(x) apply(x, 1, sum)), 2, mean), 0, 11)),
            sprintf("%.2f", c(apply(sapply(vs_read, \(x) apply(x, 1, sum)), 2, sd), 0, 0)),
            sprintf("%.2f", c(apply(sapply(vs_scie, \(x) apply(x, 1, sum)), 2, mean), 0, 11)),
            sprintf("%.2f", c(apply(sapply(vs_scie, \(x) apply(x, 1, sum)), 2, sd), 0, 0)))
colnames(vs) <- rep(c("M", "SD"), times = 3)
rownames(vs) <- c("LM tree$^c$", "LMM tree$^i$", "LMM tree$^{i,c}$", 
                  "LMM tree$^{i,r}$", "LMM tree$^{i,s,c}$", "LMM$^{t}$", "LMM$^{a}$")
kable_styling(
  add_header_above(
  add_footnote(
  row_spec(
  kable(vs, format = "latex", booktabs = TRUE, label = "application_vars_selected", align = "c", 
        escape=FALSE, linesep = "", # linesep command suppresses addlinesep every 5 rows
        caption = "Mean number of covariates used in the final models."),
  5, hline_after = TRUE),
  "\\footnotesize \\\\ \\textit{Note.} Means and standard deviations computed over 100 cross-validation repetitions. $^c$~cluster-level covariances; $^r$~estimation initialized with random effects; $^i$~random-intercept variance freely estimated; $^s$~random-slope variance freely estimated; $^t$~LMM with fixed effect of time; $^a$~LMM with fixed effects of time, all covariates and all time-by-covariate interactions.", 
  notation="none", threeparttable = TRUE, escape = FALSE), 
  c(" ", "Math" = 2, "Reading" = 2, "Science" = 2), align = "c"),
  font_size = 11, full_width=FALSE)


## Tree size
load("sizesMath.Rda")
math_sizes <- sizes
#sapply(math_sizes, function(x) table(is.na(x)))
#sapply(math_sizes, function(x) tapply(x, math_sizes$N, mean, na.rm =TRUE))

load("sizesReading.Rda")
read_sizes <- sizes
#sapply(read_sizes, function(x) table(is.na(x)))
#sapply(read_sizes, function(x) tapply(x, read_sizes$N, mean, na.rm =TRUE))

load("sizesScience.Rda")
scie_sizes <- sizes
#sapply(scie_sizes, function(x) table(is.na(x)))
#sapply(scie_sizes, function(x) tapply(x, scie_sizes$N, mean, na.rm =TRUE))

## tmp6 is math
tmp6 <- data.frame(stack(math_sizes[ , 1:5]),
                    dataset_id = factor(rep(1:nrow(math_sizes), times = 5)))
names(tmp6)[1:2] <- c("splits", "method") 
tmp6$method2 <- as.character(tmp6$method)
tmp6$method2[tmp6$method2 == "LMtree_c"] <- "LM tree\n\n(cl.cov.)"
tmp6$method2[tmp6$method2 == "LMMtree"] <- "LMM tree"
tmp6$method2[tmp6$method2 %in% c("LMMtree_c", "LMMtree_sc")] <- "LMM tree\n\n(cl.cov.)"
tmp6$method2[tmp6$method2 == "LMMtree_r"] <- "LMM tree\n\n(ran.eff.)"
tmp6$random <- ifelse(tmp6$method == "LMtree_c", "atop(hat(sigma)[b[0]] == 0, hat(sigma)[b[1]] == 0)", 
                      ifelse(tmp6$method %in% c("LMMtree", "LMMtree_c", "LMMtree_r"), 
                             "atop(hat(sigma)[b[0]] > 0, hat(sigma)[b[1]] == 0)", 
                             ifelse(tmp6$method %in% c("LMMtree_sc"), 
                                "atop(hat(sigma)[b[0]] > 0, hat(sigma)[b[1]] > 0)", NA)))
ggplot(tmp6) +
  geom_boxplot(aes(x=random, y=splits), 
               position=position_dodge(1), alpha = .5, width = .45, fill = "gray") +
  facet_grid(~method2, scales = "free", space = "free") +
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(),
        axis.text.y = element_text(size = 8), axis.title.y = element_text(size = 10)) +
  labs(x = "", y = "no. of splits")


## tmp7 is reading
tmp7 <- data.frame(stack(read_sizes[ , 1:5]),
                    dataset_id = factor(rep(1:nrow(read_sizes), times = 5)))
names(tmp7)[1:2] <- c("splits", "method") 
tmp7$method2 <- as.character(tmp7$method)
tmp7$method2[tmp7$method2 == "LMtree_c"] <- "LM tree\n\n(cl.cov.)"
tmp7$method2[tmp7$method2 == "LMMtree"] <- "LMM tree"
tmp7$method2[tmp7$method2 %in% c("LMMtree_c", "LMMtree_sc")] <- "LMM tree\n\n(cl.cov.)"
tmp7$method2[tmp7$method2 == "LMMtree_r"] <- "LMM tree\n\n(ran.eff.)"
tmp7$random <- ifelse(tmp7$method == "LMtree_c", "atop(hat(sigma)[b[0]] == 0, hat(sigma)[b[1]] == 0)", 
                      ifelse(tmp7$method %in% c("LMMtree", "LMMtree_c", "LMMtree_r"), 
                             "atop(hat(sigma)[b[0]] > 0, hat(sigma)[b[1]] == 0)", 
                             ifelse(tmp7$method %in% c("LMMtree_sc"), 
                                "atop(hat(sigma)[b[0]] > 0, hat(sigma)[b[1]] > 0)", NA)))
ggplot(tmp7) +
  geom_boxplot(aes(x=random, y=splits), 
               position=position_dodge(1), alpha = .5, width = .45, fill = "gray") +
  facet_grid(~method2, scales = "free", space = "free") +
  labs(x = "", y = "no. of splits") +
  theme(strip.background = element_blank(), strip.text.x = element_blank(),
        axis.title.x=element_blank(), axis.text.x=element_blank(),
        axis.text.y = element_text(size = 8), axis.title.y = element_text(size = 10))


## tmp8 is science
tmp8 <- data.frame(stack(scie_sizes[ , 1:5]),
                    dataset_id = factor(rep(1:nrow(scie_sizes), times = 5)))
names(tmp8)[1:2] <- c("splits", "method") 
tmp8$method2 <- as.character(tmp8$method)
tmp8$method2[tmp8$method2 == "LMtree_c"] <- "LM tree\n\n(cl.cov.)"
tmp8$method2[tmp8$method2 == "LMMtree"] <- "LMM tree"
tmp8$method2[tmp8$method2 %in% c("LMMtree_c", "LMMtree_sc")] <- "LMM tree\n\n(cl.cov.)"
tmp8$method2[tmp8$method2 == "LMMtree_r"] <- "LMM tree\n\n(ran.eff.)"
tmp8$random <- ifelse(tmp8$method == "LMtree_c", "atop(hat(sigma)[b[0]] == 0, hat(sigma)[b[1]] == 0)", 
                      ifelse(tmp8$method %in% c("LMMtree", "LMMtree_c", "LMMtree_r"), 
                             "atop(hat(sigma)[b[0]] > 0, hat(sigma)[b[1]] == 0)", 
                             ifelse(tmp8$method %in% c("LMMtree_sc"), 
                                "atop(hat(sigma)[b[0]] > 0, hat(sigma)[b[1]] > 0)", NA)))
ggplot(tmp8) +
  geom_boxplot(aes(x=random, y=splits), 
               position=position_dodge(1), alpha = .5, width = .45, fill = "gray") +
  facet_grid(~method2, scales = "free", space = "free") +
  theme(strip.background = element_blank(), strip.text.x = element_blank(),
        axis.text.y = element_text(size = 8), 
        axis.text.x = element_text(size = 9), axis.title.y = element_text(size = 10)) +
  labs(x = "", y = "no. of splits") +
  scale_x_discrete(labels = scales::parse_format()) +
  coord_cartesian(clip = 'off')


theme_set(theme_bw(base_size = 8))
cols <- rep(rainbow_hcl(10)[c(3,5)], times = 5)
ggplot(tmp[LMM_ids, ], aes(x=setting, y=tree_size)) +
  geom_boxplot(aes(fill = sigma_int), 
               position=position_dodge(1), alpha = .5, width = .6, 
               outlier.shape = NA, coef = NULL) +
  geom_count(aes(group=sigma_int), position=position_dodge(1), colour = "black", 
             fill=NA, show.legend = FALSE, alpha = .15) +
  scale_y_continuous(trans = pseudo_log_trans(), breaks = breaks, lim = lim) +
  facet_grid(~ranef, scales = "free", space = "free", labeller = label_parsed) +
  labs(x = "", y = "no. of splits", col=expression(sigma^2~(b[0]))) + 
  geom_hline(yintercept=3, col = "darkgray") + 
  theme(axis.title.x=element_blank(), axis.text.x=element_blank()) +
  labs(fill=expression(sigma[b[0]]^2)) +
  stat_summary(aes(group=sigma_int), position=position_dodge(1), fun="mean", 
               col = "black", shape=16, size = .4) 


ggplot(tmp[LMM_ids, ], aes(x=setting, y=tree_size)) +
  geom_boxplot(aes(fill = N), 
               position=position_dodge(1), alpha = .5, width = .6, 
               outlier.shape = NA, coef=NULL) + 
  geom_count(aes(group=N), 
             position=position_dodge(1), colour = "black", fill=NA, show.legend = FALSE, alpha = .15) +
  scale_y_continuous(trans = pseudo_log_trans(), breaks = breaks, lim = lim) +
  facet_grid(~ranef, scales = "free", space = "free", labeller = label_parsed) +
  labs(x = "", y = "no. of splits") +
  geom_hline(yintercept=3, col = "darkgray") + 
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(),
        strip.background = element_blank(), strip.text.x = element_blank()) +
  stat_summary(aes(group=N), position=position_dodge(1), fun="mean", 
               col = "black", shape=16, size = .4) 


ggplot(tmp[LMM_ids, ], aes(x=setting, y=tree_size)) +
  geom_boxplot(aes(fill = p_noise), position=position_dodge(1), 
               alpha = .5, width = .6, outlier.shape = NA, coef=NULL) + 
  geom_count(aes(group=p_noise), position=position_dodge(1), 
             colour = "black", fill=NA, show.legend = FALSE, alpha = .15) +
  scale_y_continuous(trans = pseudo_log_trans(), breaks = breaks, lim = lim) +
  facet_grid(~ranef, scales = "free", space = "free", labeller = label_parsed) +
  labs(x = "", y = "no. of splits", col = expression(p[noise])) +
  geom_hline(yintercept=3, col = "darkgray") +
  theme(axis.title.x=element_blank(),axis.text.x=element_blank(),
        strip.background = element_blank(), strip.text.x = element_blank()) +
  labs(fill=expression(p[noise])) +
  stat_summary(aes(group=p_noise), position=position_dodge(1), fun="mean", 
               col = "black", shape=16, size = .4) 


ggplot(tmp[LMM_ids, ], aes(x=setting, y=tree_size)) +
  geom_boxplot(aes(fill = sigma_slope), position=position_dodge(1), 
               alpha = .5, width = .6, outlier.shape = NA, coef=NULL) + 
  geom_count(aes(group=sigma_slope), position=position_dodge(1), 
             colour = "black", fill=NA, show.legend = FALSE, alpha = .15) +
  scale_y_continuous(trans = pseudo_log_trans(), breaks = breaks, lim = lim) +
  facet_grid(~ranef, scales = "free", space = "free", labeller = label_parsed) +
  labs(x = "", y = "no. of splits", col=expression(sigma^2~(b[1]))) +
  geom_hline(yintercept=3, col = "darkgray") + 
  theme(axis.title.x=element_blank(),axis.text.x=element_blank(),
        strip.background = element_blank(), strip.text.x = element_blank()) +
  labs(fill=expression(sigma[b[1]]^2)) +
  stat_summary(aes(group=sigma_slope), position=position_dodge(1), fun="mean", 
               col = "black", shape=16, size = .4) 
  


ggplot(tmp[LMM_ids, ], aes(x=setting, y=tree_size)) +
  geom_boxplot(aes(fill = rho), position=position_dodge(1), 
               alpha = .5, width = .6, outlier.shape = NA, coef=NULL) + 
  geom_count(aes(group=rho), position=position_dodge(1), 
             colour = "black", fill=NA, show.legend = FALSE, alpha = .15) +
  scale_y_continuous(trans = pseudo_log_trans(), breaks = breaks, lim = lim) +
  facet_grid(~ranef, scales = "free", space = "free", labeller = label_parsed) +
  labs(x = "", y = "no. of splits", col=expression(sigma^2~(b[1]))) +
  geom_hline(yintercept=3, col = "darkgray") + 
  theme(strip.background = element_blank(), strip.text.x = element_blank()) +
  labs(fill=expression(rho)) +
  stat_summary(aes(group=rho), position=position_dodge(1), fun="mean", 
               col = "black", shape=16, size = .4) 
theme_set(theme_bw(base_size = 11))


lim <- c(0, 18)
theme_set(theme_bw(base_size = 8))
ggplot(tmp2, aes(x=random, y=tree_size)) +
  geom_boxplot(aes(fill = N), position=position_dodge(1), alpha = .5, width = .6,
               outlier.shape = NA, coef=NULL) + 
  geom_count(aes(group=N), 
             position=position_dodge(1), colour = "black", fill=NA, show.legend = FALSE, alpha = .15) +
  scale_y_continuous(trans = pseudo_log_trans(), breaks = breaks, lim = lim) +
  facet_grid(~method, scales = "free", space = "free") +
  labs(x = "", y = "no. of splits", col=expression(N)) + 
  geom_hline(yintercept=3, col = "darkgray") + 
  theme(axis.title.x=element_blank(), axis.text.x=element_blank()) +
  labs(fill=expression(N)) +
  stat_summary(aes(group=N), position=position_dodge(1), fun="mean", 
               col = "black", shape=16, size = .4) 


ggplot(tmp2, aes(x=random, y=tree_size)) +
  geom_boxplot(aes(fill = sigma_int), position=position_dodge(1), 
               alpha = .5, width = .6, outlier.shape = NA, coef=NULL) + 
  geom_count(aes(group=sigma_int), position=position_dodge(1), 
             colour = "black", fill=NA, show.legend = FALSE, alpha = .15) +
  scale_y_continuous(trans = pseudo_log_trans(), 
                     breaks = breaks, limits = lim) +
  facet_grid(~method, scales = "free", space = "free") +
  labs(x = "", y = "no. of splits", col=expression(sigma^2~(b[0]))) +
  geom_hline(yintercept=3, col = "darkgray") + 
  theme(axis.title.x=element_blank(),axis.text.x=element_blank(),
        strip.background = element_blank(), strip.text.x = element_blank()) +
  labs(fill=expression(sigma[b[0]]^2)) +
  stat_summary(aes(group=sigma_int), position=position_dodge(1), fun="mean", 
               col = "black", shape=16, size = .4)


ggplot(tmp2, aes(x=random, y=tree_size)) +
  geom_boxplot(aes(fill = p_noise), position=position_dodge(1), 
               alpha = .5, width = .6, outlier.shape = NA, coef=NULL) + 
  geom_count(aes(group=p_noise), position=position_dodge(1), 
             colour = "black", fill=NA, show.legend = FALSE, alpha = .15) +
  facet_grid(~method, scales = "free", space = "free") +
  labs(x = "", y = "no. of splits", col = expression(p[noise])) +
  geom_hline(yintercept=3, col = "darkgray") +
  theme(axis.title.x=element_blank(),axis.text.x=element_blank(),
        strip.background = element_blank(), strip.text.x = element_blank()) +
  labs(fill=expression(p[noise])) +
  scale_y_continuous(trans = pseudo_log_trans(), breaks = breaks, lim = lim) +
  stat_summary(aes(group=p_noise), position=position_dodge(1), fun="mean", 
               col = "black", shape=16, size = .4)


ggplot(tmp2, aes(x=random, y=tree_size)) +
  geom_boxplot(aes(fill = rho), position=position_dodge(1), 
               alpha = .5, width = .6, outlier.shape = NA, coef=NULL) + 
  geom_count(aes(group=rho), position=position_dodge(1), colour = "black", 
             fill=NA, show.legend = FALSE, alpha = .15) +
  scale_y_continuous(trans = pseudo_log_trans(), breaks = breaks, lim = lim) +
  facet_grid(~method, scales = "free", space = "free") +
  labs(x = "", y = "no. of splits", col=expression(rho)) +
  geom_hline(yintercept=3, col = "darkgray") + 
  theme(axis.title.x=element_blank(),axis.text.x=element_blank(),
        strip.background = element_blank(), strip.text.x = element_blank()) +
  labs(fill=expression(rho)) + 
  stat_summary(aes(group=rho), position=position_dodge(1), fun="mean", 
               col = "black", shape=16, size = .4)


ggplot(tmp2, aes(x=random, y=tree_size)) +
  geom_boxplot(aes(fill = sigma_slope), 
               position=position_dodge(1), alpha = .5, width = .6,
               outlier.shape = NA, coef=NULL) + 
  geom_count(aes(group=sigma_slope), 
             position=position_dodge(1), colour = "black", fill=NA, show.legend = FALSE, alpha = .15) +
  scale_y_continuous(trans = pseudo_log_trans(), 
                     breaks = breaks, lim = lim) +
  facet_grid(~method, scales = "free", space = "free") +
  labs(x = "", y = "no. of splits", col=expression(sigma[b[1]]^2)) +
  geom_hline(yintercept=3, col = "darkgray") +
  theme(strip.background = element_blank(), strip.text.x = element_blank()) +
  labs(fill=expression(sigma[b[1]]^2)) + 
  scale_x_discrete(labels = scales::parse_format()) +
  stat_summary(aes(group=sigma_slope), position=position_dodge(1), fun="mean", 
               col = "black", shape=16, size = .4)
theme_set(theme_bw(base_size = 11))


resp <- "Reading"
load(paste(resp, "ability data.Rdata"))  
tmp_data <- readdata
labels <- c("kindergarten", "1st grade", "3rd grade", "5th grade", "8th grade")
tmp_data$asmtmm <- factor(tmp_data$asmtmm, labels = labels)
cols <- colorspace::rainbow_hcl(5)
ggplot(tmp_data, aes(x = score, colour = asmtmm, fill = asmtmm)) + 
  geom_density(linewidth = .75, bw = .1, alpha = .1) +
  # geom_histogram(alpha = .2) +
  labs(fill=NULL, color = NULL, x = paste(resp, "ability score")) +
  scale_colour_manual(values = cols) + 
  scale_fill_manual(values = cols)


resp <- "Math"
load(paste(resp, "ability data.Rdata"))  
tmp_data <- mathdata
labels <- c("kindergarten", "1st grade", "3rd grade", "5th grade", "8th grade")
tmp_data$asmtmm <- factor(tmp_data$asmtmm, labels = labels)
cols <- colorspace::rainbow_hcl(5)
ggplot(tmp_data, aes(x = score, colour = asmtmm, fill = asmtmm)) + 
        geom_density(linewidth = .75, bw = .1, alpha = .1) +
      # geom_histogram(alpha = .2) +
        labs(fill=NULL, color = NULL, x = paste(resp, "ability score")) +
        scale_colour_manual(values = cols) + 
        scale_fill_manual(values = cols)


resp <- "Science"
load(paste(resp, "ability data.Rdata"))  
tmp_data <- sciedata
labels <- c("3rd grade", "5th grade", "8th grade")
tmp_data$asmtmm <- factor(tmp_data$asmtmm, labels = labels)
cols <- colorspace::rainbow_hcl(5)[3:5]
ggplot(tmp_data, aes(x = score, colour = asmtmm, fill = asmtmm)) + 
  geom_density(linewidth = .75, bw = .1, alpha = .1) +
  # geom_histogram(alpha = .2) +
  labs(fill=NULL, color = NULL, x = paste(resp, "ability score")) +
  scale_colour_manual(values = cols) + 
  scale_fill_manual(values = cols)


## Create plots illustrating scaling of time
par(mfrow = c(3, 2))
for (resp in c("Reading", "Math", "Science")) {
  
  ## Load data and apply transformations if necessary
  load(paste(resp, "ability data.Rdata"))
  tmp_data <- get(ifelse(resp == "Reading", "readdata", 
                         ifelse(resp == "Math", "mathdata", "sciedata")))
  frac <- ifelse(resp == "Science", 2/3, 1/2)
  tmp_data$months <- tmp_data$months^(1/frac) 
  xlab <- switch(resp, Reading = c("time", expression(time^{1/2})), 
                 Math = c("time", expression(time^{1/2})),
                 Science = c("time", expression(time^{2/3})))
  
  ## Create plots
  for (i in 1:2) {
    y <- t(sapply(unique(tmp_data$CHILDID), function(x) 
      tmp_data$score[tmp_data$CHILDID == x]))
    x <- t(sapply(unique(tmp_data$CHILDID), function(x) 
      tmp_data$months[tmp_data$CHILDID == x]))
    plot(tmp_data$month, tmp_data$score, type = "n", main = resp, xlab = xlab[i], 
         ylab = "Score", cex.axis = .7, cex.lab = .7, cex.main = .7)
    sapply(1:nrow(x), function(row) 
      lines(x[row, ], y[row, ], col = gray(0.5, alpha = 1/20), lwd = .05))
    tmp_data$months <- tmp_data$months^if (i == 1L) frac else (1/frac)
  }
}