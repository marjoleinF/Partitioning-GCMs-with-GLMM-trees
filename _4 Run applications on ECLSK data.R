## Load libraries & save session info
library("glmertree")
library("partykit")
sInfo <- sessionInfo()
save(sInfo, file = "sessionInfo.Rda")


#################################################
##
## Cross-validate performance of LM(M) trees
##

## Set up design
methods <- c("LMtree_c", "LMMtree", "LMMtree_c", "LMMtree_r", "LMMtree_sc", 
             "LMM", "LMM_s", "LMM2", "LMM_s2", "LMM3", "LMM_s3")
nreps <- 100L
part_vars <- c("months", "GENDER", "RACE", "WKSESL", "C1GMOTOR", "C1FMOTOR", 
               "T1INTERN", "T1EXTERN", "T1INTERP", "T1CONTRO", "P1FIRKDG", "AGEBASELINE")

## Prepare formulas
LMt_form <- score ~ months | GENDER + RACE + WKSESL + C1GMOTOR + 
  C1FMOTOR + T1INTERN + T1EXTERN + T1INTERP + T1CONTRO + P1FIRKDG + AGEBASELINE

LMM_form <- score ~ (1|CHILDID) + months
LMM_s_form <- score ~ (1 + months|CHILDID) + months

LMM_form2 <- score ~ (1|CHILDID) + months + GENDER + RACE + WKSESL + C1GMOTOR + 
  C1FMOTOR + T1INTERN + T1EXTERN + T1INTERP + T1CONTRO + P1FIRKDG + AGEBASELINE
LMM_s_form2 <- score ~ (1 + months|CHILDID) + months + GENDER + RACE + WKSESL + 
  C1GMOTOR + C1FMOTOR + T1INTERN + T1EXTERN + T1INTERP + T1CONTRO + P1FIRKDG + AGEBASELINE

LMM_form3 <- score ~ (1|CHILDID) + months*(GENDER + RACE + WKSESL + C1GMOTOR + 
  C1FMOTOR + T1INTERN + T1EXTERN + T1INTERP + T1CONTRO + P1FIRKDG + AGEBASELINE)
LMM_s_form3 <- score ~ (1 + months|CHILDID) + months*(GENDER + RACE + WKSESL + 
  C1GMOTOR + C1FMOTOR + T1INTERN + T1EXTERN + T1INTERP + T1CONTRO + P1FIRKDG + AGEBASELINE)

LMMt_form <- score ~ months | CHILDID | GENDER + RACE + WKSESL + C1GMOTOR + 
  C1FMOTOR + T1INTERN + T1EXTERN + T1INTERP + T1CONTRO + P1FIRKDG + AGEBASELINE
LMMt_s_form <- score ~ months | ((1 + months) | CHILDID) | GENDER + RACE + WKSESL + 
  C1GMOTOR + C1FMOTOR + T1INTERN + T1EXTERN + T1INTERP + T1CONTRO + P1FIRKDG + 
  AGEBASELINE

set.seed(42)
for (resp in c("Math", "Reading", "Science")) {
  
  ## Read in data and prepare objects for saving results
  load(paste0(resp, " ability data.Rdata"))
  resp_short <- ifelse(resp == "Math", "math", ifelse(resp == "Reading", "read", "scie"))
  data <- get(paste0(resp_short, "data"))

  data$CHILDID <- factor(data$CHILDID)
  data$GENDER <- factor(data$GENDER)
  data$RACE <- factor(data$RACE)
  data$asmtmm <- factor(data$asmtmm)
  notes <- MSEs <- times <- sizes <- data.frame(matrix(
    NA, nrow = nreps, ncol = length(methods), dimnames = list(1:nreps, methods)))
  
  ## Keep track of selected variables
  ## For every method, create a data.frame with a column per part_var and row per repetition
  vars_selected <- list()
  for (method in methods[1:5]) {
    vars_selected[[method]] <- matrix(NA, nrow = 100, ncol = length(part_vars),
                                      dimnames = list(NULL, part_vars))
  }

  ## Generate samples for each repetition
  bag_ids <- sapply(1:nreps, function(x) sample(unique(data$CHILDID), size = 250))
  save(bag_ids, file = paste0("bag_ids ", resp, ".Rda"))
    
  for (i in 1:nreps) {
      
    print(paste0("Fold ", i, " for response: ", resp, "."))
    
    ## Set up train and test data  
    train <- data[data$CHILDID %in% bag_ids[ , i], ]
    test <- data[!data$CHILDID %in% bag_ids[ , i], ]
    ## Remove obs with levels of race in test that are not in train
    levs <- unique(test$RACE) %in% unique(train$RACE)
    if (!all(levs)) {
      test <- test[-which(!test$RACE %in% unique(test$RACE)[levs]), ]
      notes[i, ] <- paste0(notes[i, ], 
                           "Levels or race omitted from test data: ", 
                           unique(test$race)[!levs])
    }
    
    
    ## Fit LM tree(s)
    times$LMtree_c[i] <- system.time(
      tree <- lmtree(LMt_form, data = train, cluster = CHILDID, ))["elapsed"]
    preds <- predict(tree, newdata = test)
    MSEs$LMtree_c[i] <- mean((test$score - preds)^2) 
    sizes$LMtree_c[i] <- (length(tree)-1)/2
    vars_selected$LMtree_c[i, ] <- sapply(part_vars, \(x) any(grepl(x, pre:::list.rules(tree))))
    
    ## Fit LMM trees
    times$LMMtree[i] <- system.time(
      tree <- lmertree(LMMt_form, data = train))["elapsed"]
    preds <- predict(tree, newdata = test, re.form = NA)
    MSEs$LMMtree[i] <- mean((test$score - preds)^2) 
    sizes$LMMtree[i] <- (length(tree$tree)-1)/2
    vars_selected$LMMtree[i, ] <- sapply(part_vars, \(x) any(grepl(x, pre:::list.rules(tree$tree))))
                                          
    times$LMMtree_c[i] <- system.time(
      tree <- lmertree(LMMt_form, data = train, cluster = CHILDID, ))["elapsed"]
    preds <- predict(tree, newdata = test, re.form = NA)
    MSEs$LMMtree_c[i] <- mean((test$score - preds)^2) 
    sizes$LMMtree_c[i] <- (length(tree$tree)-1)/2
    vars_selected$LMMtree_c[i, ] <- sapply(part_vars, \(x) any(grepl(x, pre:::list.rules(tree$tree))))
    
    times$LMMtree_r[i] <- system.time(
      tree <- lmertree(LMMt_form, data = train, ranefstart = TRUE))["elapsed"]
    preds <- predict(tree, newdata = test, re.form = NA)
    MSEs$LMMtree_r[i] <- mean((test$score - preds)^2) 
    sizes$LMMtree_r[i] <- (length(tree$tree)-1)/2
    vars_selected$LMMtree_r[i, ] <- sapply(part_vars, \(x) any(grepl(x, pre:::list.rules(tree$tree))))
    
    times$LMMtree_sc[i] <- system.time(
      tree <- lmertree(LMMt_s_form, data = train, cluster = CHILDID))["elapsed"]
    preds <- predict(tree, newdata = test, re.form = NA)
    MSEs$LMMtree_sc[i] <- mean((test$score - preds)^2) 
    sizes$LMMtree_sc[i] <- (length(tree$tree)-1)/2
    vars_selected$LMMtree_sc[i, ] <- sapply(part_vars, \(x) any(grepl(x, pre:::list.rules(tree$tree))))
    
    
    ## Fit linear mixed models
    times$LMM[i] <- system.time(lmm <- lmer(LMM_form, data = train, ))["elapsed"]
    preds <- predict(lmm, newdata = test, re.form = NA)
    MSEs$LMM[i] <- mean((test$score - preds)^2) 
    
    times$LMM_s[i] <- system.time(lmm <- lmer(LMM_s_form, data = train, ))["elapsed"]
    preds <- predict(lmm, newdata = test, re.form = NA)
    MSEs$LMM_s[i] <- mean((test$score - preds)^2) 
    
    times$LMM2[i] <- system.time(lmm <- lmer(LMM_form2, data = train, ))["elapsed"]
    preds <- predict(lmm, newdata = test, re.form = NA)
    MSEs$LMM2[i] <- mean((test$score - preds)^2) 
    
    times$LMM_s2[i] <- system.time(lmm <- lmer(LMM_s_form2, data = train, ))["elapsed"]
    preds <- predict(lmm, newdata = test, re.form = NA)
    MSEs$LMM_s2[i] <- mean((test$score - preds)^2) 
    
    times$LMM3[i] <- system.time(lmm <- lmer(LMM_form3, data = train, ))["elapsed"]
    preds <- predict(lmm, newdata = test, re.form = NA)
    MSEs$LMM3[i] <- mean((test$score - preds)^2) 
    
    times$LMM_s3[i] <- system.time(lmm <- lmer(LMM_s_form3, data = train, ))["elapsed"]
    preds <- predict(lmm, newdata = test, re.form = NA)
    MSEs$LMM_s3[i] <- mean((test$score - preds)^2) 
    
    
    ## Save results (overwrites in each replication to prevent loss of results)
    save(times, file = paste0("times", resp, ".Rda"))
    save(MSEs, file = paste0("MSEs", resp, ".Rda"))
    save(sizes, file = paste0("sizes", resp, ".Rda"))
    save(notes, file = paste0("notes", resp, ".Rda"))
    save(vars_selected, file = paste0("vars_selected", resp, ".Rda"))

  }
}

