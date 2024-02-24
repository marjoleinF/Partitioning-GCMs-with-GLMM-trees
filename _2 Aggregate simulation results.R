########################
##
## Get tree sizes
##
##

LMM_sizes <- SEM_sizes <- Lcrt_sizes <- list()

for (i in 1:100) {
  load(paste0("Lcrt_sizes", i))
  Lcrt_sizes[[i]] <- unlist(tree_sizes) 
  load(paste0("LMM_sizes", i))
  LMM_sizes[[i]] <- tree_sizes   
  load(paste0("SEM_sizes", i))
  SEM_sizes[[i]] <- tree_sizes 
}

Lcrt_sizes <- unlist(Lcrt_sizes)
dim(Lcrt_sizes) ## should be a vector of nrep x 32
length(Lcrt_sizes)
head(Lcrt_sizes)

SEM_sizes <- do.call(rbind, SEM_sizes)
dim(SEM_sizes)
head(SEM_sizes)

LMM_sizes <- do.call(rbind, LMM_sizes)
dim(LMM_sizes)
head(LMM_sizes)

sizes <- cbind(LMM_sizes, SEM_sizes, Lcrt = Lcrt_sizes)
sizes <- apply(sizes, 2, unlist)
dim(sizes)
boxplot(sizes)
table(sizes)
save(sizes, file = "tree_sizes.Rda")







################################
##
## Get splitting variables
##
##

LMM_splits <- SEM_splits <- Lcrt_splits <- list()

for (i in 1:100) {
  load(paste0("Lcrt_splitvars", i))
  Lcrt_splits[[i]] <- unlist(splitvars) 
  load(paste0("LMM_splitvars", i))
  LMM_splits[[i]] <- splitvars   
  load(paste0("SEM_splitvars", i))
  SEM_splits[[i]] <- splitvars 
}

Lcrt_splits <- unlist(Lcrt_splits)
dim(Lcrt_splits) ## should be a vector of nrep x 32
length(Lcrt_splits)
head(Lcrt_splits)

SEM_splits <- do.call(rbind, SEM_splits)
dim(SEM_splits)
head(SEM_splits)
colnames(SEM_splits) <- colnames(SEM_sizes)

LMM_splits <- do.call(rbind, LMM_splits)
dim(LMM_splits)
head(LMM_splits)
colnames(LMM_splits) <- colnames(LMM_sizes)

splits <- cbind(LMM_splits, SEM_splits, Lcrt = Lcrt_splits)
splits <- apply(splits, 2, unlist)
dim(splits)
apply(splits, 2, table, useNA = "ifany")

save(splits, file = "splitvars.Rda")



###################################
##
## Check splits against tree size
##
##
apply(splits, 2, function(x) table(is.na(x)))
apply(sizes, 2, function(x) table(x == 1))
## NOTE: With LongCART, if first split is not implemented, tree$Treeout returns 
##   stats on the best splitting variable anyway. So the name of a splitting
##   variable will always be returned, but split may not have been implemented.





###########################################
##
## Get computation times
##
##

LMM_times <- SEM_times <- Lcrt_times <- list()

for (i in 1:100) {
  load(paste0("Lcrt_times", i))
  Lcrt_times[[i]] <- unlist(times) 
  load(paste0("LMM_times", i))
  LMM_times[[i]] <- times   
  load(paste0("SEM_times", i))
  SEM_times[[i]] <- times 
}

Lcrt_times <- unlist(Lcrt_times)
dim(Lcrt_times) ## should be nrep x 32
length(Lcrt_times)

SEM_times <- do.call(rbind, SEM_times)
dim(SEM_times)

LMM_times <- do.call(rbind, LMM_times)
dim(LMM_times)

times <- cbind(LMM_times, SEM_times, Lcrt = Lcrt_times)
times <- apply(times, 2, unlist)
head(times)
dim(times)
boxplot(log(times))

save(times, file = "comp_times.Rda")