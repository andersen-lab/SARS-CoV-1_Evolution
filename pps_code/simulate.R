library(phangorn)
set.seed(42)
source("~/sars-1/human_2030/piBUSS/src/alignment_helpers.R")

parent_dir <- "~/sars-1/human_2030/piBUSS"

pibuss_script_refx1 <- scan(paste0(parent_dir,"/src/piBUSS_randomEffects.xml"),what=character(),sep="\n",strip.white=FALSE)

#######
# log files
#######

# This section will need to be amended if runs aren't stored in parent_dir/output/job_<index>

nrun <- 1

burnin_frac <- 0.25

# Loops over multiple runs, reads in logs and removes burnin
param_logs <- vector("list",nrun)
tree_logs <- vector("list",nrun)
for (i in 1:nrun) {
  param_logs[[i]] <- read.table("~/sars-1/human_2030/piBUSS/data/Human_2003.log",header=TRUE,stringsAsFactors=FALSE)
  tree_logs[[i]] <- read.nexus("~/sars-1/human_2030/piBUSS/data/Human_2003.trees")
  
  nlog <- dim(param_logs[[i]])[1]
  keep <- (round(burnin_frac*nlog)+1):nlog
  param_logs[[i]] <- param_logs[[i]][keep,]
  tree_logs[[i]] <- tree_logs[[i]][keep]
}

# If we want to thin, now is the time, otherwise we'll simulate one alignment per sample
# And if tree sampling frequency does not match parameter frequency, we need to align those
param_log <- do.call(rbind,param_logs)
trees <- do.call(c,tree_logs)

nsim <- length(trees)

#nsim<-1


######
# Real data (for ambiguities and such)
######
real_aln <- read.phyDat("~/sars-1/human_2030/piBUSS/data/Human_2003.fasta","fasta")
real_aln_taxa <- names(real_aln)
real_aln_ambiguities <- getAmbiguities(real_aln)
nsites <- dim(as.character(real_aln))[2]

aln_taxa <- names(real_aln)

target_dir <- "SC1_human_only_sims"

dir.create(paste0(parent_dir,"/",target_dir),showWarnings=FALSE)

######
# Simulate
######

for (i in 1:nsim) {
  suppressWarnings(rm(this_script))
  
  this_script <- pibuss_script_refx1

  # Tree
  this_script <- gsub("NEWICKTREE",write.tree(trees[[i]]),this_script)
  
  # Clock
  pinv <- param_log$pInv[i]
  clock.rate <- param_log$clock.rate[i]
  clock1 <- clock.rate/(1 - pinv)

  this_script <- gsub("CLOCKRATE",clock1,this_script)
  
  # Q matrix
  log_kappa <- param_log$log.kappa[i]
  this_script <- gsub("LOGKAPPA",log_kappa,this_script)
  refx <- as.numeric(param_log[i,grepl("glmRandCoefficients",names(param_log))])
  this_script <- gsub("LOGRANDOMEFFECTS",paste0(refx,collapse=" "),this_script)
  bf <- c(param_log$frequencies1[i],param_log$frequencies2[i],param_log$frequencies3[i],param_log$frequencies4[i])
  bf <- bf/sum(bf)
  this_script <- gsub("FREQUENCIES",paste0(sprintf("%.16f",bf),collapse=" "),this_script)

  # Site counts in variable and invariant partitions
  nsites_var <- (1 - pinv) * nsites
  remainder <- nsites_var %% 1
  nsites_var <- nsites_var - remainder + rbinom(1,1,remainder)
  nsites_inv <- nsites - nsites_var
  
  this_script <- gsub("FROMPARTITION1",1,this_script)
  this_script <- gsub("TOPARTITION1",nsites_var,this_script)
  
  this_script <- gsub("FROMPARTITION2",nsites_var+1,this_script)
  this_script <- gsub("TOPARTITION2",nsites,this_script)
  
  
  tmpfile <- paste0(parent_dir,"/",target_dir,"/tmp_pibuss.xml")
  
  beastjar <- "~/beast-mcmc-hmc-clock/build/dist/beast.jar"
  
  beaglepath <- "/usr/local/lib"
  tmpaln <- paste0(dirname(tmpfile),"/sequences.fasta")
  
  # DO THE EVOLUTION, BABY!
  cat(this_script,file=tmpfile,sep="\n")
  system(paste0("java  -Djava.library.path=",beaglepath," -jar ",beastjar," -working -overwrite ",tmpfile))
  sim_aln <- read.phyDat(tmpaln,"fasta")



  sim_aln <- addAmbiguities(sim_aln,real_aln_ambiguities)
  
  write.phyDat(sim_aln,paste0(parent_dir,"/",target_dir,"/sim_",i,".fasta"),format="fasta")
  # Extra files but easier for parsimony counting
  write.tree(trees[[i]],file=paste0(parent_dir,"/",target_dir,"/sim_",i,".tre"))
  
  rm(tmpaln)
  rm(tmpfile)
}
