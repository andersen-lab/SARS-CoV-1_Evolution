library(phangorn)
source("~/Documents/SC1/piBUSS/src/test_stats.R")

parent_dir <- "~/Documents/SC1/piBUSS"
target_dir <- "SC1_142_sims"

######
# Real data (for ambiguities and such)
######
real_aln <- read.phyDat("~/Documents/SC1/piBUSS/data/HuCi_spike_fullset.fasta","fasta")
real_aln_taxa <- names(real_aln)
nsites <- dim(as.character(real_aln))[2]


aln_taxa <- names(real_aln)

sim_files <- list.files(paste0(parent_dir,"/",target_dir),full.names=TRUE)
sim_aln <- sim_files[grepl("fasta",sim_files)]
sim_aln <- sim_aln[!grepl("sequences.fasta",sim_aln,fixed=TRUE)]
nsim <- length(sim_aln)

sim_stats <- parallel::mclapply(1:nsim,function(i) {
			        sim_aln <- read.phyDat(paste0(parent_dir,"/",target_dir,"/sim_",i,".fasta"),"fasta")
			        sim_tree <- read.tree(paste0(parent_dir,"/",target_dir,"/sim_",i,".tre"))
			        return(computeStats(sim_aln,sim_tree))
},mc.preschedule=TRUE,mc.cores=10)
sim_stats <- do.call(rbind,sim_stats)

obs_stats <- parallel::mclapply(1:nsim,function(i) {
  sim_tree <- read.tree(paste0(parent_dir,"/",target_dir,"/sim_",i,".tre"))
  return(computeStats(real_aln,sim_tree))
},mc.preschedule=TRUE,mc.cores=10)
obs_stats <- do.call(rbind,obs_stats)


# These take awhile to compute, so we store them
save(obs_stats,sim_stats,file=paste0(parent_dir,"/",target_dir,"/pps_vals.Rdata"),version=2)

# # pdf("~/Downloads/SC1_human_only_posterior_predictive_analyses.pdf")
# #   for (i in 1:dim(obs_stats)[2]) {
# #   # for (i in 30:35) {
# #     toplot <- cbind(obs_stats[,i],sim_stats[,i])
# #     toplot <- toplot[!apply(toplot,1,function(x){any(is.na(x))}),]
# #     if (dim(toplot)[1] == 0) {
# #       frame()
# #       text(x=0.5,y=0.5,labels=paste0("no non-NA values for ",colnames(obs_stats)[i]))
# #     } else {
# #       plot(toplot,pch=16,col="#66666690",xlab="real value",ylab="simulated value",main=colnames(obs_stats)[i])
# #       abline(a=0,b=1,col="red")
# #       pppval <- sum(toplot[,1] < toplot[,2])/dim(toplot)[1] + 0.5 * sum(toplot[,1] == toplot[,2])/dim(toplot)[1]
# #       legend("topleft",legend=paste0("PPP-value = ",round(pppval,3)),bty="n",border=NA)
# #     }
# #   }
# # dev.off()
# 
# pdf("~/Downloads/SC1_142_posterior_predictive_analyses_take_2.pdf")
# for (i in 1:dim(obs_stats)[2]) {
#   # for (i in 30:35) {
#   toplot <- cbind(obs_stats[,i],sim_stats[,i])
#   toplot <- toplot[!apply(toplot,1,function(x){any(is.na(x))}),]
#   if (dim(toplot)[1] == 0) {
#     frame()
#     text(x=0.5,y=0.5,labels=paste0("no non-NA values for ",colnames(obs_stats)[i]))
#   } else {
#     hist(toplot[,2] - toplot[,1],breaks=50,col="#66666690",border=NA,freq=FALSE,xlab="real - simulated",main=colnames(obs_stats)[i])
#     # plot(toplot,pch=16,col="#66666690",xlab="real value",ylab="simulated value",main=colnames(obs_stats)[i])
#     abline(v=0,lwd=2,lty=2)
#     pppval <- sum(toplot[,1] < toplot[,2])/dim(toplot)[1] + 0.5 * sum(toplot[,1] == toplot[,2])/dim(toplot)[1]
#     legend("topleft",legend=paste0("PPP-value = ",round(pppval,3)),bty="n",border=NA)
#   }
# }
# dev.off()


# load("~/Downloads/piBUSS/SC1_human_only_sims/pps_vals.Rdata")
# 
# sum(obs_stats[,1] < sim_stats[,1])/dim(sim_stats)[1]
# 
#   
# 
