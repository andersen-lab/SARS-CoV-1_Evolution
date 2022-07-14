library(phangorn)

# Uses parsimony to count changes between seq s1 and seq s2
pairwiseSubstitutionCounts <- function(s1,s2,reversible=FALSE) {
  nucs <- c("a","c","g","t")
  
  # we are assuming that both sequences have all 4 nucleotides
  m <- table(s1,s2)
  from <- row.names(m)
  to <- colnames(m)
  
  k <- 0
  if (reversible) {
    add <- 0
    counts <- rep(0,6)
  } else {
    add <- 6
    counts <- rep(0,12)
  }
  for (i in 1:3) {
    for (j in (i+1):4) {
      k <- k + 1
      
      # print(paste("order", nucs[i], nucs[j], sep=" "))
      
      # ignore ambiguities
      if (nucs[i] %in% from && nucs[j] %in% to) {
        counts[k] <- counts[k] + m[nucs[i],nucs[j]]
      }
      if (nucs[i] %in% to && nucs[j] %in% from) {
        counts[k+add] <- counts[k+add] + m[nucs[j],nucs[i]]
      }
    }
  }
  return(counts)
}

# Uses parsimony to get an estimate of the number of substitutions of each type on each branch
# If reversible == TRUE, the types are A<->C, ..., G<->T
# If reversible == FALSE, the types are A->C,...,G->T,C->A,...,T->A
computeParsimonySubstitutionTable <- function(aln,phy,reversible=FALSE) {
  # recover()
  # parsimony reconstruction
  asr <- ancestral.pars(phy,aln,return="phyDat")
  char_aln <- as.character(asr)
  
  # Make sure we can map sequences to nodes
  which_node <- sapply(names(asr),function(x){
    asnum <- suppressWarnings(as.numeric(x))
    if ( !is.na(asnum) ) {
      return(asnum)
    } else {
      return(which(phy$tip.label == x))
    }
  })
  
  # loop over edges
  counts <- sapply(1:Nedge(phy),function(i){
    parent <- phy$edge[i,1]
    child <- phy$edge[i,2]
    pairwiseSubstitutionCounts(char_aln[parent,],char_aln[child,],reversible=reversible)
  })
  return(counts)
}

# option "use" can be either pvalue or statistic
#              "statistic" means we return the chi-squared statistic, sum((obs-exp)^2/exp)
#              "p-value" means we look up the result in a chi-squared distribution and return the (log) p-value
#              "standardized" means we return (statistic-df)/sqrt(2*df)
# "p-value" allows us to ignore substitution types which are not observed while keeping results comparable across simulations
parsSubstitutionChiSquared <- function(aln,phy,subsTable=NULL,reversible=FALSE,use="statistic") {
  # recover()
  if ( !is.null(subsTable) ) {
    counts <- subsTable
    if (reversible == FALSE && dim(subsTable)[1] == 6) {
      stop("Cannot use non-reversible counting when providing substitutions counted in reversible form.")
    } else if (reversible == TRUE && dim(subsTable) == 12 ) {
      counts[1:6,] <- counts[1:6,] + counts[7:12,]
      counts <- counts[1:6,]
    }
  } else {
    counts <- computeParsimonySubstitutionTable(aln,phy,reversible=reversible)
  }
  if ( grepl("val",tolower(use)) ) {
    use <- "pval"
  } else if (grepl("stan",tolower(use))) {
    use <- "standardize"
  } else {
    use <- "raw"
  }
  
  has_mut <- colSums(counts) > 0
  counts <- counts[,has_mut]
  
  if ( use != "raw" ) {counts <- counts[rowSums(counts) > 0,]}
  
  f_bar <- rowSums(counts)
  f_bar <- f_bar/sum(f_bar)
  summand <- apply(counts,2,function(f_e){
    f_e <- f_e/sum(f_e)
    return(sum( ((f_e - f_bar)^2)/f_bar ))
  })
  
  res <- sum(summand)
  
  if (use == "pval") {
    res <- pchisq(res,df=prod(dim(counts)-1),lower.tail=FALSE,log=TRUE)
  } else if (use == "standardize") {
    df <- prod(dim(counts)-1)
    res <- (res - df)/sqrt(2*df)
  }
  return(res)
}

# get branch-to-branch distances for all branches in the tree
# records distance for the node that each branch subtends
getBranchToBranchMidpointDistances <- function(phy) {
  # recover()
  if ( !is.rooted(phy) ) {
    stop("Tree must be rooted!")
  }
  ntip <- Ntip(phy)
  nnode <- Nnode(phy) + ntip
  node_dists <- dist.nodes(phy)
  branch_dists <- matrix(0,nnode,nnode)
  node_mrca <- mrca(phy,full=TRUE)
  root <- Ntip(phy) + 1
  root_len <- ifelse(is.null(phy$root.edge),NA,phy$root.edge)
  for (i in 1:(nnode-1)) {
    for (j in (i+1):nnode) {
      adj <- 0
      if (i == root) {
        adj <- adj + 0.5*root_len - 0.5*phy$edge.length[phy$edge[,2] == j]
      } else if (j == root) {
        adj <- adj + 0.5*root_len - 0.5*phy$edge.length[phy$edge[,2] == i]
      } else {
        if ( node_mrca[i,j] == i ) {
          adj <- adj + 0.5*phy$edge.length[phy$edge[,2] == i] - 0.5*phy$edge.length[phy$edge[,2] == j]
        } else if ( node_mrca[i,j] == j ) {
          adj <- adj + 0.5*phy$edge.length[phy$edge[,2] == j] - 0.5*phy$edge.length[phy$edge[,2] == i]
        } else {
          adj <- adj - 0.5*phy$edge.length[phy$edge[,2] == j] - 0.5*phy$edge.length[phy$edge[,2] == i]
        }
      }
      branch_dists[i,j] <- branch_dists[j,i] <- node_dists[i,j] + adj
    }
  }
  return(branch_dists)
}

# For each substitution type, maps all branches with that type, then examines the distances between those branches
# summary: a function to summarize the distances between branches, suggest min or mean
substitutionTypeColocolization <- function(aln,phy,subsTable=NULL,reversible=FALSE,summary.function=mean) {
  if ( !is.null(subsTable) ) {
    counts <- subsTable
    if (reversible == FALSE && dim(subsTable)[1] == 6) {
      stop("Cannot use non-reversible counting when providing substitutions counted in reversible form.")
    } else if (reversible == TRUE && dim(subsTable) == 12 ) {
      counts[1:6,] <- counts[1:6,] + counts[7:12,]
      counts <- counts[1:6,]
    }
  } else {
    counts <- computeParsimonySubstitutionTable(aln,phy,reversible=reversible)
  }
  # recover()
  
  nbranches_with_subs <- rowSums(counts > 0)
  branch_has_subs <- colSums(counts) > 0
  treelen <- sum(phy$edge.length)
  
  # Get distances between midpoints of all branches
  branch_dists <- getBranchToBranchMidpointDistances(phy)
  
  # For mapping the per-branch substitution counting to the per-node branch-to-branch distances
  branch_nodes <- phy$edge[,2]
  
  nedge <- dim(counts)[2]
  
  stats <- sapply(1:dim(counts)[1],function(stype){
    res <- 0
    if (nbranches_with_subs[stype] < 2) {
      res <- NA
    } else {
      branches_with_subs_type <- which(counts[stype,] > 0)
      nodes_with_subs_type <- branch_nodes[branches_with_subs_type]
      dists_between_subs <- branch_dists[nodes_with_subs_type,nodes_with_subs_type]
      res <- summary.function(dists_between_subs[upper.tri(dists_between_subs)])
    }
    return(res)
  })
  return(stats)
}

# A handy function that just computes all the summary values we're interested in
computeStats <- function(aln,phy){
  # recover()
  count_subs <- computeParsimonySubstitutionTable(aln,phy)
  
  
  
  
  real_data_totals <- apply(count_subs,1,sum)
  
  
  
  # parsimony-based substitution-profile chi-squared
  scs_norm <- parsSubstitutionChiSquared(NULL,sim_tree,count_subs,reversible=FALSE,use="standardized")
  scs_norm_rev <- parsSubstitutionChiSquared(NULL,sim_tree,count_subs,reversible=TRUE,use="standardized")

  scs_p <- parsSubstitutionChiSquared(NULL,sim_tree,count_subs,reversible=FALSE,use="p-value")
  scs_p_rev <- parsSubstitutionChiSquared(NULL,sim_tree,count_subs,reversible=TRUE,use="p-value")
  
  # counts of mutations
  exp_stats <- apply(count_subs,1,sum)
  total_subs <- sum(count_subs)
  
  # branch-wise variances
  var_stats <- apply(count_subs,1,var)
  
  # distances to nearest substitution of same type
  coloc_stats <- substitutionTypeColocolization(NULL,phy,count_subs,reversible=TRUE,summary.function=median)
  
  stats <- c(scs_norm,scs_norm_rev,scs_p,scs_p_rev,
             total_subs,
             exp_stats,
             var_stats,coloc_stats)
  names(stats) <- c("parsimonySubstitutionChiSquared","parsimonySubstitutionChiSquared(reversible)","parsimonySubstitutionChiSquared(p-value)","parsimonySubstitutionChiSquared(p-value,reversible)",
                    "nSubs",
                    "sum(A->C)","sum(A->G)","sum(A->T)","sum(C->G)","sum(C->T)","sum(G->T)",
                    "sum(C->A)","sum(G->A)","sum(T->A)","sum(G->C)","sum(T->C)","sum(T->G)",
                    "Var(A->C)","Var(A->G)","Var(A->T)","Var(C->G)","Var(C->T)","Var(G->T)",
                    "Var(C->A)","Var(G->A)","Var(T->A)","Var(G->C)","Var(T->C)","Var(T->G)",
                    "Colocolization(A<->C)","Colocolization(A<->G)","Colocolization(A<->T)","Colocolization(C<->G)","Colocolization(C<->T)","Colocolization(G<->T)"
                  )
  return(stats)
}
