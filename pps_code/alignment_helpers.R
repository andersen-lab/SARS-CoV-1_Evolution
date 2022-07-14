library(phangorn)

getAmbiguities <- function(aln) {
  # recover()
  char_aln <- as.character(aln)
  nucs <- c("a","c","g","t","A","C","G","T")
  ambiguities <- lapply(1:dim(char_aln)[2],function(j){
    ambs <- which(!(char_aln[,j] %in% nucs))
    if ( length(ambs) > 0 ) {
      return(list(column=j,rows=ambs,characters=char_aln[ambs,j]))
    } else {
      return(NULL)
    }
  })
  is_null <- unlist(lapply(ambiguities,is.null))
  return(ambiguities[!is_null])
}

addAmbiguities <- function(simAln,ambiguities) {
  sim <- as.character(simAln)
  if (length(ambiguities) > 0) {
    for (idx in 1:length(ambiguities)) {
      j <- ambiguities[[idx]]$column
      sim[ambiguities[[idx]]$rows,j] <- ambiguities[[idx]]$characters
    }
  }
  return(as.phyDat(sim))
}

