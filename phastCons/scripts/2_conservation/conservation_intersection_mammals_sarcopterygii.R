library(rphast)

source("~/2018/common/config/paths_config_modneuPhastCons100way.R")
source("~/2018/common/format_files.R") 

load.data <- function(chr.id)
{
  setwd(conservation.sarcopterygii.base.path)
  file.name <- paste("chr", chr.id, "_phastCons_mostConserved_conFiltro.bed", sep="")
  data.sarcopterygii <- read.table(file.name, sep="\t")
  colnames(data.sarcopterygii) <- c("chr", "start", "end")
  
  setwd(conservation.mammals.base.path)
  file.name <- paste("chr", chr.id, "_phastCons_mostConserved_conFiltro.bed", sep="")
  data.mammals <- read.table(file.name, sep="\t")
  colnames(data.mammals) <- c("chr", "start", "end")
  
  result <- list(sarcopterygii=data.sarcopterygii,
                 mammals = data.mammals)
  
  return(result)
    
}


#chr.id <- 22
for(chr.id in c(2:10))
{
  data <- load.data(chr.id)

  seq.name <- paste('chr', chr.id, sep="")
  
  feat.1 <- bed.to.feat(seq.name, data$sarcopterygii)
  feat.2 <- bed.to.feat(seq.name, data$mammals)
  
  intersection.regions <- coverage.feat(feat.1, feat.2, 
                                        or=FALSE, get.feats=TRUE)
  
  output <- feat.to.bed(chr.id, intersection.regions)
  
  setwd(conservation.intersection.output.base.path)
  file.name <- paste("chr", chr.id, "_intersection_sarcopterygii_mammals.bed", sep="")
  write.table(output, file.name, quote = FALSE, col.names = FALSE, row.names = FALSE)

}