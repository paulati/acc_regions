source("/home/rstudio/2018/common/config/paths_config_modneuPhastCons100way.R")
source("/home/rstudio/2018/common/format_files.R")

run.obsPhyloP <- function(align, neutral.mod, split.feats) {
  
  #align.copy <- copy.msa(align)
  
  obsPhyloP <- phyloP(neutral.mod, 
                      msa=align, 
                      mode="ACC",
                      features=split.feats, 
                      branches = "mammals")    
  
  return(obsPhyloP)
  
}


main <- function(chr.id, split.length)
{
  
  file.name.tail <- "_mammals.maf"
  remote.base.folder.path <- remote.mammals.align.base.folder.path
  local.base.folder.path <- mammals.align.base.path

  align <- load.msa(chr.id, file.name.tail, 
                    remote.base.folder.path, 
                    local.base.folder.path)
  
  # local.file.name.gz <- paste("chr", chr.id, "_mammals.maf.gz", sep="")
  # local.file.name <- paste("chr", chr.id, "_mammals.maf", sep="")
  # align <- read.msa.local(base.path, local.file.name.gz, local.file.name)
    
  setwd(neutral.model.base.path)
  neutral.model.file.name <- "hg38.phastCons100way_label.mod"
  neutral.mod <- read.tm(neutral.model.file.name)
  
  setwd(conservation.intersection.output.base.path)
  feat.file.name <- paste("chr", chr.id, "_intersection_sarcopterygii_mammals.bed", sep="")
  feats.bed <- read.table(feat.file.name, sep=" ")
  colnames(feats.bed) <- c("chr", "start", "end")
  seq.name <- paste("chr", chr.id, sep="")
  feats <- bed.to.feat(seq.name, feats.bed)
  
  #split de feats de long split.length
  setwd(split.feats.base.path)
  file.name <- paste("chr", chr.id, "_", split.feats.file.name, "_", split.length, ".csv", sep="")
  if(file.exists(file.name))
  {
    split.feats <- read.table(file.name, sep="\t", header = TRUE)
  }  else  {
    split.feats <- split.feat(feats, f=split.length, drop=TRUE)
    write.table(split.feats, file.name, sep="\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
  }
  
  phyloP.result <- run.obsPhyloP(align, neutral.mod, split.feats) 
  
  #ordeno los resultados de phyloP por score
  phyloP.result.sort <- phyloP.result[order(-phyloP.result$score), ]
  
  #escribo los resultados de obs.phyloP:
  setwd(obsphyloP.output.base.path)
  out.file.name <- paste("chr", chr.id, "_obsPhyloP", "_", split.length, ".csv",  sep="")
  write.table(phyloP.result.sort, out.file.name, col.names = TRUE, row.names = FALSE, sep="\t", quote=FALSE)
  

}

for(chr.id in c(21:22))
{
  split.length <- 25
  main(chr.id, split.length)
  gc()

  split.length <- 50
  main(chr.id, split.length)
  gc()
}





