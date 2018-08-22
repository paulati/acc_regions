source("~/2018/common/config/paths_config_modneuPhastCons100way.R")
source("/home/rstudio/2018/common/format_files.R")

# 
# sanity.check <- function()
# {
# 
#   setwd(neutral.model.base.path)
#   neutral.model.file.name <- "hg38.phastCons100way_label.mod"
#   neutral.mod <- read.tm(neutral.model.file.name)
#   
#   
#   file.name.tail <- "_mammals.maf"
#   remote.base.folder.path <- remote.mammals.align.base.folder.path
#   local.base.folder.path <- mammals.align.base.path
#   
#   align <- load.msa(chr.id, file.name.tail, 
#                     remote.base.folder.path, 
#                     local.base.folder.path)
#   
#   
#   setwd(conservation.intersection.output.base.path)
#   feat.file.name <- paste("chr", chr.id, "_intersection_sarcopterygii_mammals.bed", sep="")
#   feats.bed <- read.table(feat.file.name, sep=" ")
#   colnames(feats.bed) <- c("chr", "start", "end")
#   seq.name <- paste("chr", chr.id, sep="")
#   informative.elements <- bed.to.feat(seq.name, feats.bed)
#   
#   setwd(obsphyloP.output.base.path)
#   file.name <- paste("chr", chr.id, "_mar.gff", sep="")
#   #write.table(nonParaPhyloP, out.file.name, col.names = TRUE, row.names = FALSE, sep="\t", quote=FALSE)
#   nonParaSigFeats <- read.table(file.name, sep="\t", header=TRUE)
#   
#   #esto es para el extract features
#   informative.elements$seqname <- rep("hg38", nrow(informative.elements))
#   elementAlign <- extract.feature.msa(copy.msa(align), informative.elements, pointer.only=TRUE)
#   
#   
#   consEleModel <- phyloFit(elementAlign, init.mod=neutral.mod, no.opt=c("backgd", "ratematrix"))
#   rarAlign <- extract.feature.msa(copy.msa(align), nonParaSigFeats)
#   rarModel <- phyloFit(rarAlign, init.mod=neutral.mod, no.opt=c("backgd", "ratematrix"))
# 
#   library("ape")
#   par(mfrow=c(1,2), mar=c(5,2,4,2))
#   maxX <- depth.tree(rarModel$tree, "mammals")
#   plot.phylo(read.tree(text=consEleModel$tree), x.lim=c(0, 0.5), main="Conserved Elements")
#   plot.phylo(read.tree(text=rarModel$tree), x.lim=c(0, 0.5), main="MAR")  
#   
#   
# }


empirical.pval <- function(x, dist) 
{
  result <- sum(x <= dist)/length(dist)
  return(result)
}


non.para.sim <- function(nrep, align, informative.elements, split.length, neutral.mod)
{

  align.copy <- copy.msa(align)

  #lo necesito para que coincidan los nombres de secuencias en extract.feature.msa
  informative.elements$seqname <- rep("hg38", nrow(informative.elements))
  
  element.align <- extract.feature.msa(align.copy, informative.elements, pointer.only=TRUE)

  #nrep <- 100000
  
  sim.msa <- sample.msa(element.align, nrep * split.length, replace=TRUE)
  
  # produce features allowing long alignment to be interpreted as 
  # concatenation of shorter alignments
  
  
  startIdx <- seq(from=1, by=split.length, length.out=nrep)
  
  features <- feat(seqname=names.msa(sim.msa)[1], src="sim", feat=".", start=startIdx, end=startIdx+split.length-1)
  
  nonParaPhyloP <- phyloP(neutral.mod, msa=sim.msa, mode="ACC", features=features, branches="mammals")

  
  setwd(nonparasimphyloP.output.base.path)
  out.file.name <- paste("chr", chr.id, "_nonparaPhyloP_", split.length, ".csv", sep="")
  write.table(nonParaPhyloP, out.file.name, col.names = TRUE, row.names = FALSE, sep="\t", quote=FALSE)
  
  return(nonParaPhyloP)

}

main <-function(chr.id, split.length, nrep)
{
  
  file.name.tail <- "_mammals.maf"
  remote.base.folder.path <- remote.mammals.align.base.folder.path
  local.base.folder.path <- mammals.align.base.path
  
  align <- load.msa(chr.id, file.name.tail, 
                    remote.base.folder.path, 
                    local.base.folder.path)
  
  setwd(neutral.model.base.path)
  neutral.model.file.name <- "hg38.phastCons100way_label.mod"
  neutral.mod <- read.tm(neutral.model.file.name)
  
  
  setwd(conservation.intersection.output.base.path)
  feat.file.name <- paste("chr", chr.id, "_intersection_sarcopterygii_mammals.bed", sep="")
  feats.bed <- read.table(feat.file.name, sep=" ")
  colnames(feats.bed) <- c("chr", "start", "end")
  seq.name <- paste("chr", chr.id, sep="")
  informative.elements <- bed.to.feat(seq.name, feats.bed)
  
  
  nonParaPhyloP <- non.para.sim(nrep, align, informative.elements, split.length, neutral.mod)
  
  setwd(obsphyloP.output.base.path)
  obsPhyloP.file.name <- paste("chr", chr.id, "_obsPhyloP_", split.length, ".csv", sep="")
  obsPhyloP <- read.table(obsPhyloP.file.name, sep="\t", header=TRUE)
  
  nonParaPval <- sapply(obsPhyloP$lnlratio, empirical.pval, nonParaPhyloP$lnlratio)
  nonParaFDR <- p.adjust(nonParaPval, method="BH")
  
  setwd(split.feats.base.path)
  file.name <- paste("chr", chr.id, "_", split.feats.file.name, "_", split.length, ".csv", sep="")
  split.feats <- read.table(file.name, sep="\t", header = TRUE)
  
  nonParaSigFeats <- split.feats[nonParaFDR < 0.05,]  

  nonParaSigFeats$feature <- "mammalsAR"
  nonParaSigFeats$score <- obsPhyloP$lnlratio[nonParaFDR < 0.05]
  nonParaSigFeats$seqname <- paste("hg38.chr", chr.id, sep="")
  
  #ordeno por score:
  nonParaSigFeats.sort <- nonParaSigFeats[order(- nonParaSigFeats$score), ]
  file.name <- paste("chr", chr.id, "_mar_", split.length, ".gff", sep="")
  setwd(nonparasimphyloP.gff.output.base.path)
  #no usar write feat porque redondea el score a 0
  #write.feat(nonParaSigFeats.sort, file.name)  
  write.table(nonParaSigFeats.sort, file.name, sep="\t", 
              col.names = FALSE, row.names = FALSE, quote = FALSE)  

  q <- nrow(nonParaSigFeats.sort)
  result.bed <- data.frame(chr=character(q),
                           start=integer(q),
                           end=character(q))
  
  chr.name <- paste("chr", chr.id, sep="")
  result.bed$chr <- rep(chr.name, q)
  result.bed$start <- nonParaSigFeats.sort$start
  result.bed$end <- nonParaSigFeats.sort$end
  
  setwd(nonparasimphyloP.bed.output.base.path)
  file.name <- paste("chr", chr.id, "_acc_", split.length, ".bed", sep="")
  write.table(result.bed, file.name, sep="\t", col.names = FALSE, row.names = FALSE, quote = FALSE)
   
}


nrep <- 100000
for(chr.id in 5:9)
{
  split.length <- 50
  main(chr.id, split.length, nrep)
  gc()

  split.length <- 25
  main(chr.id, split.length, nrep)
  gc()
  
}

