source('/home/rstudio/2018/phastCons/scripts/config/paths_config.R')
source('/home/rstudio/2018/phastCons/scripts/common/aws_base.R')


library(rphast)
library(R.utils)


load.feats <- function(chr.id) {
  
  remote.base.folder.path <- remote.feats.base.folder.path
  file.name.gz <- paste("hg38_chr", chr.id, "_ncbiRefSeq.gtf.gz", sep="")
  local.base.folder.path <- feats.base.path
  feat.local.file.name <- paste("hg38_chr", chr.id, "_ncbiRefSeq.gtf", sep="")
  
  download.from.s3(account.key, account.secret,
                   remote.base.folder.path, file.name.gz,
                   local.base.folder.path, file.name.gz,
                   feat.local.file.name, bucket.name)
  
  setwd(feats.base.path)
  feats <- read.feat(feat.local.file.name)  

  return(feats)  
}

load.alignment <- function(chr.id) {
  
  remote.base.folder.path <- remote.align.base.folder.path
  file.name.gz <- paste("chr", chr.id, ".maf.gz", sep="")
  local.base.folder.path <- align.base.path
  align.local.file.name <- paste("chr", chr.id, ".maf", sep="")
  
  download.from.s3(account.key, account.secret,
                   remote.base.folder.path, file.name.gz,
                   local.base.folder.path, file.name.gz,
                   align.local.file.name, bucket.name)
  
  setwd(align.base.path)
  align <- read.msa(align.local.file.name)  
  
  return(align)  

}


chr.ids <-  c(c(1:22, "X", "Y", "M"))

chr.id <- "11"

#feat
feats <- load.feats(chr.id)

#alignment
#align <- load.alignment(chr.id)

align.local.file.name <- "chr11.maf"

setwd(align.base.path)
align <- read.msa(align.local.file.name)  




#tree
setwd(tree.base.path)
tree <- read.newick.tree(tree.file.name)

#pensar: que especies mantener para la neutralidad? las 100? 

