source("/home/rstudio/2018/common/aws_base.R")
library(rphast)

bed.to.feat <- function(seq.name, bed)
{
  q <- nrow(bed)
  df <- data.frame(seqname=character(q),
                   src=character(q),
                   feature=character(q),
                   start=integer(q),
                   end=integer(q))
  
  df$seqname <- rep(seq.name, q)
  df$src <- rep(".", q)
  df$feature <- rep(".", q)
  df$start <- bed$start
  df$end <- bed$end
  
  result <- feat(df$seqname, df$src, df$feature, df$start, df$end)
  return(result)
}


feat.to.bed <- function(chr.id, feat){
  
  q <- nrow(feat)
  result <- data.frame(chr=character(q),
                       start=character(q),
                       end=character(q))
  chr.name <- paste("chr", chr.id, sep="")
  result$chr <- rep(chr.name, q)
  result$start <- feat$start
  result$end <- feat$end
  return(result)
  
}

load.msa <- function(chr.id, file.name.tail, 
                           remote.base.folder.path, 
                           local.base.folder.path)
{
  
  file.name.gz <- paste("chr", chr.id, file.name.tail, ".gz", sep="")
  align.local.file.name <- paste("chr", chr.id, file.name.tail, sep="")
  
  setwd(local.base.folder.path)
  
  if(! file.exists(align.local.file.name)){
    if(file.exists(file.name.gz))
    {
      gunzip(local.file.name.gz, align.local.file.name, remove = FALSE)
      
    }
    else
    {
      download.from.s3(account.key, account.secret,
                       remote.base.folder.path, file.name.gz,
                       local.base.folder.path, file.name.gz,
                       align.local.file.name, bucket.name)
    }
  }
  
  align <- read.msa(align.local.file.name)  
  
  return(align)    
    
}


