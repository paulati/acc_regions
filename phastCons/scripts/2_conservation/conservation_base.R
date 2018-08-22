#source('/home/rstudio/2018/phastCons/scripts/config/paths_config.R')
source('/home/rstudio/2018/common/config/paths_config_modneuPhastCons100way.R')
source('/home/rstudio/2018/common/aws_base.R')


library(rphast)
library(R.utils)

# 
# load.alignment.base <- function(chr.id, file.name.tail, 
#                                 remote.base.folder.path, 
#                                 local.base.folder.path) {
#   
#   #remote.base.folder.path <- remote.sarcopterygii.align.base.folder.path
#   #file.name.gz <- paste("chr", chr.id, "_sarcopterygii.maf.gz", sep="")
#   file.name.gz <- paste("chr", chr.id, file.name.tail, ".gz", sep="")
#   #local.base.folder.path <- sarcopterygii.align.base.path
#   #align.local.file.name <- paste("chr", chr.id, "_sarcopterygii.maf", sep="")
#   align.local.file.name <- paste("chr", chr.id, file.name.tail, sep="")
#   
#   setwd(local.base.folder.path)
#   
#   if(! file.exists(align.local.file.name)){
#     download.from.s3(account.key, account.secret,
#                      remote.base.folder.path, file.name.gz,
#                      local.base.folder.path, file.name.gz,
#                      align.local.file.name, bucket.name)
#     
#   }
#   
#   align <- read.msa(align.local.file.name)  
#   
#   return(align)  
#   
# }
# 


main <- function()
{
  
  #chr.id <- 11
  
  for(chr.id in c(2:2))
  {
    print(chr.id)
    
    align <- load.alignment(chr.id)
    
    #libero mem, rompe con chr 1 y 2
    gc()
    
    required.feats <- required.species.feats(align)
    
    setwd(neutral.model.base.path)
    neutralMod <- read.tm(neutral.model.file.name)

    #libero mem, rompe con chr 1 y 2
    gc()
    
    elements <- phastCons(align, neutralMod, expected.length=45,
                          target.coverage=0.3, rho=0.3, viterbi=TRUE)
    
    #cons.elements es un objeto de tipo feat
    cons.elements <- elements$most.conserved
    
    intersection.cons.reqs <- coverage.feat(required.feats, cons.elements, or=FALSE, get.feats=TRUE)
    
    q <- nrow(intersection.cons.reqs)
    bed.data <- data.frame(chr=character(q),
                           start=integer(q),
                           end=integer(q))
    
    char.name <- paste("chr", chr.id, sep="")
    
    bed.data$chr <- rep(char.name, q)
    bed.data$start <- intersection.cons.reqs$start
    bed.data$end <- intersection.cons.reqs$end
    out.file.name <- paste("chr", chr.id, "_phastCons_mostConserved_conFiltro.bed", sep="")
    setwd(phastConsElements.out.base.path)
    write.table(bed.data, out.file.name, sep="\t", col.names = FALSE, row.names = FALSE, quote=FALSE)
    
    
    rm()
    gc()
    
  }  
  
}