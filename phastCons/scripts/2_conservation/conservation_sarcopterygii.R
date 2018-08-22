source("~/2018/phastCons/scripts/config/paths_config_modneuPhastCons100way.R")
source("~/2018/phastCons/scripts/2_conservation/conservation_base.R")


load.alignment <- function(chr.id) {
  
  
  file.name.tail <- "_sarcopterygii.maf"
  remote.base.folder.path <- remote.sarcopterygii.align.base.folder.path
  local.base.folder.path <- sarcopterygii.align.base.path
  align <- load.alignment.base(chr.id, file.name.tail, 
                               remote.base.folder.path, 
                               local.base.folder.path)
  return(align)  
  
}

#paula rehacer, esta todo hardcode
# save.bed.file <- function(data, chr.id)
# {
#   
#   q <- nrow(data)
#   bed.data <- data.frame(chr=character(q),
#                          start=integer(q),
#                          end=integer(q))
#   
#   char.name <- paste("chr", chr.id, sep="")
#   
#   bed.data$chr <- rep(char.name, q)
#   bed.data$start <- data$start
#   bed.data$end <- data$end
#   out.file.name <- paste("chr", chr.id, "_intersect_req.bed", sep="")
#   setwd(phastConsElements.out.base.path)
#   write.table(bed.data, out.file.name, sep="\t", col.names = FALSE, row.names = FALSE, quote=FALSE)
#   
#   
# }
  


required.species.feats <- function(align)
{
  cons.req.species <- c("allMis1", "anoCar2")
  cons.req.turtles <- c("cheMyd1", "chrPic2", "pelSin1", "apaSpi1")
  #esta no esta en el alineamiento 
  
  #informative.regions.msa: This function will not alter the value of x even if it is stored as a pointer.
  req.regions <- informative.regions.msa(align, min.numspec=2, spec=cons.req.species, 
                                         refseq="hg38", gaps.inf=FALSE)    
  
  #informative.regions.msa: This function will not alter the value of x even if it is stored as a pointer.  
  turtle.regions <- informative.regions.msa(align, min.numspec=1, spec=cons.req.turtles, 
                                            refseq="hg38", gaps.inf=FALSE)    
  
  #calculo la interseccion entre req.regions y turtle.regions
  #coverage.feat: Any features object passed into this function which is stored as a pointer 
  #to an object stored in C may be reordered (sorted) by this function.
  intersection.regions <- coverage.feat(req.regions, turtle.regions, or=FALSE, get.feats=TRUE)
  
  intersection.regions.order <- sort(intersection.regions, decreasing = FALSE)
  
  return(intersection.regions.order)
  
}


main()

# 
# 
# 
# setwd(phastConsElements.out.base.path)
# out.file.name <- paste("chr", chr.id, "_phastCons_mostConserved_sinFiltro.csv", sep="")
# write.table(cons.elements, out.file.name, sep="\t", col.names = TRUE, row.names = FALSE, quote=FALSE)
# 
# 
# q <- nrow(cons.elements)
# bed.data <- data.frame(chr=character(q),
#                        start=integer(q),
#                        end=integer(q))
# 
# char.name <- paste("chr", chr.id, sep="")
# 
# bed.data$chr <- rep(char.name, q)
# bed.data$start <- cons.elements$start
# bed.data$end <- cons.elements$end
# out.file.name <- paste("chr", chr.id, "_phastCons_mostConserved_sinFiltro.bed", sep="")
# write.table(bed.data, out.file.name, sep="\t", col.names = FALSE, row.names = FALSE, quote=FALSE)
# 



