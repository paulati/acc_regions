
source("~/2018/common/config/paths_config_modneuPhastCons100way.R")
source("~/2018/phastCons/scripts/2_conservation/conservation_base.R")
source("~/2018/common/format_files.R")


load.alignment <- function(chr.id) {
  
  
  file.name.tail <- "_mammals.maf"
  remote.base.folder.path <- remote.mammals.align.base.folder.path
  local.base.folder.path <- mammals.align.base.path
  # align <- load.alignment.base(chr.id, file.name.tail, 
  #                              remote.base.folder.path, 
  #                              local.base.folder.path)
  
  align <- load.msa(chr.id, file.name.tail, 
                    remote.base.folder.path, 
                    local.base.folder.path)
  
  return(align)  
  
}


# paula rehacer, esta todo hardcode
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
  
  # Primate subset: al menos hg38 y otoGar3. (No sé si incorporar a otro primate, como Rhesus)
  # Euarchontoglires subset: al menos mm10, oryCun2, ochPri3
  # Laurasiatheria subset: al menos susScr3, turTru2, bosTau8, felCat8, myoLuc2
  # Afrotheria subset: al menos loxAfr3 o echTel2 
  # Mammal subset: al menos dasNov3 o monDom5 o macEug2, ornAna1.
  
  cons.req.primates <- c("hg38", "otoGar3" )
  cons.req.euarchontoglires <- c("mm10", "oryCun2", "ochPri3")
  cons.req.laurasiatheria <- c("susScr3", "turTru2", "bosTau8", "felCat8", "myoLuc2")
  cons.req.mammal <- c("ornAna1")
  
  cons.req <- c(cons.req.primates, cons.req.euarchontoglires, cons.req.laurasiatheria, cons.req.mammal)
  
  #informative.regions.msa: This function will not alter the value of x even if it is stored as a pointer.
  req.regions.1 <- informative.regions.msa(align, min.numspec=11, spec=cons.req, 
                                         refseq="hg38", gaps.inf=FALSE)    
  
  cons.opc.afrotheria <- c("loxAfr3", "echTel2")
  req.regions.2 <- informative.regions.msa(align, min.numspec=1, spec=cons.opc.afrotheria, 
                                            refseq="hg38", gaps.inf=FALSE)    
  
  cons.opc.mammal <- c("dasNov3", "monDom5", "macEug2")
  req.regions.3 <- informative.regions.msa(align, min.numspec=1, spec=cons.opc.mammal, 
                                           refseq="hg38", gaps.inf=FALSE)    
  
  
  #coverage.feat: Any features object passed into this function which is stored as a pointer 
  #to an object stored in C may be reordered (sorted) by this function.
  intersection.regions <- coverage.feat(req.regions.1, 
                                        req.regions.2, 
                                        req.regions.3, 
                                        or=FALSE, get.feats=TRUE)
  
  intersection.regions.order <- sort(intersection.regions, decreasing = FALSE)
  
  return(intersection.regions.order)
  
}

#main()


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



