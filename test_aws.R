library(aws.s3)

account.key <- ""
account.secret <- ""


bucket.name <- "acc-regions-2018"

remote.base.folder.path <- "/download/alignments/multiz100way/"

local.base.folder.path <- "/home/rstudio/2018/phastCons/data/tmp"

file.name.gz <- "chr11.maf.gz"

local.file.name.gz <- "chr11.maf.gz"
  
setwd(local.base.folder.path)

object.name <- paste(remote.base.folder.path, file.name.gz, sep="")

setwd(local.base.folder.path)

aws.s3::save_object(object = object.name,
                    key = account.key,
                    secret = account.secret,
                    bucket = bucket.name,
                    file = local.file.name.gz)


gc()

raw.object <- aws.s3::get_object(object = object.name, 
                   key = account.key,
                   secret = account.secret,
                   bucket = bucket.name)




#getwd()
