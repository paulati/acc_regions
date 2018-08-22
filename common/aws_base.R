library(aws.s3)

download.from.s3 <- function(account.key, account.secret,
                             remote.base.folder.path, file.name.gz,
                             local.base.folder.path, local.file.name.gz,
                             local.file.name, bucket.name)
{
  setwd(local.base.folder.path)

  object.name <- paste(remote.base.folder.path, file.name.gz, sep="")
  aws.s3::save_object(object = object.name,
                      key = account.key,
                      secret = account.secret,
                      bucket = bucket.name,
                      file = local.file.name.gz)

 #borra el gz:
 gunzip(local.file.name.gz, local.file.name)
}