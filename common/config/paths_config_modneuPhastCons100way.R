account.key <- ""
account.secret <- ""

base.path <- "/home/rstudio/2018"


bucket.name <- "acc-regions-2018"

remote.feats.base.folder.path <- "/download/hg38_feats/"
feats.base.path <- paste(base.path, "/phastCons/data/tmp/feats", sep="")

remote.align.base.folder.path <- "/download/alignments/multiz100way/"
#align.base.path <- paste(base.path, "/phastCons/data/tmp/align", sep="")
align.base.path <- "/home/rstudio/disco_tmp"

remote.sarcopterygii.align.base.folder.path <- "/preparation/sarcopterygii/"
#align.base.path <- paste(base.path, "/phastCons/data/tmp/align", sep="")
sarcopterygii.align.base.path <- "/home/rstudio/2018"


remote.mammals.align.base.folder.path <- "/preparation/mammals/"
#align.base.path <- paste(base.path, "/phastCons/data/tmp/align", sep="")
mammals.align.base.path <- "/home/rstudio/disco_tmp"




tree.base.path <- paste(base.path, "/phastCons/data/tree", sep="")
tree.file.name <- "hg38.100way.nh"


neutral.model.base.path <- paste(base.path, "/phastCons/data/neutral_model", sep="")
neutral.model.file.name <- "hg38.phastCons100way.mod"

phastConsElements.out.base.path <- paste(base.path, "/phastCons/data/output/modneu_phastCons100way_mammals", sep="")


conservation.sarcopterygii.base.path <- paste(base.path, "/phastCons/data/output/modneu_phastCons100way_sarcopterygii", sep="")
conservation.mammals.base.path <- paste(base.path, "/phastCons/data/output/modneu_phastCons100way_mammals", sep="")

conservation.intersection.output.base.path <- paste(base.path, "/phastCons/data/output/intersection_modneu_phastCons100way_sarcopterygii_mammals", sep="")

obsphyloP.output.base.path <- paste(base.path, "/acc/data/obs_phyloP/sarcopterygii_mammals", sep="")

nonparasimphyloP.output.base.path <- paste(base.path, "/acc/data/non_param_sim_phyloP/sarcopterygii_mammals/nonParaPhyloP", sep="")
nonparasimphyloP.gff.output.base.path <- paste(base.path, "/acc/data/non_param_sim_phyloP/sarcopterygii_mammals/acc_elements_gff", sep="")
nonparasimphyloP.bed.output.base.path <- paste(base.path, "/acc/data/non_param_sim_phyloP/sarcopterygii_mammals/acc_elements_bed", sep="")

split.feats.base.path <- paste(base.path, "/acc/data/obs_phyloP/split_feats", sep="")
split.feats.file.name <- "sarcopterygii_mammals_split_feats"

