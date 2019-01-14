library(dplyr)
setwd("/media/hdd/Prabakaran_Lab/data")
novel_peptides <- read.csv("novel_peptide_region_features.csv",header=T)
peptide_alignments <- read.csv("novel_peptide_frame_evidence.csv",header=T,stringsAsFactors = F)
sense <- novel_peptides[novel_peptides$strand_relationship=="sense",]
antisense <- novel_peptides[novel_peptides$strand_relationship=="anti-sense",]
piechartlabels <- c("protein_coding","sense", "antisense","lincRNA","miRNA","Other RNA","pseudogene","bidirectional_promoter","TEC")
sensevals <- c(243,5,50,20,15,4+3+0+13+9,7+1+2+4+54,1,49)
antisensevals <- c(259,7+6+0,37,16,27,1+3+2+22+19,2+7+0+1+63,1,50)
pie(sensevals,paste0(piechartlabels,"(",sensevals,")"),main="Biotype Relationship of Genomic Features Sense to Novel Transcripts")
pie(antisensevals,paste0(piechartlabels,"(",antisensevals,")"),main="Biotype Relationship of Genomic Features Anti-sense to Novel Transcripts")

fragmentclassification <- distinct(read.table("novel_fragmennt_alignment_classification.txt",header=T))
#sapply(unique(as.character(fragmentclassification$region)), function(x){
#  temp <- fragmentclassification[fragmentclassification$region==x,]
#  print(temp)
#})

hist(table(fragmentclassification$region), main="Features per Peptide Fragment Alignment", xlab="Number of features per fragment")
pie(summary(fragmentclassification$strand_relationship), paste0(names(summary(fragmentclassification$strand_relationship)),"(",summary(fragmentclassification$strand_relationship),")"),main="Relationship of peptide fragments to neighbouring features")
peptidesense <- fragmentclassification[fragmentclassification$strand_relationship=="sense",]
peptideantisense <- fragmentclassification[fragmentclassification$strand_relationship=="anti-sense",]
intergenic <- fragmentclassification[fragmentclassification$strand_relationship=="intergenic",]
intronic <- fragmentclassification[fragmentclassification$exonic_relationship=="intronic" & !is.na(fragmentclassification$exonic_relationship),]
UTR3 <- fragmentclassification[!is.na(fragmentclassification$exonic_relationship) & fragmentclassification$exonic_relationship=="within" & fragmentclassification$cds_relationship=="3UTR" & !is.na(fragmentclassification$cds_relationship), ]
UTR5 <- fragmentclassification[!is.na(fragmentclassification$exonic_relationship) & fragmentclassification$exonic_relationship=="within" & fragmentclassification$cds_relationship=="5UTR" & !is.na(fragmentclassification$cds_relationship), ]
cds <- fragmentclassification[!is.na(fragmentclassification$cds_relationship) & (fragmentclassification$cds_relationship == "within" | fragmentclassification$cds_relationship == "upstream_overlap" | fragmentclassification$cds_relationship == "downstream_overlap"), ]
ncrna <- peptidesense[grep("RNA",peptidesense$biotype),]
pseudogene <- peptidesense[grep("pseudogene",peptidesense$biotype),]

peptidefragmentfeatureclasses <- c("antisense","ncRNA","proteincoding","pseudogene","TEC")
peptidefragmentfeatureclassescounts <- c(45,6+1,941,23+5+25,4)
pie(peptidefragmentfeatureclassescounts,paste0(peptidefragmentfeatureclasses,"(",peptidefragmentfeatureclassescounts,")"),main="Biotype of features sense to peptide alignments")

overalltags <- c("3UTR","5UTR","Intergenic","Within or Overlap CDS","Intronic","ncRNA","antisense to some feature","pseudogene")
overalltagcounts <- c(12,2,22,60,900,16+5+1+1,341,5+8)
overalltagcounts <- c(nrow(UTR3),nrow(UTR5), nrow(intergenic), nrow(cds), nrow(intronic), nrow(ncrna), nrow(peptideantisense),nrow(pseudogene))
pie(overalltagcounts,paste0(overalltags,"(",overalltagcounts,")"),main="Relationship of peptides to genomic features\n(617/632 peptide alignments) (804/835 annotations)")

