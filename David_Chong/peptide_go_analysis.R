# GO terms obtained using interpro are collated to determine enriched/depleted terms with respect to known proteins in T and B cells
# Significantly enriched/depleted terms are clustered using GOSim and a bar plot of enrichment/depletion fold to aid interpretation.
library(qvalue)
library(GO.db)
library(GOSim)

setwd("/media/hdd/Prabakaran_Lab/data")
sorfs <- read.table("sorfs_with_evidence_of_transcription_or_translation.txt",header=F)
desorfs <- read.table("de_sorfs.txt",header=F)
non_redundant_sorfs <- data.frame(read.table("non_redundant_sorfs.txt",header=F)[,1])
altorfs <- read.table("altorf_evidence_of_transcription_or_translation.txt",header=F)
known_proteins <- read.table("transcribed_and_translated_proteins_genenames.txt",header=F)
known_protein_mapping <- read.table("transcribed_and_translated_proteins_mapping.txt",header=F)

sorfgoterms <- read.table("interproanalysis/sorf_goterms_table.tsv",header=F)
desorfgoterms <- read.table("interproanalysis/desorf_goterms_table.tsv",header=F)
non_redundant_sorfsgoterms <- read.table("interproanalysis/non_redundant_sorf_goterms_table.tsv",header=F)
altorfgoterms <- read.table("interproanalysis/altorf_goterms_table.tsv",header=F)
knownproteingoterms <- read.table("interproanalysis/mus_muscularis_transcribed_and_translated_known_protein_goterms_table.tsv",header=F)

sorfgotermset <- read.table("interproanalysis/sorf_goterms_unique_goterms.tsv",header=F)
desorfgotermset <- read.table("interproanalysis/desorf_goterms_unique_goterms.tsv",header=F)
non_redundant_sorfsgotermset <- read.table("interproanalysis/non_redundant_sorf_goterms_unique_goterms.tsv",header=F)
altorfgotermset <- read.table("interproanalysis/altorf_goterms_unique_goterms.tsv",header=F)
knownproteingotermset <- read.table("interproanalysis/mus_muscularis_transcribed_and_translated_known_protein_goterms_unique_goterms.tsv",header=F)

print("sORF and altORF common GO terms:")
print(sum(!is.na(match(sorfgotermset$V1,altorfgotermset$V1))))

print("sORF and known protein common GO terms:")
print(sum(!is.na(match(sorfgotermset$V1,knownproteingotermset$V1))))

print("altORF and known protein common GO terms:")
print(sum(!is.na(match(altorfgotermset$V1,knownproteingotermset$V1))))

print("sORFs with at least 1 GO term")
print(length(unique(sorfgoterms$V1)))

print("altORFs with at least 1 GO term")
print(length(unique(altorfgoterms$V1)))

knowngenegoterms <- matrix(c(as.character(unique(known_protein_mapping$V1)), sapply(unique(known_protein_mapping$V1), function(x){
  uniprotids <- known_protein_mapping[!is.na(match(known_protein_mapping$V1,x)),2]
  goterms <- as.character(knownproteingoterms[!is.na(match(knownproteingoterms$V1,uniprotids)),2])
  if (length(goterms) > 0) {
    uniquegoterms <- c()
    for(i in 1:length(goterms)) {
      uniquegoterms <- c(uniquegoterms, unlist(strsplit(goterms[i],"|",fixed=T)))
    }
    uniquegoterms <- unique(uniquegoterms)
    return(paste(uniquegoterms,sep="|",collapse="|"))
  }
  else return("")
})),ncol=2)

goterms = unlist(Term(GOTERM))
print("Known genes with at least 1 GO term")
print(sum(knowngenegoterms[,2]!=""))

sorfknowngooverlap <- knownproteingotermset[match(sorfgotermset$V1,knownproteingotermset$V1),]
sorfknowngooverlap <- sorfknowngooverlap[!is.na(sorfknowngooverlap)]

nrsorfknowngooverlap <- knownproteingotermset[match(non_redundant_sorfsgotermset$V1,knownproteingotermset$V1),]
nrsorfknowngooverlap <- nrsorfknowngooverlap[!is.na(nrsorfknowngooverlap)]

sorfdesorfoverlap <- sorfgotermset[match(desorfgotermset$V1,sorfgotermset$V1),]
sorfdesorfoverlap <- sorfdesorfoverlap[!is.na(sorfdesorfoverlap)]

sorfgotermprops <- sapply(sorfknowngooverlap,function(x) {length(grep(x,sorfgoterms$V2))/nrow(sorfs)})
nrsorfgotermprops <- sapply(nrsorfknowngooverlap,function(x) {length(grep(x,non_redundant_sorfsgoterms$V2))/nrow(non_redundant_sorfs)})
knowngenesorfgotermprops <- sapply(sorfknowngooverlap,function(x){length(grep(x,knowngenegoterms[,2]))/nrow(knowngenegoterms)})
knowngenenrsorfgotermprops <- sapply(nrsorfknowngooverlap,function(x){length(grep(x,knowngenegoterms[,2]))/nrow(knowngenegoterms)})

sorfdesorfgotermprops <- sapply(sorfdesorfoverlap,function(x) {length(grep(x,sorfgoterms$V2))/nrow(sorfs)})
desorfdesorfgotermprops <- sapply(sorfdesorfoverlap,function(x) {length(grep(x,desorfgoterms$V2))/nrow(desorfs)})

sorfknownoverlapchisqtests <- sapply(sorfknowngooverlap,function(x) {chisq.test(x=c(length(grep(x,sorfgoterms$V2)),nrow(sorfs)-length(grep(x,sorfgoterms$V2))),p=c(length(grep(x,knowngenegoterms[,2]))/nrow(knowngenegoterms),1-length(grep(x,knowngenegoterms[,2]))/nrow(knowngenegoterms)),simulate.p.value=T)})
chisqpvals <- unlist(sorfknownoverlapchisqtests[seq(3,4410,9)])
qvals <- qvalue(chisqpvals)$qvalue

nrsorfknownoverlapchisqtests <- sapply(nrsorfknowngooverlap,function(x) {chisq.test(x=c(length(grep(x,non_redundant_sorfsgoterms$V2)),nrow(sorfs)-length(grep(x,non_redundant_sorfsgoterms$V2))),p=c(length(grep(x,knowngenegoterms[,2]))/nrow(knowngenegoterms),1-length(grep(x,knowngenegoterms[,2]))/nrow(knowngenegoterms)),simulate.p.value=T)})
nrchisqpvals <- unlist(nrsorfknownoverlapchisqtests[seq(3,4410,9)])
nrqvals <- qvalue(nrchisqpvals)$qvalue

sorfdesorfoverlapchisqtests <- sapply(sorfdesorfoverlap,function(x) {chisq.test(x=c(length(grep(x,desorfgoterms$V2)),nrow(desorfs)-length(grep(x,desorfgoterms$V2))),p=c(length(grep(x,sorfgoterms[,2]))/nrow(sorfgoterms),1-length(grep(x,sorfgoterms[,2]))/nrow(sorfgoterms)),simulate.p.value=T)})
desorfchisqpvals <- unlist(sorfdesorfoverlapchisqtests[seq(3,4410,9)])
qvals <- qvalue(desorfchisqpvals)$qvalue

print("Number of significantly different sORF GO terms:")
print(sum(qvalue(chisqpvals)$qvalue < 0.01 & chisqpvals < 0.01))

print("Number of significantly different non-redundant sORF GO terms:")
print(sum(qvalue(nrchisqpvals)$qvalue < 0.01 & nrchisqpvals < 0.01))

print("Number of significantly GO terms DE:")
print(sum(qvalue(desorfchisqpvals)$qvalue < 0.01 & desorfchisqpvals < 0.01))

significantgoterms <- sorfknowngooverlap[qvalue(chisqpvals)$qvalue < 0.01 & chisqpvals < 0.01]
significantgotermprops <- sorfgotermprops[qvalue(chisqpvals)$qvalue < 0.01 & chisqpvals < 0.01] - knowngenesorfgotermprops[qvalue(chisqpvals)$qvalue < 0.01 & chisqpvals < 0.01]
enriched <- significantgoterms[(sorfgotermprops[qvalue(chisqpvals)$qvalue < 0.01 & chisqpvals < 0.01] - knowngenesorfgotermprops[qvalue(chisqpvals)$qvalue < 0.01 & chisqpvals < 0.01]) > 0]
depleted <- significantgoterms[(sorfgotermprops[qvalue(chisqpvals)$qvalue < 0.01 & chisqpvals < 0.01] - knowngenesorfgotermprops[qvalue(chisqpvals)$qvalue < 0.01 & chisqpvals < 0.01]) < 0]
enricheddescriptions <- goterms[as.character(enriched)]
depleteddescriptions <- goterms[as.character(depleted)]
significantgotermsdescriptions <- goterms[as.character(significantgoterms)]
png("../images/sORF_significant_goterm_props.png")
par(las=1,mar=c(5.1,25.1,2.1,2.1),cex=0.6)
barplot(significantgotermprops[order(significantgotermprops)], names.arg=significantgotermsdescriptions[order(significantgotermprops)], horiz=T, main="Significantly Enriched and Depleted sORF GO Terms", xlab="Proportion compared to known proteins")
dev.off()

nrsignificantgoterms <- nrsorfknowngooverlap[qvalue(nrchisqpvals)$qvalue < 0.01 & nrchisqpvals < 0.01]
nrsignificantgotermprops <- nrsorfgotermprops[qvalue(nrchisqpvals)$qvalue < 0.01 & nrchisqpvals < 0.01] - knowngenenrsorfgotermprops[qvalue(nrchisqpvals)$qvalue < 0.01 & nrchisqpvals < 0.01]
nrenriched <- nrsignificantgoterms[(nrsorfgotermprops[qvalue(nrchisqpvals)$qvalue < 0.01 & nrchisqpvals < 0.01] - knowngenenrsorfgotermprops[qvalue(nrchisqpvals)$qvalue < 0.01 & nrchisqpvals < 0.01]) > 0]
nrdepleted <- nrsignificantgoterms[(nrsorfgotermprops[qvalue(nrchisqpvals)$qvalue < 0.01 & nrchisqpvals < 0.01] - knowngenenrsorfgotermprops[qvalue(nrchisqpvals)$qvalue < 0.01 & nrchisqpvals < 0.01]) < 0]
nrenricheddescriptions <- goterms[as.character(nrenriched)]
nrdepleteddescriptions <- goterms[as.character(nrdepleted)]
nrsignificantgotermsdescriptions <- goterms[as.character(nrsignificantgoterms)]
png("../images/nrsORF_significant_goterm_props.png")
par(las=1,mar=c(5.1,20.1,2.1,2.1),cex=0.6)
barplot(nrsignificantgotermprops[order(nrsignificantgotermprops)], names.arg=nrsignificantgotermsdescriptions[order(nrsignificantgotermprops)], horiz=T, main="Significantly Enriched and Depleted non redundant sORF GO Terms", xlab="Proportion compared to known proteins")
dev.off()

print("Enriched GO terms:")
print(length(enriched))
print("Depleted GO terms:")
print(length(depleted))

print("NR Enriched GO terms:")
print(length(nrenriched))
print("NR Depleted GO terms:")
print(length(nrdepleted))

setEvidenceLevel(organism="org.Mm.eg.db")
gotermdistmat <- getTermSim(as.character(significantgoterms))
colnames(gotermdistmat) <- goterms[colnames(gotermdistmat)]
rownames(gotermdistmat) <- goterms[rownames(gotermdistmat)]
png("../images/sORF_Overall_significant_go_terms_heatmap.png",width=800,height=800)
heatmap(gotermdistmat[apply(gotermdistmat,1,sum) != 0,apply(gotermdistmat,1,sum) != 0], main="Clustering of sORF significant GO terms within Mouse GO Set", margins=c(20,20))
dev.off()

gotermdistmat <- getTermSim(as.character(enriched))
colnames(gotermdistmat) <- goterms[colnames(gotermdistmat)]
rownames(gotermdistmat) <- goterms[rownames(gotermdistmat)]
png("../images/sORF_Enriched_significant_go_terms_heatmap.png",width=800,height=800)
heatmap(gotermdistmat[apply(gotermdistmat,1,sum) != 0,apply(gotermdistmat,1,sum) != 0], main="Clustering of Enriched sORF significant GO terms within Mouse GO Set", margins=c(20,20))
dev.off()

gotermdistmat <- getTermSim(as.character(depleted))
colnames(gotermdistmat) <- goterms[colnames(gotermdistmat)]
rownames(gotermdistmat) <- goterms[rownames(gotermdistmat)]
png("../images/sORF_Depleted_significant_go_terms_heatmap.png",width=800,height=800)
heatmap(gotermdistmat[apply(gotermdistmat,1,sum) != 0,apply(gotermdistmat,1,sum) != 0], main="Clustering of Depleted sORF significant GO terms within Mouse GO Set", margins=c(20,20))
dev.off()

gotermdistmat <- getTermSim(as.character(nrsignificantgoterms))
colnames(gotermdistmat) <- goterms[colnames(gotermdistmat)]
rownames(gotermdistmat) <- goterms[rownames(gotermdistmat)]
png("../images/nrsORF_Overall_significant_go_terms_heatmap.png",width=800,height=800)
heatmap(gotermdistmat[apply(gotermdistmat,1,sum) != 0,apply(gotermdistmat,1,sum) != 0], main="Clustering of non-redundant sORF significant GO terms within Mouse GO Set", margins=c(20,20))
dev.off()

gotermdistmat <- getTermSim(as.character(nrenriched))
colnames(gotermdistmat) <- goterms[colnames(gotermdistmat)]
rownames(gotermdistmat) <- goterms[rownames(gotermdistmat)]
png("../images/nrsORF_Enriched_significant_go_terms_heatmap.png",width=800,height=800)
heatmap(gotermdistmat[apply(gotermdistmat,1,sum) != 0,apply(gotermdistmat,1,sum) != 0], main="Clustering of Enriched non-redundant sORF significant GO terms within Mouse GO Set", margins=c(20,20))
dev.off()

gotermdistmat <- getTermSim(as.character(nrdepleted))
colnames(gotermdistmat) <- goterms[colnames(gotermdistmat)]
rownames(gotermdistmat) <- goterms[rownames(gotermdistmat)]
png("../images/nrsORF_Depleted_significant_go_terms_heatmap.png",width=800,height=800)
heatmap(gotermdistmat[apply(gotermdistmat,1,sum) != 0,apply(gotermdistmat,1,sum) != 0], main="Clustering of Depleted non-redundant sORF significant GO terms within Mouse GO Set", margins=c(20,20))
dev.off()
