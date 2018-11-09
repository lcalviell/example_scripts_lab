#LOAD PACKAGES

suppressMessages(library("GenomicFeatures"))
suppressMessages(library("GenomicAlignments"))
suppressMessages(library("GenomicFiles"))

#LOAD SOME FUNCTIONS
source("/bfd/lcalviel/Ribo-seQC/riboseqc-master/riboseqc_functions.R")


gtf_file="/bfd/shared/annotations/human/gencode25/gencode.v25.annotation.gtf"
twobit_file="/bfd/shared/genomes/human/genome_Homo_sapiens/genome.2bit"
seqinfotwob<-seqinfo(TwoBitFile(twobit_file))
annotation<-makeTxDbFromGFF(file=gtf_file,format="gtf",chrominfo = seqinfotwob)

genes<-genes(annotation)
exons_ge<-exonsBy(annotation,by="gene")
cds_tx<-cdsBy(annotation,"tx",use.names=T)
cds_gen<-cdsBy(annotation,"gene")
cds_gen[[1]][7:8]

reduce(cds_gen[[1]][7:8])

cds_ge<-reduce(cds_gen)
exons_ge<-reduce(exons_ge)

three_utrs<-threeUTRsByTranscript(annotation,use.names=T)

#get threeutr with no cds overlapping

threeutrs<-reduce(setdiff(unlist(three_utrs),unlist(cds_ge),ignore.strand=FALSE))
findOverlaps(threeutrs,cds_ge)

head(threeutrs%over%cds_ge)
table(threeutrs%over%cds_ge)


#PREPARE OBJECTS AND SAVE THEM

load_annotation(path = "/bfd/shared/annotations/human/gencode25/annotation_SaTAnn_RiboseQC/gencode.v25.annotation.gtf_Rannot")

ddx3x_txs<-GTF_annotation$trann$transcript_id[which(GTF_annotation$trann$gene_name=="DDX3X")]
ddx3_cds_tx<-cds_tx[names(cds_tx)%in%ddx3x_txs]

seqqq<-extractTranscriptSeqs(genome_seq,ddx3_cds_tx[1:4])
translate(seqqq)



GTF_annotation$cds_txs_coords

#read BAMfile in one chunk
opts <- BamFile(file=, yieldSize=5000000)
param <- ScanBamParam(flag=scanBamFlag(isDuplicate=FALSE,isSecondaryAlignment=FALSE),what=c("mapq"),tag = "MD")

x<-readGAlignments(opts,param = param)
maps_gene<-assay(summarizeOverlaps(x,features = GTF_annotation$genes,ignore.strand=F,mode="Union",inter.feature=FALSE))
maps_gene_cds<-assay(summarizeOverlaps(x,features = GTF_annotation$cds_genes,ignore.strand=F,mode="Union",inter.feature=FALSE))

#GRangesList

regions <- list(reduce(unlist(GTF_annotation$cds_genes)), GTF_annotation$fiveutrs, GTF_annotation$threeutrs,
                GTF_annotation$ncIsof, GTF_annotation$ncRNAs, GTF_annotation$introns, GTF_annotation$intergenicRegions)
names(regions) <- c("cds","fiveutrs","threeutrs",
                    "ncIsof","ncRNAs","introns","intergenic")
regions<-GRangesList(endoapply(regions,function(x){mcols(x)<-NULL;x}))

maps_all<-assay(summarizeOverlaps(reads=x,features=regions,ignore.strand=F,mode="Union",inter.feature=FALSE))
maps_uniq<-assay(summarizeOverlaps(reads=x[mcols(x)$mapq>50],features=regions,ignore.strand=F,mode="Union",inter.feature=FALSE))
barplot(t(maps_all))
barplot(t(maps_uniq))

summarizeJunctions(x)


most_cnts<-rownames(maps_gene_cds)[which.max(maps_gene_cds)]
cover<-coverage(x)
plot(cover[GTF_annotation$genes[most_cnts]][[1]],type="h")
pcd<-GTF_annotation$trann$transcript_biotype=="protein_coding"
txs_coding_most_cnts<-GTF_annotation$trann$transcript_id[pcd][GTF_annotation$trann$gene_id[pcd]==most_cnts]

#From genomic to transcript coordinates
tx_coords<-GTF_annotation$cds_txs_coords[seqnames(GTF_annotation$cds_txs_coords)=="ENST00000234590.8"]

map_tx<-mapToTranscripts(GRanges(x),transcripts = GTF_annotation$exons_txs[txs_coding_most_cnts][1])
cov_tx<-coverage(map_tx)
plot(cov_tx[["ENST00000234590.8"]],type="h")

abline(v = c(start(tx_coords),end(tx_coords)),col="red",lty=2)


load("/bfd/lcalviel/data/SaTAnn/totake_sep5/SaTAnn_Figures_script/Jun18/list_toSat_Mar5news")

names(list_tosat)
par(mfrow=c(2,3))
for(i in names(list_tosat)){
    y<-list_tosat[[i]]$P_sites_all
    y
    map_tx<-mapToTranscripts(y,transcripts = GTF_annotation$exons_txs[txs_coding_most_cnts][1])
    map_tx$score<-y$score[map_tx$xHits]
    cov_tx<-coverage(map_tx,weight = map_tx$score)
    plot(cov_tx[["ENST00000234590.8"]],type="h",col=c("red","blue","forestgreen"),main=i)
    abline(v = c(start(tx_coords),end(tx_coords)),col="black",lty=2)
    
}

#Single nt resolution!


#for things like par-clip:
load(example"mut_rl_Robj")
rds<-res_mut_rl$reads
rd<-rds[sum(rds$conversion=="T>C")>0]
plot(coverage(rds)[[1]],type="h")
plot(coverage(rd)[[1]],type="h")



