library(trackViewer) # load package
library(VariantAnnotation)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(org.Hs.eg.db)
library(rtracklayer)
fl <- system.file("extdata", "chr22.vcf.gz", package="VariantAnnotation")
## set the track range
gr <- GRanges("22", IRanges(50968014, 50970514, names="TYMP"))
## read in vcf file
tab <- TabixFile(fl)
vcf <- readVcf(fl, "hg19", param=gr)
## get GRanges from VCF object 
mutation.frequency <- rowRanges(vcf)
## keep the metadata
mcols(mutation.frequency) <- 
  cbind(mcols(mutation.frequency), 
        VariantAnnotation::info(vcf))
## set colors
mutation.frequency$border <- "gray30"
mutation.frequency$color <-
  ifelse(grepl("^rs", names(mutation.frequency)), 
         "lightcyan", "lavender")
## plot Global Allele Frequency based on AC/AN
mutation.frequency$score <- round(mutation.frequency$AF*100)
## change the color for score==0;
mutation.frequency$color[mutation.frequency$score==0] <- "#FFFFFF"
#change it to 1 to make the color changable.
mutation.frequency$score[mutation.frequency$score==0] <- 1
## keep sequence level style same
seqlevelsStyle(gr) <- seqlevelsStyle(mutation.frequency) <- "UCSC"
## extract transcripts in the range
trs <- geneModelFromTxdb(TxDb.Hsapiens.UCSC.hg19.knownGene, 
                         org.Hs.eg.db, gr=gr)
## subset the features to show the interested transcripts only
features <- GRangesList(trs[[1]]$dat, trs[[5]]$dat, trs[[6]]$dat)
flen <- elementNROWS(features)
features <- unlist(features)
## define the feature track layers
features$featureLayerID <- rep(1:2, c(sum(flen[-3]), flen[3]))
## define the feature labels
names(features) <- features$symbol
## define the feature colors
features$fill <- rep(c("lightblue", "mistyrose", "mistyrose"), flen)
## define the feature heights
features$height <- ifelse(features$feature=="CDS", .04, .02)
## import methylation data from a bed file
methy <- import(system.file("extdata", "methy.bed", package="trackViewer"), "BED")
## subset the data
methy <- methy[methy$score > 20]
## simulate multiple patients
rand.id <- sample.int(length(methy), 3*length(methy), replace=TRUE)
rand.id <- sort(rand.id)
methy.mul.patient <- methy[rand.id]
## pie.stack require metadata "stack.factor", and the metadata can not be 
## stack.factor.order or stack.factor.first
len.max <- max(table(rand.id))
stack.factors <- paste0("patient", 
                        formatC(1:len.max, 
                                width=nchar(as.character(len.max)), 
                                flag="0"))
methy.mul.patient$stack.factor <- 
  unlist(lapply(table(rand.id), sample, x=stack.factors))
methy.mul.patient$score <- 
  sample.int(100, length(methy.mul.patient), replace=TRUE)
## for a pie plot, two or more numeric meta-columns are required.
methy.mul.patient$score2 <- 100 - methy.mul.patient$score
## set different color set for different patient
patient.color.set <- as.list(as.data.frame(rbind(rainbow(length(stack.factors)), 
                                                 "#FFFFFFFF"), 
                                           stringsAsFactors=FALSE))
names(patient.color.set) <- stack.factors
methy.mul.patient$color <- 
  patient.color.set[methy.mul.patient$stack.factor]
## lollipop plot
mutations = new("track", dat=features, dat2=mutation.frequency, type="transcript")
methylations = new("track", dat=features, dat2=methy.mul.patient, type="transcript")
browseTracks(trackList(mutations, methylations), 
          gr=gr)
