library(trackViewer) #load package
## loading data.
data <- read.delim("IARC-TP53/datasets/somaticMutationDataIARC TP53 Database, R19.txt", 
                   stringsAsFactors = FALSE)
data <- data[data$Morphology %in% "Small cell carcinoma, NOS" & 
               data$Effect %in% c("nonsense", "missense", "silent", "intronic"), 
             c("hg38_Chr17_coordinates", "g_description_GRCh38", "Effect")]
counts <- table(data$g_description_GRCh38)
data$counts <- as.numeric(counts[data$g_description_GRCh38])
data <- unique(data)
## prepare a GRanges object for mutations
snps <- with(data, GRanges("chr17", IRanges(hg38_Chr17_coordinates, width=1), 
                           effect=factor(Effect), score=counts))
## set the bristles head colors of the pappus by mutation types
snps$color <- as.numeric(snps$effect)+1
## set the beak color of dandelion seeds.
snps$border <- "gray"
## set plotting region
gr <- GRanges("chr17", IRanges(7669000, 7677000))
## extract transcripts in the range
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)
trs <- geneModelFromTxdb(TxDb.Hsapiens.UCSC.hg38.knownGene, 
                         org.Hs.eg.db, gr=gr)
## subset the features to show the interested transcripts only
features <- c(trs[[1]]$dat, trs[[3]]$dat, trs[[4]]$dat)
lens <- sapply(trs[c(1, 3, 4)], function(.e) length(.e$dat))
## define the feature legend name
names(features) <- rep(names(trs)[c(1, 3, 4)], lens)
## define the feature track layers
features$featureLayerID <- rep(seq.int(3), lens)
## define the feature colors
features$fill <- rep(c("lightblue", "mistyrose", "orange"), lens)
## define the feature heights
features$height <- ifelse(features$feature=="CDS", 0.02, 0.01)
snps$type="dandelion"
SNP <- new("track", dat=features, dat2=snps, type="transcript")
## interactive plot
browseTracks(trackList(SNP), gr=gr)
