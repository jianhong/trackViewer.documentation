library(trackViewer) #load package
library(Biostrings)
library(motifStack)
## motif was downloaded from cis-bp
motif <- importMatrix("TET1.PWM.txt", format = "cisbp") 
pwm <- pfm2pwm(motif[[1]])
fa <- readDNAStringSet("FMR1.ups3K.dws3K.fasta")
names(fa) <- "FMR1"
TET1.binding.sites.v <- matchPWM(pwm, subject = fa[[1]], 
                                 min.score = "95%", with.score = TRUE)
TET1.binding.sites <- shift(ranges(TET1.binding.sites.v), shift = 147908950)
TET1.binding.sites <- GRanges("chrX", TET1.binding.sites, strand = "+", 
                              score=mcols(TET1.binding.sites.v)$score)
width(TET1.binding.sites) <- 1
TET1.binding.sites$border <- "gray80"
TET1.binding.sites$color <- 3
## set lollipop plot type to pin.
TET1.binding.sites$type <- "pin"
TET1.binding.sites <- new("track", dat=TET1.binding.sites, 
                          type="lollipopData")
gr <- GRanges("chrX", IRanges(147910500, 147914000))
## extract transcripts in the range
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)
FMR1 <- geneTrack(get("FMR1", org.Hs.egSYMBOL2EG), 
                  TxDb.Hsapiens.UCSC.hg38.knownGene)[[1]]
FMR1$dat2 <- 
  GRanges("chrX", 
          IRanges(c(147911604, 147911617, 147911727, 147911743, 147911758, 
                    147911768, 147911810, 147911821, 147911854, 147911877,
                    147911882, 147911902, 147911963), width = 1, 
                  names = c("AP2", "UBP1", "Sp1", "Sp1", "NRF1", 
                            "Sp1", "AGP", "NRF1", "Sp1", "AP2", 
                            "Sp1-like", "Myc", "Zeste")))
FMR1$dat2$color <- as.numeric(factor(names(FMR1$dat2)))
FMR1$dat2$border <- "gray"
## set lollipop label parameter.
FMR1$dat2$label.parameter.rot <- 45
## add methylation counts
maxX <- GRanges("chrX", IRanges(147911550, width=1), score=9, 
                color="white", border="white")
FX52_mock_methy <- 
  GRanges("chrX", IRanges(147911556+seq.int(35)*4, width=1),
          score=c(8, 8, 8, 7, 8, 9, 9, 9, 8, 9, 9, 9, 9, 8, 9,
                  9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9,
                  9, 9, 9, 8, 9), color=5, border="gray")
FX52_mock_methy <- new("track", dat=c(FX52_mock_methy, maxX), 
                       type="lollipopData")
FX52_dC_T_methy <- 
  GRanges("chrX", IRanges(147911556+seq.int(35)*4, width=1),
          score=c(3, 5, 3, 2, 0, 1, 1, 2, 0, 1, 1, 0, 0, 1, 0,
                  0, 2, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,
                  0, 1, 0, 0, 1), color=4, border="gray")
FX52_dC_T_methy$color[FX52_dC_T_methy$score==0] <- "white"
FX52_dC_T_methy <- new("track", dat=c(FX52_dC_T_methy, maxX), 
                       type="lollipopData")
## import RNA-seq tracks
iPSC_dC_T.RNAseq <- 
  importScore("iPSC_dC-T.bw", ranges = gr, format = "BigWig") 
iPSC_mock.RNAseq <- 
  importScore("iPSC_mock.bw", ranges = gr, format = "BigWig")
## import ChIP-Bisulfite-Seq tracks, two data in one track
iPSC_dC_T.BSseq <- 
  importScore("IPSC_dC-T.methy.bedgraph",
              "IPSC_dC-T.demethy.bedgraph",
              ranges = gr, format = "bedGraph")
iPSC_dC_dT.BSseq <- 
  importScore("IPSC_dC-dT.methy.bedgraph", 
              "IPSC_dC-dT.demethy.bedgraph", 
              ranges = gr, format = "bedGraph")
##stronger the signals
width(iPSC_dC_T.BSseq$dat) <- 
  width(iPSC_dC_T.BSseq$dat) + 1
width(iPSC_dC_dT.BSseq$dat) <- 
  width(iPSC_dC_dT.BSseq$dat) + 1
## interactive plot
browseTracks(trackList(FMR1, TET1.binding.sites,
                       FX52_dC_T_methy, FX52_mock_methy,
                       iPSC_dC_dT.BSseq, iPSC_dC_T.BSseq, 
                       iPSC_mock.RNAseq, iPSC_dC_T.RNAseq, 
                       heightDist=c(2, 1, 1, 1, 1, 1, 1, 1)), 
             gr)
## guildline at 147911556, 147911695, 147912052, 147912111