#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(rtracklayer)
library(VariantAnnotation)
library(trackViewer)
# trackViewer version must be no less than 1.19.10
stopifnot(packageVersion("trackViewer")>="1.19.10")

## data folder is used to save the data for ploting
datafolder <- "data"
## installed.packages
ip <- as.data.frame(installed.packages()[, c(1, 3:4)])
ip <- unique(as.character(ip$Package))
txdbs <- ip[grepl("^TxDb", ip)]
orgs <- ip[grepl("^org", ip)]

# Define UI for trackViewer
ui <- fluidPage(
   
   # Application title
   titlePanel("trackViewer"),
   
   # Sidebar with a slider input for number of bins 
   sidebarLayout(
      sidebarPanel(
         selectInput("TxDb", label="select TxDb package",
                     choices = txdbs),
         selectInput("org", label="select Org package", 
                     choices = orgs),
         
         textInput("chr", label = "chromosome", value = "chr3"),
         numericInput("start", label = "start", value = 108476000),
         numericInput("end", label = "end", value = 108485000),
         
         actionButton("add", "add file"),
         actionButton("lolli", "add lollipop plot track"),
         tags$hr(),
         actionButton("refresh", label="apply change")
      ),
      
      # Show a plot of the generated distribution
      mainPanel(
        browseTracksOutput("trackViewer")
      )
   )
)

# Define server logic required to draw tracks
server <- function(input, output) {
   global <- reactiveValues(refresh = FALSE, fileIndex=0, lolliIndex=0)
   observe({
     if(input$refresh){
       isolate(global$refresh <- TRUE)
     } 
   })
   plot <- renderbrowseTracks({
     if(!global$refresh) return()
     gr <- GRanges(input$chr, IRanges(as.numeric(input$start), as.numeric(input$end)))
     require(input$TxDb, character.only = TRUE)
     require(input$org, character.only = TRUE)
     trs <- tryCatch(geneModelFromTxdb(get(input$TxDb), get(input$org), gr=gr), 
                     error = function(e){ NULL })
     tks <- list()
     if(global$fileIndex>0){
       for(i in seq.int(global$fileIndex)){
         tks[[input[[paste0("sample", i)]]]] <- 
           tryCatch( importScore(file = file.path(datafolder, input[[paste0("file", i)]]),
                                 format = input[[paste0("format", i)]],
                                 ranges = gr),
                     error = function(e){ NULL})
       }
     }
     if(global$lolliIndex>0){
       for(i in seq.int(global$lolliIndex)){
         if(input[[paste0("lolliformat", i)]]=="VCF"){
           tab <- TabixFile(file.path(datafolder, input[[paste0("lollifile", i)]]))
           gen <- strsplit(input$TxDb, "\\.")[[1]]
           gen <- gen[length(gen)-1]
           vcf <- readVcf(fl, genome=gen, param = gr)
           mutation.frequency <- rowRanges(vcf)
           mcols(mutation.frequency) <- cbind(mcols(mutation.frequency), 
                                              VariantAnnotation::info(vcf))
           if(length(mutation.frequency$score)==0){
             mutation.frequency$score <- mutation.frequency$AF*100
           }
         }else{
           mutation.frequency <- import(con = file.path(datafolder, input[[paste0("lollifile", i)]]),
                                        format = input[[paste0("lolliformat", i)]],
                                        which = gr)
           if(length(mutation.frequency$itemRgb)>0) mutation.frequency$color <- mutation.frequency$itemRgb
           if(length(mutation.frequency$score)==0){
             mutation.frequency$score <- 1
           }
         }
         mutation.frequency <- promoters(mutation.frequency, upstream = 0, downstream = 1)
         seqlevelsStyle(mutation.frequency) <- "UCSC"
         thislolli <- new("track", dat=mutation.frequency, 
                          name=input[[paste0("lollisample", i)]], 
                          type="lollipopData")
         if(input[[paste0("lolliradio", i)]]!="none"){
           if(input[[paste0("lolliradio", i)]]=="default"){
             if(length(trs)>0){
               thislolli$type <- "transcript"
               thislolli$dat2 <- thislolli$dat
               ## get the features granges
               features <- unlist(GRangesList(lapply(trs, function(.e) .e$dat)))
               lens <- sapply(trs, function(.e) length(.e$dat))
               ## define the feature legend name
               names(features) <- rep(names(trs), lens)
               ## define the feature track layers
               features$featureLayerID <- rep(seq_along(trs), lens)
               ## define the feature colors
               features$fill <- rep(seq_along(trs), lens)
               ## define the feature heights
               features$height <- ifelse(features$feature=="CDS", 0.02, 0.01)
               thislolli$dat <- features
             }
           }else{
             thislolli$type <- "transcript"
             thislolli$dat2 <- thislolli$dat
             features <- import(file.path(datafolder, input[[paste0("lollitxfile", i)]]),
                                format = input[[paste0("lollitxformat", i)]],
                                which = gr)
             if(length(features$feature)<1){
               features$feature <- "CDS"
             }
             features$featureLayerID <- as.numeric(levels(features$gene))
             thislolli$dat <- features
           }
         }
       }
       tks[[input[[paste0("lollisample", i)]]]] <- thislolli
     }
     trackList <- list()
     if(length(trs)>0 && length(tks)==0){
       trackList <- trackList(trs)
     }
     if(length(trs)==0 && length(tks)>=0){
       trackList <- trackList(tks)
     }
     if(length(trs)>0 && length(tks)>0){
       trackList <- trackList(trs, tks, heightDist = c(length(trs), length(tks)))
     }
     
     if(length(trackList)>0){
       optSty <- optimizeStyle(trackList, theme="col")
       trackList <- optSty$tracks
       viewerStyle <- optSty$style
       browseTracks(trackList, gr=gr)
     }
   })
   output$trackViewer <- plot
   
   observeEvent(input$add, {
     isolate(global$fileIndex <- global$fileIndex+1)
     insertUI(selector = "#add",
              where = "beforeBegin",
              ui = tagList(
                tags$h4("Add data track from file"),
                selectInput(paste0("file", global$fileIndex), label="select file",
                            choices = dir(datafolder), multiple = FALSE),
                selectInput(paste0("format", global$fileIndex), label="file format", choices = c("BED", "bedGraph", "BigWig")),
                textInput(paste0("sample", global$fileIndex), label = "sample name", value = ""),
                tags$hr()
              ))
   })
   observeEvent(input$lolli, {
     isolate(global$lolliIndex <- global$lolliIndex+1)
     insertUI(selector = "#add",
              where = "beforeBegin",
              ui = tagList(
                tags$h4("Add lollipop plot track from file"),
                selectInput(paste0("lollifile", global$lolliIndex), label="select file",
                            choices = dir(datafolder), multiple = FALSE),
                selectInput(paste0("lolliformat", global$lolliIndex), label="file format", choices = c("BED", "bedGraph", "VCF")),
                textInput(paste0("lollisample", global$lolliIndex), label = "sample name", value = ""),
                radioButtons(paste0("lolliradio", global$lolliIndex), label = "gene model",
                             choices = c("none"="none",
                                         "use TxDb package"="default",
                                         "load from following file"="file")),
                selectInput(paste0("lollitxfile", global$lolliIndex), label="select transcript file",
                            choices = dir(datafolder), multiple = FALSE),
                selectInput(paste0("lollitxformat", global$lolliIndex), label="transcript file format", choices = c("GFF", "BED", "bedGraph")),
                tags$hr()
              ))
   })
}

# Run the application 
shinyApp(ui = ui, server = server)

