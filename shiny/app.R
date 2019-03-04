#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(trackViewer)

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
   global <- reactiveValues(refresh = FALSE, fileIndex=0)
   observe({
     if(input$refresh) isolate(global$refresh <- TRUE)
   })
   output$trackViewer <- renderbrowseTracks({
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
   
   observeEvent(input$add, {
     isolate(global$fileIndex <- global$fileIndex+1)
     insertUI(selector = "#add",
              where = "beforeBegin",
              ui = tagList(
                selectInput(paste0("file", global$fileIndex), label="select file",
                            choices = dir(datafolder), multiple = FALSE),
                selectInput(paste0("format", global$fileIndex), label="file format", choices = c("BED", "bedGraph", "BigWig")),
                textInput(paste0("sample", global$fileIndex), label = "sample name", value = ""),
                tags$hr()
              ))
   })
}

# Run the application 
shinyApp(ui = ui, server = server)

