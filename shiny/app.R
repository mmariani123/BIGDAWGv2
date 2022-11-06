#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(shinydashboard)
library(shinyWidgets)
library(shinyFiles)
library(stringr)
#library(devtools)
library(shinyBS) #Tooltip doesn't seem to work with bsButton or action button
library(BIGDAWG)
library(BIGDAWGv2)
#library(shiny.fluent)
#devtools::load_all()

#human
#bola (bovine)
#chicken
#dog

ui <- dashboardPage(
  dashboardHeader(title = "BIGDAWGv2"),
  dashboardSidebar(
    #tags$style('
    #    body{
    #      background-color: white;
    #    }
    #   div{
    #   background-color: black;
    #      color: white;
    #      border-radius: 5px;
    #      margin: 0px;
    #      padding: 5px;
    #   }'
    #),
    #tags$head(tags$script(src = "message-handler.js")),

    #ui <- fluidPage(
    #  tags$style(
    #    ".first-p {
    #  color: red;
    #}
    ##element {
    #  color: red;
    #}
    #"
    #  ),
    #  p(class = "first-p", "Hello World"),
    #  p("Another text"),
    #  div(id = "element", "A block")
    #)

    tags$style(".span_1 {margin-left:15px}"),
    tags$style(".d_button {margin-left:15px; margin-bottom:10px;}"),
    tags$style(".text_entry {margin-top:-10px;}"),
    tags$style("#select {margin-top:-40px;}"),
    tags$style("#dropdown {margin-top:-15px;}"),
    tags$style("#run_button {margin-top:20px;}"),

    div(id="run_button",
      #actionButton("Run.Bigdawg", "Run BIGDAWG"),
      tags$div(
      title = paste0("Click this button to run BIGDAWG ",
                     "after you have uploaded your desired input file ",
                     "and selected your parameters below"),
      shinyBS::bsButton(inputId="Run.Bigdawg", label="Run BIGDAWG"),
      bsTooltip(id="div",
                title=paste0("Click this button to run BIGDAWG ",
                             "after you have uploaded your desired input file ",
                             "and selected your parameters below"),
                             placement = "bottom",
                             trigger = "hover",
                             options = NULL)
      )
    ),
    div(
      fileInput("Data",
                "Choose data file",
                accept = c(
                  "text",
                  "text/csv",
                  "text/comma-separated-values,text/plain",
                  ".csv")),
                bsTooltip(id="Data",
                          title=paste0("Choose your input data file ",
                          "(make sure it is in the proper format - ",
                          "see documentation for details)"),
                          placement = "bottom",
                          trigger = "hover",
                          options = NULL)
    ),
    div(id="select",
      selectInput("speciesSelect", "Select Species",
                  choices = c("HUMAN", "DOG", "CHICKEN", "COW"),
                  multiple=F,
                  selectize=T,
                  selected="HUMAN",
                  width = '98%'),
      bsTooltip(id="speciesSelect",
                title="Select your preferred species",
                placement = "bottom",
                trigger = "hover",
                options = NULL)
    ),
    div(
      span("HLA",class="span_1"),
      shinyWidgets::prettyToggle(
        inputId = "HLA",
        label_on = "TRUE",
        label_off = "FALSE",
        value = FALSE,
        outline = FALSE,
        plain = FALSE,
        #icon_on = icon("thumbs-up"),
        #icon_off = icon("thumbs-down"),
        animation="smooth"
      ),
      bsTooltip(id="HLA",
                title=paste0("Click this button to run BIGDAWG ",
                             "after you have uploaded your desired input file ",
                             "and selected your parameters below"),
                placement = "bottom",
                trigger = "hover",
                options = NULL)
    ),
    div(id="dropdown",
      selectInput("Run.Tests",
                "Select Test(s)",
                choices = c("HWE","H","L","A"),
                multiple=T,
                selectize=T,
                width = '98%'),
      bsTooltip(id="Run.Tests",
                title=paste0("Click this button to run BIGDAWG ",
                             "after you have uploaded your desired input file ",
                             "and selected your parameters below"),
                placement = "bottom",
                trigger = "hover",
                options = NULL)
    ),
    div(class="text_entry",
              textInput("Loci.Set",
              label = "Loci.Set",
              value = 'A;DRB1',
              #value = list("A"),
              width = NULL,
              placeholder = NULL),
              bsTooltip(id="Loci.Set",
                        title=paste0("Input list defining which loci to use ",
                                     "for analyses"),
                        placement = "bottom",
                        trigger = "hover",
                        options = NULL)
              #verbatimTextOutput("Loci.Set.value")
    ),
    div(class="text_entry",
              textInput("Exon",
              label="Exon",
              value="3,5,6",
              width = NULL,
              placeholder = NULL),
              #verbatimTextOutput("Exon.value")
              bsTooltip(id="Exon",
                        title=paste0("A single numeric or numeric vector that ",
                              "defines exons to target in the amino acid ",
                              "analysis."),
                        placement = "bottom",
                        trigger = "hover",
                        options = NULL)
    ),
    div(
      span("All.Pairwise", class="span_1"),
      shinyWidgets::prettyToggle(
        inputId = "All.Pairwise",
        label_on = "TRUE",
        label_off = "FALSE",
        value = FALSE,
        outline = FALSE,
        plain = FALSE,
        animation="smooth"),
      bsTooltip(id="All.Pairwise",
                title=paste0("Should pairwise combinations of loci be run in ",
                             "the haplotype analysis? Only relevant to ",
                             "haplotype analysis."),
                placement = "bottom",
                trigger = "hover",
                options = NULL)
    ),
    div(
      span("Trim", class="span_1"),
      shinyWidgets::prettyToggle(
        inputId = "Trim",
        label_on = "TRUE",
        label_off = "FALSE",
        value = FALSE,
        outline = FALSE,
        plain = FALSE,
        animation="smooth"),
      bsTooltip(id="Trim",
                title=paste0("Flags whether or not to Trim HLA alleles ",
                             "to a specified resolution. Should not be ",
                             "optioned for data that does not conform to ",
                             "HLA naming conventions."),
                placement = "bottom",
                trigger = "hover",
                options = NULL)
    ),
    div(class="text_entry",
              textInput("Res",
              label="Res",
              value = "2",
              width = NULL,
              placeholder = NULL),
              bsTooltip(id="Res",
                        title=paste0("Sets the desired resolution when ",
                                     "trimming HLA alleles."),
                        placement = "bottom",
                        trigger = "hover",
                        options = NULL),
              #verbatimTextOutput("Res.value")
    ),
    div(
      span("EVS.rm", class="span_1"),
      shinyWidgets::prettyToggle(
        inputId = "EVS.rm",
        label_on = "TRUE",
        label_off = "FALSE",
        value = FALSE,
        outline = FALSE,
        plain = FALSE,
        animation="smooth"),
      bsTooltip(id="EVS.rm",
                title=paste0("Flags whether or not to strip expression ",
                             "variant suffixes from HLA alleles. ",
                             "Example: A*01:11N will convert to A*01:11."),
                placement = "bottom",
                trigger = "hover",
                options = NULL)
    ),
    div(class="text_entry",
              textInput("Missing",
              label="Missing",
              value = "0",
              width = NULL,
              placeholder = NULL),
              bsTooltip(id="Missing",
                        title=paste0("Sets the allowable per subject threshold",
                                     " for missing alleles. Relevant to ",
                                     "running the haplotype analysis."),
                        placement = "bottom",
                        trigger = "hover",
                        options = NULL)
              #verbatimTextOutput("Missing.value")
    ),
    div(
      span("Strict.Bin", class="span_1"),
      shinyWidgets::prettyToggle(
        inputId = "Strict.Bin",
        label_on = "TRUE",
        label_off = "FALSE",
        value = FALSE,
        outline = FALSE,
        plain = FALSE,
        animation="smooth"),
      bsTooltip(id="Strict.Bin",
                title=paste0("Sets whether strict binning should be used ",
                             "during Chi Square testing."),
                placement = "bottom",
                trigger = "hover",
                options = NULL)
    ),
    div(class="text_entry",
              textInput("Cores.Lim",
              label="Cores.Lim",
              value = "1L",
              width = NULL,
              placeholder = NULL),
              bsTooltip(id="Cores.Lim",
                        title=paste0("Specifies the number of cores ",
                                     "accessible by BIGDAWG in amino acid ",
                                     "analysis. Not relevant to Windows ",
                                     "operating systems which will use only ",
                                     "a single core"),
                        placement = "bottom",
                        trigger = "hover",
                        options = NULL)
              #verbatimTextOutput("Cores.Lim.value")
    ),
    #div(
    #               shinyFiles::shinyDirButton(
    #               id="Results.Dir",
    #               label="Select a folder",
    #               title="Please select a folder",
    #               multiple=FALSE)
    #),
    div(
      downloadButton("Results.Dir",
                     "Download Results",
                     class="d_button"),
      bsTooltip(id="Results.Dir",
                title=paste0("String name of a folder for BIGDAWG output. ",
                             "Subfolder for each locus set will be generated ",
                             "within any output folder specified."),
                placement = "bottom",
                trigger = "hover",
                options = NULL)
    ),
    div(
      span("Return", class="span_1"),
      shinyWidgets::prettyToggle(
        inputId = "Return",
        label_on = "TRUE",
        label_off = "FALSE",
        value = TRUE,
        outline = FALSE,
        plain = FALSE,
        animation="smooth"),
      bsTooltip(id="Return",
                title=paste0("Specifies if BIGDAWG should output analysis ",
                             "results to a specified object."),
                placement = "bottom",
                trigger = "hover",
                options = NULL)
    ),
    div(
      span("Output", class="span_1"),
      shinyWidgets::prettyToggle(
        inputId = "Output",
        label_on = "TRUE",
        label_off = "FALSE",
        value = TRUE,
        outline = FALSE,
        plain = FALSE,
        animation="smooth"),
      bsTooltip(id="Output",
                title=paste0("Turns on or off the writing of results to files.",
                             " The default will write all results to text ",
                             "files in the output directory."),
                placement = "bottom",
                trigger = "hover",
                options = NULL)
    ),
    div(
      span("Merge.Output", class="span_1"),
      shinyWidgets::prettyToggle(
        inputId = "Merge.Output",
        label_on = "TRUE",
        label_off = "FALSE",
        value = FALSE,
        outline = FALSE,
        plain = FALSE,
        animation="smooth"),
      bsTooltip(id="Merge.Output",
                title=paste0("Turns on or off the merging of all analysis ",
                             "results into single files labeled ",
                             "’Merged_xxxx.txt'."),
                placement = "bottom",
                trigger = "hover",
                options = NULL)
    ),
    div(
      span("Verbose", class="span_1"),
      shinyWidgets::prettyToggle(
        inputId = "Verbose",
        label_on = "TRUE",
        label_off = "FALSE",
        value = TRUE,
        outline = FALSE,
        plain = FALSE,
        animation="smooth"),
      bsTooltip(id="Verbose",
                title=paste0("Sets the levels of detail that should be ",
                             "displayed on the console. The default will ",
                             "display summaries of the analysis from each ",
                             "specified test. When turned off, only the ",
                             "completion status of each test is displayed."),
                placement = "bottom",
                trigger = "hover",
                options = NULL)
    ),
    div(
      downloadButton("Download.Test.Data",
                    "Download Test Data",
                    class="d_button"),
        bsTooltip(id="Download.Test.Data",
                  title=paste0("Download the HLA test data file (.txt)"),
                  placement = "bottom",
                  trigger = "hover",
                  options = NULL)
    )
    #HLA
    #Run.Tests,
    #Loci.Set,
    #Exon,
    #All.Pairwise     = FALSE,
    #Trim             = FALSE,
    #Res              = 2,
    #EVS.rm           = FALSE,
    #Missing          = 2,
    #Strict.Bin       = FALSE,
    #Cores.Lim        = 1L,
    #Results.Dir,
    #Return           = FALSE,
    #Output           = TRUE,
    #Merge.Output     = FALSE,
    #Verbose          = TRUE
  ),

  dashboardBody(
    # Listen for messages
    tags$script(
      type="text/javascript",
      "$(document).on('shiny:connected',
                     function(event){
        alert('Welcome to BIGDAWGv2!');
      });"
    ),

    #  "$(document).ready(function(){
    #    Shiny.addCustomMessageHandler(\'alert\',
    #      function(data){
    #            new Notification(data.msg);
    #      })
    #  });"

    tabBox(title=NULL,
           width=12,
           id = "tabset1",
           height = "1300px",
    tabPanel(tabName = "run.panel","Run",
    # Boxes need to be put in a row (or column)
    fluidRow(
      div(
      HTML("Welcome to BIGDAWGv2. Upload your input LA data file, fill
            out the parameters on the left, and hit the 'Run BIGDAWG' button.
            results will be available via download using the download button and
            output to the screen as well. The test data can also be downloaded
            (tab-delimited .txt format) by clicking the 'download test data'
            button in the parameters sections. For further information,
            consult the documentation. We hope that you enjoy using BIGDDAWG."),
            style="text-align:justify;
                   color:black;
                   margin-top:5px;
                   margin-right:20px;
                   margin-left:20px;"),
      tableOutput("contents")
    )
    ),
    tabPanel(tabName = "citation.panel", "Citation",
    fluidRow(
      div(
        HTML('Citation:<br><br>"Data sets and functions for chi-squared Hardy-Weinberg and case-control
        association tests of highly polymorphic genetic data [e.g., human leukocyte antigen
        (HLA) data]. Performs association tests at multiple levels of polymorphism
        (haplotype, locus and HLA amino-acids) as described in Pappas DJ, Marin W, Hollenbach
        JA, Mack SJ (2016) <doi:10.1016/j.humimm.2015.12.006>. Combines rare variants to a
        common class to account for sparse cells in tables as described by Hollenbach JA,
        Mack SJ, Thomson G, Gourraud PA (2012) <doi:10.1007/978-1-61779-842-9_14>."
             <br><br> *Additional updates in progress for BIGDAWGv2 by Michael P. Mariani and Richard Single'),
        style="text-align:justify;
               color:black;
               margin-left:20px;
               margin-right:20px;
               margin-top:10px;"
      )
    )
    ),
    tabPanel(tabName = "debug.panel", "Debug",
      fluidRow(infoBoxOutput("info.out"))
      )
    )
  )
)

# Define server logic required to draw a histogram
server <- function(input, output, session) {

#Javascripting:

  #observeEvent(eventExpr=input$Run.Bigdawg,
  #             handlerExpr={
  #              Sys.sleep(2);
  #              session$sendCustomMessage(
  #                type = 'alert',
  #                message = list(msg = "my message")
  #              )
  #            }
  #           )

    #observeEvent(input$Run.Bigdawg, {
    #  session$sendCustomMessage(type = 'testmessage',
    #                            message = 'Thank you for trying BIGDAWG')
    #})

    #observe({
    #
    #  loci.vals <- switch(input$speciesSelect,
    #         "HUMAN" = "A",
    #         "DOG" = "B",
    #         "CHICKEN" = "In Progress",
    #         "COW" = "In Progress")
    #
    #  updateSelectInput(
    #    session = session,
    #    inputId = "Loci",
    #    choices = loci.vals,
    #    selected = head(loci.vals, 1)
    #  )
    #
    #})

    output$info.out <- renderInfoBox({
      infoBox("Debug output below",
              input$debug.panel,
              icon = icon("info-circle")
              )
    })

    bigdawg.output <- eventReactive(input$Run.Bigdawg, {
      #if(input$HLA           &
      #   input$Run.Tests     &
      #   input$Loci.Set      &
      #   input$Exon          &
      #   input$All.Pairwise  &
      #   input$Trim          &
      #   input$Res           &
      #   input$EVS.rm        &
      #   input$Missing       &
      #   input$Strict.Bin    &
      #   input$Cores.Lim     &
      #   input$Results.Dir   &
      #   input$Return        &
      #   input$Output        &
      #   input$Merge.Output  &
      #   input$Verbose){
      print(input$Data$datapath)
      print(input$Loci.Set)
      bigdawg.results <- BIGDAWGv2::BIGDAWGv2(
        Data         = input$Data$datapath,
        HLA          = input$HLA,
        Run.Tests    = input$Run.Tests,
        Loci.Set     = strsplit(input$Loci.Set, split=";"),
        Exon         = {if(grepl(":",input$Exon)){
                       ends=as.numeric(unlist(strsplit(input$Exon, split=":")))
                       seq(from=ends[1],to=ends[length(ends)])
                       }else if(grepl(",",input$Exon)){
                       as.numeric(unlist(strsplit(input$Exon, split=",")))
                       }else{
                       input$Exon
                       }},
        All.Pairwise = input$All.Pairwise,
        Trim         = input$Trim,
        Res          = as.numeric(input$Res),
        EVS.rm       = input$EVS.rm,
        Missing      = as.numeric(input$Missing),
        Strict.Bin   = input$Strict.Bin,
        Cores.Lim    = as.integer(str_extract(input$Cores.Lim,"(\\d+)")),
        #Results.Dir  = input$Results.Dir,
        Return       = input$Return,
        Output       = input$Output,
        Merge.Output = input$Merge.Output,
        Verbose      = input$Verbose)
        return(bigdawg.results)

        #}
    })

    output$contents <- renderTable(bigdawg.output())

    df.in <- reactive({

      fileIn <- input$fileIn$datapath
      if(is.null(fileIn)){
        return(NULL)
      }else{
        #read.table(fileIn$datapath,
        #           header = TRUE,
        #           sep = "\t",
        #           stringsAsFactors = FALSE)
        return(fileIn)
      }

    })

    output$HLA.value <- renderText({
      sprintf("Value: %s", input$HLA)
    })

    #output$Loci.Set.value <- renderText({input$Loci.Set})

    #output$Exon.value <- renderText({input$Exon})

    output$All.Pairwise.value <- renderText({
      sprintf("Value: %s", input$All.Pairwise)
    })

    output$Trim.value <- renderText({
      sprintf("Value: %s", input$Trim)
    })

    output$Res.value <- renderText({input$Res})

    output$EVS.Remove.value <- renderText({
      sprintf("Value: %s", input$EVS.Remove)
    })

    output$Missing.value <- renderText({input$Missing})

    output$Strict.Bin.value <- renderText({
      sprintf("Value: %s", input$Strict.Bin)
    })

    output$Cores.Lim.value <- renderText({input$Cores.Lim})

    #shinyFiles::shinyDirChoose(input=input,
    #               id='Results.Dir',
    #               roots=c(wd='.'),
    #               filetypes=c('', 'txt')
    #               )
    #observe({
    #  print(input$Results.Dir.folder)
    #})

    #output$Results.Dir <- downloadHandler(
    #  filename = function(){
    #    paste0("bigdawg_output_", Sys.Date(), ".txt")
    #  },
    #  content = function(file){
    #    #writeLines(paste(text,
    #    #                 collapse = ", "),
    #    #           file)
    #    write.table(bigdawg.output(),
    #                file,
    #                col.names = FALSE,
    #                row.names = FALSE,
    #                quote = FALSE,
    #                sep="\t")
    #  }
    #)

    #Can be zipped as well:
    output$Results.Dir <- downloadHandler(
      filename <- function(){
        return("BIGDAWGv2_results.zip")
      },
      content <- function(file){
        output.folder <- grep("output ",
                              list.dirs(full.names = TRUE),
                              value=TRUE)
        print(output.folder)
        if(.Platform$OS.type=="windows"){
          system(paste0("powershell Compress-Archive ",
                        "'",
                        output.folder,
                        "'",
                        " ",
                        "'",
                        output.folder,
                        ".zip",
                        "'"))
        }else if(.Platform$OS.type=="linux"){
          system("zip -r ",
                 output.folder,
                 ".zip ",
                 output.folder)
        }
        file.copy(paste0(output.folder,".zip"), file)
        unlink(output.folder, recursive=TRUE)
        unlink(paste0(output.folder,".zip"), recursive=TRUE)
      },
      contentType = "application/zip"
    )

    output$Return.value <- renderText({
      sprintf("Value: %s", input$Return)
    })

    output$Output.value <- renderText({
      sprintf("Value: %s", input$Output)
    })

    output$Merge.Output.value <- renderText({
      sprintf("Value: %s", input$Merge.Output)
    })

    output$Verbose.value <- renderText({
      sprintf("Value: %s", input$Verbose)
    })

    output$Download.Test.Data <- downloadHandler(
      filename = function(){
        paste0("HLA_data.txt")
      },
      content = function(file){
        file.copy("www/HLA_data.txt", file)
      }
    )

    #Can be zipped as well:
    output$Download.Test.Data.Zip <- downloadHandler(
      filename <- function() {
        return("HLA_data.txt")
      },
      content <- function(file) {
        file.copy("HLA_data.txt.zip", file)
      },
      contentType = "application/zip"
    )

}

# Run the application
shinyApp(ui = ui, server = server)
