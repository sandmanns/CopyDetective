shinyUI(fluidPage(
    theme = shinytheme("flatly"),
    titlePanel(div("CopyDetective",
               img(height = 55, width = 63, src = "CD.png"))
               ),
    sidebarPanel(
        shinyjs::useShinyjs(),
        tabsetPanel(
            tabPanel("Perform analysis",
                h4("Input:"),
                textInput('input_folder','Define input folder containing polymorphism information',
                          "/path/to/input/data/"),
                fileInput('sampleFile',label = "Upload sample names file"),
                h5("Note: a 2-column txt-file is required. The germline sample 
           names are defined in the first column, the matching tumor 
           sample names are defined in the second column (no header)."),
                hr(),
                h4("Output:"),
                textInput('output_folder','Define output folder',
                          "/path/to/output/"),
                checkboxGroupInput('output_files',"Select output files",
                                   choices = c("Raw CNV calls","Merged CNV calls","Filtered CNV calls"),
                                   selected = c("Raw CNV calls","Merged CNV calls","Filtered CNV calls"),inline = F),
                checkboxGroupInput('output_plots',"Select output plots",
                                   choices = c("Raw CNV calls (all)","Raw CNV calls (sig)","Merged CNV calls","Filtered CNV calls"),
                                   selected = c("Raw CNV calls (all)","Merged CNV calls"),inline=F),
                
                hr(),
                h4("1. Quality Analysis:"),
                radioButtons('simulation',label="Are detection thresholds already available?",
                             choices = c("No","Yes"),selected = "No",inline = T),
                uiOutput("simulationUI"),
                uiOutput("simulationUI2"),
                uiOutput("simulationUI3"),
                uiOutput("simulationUI4"),
                hr(),
                h4("2. CNV calling:"),
                sliderInput('aimSens',"Desired sensitivity",min=0,max=1,value=0.95),
                hr(),
                h4("3. Merging:"),
                numericInput('maxDist',"Maximum allowed distance between raw calls [bp]",
                             min=100000,max=59128983,value = 20000000),
                hr(),
                h4("4. Filtration:"),
                radioButtons('final_filter',label = "Perform final filtration?",
                             choices = c("No","Yes"),selected = "Yes",inline = T),
                uiOutput("filtration_thresholdUI"),
                hr(),
                actionButton('do',"Start analysis",class = "btn-primary")
            ),
            tabPanel("Display results",
                     shinyjs::useShinyjs(),
                radioButtons("select_samples",label="Available samples",
                            choices = c(NA,rep(NULL,100))),
                radioButtons('output_files2',"Select CNV calls to display",
                                   choices = c("Merged CNV calls","Filtered CNV calls"),
                                   selected = c("Merged CNV calls"),inline = F),
                checkboxGroupInput('output_plots2',"Select plots to display",
                                   choices = c("Raw CNV calls (all)","Raw CNV calls (sig)","Merged CNV calls","Filtered CNV calls"),
                                   selected = c("Raw CNV calls (all)","Merged CNV calls"),inline=F),
                
                hr(),
                actionButton('do2',"Display results",class = "btn-primary")
            )
        )

        ),
    mainPanel(
        shinyjs::useShinyjs(),
        tabsetPanel(
            tabPanel("Log",
                     div(id = "text")
                     ),
            tabPanel("Detection thresholds",
                     DT::dataTableOutput('table_dt')
                     ),
            tabPanel("CNV calls",
                     h3(textOutput("sample")),
                     hr(),
                     DT::dataTableOutput('table_cnvs')
                     ),
            tabPanel("Plots",
                     h3(textOutput("sample2")),
                     h4(textOutput("text_plot1")),
                     imageOutput("plot1"),
                     hr(),
                     h4(textOutput("text_plot2")),
                     imageOutput("plot2"),
                     hr(),
                     h4(textOutput("text_plot3")),
                     imageOutput("plot3"),
                     hr(),
                     h4(textOutput("text_plot4")),
                     imageOutput("plot4")
            )
        )
    )
    )
)




