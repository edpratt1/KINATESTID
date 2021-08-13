get_manualcandidateaa <- function(fisher_candidates){
  subset_candidates <- runApp(manualcandidateaa_app(fisher_candidates))
  return(subset_candidates)
}

manualcandidateaa_app <- function(fisher_candidates){
  ui <- fluidPage(
    title = 'Select Table Rows',
    h1('A Server-side Table'),
    fluidRow(
      column(9, DT::dataTableOutput('x3'),
             h4("Selected Amino Acids"),
             verbatimTextOutput('y3'),
             h4("No. Permutations (Max 200,000)"),
             verbatimTextOutput('x4', placeholder = TRUE),
             actionButton('done', "Done"),
             tags$style(type="text/css",
                        ".shiny-output-error { visibility: hidden; }",
                        ".shiny-output-error:before { visibility: hidden; }"
             ),
      )
    )
  )
  server <- function(input, output, session) {
    options(shiny.sanitize.errors = TRUE)
    dt <- fisher_candidates[order(match(flank_pos, aa_cols))]
    perm_log <- reactiveVal("")
    
    output$x3 = DT::renderDataTable(
      dt, server = TRUE,
      selection = list(target = "row", selected = which(dt[, flank_pos] == 0,
                                                        arr.ind = TRUE)))

    output$y3 = renderPrint(dt[input$x3_rows_selected, ]$barcode)

    observeEvent(input$x3_rows_selected,{
      perm <- prod(dt[input$x3_rows_selected][, length(amino_acid), 
                                              by = flank_pos][, V1])
      perm_log(perm)
    })
    
    output$x4 = renderPrint(perm_log())

    return_values <- reactive({dt[input$x3_rows_selected, ]})
    observe({
      if (input$done > 0){
        stopApp(return_values())
      }
    })
  }
  shinyApp(ui, server)
}