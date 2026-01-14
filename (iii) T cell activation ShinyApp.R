library(deSolve)              #Version 1.40
library(ggplot2)              #Version 3.5.1
library(scales)               #Version 1.3.0
library(viridis)              #Version 0.6.5
library(shiny)                #Version 1.10.0
library(shinyWidgets)         #Version 0.9.0
library(plotly)               #Version 4.10.4
library(shinythemes)          #Version 1.2.0
library(digest)               #Version 0.6.37
library(shinyjs)              #Version 2.1.0
library(openxlsx)             #Version 4.2.8
library(png)                  #Version 4.4.2
library(shinyscreenshot)      #Version 0.2.1



# Press Run App in the top right for the embedded picture to be displayed 

# The parameters displayed in the app differ from the code in their naming:
                           #"kon" = "k<sub>on</sub>",
                           #"koff" = "k<sub>off</sub>",
                           #"kp" = "k<sub>p</sub>",
                           #"phi" = "k<sub>lim</sub>",
                           #"Ymax" = "M<sub>max</sub>",
                           #"lambda" = "k<sub>amp</sub>",
                           #"delta" = "k<sub>act</sub>",
                           #"my" = "k<sub>inh</sub>")

setwd(dirname(rstudioapi::getSourceEditorContext()$path))


# Define initial values for sliders
initial_values <- list(
  phi = 0.01,
  lambda = 100000,
  delta = 50000,
  my = 50000,
  kon = 1,
  koff = 1,
  kp = 50,
  Pmax = 10,
  Ymax = 10,
  R = 50
)



#1 UI ----
ui <- {fluidPage(
  useShinyjs(),
  theme = shinytheme("simplex"),
  sidebarLayout(
    sidebarPanel(
      chooseSliderSkin(skin = "Flat"),
      fluidRow( 
        column(9, sliderTextInput(
          inputId = "phi",
          label = HTML('<span style="font-size: 14px;">k<sub>lim</sub>:</span>'),
          grid = TRUE,
          choices = c("N/A",seq(0.0001,0.001,0.0001),seq(0.001,0.01,0.001),seq(0.01,0.1,0.01),seq(0.1,1,0.1),seq(1,10,1),seq(10,100,10)),
          selected = initial_values$phi
        )),
        column(3, numericInput(
          inputId = "phi_value",
          label = HTML('<span style="font-size: 14px;">k<sub>lim</sub>:</span>'),
          value = initial_values$phi
        )),
      ),
      fluidRow(
        column(9, sliderTextInput(
          inputId = "lambda",
          label = HTML('<span style="font-size: 14px;">k<sub>amp</sub>:</span>'),
          grid = TRUE,
          choices = c("N/A",seq(1,10,1),seq(10,100,10),seq(100,1000,100),seq(1000,10000,1000),seq(10000,100000,10000),seq(100000,1000000,100000)),
          selected = initial_values$lambda
        )),
        column(3, numericInput(
          inputId = "lambda_value",
          label = HTML('<span style="font-size: 14px;">k<sub>amp</sub>:</span>'),
          value = initial_values$lambda
        ))
      ),
      fluidRow(
        column(9, sliderTextInput(
          inputId = "delta",
          label = HTML('<span style="font-size: 14px;">k<sub>act</sub>:</span>'),
          grid = TRUE,
          choices = c("N/A",seq(1,10,1),seq(10,100,10),seq(100,1000,100),seq(1000,10000,1000),seq(10000,100000,10000),seq(100000,1000000,100000)),
          selected = initial_values$delta
        )),
        column(3, numericInput(
          inputId = "delta_value",
          label = HTML('<span style="font-size: 14px;">k<sub>act</sub>:</span>'),
          value = initial_values$delta
        ))
      ),
      fluidRow(
        column(9, sliderTextInput(
          inputId = "my",
          label = HTML('<span style="font-size: 14px;">k<sub>inh</sub>:</span>'),
          grid = TRUE,
          choices = c("N/A",seq(1,10,1),seq(10,100,10),seq(100,1000,100),seq(1000,10000,1000),seq(10000,100000,10000),seq(100000,1000000,100000)),
          selected = initial_values$my
        )),
        column(3, numericInput(
          inputId = "my_value",
          label = HTML('<span style="font-size: 14px;">k<sub>inh</sub>:</span>'),
          value = initial_values$my
        ))
      ),
      fluidRow(
        column(9, sliderTextInput(
          inputId = "kon",
          label = HTML('<span style="font-size: 14px;">k<sub>on</sub>:</span>'),
          grid = TRUE,
          choices = c("N/A",seq(0.0001,0.001,0.0001),seq(0.001,0.01,0.001),seq(0.01,0.1,0.01),seq(0.1,1,0.1),seq(1,10,1),seq(10,100,10)),
          selected = initial_values$kon
        )),
        column(3, numericInput(
          inputId = "kon_value",
          label = HTML('<span style="font-size: 14px;">k<sub>on</sub>:</span>'),
          value = 10^initial_values$kon
        ))
      ),
      fluidRow(
        column(9, sliderTextInput(
          inputId = "koff",
          label = HTML('<span style="font-size: 14px;">k<sub>off</sub>:</span>'),
          grid = TRUE,
          choices = c("N/A",seq(0.0001,0.001,0.0001),seq(0.001,0.01,0.001),seq(0.01,0.1,0.01),seq(0.1,1,0.1),seq(1,10,1),seq(10,100,10)),
          selected = initial_values$koff
        )),
        column(3, numericInput(
          inputId = "koff_value",
          label = HTML('<span style="font-size: 14px;">k<sub>off</sub>:</span>'),
          value = 10^initial_values$koff
        ))
      ),
      fluidRow(
        column(9, sliderTextInput(
          inputId = "kp",
          label = HTML('<span style="font-size: 14px;">k<sub>p</sub>:</span>'),
          grid = TRUE,
          choices = c("N/A",seq(0.0001,0.001,0.0001),seq(0.001,0.01,0.001),seq(0.01,0.1,0.01),seq(0.1,1,0.1),seq(1,10,1),seq(10,100,10)),
          selected = initial_values$kp
        )),
        column(3, numericInput(
          inputId = "kp_value",
          label = HTML('<span style="font-size: 14px;">k<sub>p</sub>:</span>'),
          value = initial_values$kp
        ))
      ),
      fluidRow(
        column(9, sliderTextInput(
          inputId = "Pmax",
          label = HTML('<span style="font-size: 14px;">E<sub>max</sub>:</span>'),
          grid = TRUE,
          choices = c("N/A",seq(0.1, 0.9, 0.1), seq(1, 9), seq(10, 99, 10), seq(100, 1000, 100)),
          selected = initial_values$Pmax
        )),
        column(3, numericInput(
          inputId = "Pmax_value",
          label = HTML('<span style="font-size: 14px;">E<sub>max</sub>:</span>'),
          value = initial_values$Pmax
        ))
      ),
      fluidRow(
        column(9, sliderTextInput(
          inputId = "Ymax",
          label = HTML('<span style="font-size: 14px;">M<sub>max</sub>:</span>'),
          grid = TRUE,
          choices = c("N/A",seq(0.1, 0.9, 0.1), seq(1, 9), seq(10, 99, 10), seq(100, 1000, 100)),
          selected = initial_values$Ymax
        )),
        column(3, numericInput(
          inputId = "Ymax_value",
          label = HTML('<span style="font-size: 14px;">M<sub>max</sub>:</span>'),
          value = initial_values$Ymax
        ))
      ),
      fluidRow(
        column(9, sliderTextInput(
          inputId = "R",
          label = HTML('<span style="font-size: 14px;">TCR:</span>'),
          grid = TRUE,
          choices = c("N/A",seq(1,10,1),seq(10,100,10),seq(100,1000,100)),
          selected = initial_values$R
        )),
        column(3, numericInput(
          inputId = "R_value",
          label = HTML('<span style="font-size: 14px;">TCR:</span>'),
          value = initial_values$R
        ))
      ),
      div(
        icon("exclamation-triangle", style = "color: orange;"), 
        " For calculations the numeric inputs are used!", 
        style = "margin-top: 5px; margin-bottom: 5px; font-weight: bold; text-align: center; display: block;"),
      fluidRow(
        column(width = 4, class = "col-sm-4 col-md-4 col-lg-4", 
               actionButton("resetsliders", "Reset Sliders", icon = icon("rotate-left"), style = "width: 100%; font-size: 0.7vw;")),
        column(width = 4, class = "col-sm-4 col-md-4 col-lg-4", 
               actionButton("deleteruns", "Delete All Runs", icon = icon("trash-can"), style = "width: 100%; font-size: 0.7vw;")),
        column(width = 4, class = "col-sm-4 col-md-4 col-lg-4", 
               actionButton("recalculate", "Calculate Graphs", icon = icon("gears"), style = "width: 100%; font-size: 0.7vw;")),
      ),
      
    ),
    
    mainPanel(
      div(style = "border: 4px solid #3D3D3D; padding: 20px;",
          fluidRow(
            column(8, tags$div(style = "color: #3D3D3D; font-size: 28px; font-weight: bold;", "T-cell activation model")),
            column(4, tags$img(src = "BayerLogo.png", class = "img-responsive", style = "float: right; max-width: 100px; height: auto;"))),
          plotlyOutput("graph"),
          div(style = "height: 20px;"),  #adding some space between plot and model pic etc.
          fluidRow(
            column(width = 4, class = "col-sm-4 col-md-4 col-lg-4", tags$img(src = "ModelTCA.png", class = "img-responsive", style = "border: 3px white; padding: 0px; max-width: 155%; height: auto;"),
                   div(style = "height: 10px;"),
                   actionButton("Screenshot", "Screenshot",
                                style = "width: 155%; font-size: 0.7vw;")
            ),
            column(width = 3, offset = 2,
                   p("In vitro activation of T cells can be described by the model of kinetic proofreading with limited signaling and an incoherent feed-forward loop (KPL-IFF) [1]."),
                   p(HTML("k<sub>lim</sub>: limited signaling parameter")),
                   p(HTML("k<sub>amp</sub>, k<sub>act</sub>, k<sub>inh</sub>: loop variables")),
                   p(HTML("k<sub>on</sub>: on-rate of ligand")),
                   p(HTML("k<sub>off</sub>: off-rate of ligand")),
                   p(HTML("k<sub>p</sub>: proofreading variable")),
                   p(HTML("E<sub>max</sub>: maximal effect")),
                   p(HTML("M<sub>max</sub>: maximal mediator")),
                   p(HTML("TCR: receptor amount")),
                   p(HTML("pMHC: ligand amount")),
                   div(style = "height: 20px;"),  #adding some space 
                   p("[1] M. Lever et al., Architecture of a minimal signaling pathway explains the T-cell response to a 1 million-fold variation in antigen affinity and dose, 2016, DOI: ",
                     tags$a(href = "https://doi.org/10.1073/pnas.1608820113", "https://doi.org/10.1073/pnas.1608820113"))),
            column(width = 3, class = "col-sm-3 col-md-3 col-lg-3",actionButton("saveruns", "Save Selected Runs", style = "width: 100%; font-size: 0.7vw; margin-bottom: 10px;"), uiOutput("file_list"))
          )
      )
    )
    
  ))} 

#2 Server Function ----

server <- function(input,output,session) {
  #Initializing
  click_count <- reactiveVal(0)
  selected_runs <- reactiveVal(NULL) 
  run_colors <- reactiveVal(list())
  
  # Update selected_runs based on user interaction or recalculation
  observeEvent(input$selected_files, {
    selected_runs(input$selected_files)
  }, ignoreNULL = FALSE)
  
  
  observeEvent(input$recalculate, {      # this observeEvent part is used to connect the observation (clicking the actionbutton) to the calculating of the plot. 
    
    # Increment the counter each time the button is clicked
    click_count(click_count() + 1)
    
    #Disable the Calculate button 
    disable(id = "recalculate")
    
    
    # Extract the input values and store them in variables
    R_value <- input$R_value
    phi_value <- input$phi_value
    lambda_value <- input$lambda_value
    delta_value <- input$delta_value
    my_value <- input$my_value
    kon_value <- input$kon_value
    koff_value <- input$koff_value
    kp_value <- input$kp_value
    Pmax_value <- input$Pmax_value
    Ymax_value <- input$Ymax_value
    L_values <- 10^seq(-4, 4, 0.1)
    
    
    # Create parameter sets 
    parameters <- list(
      yYgen = 100,
      yYdegr = 50000,
      yPgen = 100,
      yPdegr = 50000,
      kp = kp_value,
      phi = phi_value,
      lambda = lambda_value,
      delta = delta_value,
      my = my_value,
      kon = kon_value,
      koff = koff_value,
      Ymax = Ymax_value,
      Pmax = Pmax_value,
      R = R_value
    )
    
    # Define the ODE system
    ode_system <- function(t, state, parameters) {
      with(as.list(c(state, parameters)), {
        dL <- - L * R * kon + koff * (C0 + C1 + C2)
        dR <- - L * R * kon + koff * (C0 + C1 + C2)
        dC0 <- kon * L * R - (koff + kp) * C0
        dC1 <- kp * C0 - (koff + phi * kp) * C1
        dC2 <- phi * kp * C1 - koff * C2
        dY <- yYgen * (Ymax - Y) - yYdegr * Y + lambda * C1 * (Ymax - Y)
        dP <- yPgen * (Pmax - P) - yPdegr * P + delta * Y * (Pmax - P) - my * C1 * P
        return(list(c(dL, dR, dC0, dC1, dC2, dY, dP)))
      })
    }
    
    
    # Function to solve ODEs 
    solve_ode <- function(L_value, params) {
      
      initial_state <- c(L = L_value, R = params$R, C0 = 0, C1 = 0, C2 = 0, Y = 0, P = 0)
      
      ode_output <- ode(y = initial_state, times = 0:30, func = ode_system, parms = params)
      last_row <- ode_output[nrow(ode_output), ]
      return(data.frame(L = L_value, P = last_row["P"]))
    }
    
    # Solve the ODE system for each L value
    results <- lapply(L_values, solve_ode, params = parameters)
    
    # Combine all results into a single data frame
    results <- do.call(rbind, results)
    
    # Add columns for each numeric input, repeating the values for the length of the results dataframe
    results$phi <- rep(input$phi_value, nrow(results))
    results$lambda <- rep(input$lambda_value, nrow(results))
    results$delta <- rep(input$delta_value, nrow(results))
    results$my <- rep(input$my_value, nrow(results))
    results$kon <- rep(input$kon_value, nrow(results))  
    results$koff <- rep(input$koff_value, nrow(results))
    results$kp <- rep(input$kp_value, nrow(results))
    results$Pmax <- rep(input$Pmax_value, nrow(results))
    results$Ymax <- rep(input$Ymax_value, nrow(results))
    results$R <- rep(input$R_value, nrow(results))
    
    # Save files as rds
    file_name <- paste0("results_run_", click_count(),".rds")
    saveRDS(results, file = file_name)
    
    
    
    # This block auto-selects a newly generated run 
    current_selected_runs <- selected_runs()
    if (is.null(current_selected_runs)) {
      current_selected_runs <- character()
    }
    new_selected_runs <- c(current_selected_runs, file_name)
    selected_runs(new_selected_runs)  # Update the reactive value
    
    
    # Assign a unique run_id
    results$run_id <- paste0("run_", click_count()) 
    
    
    # Reactive expression to list the.rds files
    rds_files <- reactive({
      files <- list.files(pattern = "results_run_.*\\.rds$")
      run_numbers <- as.numeric(gsub("results_run_|\\.rds$", "", files))
      files_ordered <- files[order(run_numbers)]
      files_ordered
    })
    
    
    # Render the checkbox group input for the list of .rds files (this is what you see in the shinyApp)
    output$file_list <- renderUI({
      checkboxGroupInput("selected_files", "Select results to display:",
                         choices = rds_files(), selected = selected_runs())
    })
    
    # Debounce (delay) plotting because the checkboxlist bugs out if you click too fast
    debounced_selected_files <- reactive({
      input$selected_files
    }) %>% throttle(1000) #1000 ms debounce 
    
    # PLOTTING         
    #Generate colors for plotting 
    generate_colors <- function(num_runs) {
      # Load turbo-palette
      all_colors <- viridis::viridis_pal(option = "turbo")(256)
      selected_colors <- character(num_runs)
      
      # First color always red
      selected_colors[1] <- "#ff0000" 
      
      if (num_runs > 1) {
        for (i in 2:num_runs) {
          remaining_colors <- all_colors[!all_colors %in% selected_colors]
          max_dist <- 0
          for (color in remaining_colors) {
            # Calculate disctance between current color and all previously selected colors
            dists <- sapply(selected_colors[1:(i-1)], function(selected_color) {
              sum((col2rgb(color) - col2rgb(selected_color))^2)
            })
            min_dist <- min(dists)
            if (min_dist > max_dist) {
              max_dist <- min_dist
              selected_colors[i] <- color
            }
          }
        }
      }
      
      return(selected_colors)
    }
    
    # run_colors saves the colors for each run
    current_colors <- run_colors() 
    
    # Generating colors based on run amount +1 for a new run
    num_runs <- length(current_colors) + 1 
    new_colors <- generate_colors(num_runs)
    
    # Take the last color for the new run
    new_color <- tail(new_colors, 1)
    
    # Adding the new color to the new run
    current_colors[[paste0("run_", click_count())]] <- new_color 
    
    # Saving the color choice
    run_colors(current_colors)
    
    # Observe changes in the checkbox input and render the selected graphs
    output$graph <- renderPlotly({
      # Disable the selection inputs
      shinyjs::disable("selected_files")
      
      # Ensure that the inputs are re-enabled when the function exits
      on.exit(shinyjs::enable("selected_files"), add = TRUE)
      
      # Use debounced input
      req(debounced_selected_files())  # Require that at least one file is selected
      selected_files <- debounced_selected_files()
      
      
      
      # Combine the selected.rds files into one data frame
      combined_results <- do.call(rbind, lapply(selected_files, function(file) {
        result <- readRDS(file)
        # Extract just the numeric part of the filename (assuming it ends with a number)
        run_number <- gsub("results_run_|\\.rds$", "", file)
        # Construct the run_id using the naming convention
        run_id <- paste0("run_", run_number)
        result$run_id <- run_id
        return(result)
      }))
      
      # Sort combined_results by run_id
      combined_results$run_id <- factor(combined_results$run_id, levels = unique(combined_results$run_id))
      
      # Labels for x axis formatting
      custom_labels <- function(breaks) {
        sapply(breaks, function(x) {
          if (!is.na(x)) {
            if(x == 0) {
              "10^0"
            } else {
              paste0("10<sup>", log10(x), "</sup>")
            }
          } else {
            NA 
          }
        })
      }
      
      # Create list for colors based on the reactiveValue
      color_mapping <- run_colors()
      
      # Create the ggplot object
      graph <- ggplot(combined_results, aes(x = L, y = P, color = as.factor(run_id))) +
        geom_line(aes(group = run_id), linewidth = 1.3) +
        labs(title = "Effect versus ligand amount",
             x = "L (ligand amount)",
             y = "E (effect)",
             color = "run ID") +
        scale_color_manual(values = color_mapping) +
        scale_x_log10(breaks = log_breaks(base = 10, n = 9), labels = custom_labels) +
        theme_classic() +
        guides(color = guide_legend(title = "Run ID", order = 1)) + # Ensure legend is ordered
        theme(
          axis.title.x = element_text(size = 14),  
          axis.title.y = element_text(size = 14),
          axis.text.x = element_text(size = 12), 
          axis.text.y = element_text(size = 12),
          panel.background = element_rect(fill = NA), # Make background transparent
          plot.background = element_rect(fill = NA, color = NA),# Make plot background transparent
          legend.background = element_rect(fill = NA)) 
      
      # Convert the ggplot object to a Plotly object with custom tooltip
      ggplotly(graph, tooltip = c("x", "y", "color"))
      
      graph
      
      
    }) 
    
    # Enable the Calculate button again
    enable(id = "recalculate")
  })
  
  
  # Reset slider inputs to initial value when button is clicked
  observeEvent(input$resetsliders, {
    updateSliderTextInput(session, "phi", selected = initial_values$phi)
    updateSliderTextInput(session, "lambda", selected = initial_values$lambda)
    updateSliderTextInput(session, "delta", selected = initial_values$delta)
    updateSliderTextInput(session, "my", selected = initial_values$my)
    updateSliderTextInput(session, "kon", selected = initial_values$kon)
    updateSliderTextInput(session, "koff", selected = initial_values$koff)
    updateSliderTextInput(session, "kp", selected = initial_values$kp)
    updateSliderTextInput(session, "Pmax", selected = initial_values$Pmax)
    updateSliderTextInput(session, "Ymax", selected = initial_values$Ymax)
    updateSliderTextInput(session, "R", selected = initial_values$R)
  })
  
  # This block assures that the numeric input is updated when the slider is moved (excluding NA values)
  observeEvent(input$phi, {
    if (!is.na(as.numeric(input$phi))) {
      updateNumericInput(session, "phi_value", value = as.numeric(input$phi))}})
  observeEvent(input$lambda, {
    if (!is.na(as.numeric(input$lambda))) {
      updateNumericInput(session, "lambda_value", value = as.numeric(input$lambda))}})
  observeEvent(input$delta, {
    if (!is.na(as.numeric(input$delta))) {
      updateNumericInput(session, "delta_value", value = as.numeric(input$delta))}})
  observeEvent(input$my, {
    if (!is.na(as.numeric(input$my))) {
      updateNumericInput(session, "my_value", value = as.numeric(input$my))}})
  observeEvent(input$kon, {
    if (!is.na(as.numeric(input$kon))) {
      updateNumericInput(session, "kon_value", value = as.numeric(input$kon))}})
  observeEvent(input$koff, {
    if (!is.na(as.numeric(input$koff))) {
      updateNumericInput(session, "koff_value", value = as.numeric(input$koff))}})
  observeEvent(input$kp, {
    if (!is.na(as.numeric(input$kp))) {
      updateNumericInput(session, "kp_value", value = as.numeric(input$kp))}})
  observeEvent(input$Pmax, {
    if (!is.na(as.numeric(input$Pmax))) {
      updateNumericInput(session, "Pmax_value", value = as.numeric(input$Pmax))}})
  observeEvent(input$Ymax, {
    if (!is.na(as.numeric(input$Ymax))) {
      updateNumericInput(session, "Ymax_value", value = as.numeric(input$Ymax))}})
  observeEvent(input$R, {
    if (!is.na(as.numeric(input$R))) {
      updateNumericInput(session, "R_value", value = as.numeric(input$R))}})
  
  
  # Observe clicks on the deleteruns button
  observeEvent(input$deleteruns, {
    # Delete all.rds files related to the runs
    results_files <- list.files(pattern = "results_run_.*\\.rds$")
    if (length(results_files) > 0) {
      file.remove(results_files)
    }
    
    # Reset click count to 0, reset tickbox, reset colors
    click_count(0)
    selected_runs(NULL)
    run_colors(list()) 
    
    
    # Reactive expression to list the.rds files
    rds_files <- reactive({
      files <- list.files(pattern = "results_run_.*\\.rds$")
      run_numbers <- as.numeric(gsub("results_run_|\\.rds$", "", files))
      files_ordered <- files[order(run_numbers)]
      files_ordered
    })
    
    
    output$file_list <- renderUI({
      checkboxGroupInput("selected_files", "Select results to display:",
                         choices = rds_files(),
                         selected = selected_runs())
    })
    
    # Optionally, clear the plot by rendering an empty plot or hiding it
    output$graph <- renderPlotly({})
  })
  
  # Checking if numeric and slider are the same value (with delay), if not -> set slider to NA
  observe({
    # Use debounce to introduce a delay, here it's set to 500 milliseconds
    # Adjust the delay time as needed
    debounced_phi_value <- debounce(reactive(input$phi_value), 500)
    
    observeEvent(debounced_phi_value(), {
      # Ensure that input$phi is not NULL or NA before proceeding
      if (!is.null(input$phi) && !is.na(input$phi) && input$phi != "N/A") {
        # Convert input$phi to a numeric value for comparison, handling potential errors
        phi_numeric <- tryCatch(as.numeric(input$phi), error = function(e) NA)
        
        # Now check if the 'phi' slider text and 'phi_value' numeric input are the same
        # Ensure that phi_numeric is not NA before comparing
        if (!is.na(phi_numeric) && phi_numeric != input$phi_value) {
          # If they are not the same, set the slider input to "N/A"
          updateSliderTextInput(session, "phi", selected = "N/A")
        }
      } else {
        
      }
    }, ignoreInit = TRUE)
  })
  
  # For lambda
  observe({
    debounced_lambda_value <- debounce(reactive(input$lambda_value), 500)
    observeEvent(debounced_lambda_value(), {
      if (!is.null(input$lambda) && !is.na(input$lambda) && input$lambda != "N/A") {
        lambda_numeric <- tryCatch(as.numeric(input$lambda), error = function(e) NA)
        if (!is.na(lambda_numeric) && lambda_numeric != input$lambda_value) {
          updateSliderTextInput(session, "lambda", selected = "N/A")
        }
      }
    }, ignoreInit = TRUE)
  })
  
  # For delta
  observe({
    debounced_delta_value <- debounce(reactive(input$delta_value), 500)
    observeEvent(debounced_delta_value(), {
      if (!is.null(input$delta) && !is.na(input$delta) && input$delta != "N/A") {
        delta_numeric <- tryCatch(as.numeric(input$delta), error = function(e) NA)
        if (!is.na(delta_numeric) && delta_numeric != input$delta_value) {
          updateSliderTextInput(session, "delta", selected = "N/A")
        }
      }
    }, ignoreInit = TRUE)
  })
  
  # For my 
  observe({
    debounced_my_value <- debounce(reactive(input$my_value), 500)
    observeEvent(debounced_my_value(), {
      if (!is.null(input$my) && !is.na(input$my) && input$my != "N/A") {
        my_numeric <- tryCatch(as.numeric(input$my), error = function(e) NA)
        if (!is.na(my_numeric) && my_numeric != input$my_value) {
          updateSliderTextInput(session, "my", selected = "N/A")
        }
      }
    }, ignoreInit = TRUE)
  })
  
  # For kon
  observe({
    debounced_kon_value <- debounce(reactive(input$kon_value), 500)
    observeEvent(debounced_kon_value(), {
      if (!is.null(input$kon) && !is.na(input$kon) && input$kon != "N/A") {
        kon_numeric <- tryCatch(as.numeric(input$kon), error = function(e) NA)
        if (!is.na(kon_numeric) && kon_numeric != input$kon_value) {
          updateSliderTextInput(session, "kon", selected = "N/A")
        }
      }
    }, ignoreInit = TRUE)
  })
  
  # For koff
  observe({
    debounced_koff_value <- debounce(reactive(input$koff_value), 500)
    observeEvent(debounced_koff_value(), {
      if (!is.null(input$koff) && !is.na(input$koff) && input$koff != "N/A") {
        koff_numeric <- tryCatch(as.numeric(input$koff), error = function(e) NA)
        if (!is.na(koff_numeric) && koff_numeric != input$koff_value) {
          updateSliderTextInput(session, "koff", selected = "N/A")
        }
      }
    }, ignoreInit = TRUE)
  })
  
  # For kp
  observe({
    debounced_kp_value <- debounce(reactive(input$kp_value), 500)
    observeEvent(debounced_kp_value(), {
      if (!is.null(input$kp) && !is.na(input$kp) && input$kp != "N/A") {
        kp_numeric <- tryCatch(as.numeric(input$kp), error = function(e) NA)
        if (!is.na(kp_numeric) && kp_numeric != input$kp_value) {
          updateSliderTextInput(session, "kp", selected = "N/A")
        }
      }
    }, ignoreInit = TRUE)
  })
  
  # For Pmax
  observe({
    debounced_Pmax_value <- debounce(reactive(input$Pmax_value), 500)
    observeEvent(debounced_Pmax_value(), {
      if (!is.null(input$Pmax) && !is.na(input$Pmax) && input$Pmax != "N/A") {
        Pmax_numeric <- tryCatch(as.numeric(input$Pmax), error = function(e) NA)
        if (!is.na(Pmax_numeric) && Pmax_numeric != input$Pmax_value) {
          updateSliderTextInput(session, "Pmax", selected = "N/A")
        }
      }
    }, ignoreInit = TRUE)
  })
  
  # For Ymax
  observe({
    debounced_Ymax_value <- debounce(reactive(input$Ymax_value), 500)
    observeEvent(debounced_Ymax_value(), {
      if (!is.null(input$Ymax) && !is.na(input$Ymax) && input$Ymax != "N/A") {
        Ymax_numeric <- tryCatch(as.numeric(input$Ymax), error = function(e) NA)
        if (!is.na(Ymax_numeric) && Ymax_numeric != input$Ymax_value) {
          updateSliderTextInput(session, "Ymax", selected = "N/A")
        }
      }
    }, ignoreInit = TRUE)
  })
  
  # For R
  observe({
    debounced_R_value <- debounce(reactive(input$R_value), 500)
    observeEvent(debounced_R_value(), {
      if (!is.null(input$R) && !is.na(input$R) && input$R != "N/A") {
        R_numeric <- tryCatch(as.numeric(input$R), error = function(e) NA)
        if (!is.na(R_numeric) && R_numeric != input$R_value) {
          updateSliderTextInput(session, "R", selected = "N/A")
        }
      }
    }, ignoreInit = TRUE)
  })
  
  generate_and_save_plots <- function(selected_files) {
    # Create a new workbook
    wb <- createWorkbook()
    
    # Create a new worksheet for plots
    addWorksheet(wb, "Plots")
    # Also create a new worksheet for raw data
    addWorksheet(wb, "Raw Data")
    raw_data_start_row <- 1  
    raw_data_start_col <- 1
    

    # Loop through selected rds files
    for (i in seq_along(selected_files)) {
      file <- selected_files[i]
      # Read the rds file
      results <- readRDS(file)
      
      # Generate a title for the plot and data, incorporating the Run ID
      # Assuming the Run ID can be derived from the file name or its index
      run_id <- gsub("results_run_|\\.rds$", "", basename(file))
      
      # Extract L and P columns for the raw data sheet
      lp_data <- results[, c("L", "P")]
      
      # Exclude columns 'L' and 'P' for the data frame to be added to Excel
      results_excluded <- results[1, !(names(results) %in% c("L", "P"))]
      
      # Check if there are columns to exclude and process
      if (ncol(results_excluded) > 0) {
        # Transform to a two-column format: Parameter and Value
        params <- names(results_excluded)
        values <- as.character(unlist(results_excluded))
        data_for_excel <- data.frame(Parameter = params, Value = values, stringsAsFactors = FALSE)
      } else {
        # Create a placeholder to avoid errors
        data_for_excel <- data.frame(Parameter = NA, Value = NA)
      }
      
      
      
      # Create a plot for the current file
      plot <- ggplot(results, aes(x = L, y = P)) +
        geom_line(linewidth = 1.1) +
        labs(title =  paste("Run", run_id),
             x = "L (ligand amount)",
             y = "P (effect)") +
        scale_x_log10(breaks = log_breaks(base = 10, n = 9)) +
        theme_classic() +
        theme(
          axis.title.x = element_text(size = 14),  
          axis.title.y = element_text(size = 14),
          axis.text.x = element_text(size = 12), 
          axis.text.y = element_text(size = 12))
      
      # Save the plot as an image
      plot_filename <- paste0(tempdir(), "/plot_", i, ".png")
      ggsave(plot_filename, plot, width = 10, height = 5, dpi = 300)
      
      # Add the plot image to the worksheet
      insertImage(wb, sheet = "Plots", file = plot_filename, startRow = (i-1)*22 + 1, startCol = 1, width = 8, height = 4)
      
      # Add the parameter names and their values next to the plot image
      writeData(wb, sheet = "Plots", x = data_for_excel, startRow = (i-1)*22 + 1, startCol = 12, colNames = TRUE)
      
      # Write a header for the raw data indicating the run number
      run_header <- data.frame(Run_Header = paste("Run", run_id))
      writeData(wb, sheet = "Raw Data", x = run_header, startRow = raw_data_start_row, startCol = raw_data_start_col, colNames = FALSE)
      raw_data_start_row <- raw_data_start_row + 1  # Increment the row for the actual data
      
      # Add the raw data (L and P) to the "Raw Data" sheet, now with a header for the run
      writeData(wb, sheet = "Raw Data", x = lp_data, startRow = raw_data_start_row, startCol = raw_data_start_col, colNames = TRUE)
      
      
      
      # Reset the starting row for the next run's header and data
      raw_data_start_row <- 1  # Reset row to 1 for the next run's data
      raw_data_start_col <- raw_data_start_col + 3  # Move to the next set of columns, leaving a space for separation
    }
    
    # Save the workbook
    saveWorkbook(wb, "selected_runs_plots_and_data.xlsx", overwrite = TRUE)
  }
  
  
  
  # Observe clicks on the 'saveruns' action button
  observeEvent(input$saveruns, {
    # Ensure that there are selected runs to save
    if (!is.null(input$selected_files) && length(input$selected_files) > 0) {
      # Retrieve the full paths of the selected RDS files
      selected_files_paths <- sapply(input$selected_files, function(file) {
        # Construct the full path to the file
        # Assuming the files are stored in a specific directory, adjust as necessary
        file_path <- file.path(getwd(), file)
        return(file_path)
      })
      
      # Call the function to generate and save plots
      # Pass the full paths of the selected files
      generate_and_save_plots(selected_files_paths)
      
      # Optional: Display a message to indicate completion
      showNotification("Selected runs have been saved with plots.", type = "message")
    } else {
      # No files selected, show a warning
      showNotification("No runs selected to save.", type = "warning")
    }
  })  
  
  
  
  # Stop the parallel backend when the app is closed
  onStop(function() {
    # Delete all results files when the app is closed
    results_files <- list.files(pattern = "results_run_.*\\.rds$")
    file.remove(results_files)
  })
  
  
  
  
  # slider color
  css_color <- "#20B2AA"
  
  session$onFlushed(function() {
    insertUI(
      selector = "head", where = "beforeEnd",
      ui = tags$style(HTML(sprintf("
      /* Make handle fully invisible, but still clickable/draggable */
      .irs--flat .irs-handle {
        background-color: transparent !important;
        border: none !important;
        box-shadow: none !important;
        opacity: 0 !important;           /* invisible */
      }
      .irs--flat .irs-handle.state_hover,
      .irs--flat .irs-handle:hover {
        opacity: 0 !important;           /* remain invisible on hover */
      }
      .irs--flat .irs-handle > i {
        display: none !important;        /* hide stripes */
      }

      /* Bar & edge in turquoise */
      .irs--flat .irs-bar,
      .irs--flat .irs-bar-edge {
        background-color: %s !important;
        border-color: %s !important;
      }

      /* Tooltip bubbles in turquoise with white text */
      .irs--flat .irs-single,
      .irs--flat .irs-from,
      .irs--flat .irs-to {
        background-color: %s !important;
        border-color: %s !important;
        color: #fff !important;
      }

      /* Tooltip pointer triangle (make it turquoise) */
      .irs--flat .irs-single:before,
      .irs--flat .irs-from:before,
      .irs--flat .irs-to:before {
        border-top-color: %s !important;
      }

      /* Base line neutral */
      .irs--flat .irs-line { background-color: #eaeaea !important; }
    ", css_color, css_color, css_color, css_color, css_color)))
    )
  }, once = TRUE)
  
  
  observeEvent(input$Screenshot, {
    screenshot(scale = 10)
  })
}



#3 server call to action ----
shinyApp(ui = ui, server = server)

