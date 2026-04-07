# ELISA Plate View Shiny App
# Displays 96-well plate heatmaps for:
#   1. Sample dilution factors (from mapping file)
#   2. Assay performance / Within_Range (from Output CSV)

library(shiny)
library(ggplot2)
library(dplyr)
library(tidyr)

# ---- Path helpers ----

# Determine the project root directory.
# Supports running from the shiny_plate_view/ subdirectory or from root.
get_root_dir <- function() {
  wd <- normalizePath(getwd())
  if (basename(wd) == "shiny_plate_view") dirname(wd) else wd
}

# ---- Data loading helpers ----

get_experiment_ids <- function(root_dir) {
  pipeline_dir <- file.path(root_dir, "Analysis_Pipeline")
  if (!dir.exists(pipeline_dir)) return(character(0))
  ids <- list.dirs(pipeline_dir, full.names = FALSE, recursive = FALSE)
  ids[nchar(ids) > 0]
}

load_mapping <- function(root_dir, exp_id) {
  path <- file.path(root_dir, "mapping files", paste0(exp_id, "_MappingFile.csv"))
  if (!file.exists(path)) return(NULL)
  read.csv(path, stringsAsFactors = FALSE, na.strings = c("", "NA"))
}

load_output <- function(root_dir, exp_id) {
  path <- file.path(root_dir, "Analysis_Pipeline", exp_id,
                    paste0(exp_id, "Output.csv"))
  if (!file.exists(path)) return(NULL)
  read.csv(path, stringsAsFactors = FALSE, na.strings = c("", "NA"))
}

# ---- Plate layout helpers ----

# Parse well ID (e.g. "A1", "H12") into separate Row / Col columns.
add_well_coords <- function(df) {
  df %>%
    mutate(
      Row = factor(substr(Well, 1, 1),
                   # Reversed: ggplot2 plots lowest level at bottom,
                   # so rev() ensures Row A appears at the top of the plate.
                   levels = rev(LETTERS[1:8])),
      Col = as.integer(substr(Well, 2, nchar(Well)))
    )
}

# ---- Plot helpers ----

WITHIN_RANGE_COLOURS <- c(
  "Below LLOD"           = "#FDE725",   # yellow (viridis end)
  "Between LLOD and LLOQ" = "#B0B0B0",  # grey
  "Within Range"          = "#21908C",  # teal (viridis mid)
  "Above ULOQ"            = "#E31A1C"   # red
)

NA_FILL <- "grey92"

# Build a 96-well plate heatmap for dilution factor (continuous numeric fill).
make_dilution_heatmap <- function(data, plate_id) {
  plate_data <- data %>%
    filter(PlateID == plate_id) %>%
    add_well_coords()

  ggplot(plate_data,
         aes(x = Col, y = Row, fill = Dilution_Factor)) +
    geom_tile(colour = "white", linewidth = 0.6) +
    scale_fill_viridis_c(
      name = "Dilution\nfactor",
      option = "plasma",
      na.value = NA_FILL,
      direction = 1
    ) +
    scale_x_continuous(
      breaks = 1:12, limits = c(0.5, 12.5), expand = c(0, 0)
    ) +
    scale_y_discrete(expand = c(0, 0)) +
    coord_fixed() +
    labs(
      title = paste("Plate", plate_id, "\u2013 Sample Dilution Factor"),
      x = "Column", y = "Row"
    ) +
    theme_minimal(base_size = 13) +
    theme(
      panel.grid  = element_blank(),
      plot.title  = element_text(face = "bold", hjust = 0.5),
      axis.ticks  = element_blank(),
      legend.position = "right"
    )
}

# Build a 96-well plate heatmap for Within_Range (categorical fill).
# Joins the full well list from the mapping file so standards/blanks are shown.
make_range_heatmap <- function(mapping_data, output_data, plate_id, cytokine) {
  # All wells for this plate (from mapping)
  all_wells <- mapping_data %>%
    filter(PlateID == plate_id) %>%
    select(Well, PlateID, Content, Sample_ID)

  # Sample results for this plate + cytokine
  results <- output_data %>%
    filter(PlateID == plate_id, Cytokine == cytokine) %>%
    select(Well, PlateID, Within_Range)

  # Left-join: all wells get results where available, NA elsewhere
  plot_data <- all_wells %>%
    left_join(results, by = c("Well", "PlateID")) %>%
    add_well_coords() %>%
    mutate(
      Within_Range = factor(Within_Range,
                            levels = names(WITHIN_RANGE_COLOURS))
    )

  ggplot(plot_data,
         aes(x = Col, y = Row, fill = Within_Range)) +
    geom_tile(colour = "white", linewidth = 0.6) +
    scale_fill_manual(
      name   = "Within Range",
      values = WITHIN_RANGE_COLOURS,
      na.value = NA_FILL,
      drop = FALSE
    ) +
    scale_x_continuous(
      breaks = 1:12, limits = c(0.5, 12.5), expand = c(0, 0)
    ) +
    scale_y_discrete(expand = c(0, 0)) +
    coord_fixed() +
    labs(
      title = paste("Plate", plate_id, "\u2013", cytokine, "Assay Performance"),
      x = "Column", y = "Row"
    ) +
    theme_minimal(base_size = 13) +
    theme(
      panel.grid  = element_blank(),
      plot.title  = element_text(face = "bold", hjust = 0.5),
      axis.ticks  = element_blank(),
      legend.position = "right"
    )
}

# ---- UI ----

ui <- fluidPage(

  titlePanel("ELISA Plate View"),

  sidebarLayout(
    sidebarPanel(
      width = 3,

      selectInput(
        "exp_id",
        label    = "Experiment ID:",
        choices  = NULL,
        selected = NULL
      ),

      selectInput(
        "plate_id",
        label    = "Plate:",
        choices  = c("P1", "P2"),
        selected = "P1"
      ),

      selectInput(
        "cytokine",
        label    = "Cytokine (assay performance):",
        choices  = NULL,
        selected = NULL
      ),

      hr(),

      tags$p(tags$b("Dilution Factor"),
             "– colour shows the fold dilution of each sample well.",
             "Grey wells are standards, blanks, or controls."),

      tags$p(tags$b("Assay Performance"),
             "– colour shows whether each sample measurement fell",
             "within the assay dynamic range."),

      hr(),

      tags$p(tags$b("Legend:")),
      tags$ul(
        tags$li(tags$span(
          style = paste0("background:", WITHIN_RANGE_COLOURS["Within Range"],
                         "; color:white; padding:1px 6px; border-radius:3px;"),
          "Within Range")),
        tags$li(tags$span(
          style = paste0("background:", WITHIN_RANGE_COLOURS["Between LLOD and LLOQ"],
                         "; padding:1px 6px; border-radius:3px;"),
          "Between LLOD and LLOQ")),
        tags$li(tags$span(
          style = paste0("background:", WITHIN_RANGE_COLOURS["Below LLOD"],
                         "; padding:1px 6px; border-radius:3px;"),
          "Below LLOD")),
        tags$li(tags$span(
          style = paste0("background:", WITHIN_RANGE_COLOURS["Above ULOQ"],
                         "; color:white; padding:1px 6px; border-radius:3px;"),
          "Above ULOQ"))
      )
    ),

    mainPanel(
      width = 9,
      fluidRow(
        column(
          6,
          h4("Dilution Factor", style = "text-align:center;"),
          plotOutput("dilution_heatmap", height = "380px")
        ),
        column(
          6,
          h4("Assay Performance", style = "text-align:center;"),
          plotOutput("range_heatmap", height = "380px")
        )
      )
    )
  )
)

# ---- Server ----

server <- function(input, output, session) {

  root_dir <- get_root_dir()

  # Populate experiment dropdown on startup
  observe({
    ids <- get_experiment_ids(root_dir)
    updateSelectInput(session, "exp_id",
                      choices  = ids,
                      selected = if (length(ids) > 0) ids[1] else NULL)
  })

  # Reactive: load mapping file for the selected experiment
  mapping_data <- reactive({
    req(input$exp_id)
    load_mapping(root_dir, input$exp_id)
  })

  # Reactive: load output CSV for the selected experiment
  output_data <- reactive({
    req(input$exp_id)
    load_output(root_dir, input$exp_id)
  })

  # Populate cytokine dropdown from the output CSV
  observe({
    df <- output_data()
    if (is.null(df)) return()
    cytokines <- sort(unique(df$Cytokine))
    updateSelectInput(session, "cytokine",
                      choices  = cytokines,
                      selected = if (length(cytokines) > 0) cytokines[1] else NULL)
  })

  # Dilution factor heatmap
  output$dilution_heatmap <- renderPlot({
    mapping <- mapping_data()
    if (is.null(mapping)) {
      return(ggplot() +
               annotate("text", x = 0.5, y = 0.5,
                        label = "Mapping file not found", size = 5) +
               theme_void())
    }
    make_dilution_heatmap(mapping, input$plate_id)
  })

  # Assay performance (Within_Range) heatmap
  output$range_heatmap <- renderPlot({
    mapping <- mapping_data()
    output_df <- output_data()
    if (is.null(mapping) || is.null(output_df)) {
      return(ggplot() +
               annotate("text", x = 0.5, y = 0.5,
                        label = "Output file not found", size = 5) +
               theme_void())
    }
    req(input$cytokine)
    make_range_heatmap(mapping, output_df, input$plate_id, input$cytokine)
  })
}

shinyApp(ui, server)
