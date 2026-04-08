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
  data_dir <- file.path(root_dir, "data_files")
  if (!dir.exists(data_dir)) return(character(0))
  files <- list.files(data_dir, pattern = "Output\\.csv$", full.names = FALSE)
  sub("Output\\.csv$", "", files)
}

load_mapping <- function(root_dir, exp_id) {
  path <- file.path(root_dir, "map_files", paste0(exp_id, "_MappingFile.csv"))
  if (!file.exists(path)) return(NULL)

  # Detect delimiter from the first line
  first_line <- readLines(path, n = 1)
  sep <- if (grepl(";", first_line)) ";" else ","

  df <- read.csv(path, sep = sep, stringsAsFactors = FALSE,
                 na.strings = c("", "NA"))

  # Flexible column name: rename Sample -> Sample_ID if needed
  if ("Sample" %in% names(df) && !"Sample_ID" %in% names(df)) {
    names(df)[names(df) == "Sample"] <- "Sample_ID"
  }

  # Add PlateID if absent (single-plate file: assume P1)
  if (!"PlateID" %in% names(df)) {
    df$PlateID <- "P1"
  }

  # Normalise well IDs: uppercase row letter, strip leading zeros (A01 -> A1)
  df$Well <- toupper(gsub("^([A-Ha-h])0*(\\d+)$", "\\1\\2", df$Well))

  # Fix decimal-comma values in DF columns (e.g. "1,5" -> 1.5)
  df_cols <- grep("^DF_|^Dilution_Factor$", names(df), value = TRUE)
  for (col in df_cols) {
    if (is.character(df[[col]])) {
      df[[col]] <- as.numeric(gsub(",", ".", df[[col]]))
    }
  }

  df
}

load_output <- function(root_dir, exp_id) {
  path <- file.path(root_dir, "data_files", paste0(exp_id, "Output.csv"))
  if (!file.exists(path)) return(NULL)

  # Detect delimiter from the first line
  first_line <- readLines(path, n = 1)
  sep <- if (grepl(";", first_line)) ";" else ","

  df <- read.csv(path, sep = sep, stringsAsFactors = FALSE,
                 na.strings = c("", "NA"))

  # Normalise well IDs: uppercase row letter, strip leading zeros (B01 -> B1)
  df$Well <- toupper(gsub("^([A-Ha-h])0*(\\d+)$", "\\1\\2", df$Well))

  df
}

# Find the dilution-factor column that matches a given cytokine.
# Tries DF_{cytokine} (case-insensitive), then "Dilution_Factor",
# then the first DF_* column found.
get_df_column <- function(mapping_df, cytokine) {
  cols <- names(mapping_df)
  candidate <- paste0("DF_", cytokine)
  idx <- which(tolower(cols) == tolower(candidate))
  if (length(idx) > 0) return(cols[idx[1]])
  if ("Dilution_Factor" %in% cols) return("Dilution_Factor")
  df_cols <- grep("^DF_", cols, value = TRUE)
  if (length(df_cols) > 0) return(df_cols[1])
  NULL
}

# Return "B" for blank wells, "S" for standard wells, NA otherwise.
well_label <- function(content) {
  dplyr::case_when(
    grepl("^blank$",   content, ignore.case = TRUE) ~ "B",
    grepl("standard",  content, ignore.case = TRUE) ~ "S",
    TRUE ~ NA_character_
  )
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
# df_col: name of the column to use as dilution factor (e.g. "DF_IFNa").
make_dilution_heatmap <- function(data, plate_id, df_col) {
  if (is.null(df_col) || !df_col %in% names(data)) {
    return(ggplot() +
             annotate("text", x = 0.5, y = 0.5,
                      label = "Dilution factor column not found", size = 5) +
             theme_void())
  }

  plate_data <- data %>%
    filter(PlateID == plate_id) %>%
    add_well_coords()

  # Mapping file may only define one layout (all plates share the same wells).
  # If no rows found for this plate, fall back to the full mapping.
  if (nrow(plate_data) == 0) {
    plate_data <- data %>% add_well_coords()
  }

  if ("Content" %in% names(plate_data)) {
    plate_data <- plate_data %>% mutate(Well_Label = well_label(Content))
  } else {
    plate_data$Well_Label <- NA_character_
  }

  ggplot(plate_data,
         aes(x = Col, y = Row, fill = as.numeric(.data[[df_col]]))) +
    geom_tile(colour = "white", linewidth = 0.6) +
    geom_text(aes(label = Well_Label), size = 3, fontface = "bold",
              colour = "gray30", na.rm = TRUE, show.legend = FALSE) +
    scale_fill_viridis_c(
      name = "Dilution factor",
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
      panel.grid      = element_blank(),
      plot.title      = element_text(face = "bold", hjust = 0.5),
      axis.ticks      = element_blank(),
      legend.position = "none"
    )
}

# Build a 96-well plate heatmap for Within_Range (categorical fill).
# Joins the full well list from the mapping file so standards/blanks are shown.
make_range_heatmap <- function(mapping_data, output_data, plate_id, cytokine) {
  # All wells for this plate (from mapping).
  # If the mapping has no row for plate_id (single-layout file), use any plate's
  # rows and re-label them so the join with output data works correctly.
  mapping_plates <- unique(mapping_data$PlateID)
  map_plate <- if (plate_id %in% mapping_plates) plate_id else mapping_plates[1]

  all_wells <- mapping_data %>%
    filter(PlateID == map_plate) %>%
    select(any_of(c("Well", "Content", "Sample_ID"))) %>%
    mutate(PlateID = plate_id)

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
                            levels = names(WITHIN_RANGE_COLOURS)),
      Well_Label = if ("Content" %in% names(.)) well_label(Content)
                  else NA_character_
    )

  ggplot(plot_data,
         aes(x = Col, y = Row, fill = Within_Range)) +
    geom_tile(colour = "white", linewidth = 0.6) +
    geom_text(aes(label = Well_Label), size = 3, fontface = "bold",
              colour = "gray30", na.rm = TRUE, show.legend = FALSE) +
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
      panel.grid      = element_blank(),
      plot.title      = element_text(face = "bold", hjust = 0.5),
      axis.ticks      = element_blank(),
      legend.position = "none"
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
        choices  = NULL,
        selected = NULL
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

      uiOutput("sidebar_legends")
    ),

    mainPanel(
      width = 9,
      uiOutput("plates_ui")
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

  # Populate plate dropdown from the output CSV (plus "All" option)
  observe({
    df <- output_data()
    if (is.null(df)) return()
    plates <- sort(unique(df$PlateID))
    updateSelectInput(session, "plate_id",
                      choices  = c("All", plates),
                      selected = if (length(plates) > 0) plates[1] else NULL)
  })

  # Dilution factor heatmap (single-plate view)
  output$dilution_heatmap <- renderPlot({
    req(input$plate_id != "All")
    mapping <- mapping_data()
    if (is.null(mapping)) {
      return(ggplot() +
               annotate("text", x = 0.5, y = 0.5,
                        label = "Mapping file not found", size = 5) +
               theme_void())
    }
    req(input$cytokine)
    df_col <- get_df_column(mapping, input$cytokine)
    make_dilution_heatmap(mapping, input$plate_id, df_col)
  })

  # Assay performance (Within_Range) heatmap (single-plate view)
  output$range_heatmap <- renderPlot({
    req(input$plate_id != "All")
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

  # Dynamic per-plate plots used in the "All" view
  observe({
    mapping <- mapping_data()
    output_df <- output_data()
    req(mapping, output_df, input$cytokine)

    plates <- sort(unique(output_df$PlateID))
    lapply(plates, function(pid) {
      local({
        p_id  <- pid
        cyt   <- input$cytokine
        m_df  <- mapping
        o_df  <- output_df

        output[[paste0("plot_df_",    p_id)]] <- renderPlot({
          df_col <- get_df_column(m_df, cyt)
          make_dilution_heatmap(m_df, p_id, df_col)
        })

        output[[paste0("plot_range_", p_id)]] <- renderPlot({
          make_range_heatmap(m_df, o_df, p_id, cyt)
        })
      })
    })
  })

  # Sidebar legends: dilution factor (discrete swatches) + Within_Range
  output$sidebar_legends <- renderUI({
    mapping <- mapping_data()

    # ---- Dilution Factor discrete legend ----
    df_legend <- NULL
    if (!is.null(mapping) && !is.null(input$cytokine)) {
      df_col <- get_df_column(mapping, input$cytokine)
      if (!is.null(df_col) && df_col %in% names(mapping)) {
        vals <- sort(unique(as.numeric(mapping[[df_col]])), na.last = NA)
        if (length(vals) > 0) {
          # Map each value to a plasma colour matching the continuous ggplot2 scale
          rng  <- range(vals)
          norm <- if (diff(rng) == 0) rep(0.5, length(vals))
                  else (vals - rng[1]) / diff(rng)
          pal    <- viridisLite::viridis(256, option = "plasma")
          colors <- pal[pmax(1L, round(norm * 255) + 1L)]
          # Use white text on dark swatches, black on bright yellow end
          txt_colors <- ifelse(norm > 0.8, "black", "white")

          swatches <- lapply(seq_along(vals), function(i) {
            tags$li(tags$span(
              style = paste0("background:", colors[i],
                             "; color:", txt_colors[i],
                             "; padding:1px 8px; border-radius:3px;",
                             " font-family:monospace;"),
              as.character(vals[i])
            ))
          })
          df_legend <- tagList(
            tags$p(tags$b("Dilution factor:")),
            tags$ul(style = "padding-left:18px;", swatches)
          )
        }
      }
    }

    # ---- Within_Range categorical legend ----
    range_legend <- tagList(
      tags$p(tags$b("Assay performance:")),
      tags$ul(
        style = "padding-left:18px;",
        tags$li(tags$span(
          style = paste0("background:", WITHIN_RANGE_COLOURS["Within Range"],
                         "; color:white; padding:1px 8px; border-radius:3px;"),
          "Within Range")),
        tags$li(tags$span(
          style = paste0("background:", WITHIN_RANGE_COLOURS["Between LLOD and LLOQ"],
                         "; padding:1px 8px; border-radius:3px;"),
          "Between LLOD and LLOQ")),
        tags$li(tags$span(
          style = paste0("background:", WITHIN_RANGE_COLOURS["Below LLOD"],
                         "; padding:1px 8px; border-radius:3px;"),
          "Below LLOD")),
        tags$li(tags$span(
          style = paste0("background:", WITHIN_RANGE_COLOURS["Above ULOQ"],
                         "; color:white; padding:1px 8px; border-radius:3px;"),
          "Above ULOQ"))
      )
    )

    tagList(hr(), df_legend, range_legend)
  })

  # Main panel UI: single plate or all plates
  output$plates_ui <- renderUI({
    req(input$plate_id)

    if (input$plate_id == "All") {
      df     <- output_data()
      plates <- if (!is.null(df)) sort(unique(df$PlateID)) else character(0)

      tagList(
        fluidRow(
          column(6, h4("Dilution Factor",    style = "text-align:center;")),
          column(6, h4("Assay Performance",  style = "text-align:center;"))
        ),
        lapply(plates, function(pid) {
          fluidRow(
            column(6, plotOutput(paste0("plot_df_",    pid), height = "430px")),
            column(6, plotOutput(paste0("plot_range_", pid), height = "430px"))
          )
        })
      )
    } else {
      fluidRow(
        column(6,
               h4("Dilution Factor",   style = "text-align:center;"),
               plotOutput("dilution_heatmap", height = "430px")),
        column(6,
               h4("Assay Performance", style = "text-align:center;"),
               plotOutput("range_heatmap",    height = "430px"))
      )
    }
  })
}

shinyApp(ui, server)
