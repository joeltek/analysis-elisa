# analysis-elisa
ELISA data analysis pipeline based on R

## Plate View Shiny App

An interactive Shiny app (`shiny_plate_view/app.R`) displays 96-well plate
heatmaps for any experiment that has been processed by the analysis pipeline.

### What it shows

| Panel | Source | Description |
|---|---|---|
| **Dilution Factor** | `mapping files/<ID>_MappingFile.csv` | Colour of each well shows the sample dilution factor used. Standards, blanks, and controls are shown in grey (no dilution factor). |
| **Assay Performance** | `Analysis_Pipeline/<ID>/<ID>Output.csv` | Colour of each well shows whether the measurement fell within the assay dynamic range (`Within Range`, `Below LLOD`, `Between LLOD and LLOQ`, `Above ULOQ`). |

### How to run

```r
# From the project root in R:
shiny::runApp("shiny_plate_view")
```

Use the **Experiment ID** dropdown to switch between experiments.  
Use the **Plate** dropdown (P1/P2) to view individual plates.  
Use the **Cytokine** dropdown to select which ELISA to display in the
Assay Performance panel.

### Adding a new experiment

1. Place the mapping file at `mapping files/<EXP_ID>_MappingFile.csv`
   with columns: `PlateID, Well, Content, Sample_ID, Dilution_Factor`
2. Place the output CSV at `Analysis_Pipeline/<EXP_ID>/<EXP_ID>Output.csv`
   with columns: `Well, PlateID, Cytokine, Donor_ID, Sample, PreTreatment,
   Stim, Stim_Conc, Timepoint, Final_Concentration_pg_mL, Within_Range`

The new experiment ID will appear automatically in the dropdown.

## Project structure

```
├── Analysis_Pipeline/          # Per-experiment analysis outputs
│   └── EXP-26-EI5039/
│       └── EXP-26-EI5039Output.csv
├── mapping files/              # Per-experiment plate mapping files
│   └── EXP-26-EI5039_MappingFile.csv
├── scripts/                    # Analysis scripts (Quarto / R Markdown)
│   └── EXP-26-EI5039-40_Comparison_1hvs20hStim.qmd
└── shiny_plate_view/           # Interactive plate view Shiny app
    └── app.R
```
