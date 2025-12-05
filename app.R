## ======================================================================
##  Immune Relevance Shiny App
##  Upload any DEG CSV → compute immune scores → rank genes → volcano plot
## ======================================================================

library(shiny)
library(readr)
library(dplyr)
library(janitor)
library(stringr)
library(msigdbr)
library(ggplot2)
library(ggrepel)

## -------------------------------------------------------------
## Helper: choose first matching column by name or regex
## -------------------------------------------------------------
pick <- function(nms, choices, rx = NULL) {
  hit <- intersect(choices, nms)
  if (length(hit)) return(hit[1])
  if (!is.null(rx)) {
    j <- which(stringr::str_detect(nms, rx))
    if (length(j)) return(nms[j[1]])
  }
  return(NA_character_)
}

## -------------------------------------------------------------
## UI
## -------------------------------------------------------------
ui <- fluidPage(
  titlePanel("Immune Relevance Scoring & Volcano Plot"),
  
  sidebarLayout(
    sidebarPanel(
      fileInput("file", "Upload DEG CSV", accept = c(".csv")),
      
      numericInput(
        "topn",
        "Number of up/down genes to label:",
        value = 10, min = 1, max = 50, step = 1
      ),
      
      hr(),
      downloadButton("download_table", "Download Relevance Table"),
      br(), br(),
      downloadButton("download_plot", "Download Volcano Plot")
    ),
    
    mainPanel(
      tabsetPanel(
        tabPanel("Relevance Table", tableOutput("table")),
        tabPanel("Volcano Plot", plotOutput("volcano", height = "600px"))
      )
    )
  )
)

## -------------------------------------------------------------
## SERVER
## -------------------------------------------------------------
server <- function(input, output, session) {
  
  ## --------- 1) Load DEG data from uploaded CSV ----------
  deg_data <- reactive({
    req(input$file)
    
    df <- read_csv(
      input$file$datapath,
      guess_max = 100000,
      show_col_types = FALSE
    ) |>
      clean_names()
    
    nms <- names(df)
    
    gene_col <- pick(
      nms,
      c("gene", "gene_name", "symbol", "hgnc_symbol"),
      "(^|_)gene"
    )
    
    lfc_col <- pick(
      nms,
      c("log2fc", "log2_fold_change", "log2foldchange", "lfc"),
      "log2.*fold"
    )
    
    padj_col <- pick(
      nms,
      c("padj", "p_adj", "qvalue", "fdr"),
      "padj|qvalue|fdr"
    )
    
    if (is.na(gene_col) || is.na(lfc_col)) {
      stop(
        "Could not identify gene or log2FC columns.\n",
        "Column names are: ", paste(nms, collapse = ", ")
      )
    }
    
    df |>
      transmute(
        gene   = .data[[gene_col]],
        log2fc = suppressWarnings(as.numeric(.data[[lfc_col]])),
        padj   = if (!is.na(padj_col))
          suppressWarnings(as.numeric(.data[[padj_col]]))
        else
          NA_real_
      ) |>
      distinct(gene, .keep_all = TRUE) |>
      filter(!is.na(gene), !is.na(log2fc))
  })
  
  ## --------- 2) Compute immune scores + relevance ----------
  scored <- reactive({
    df <- deg_data()
    
    # Immune-related GO:BP terms
    immune_rx <- "(IMMUNE|INFLAMM|LEUKOCYTE|CYTOKINE|INTERFERON|CHEMOKINE|ANTIGEN|T_CELL|B_CELL)"
    
    m_df <- msigdbr(
      species    = "Mus musculus",
      category   = "C5",
      subcategory = "GO:BP"
    ) |>
      select(gs_name, gene_symbol) |>
      filter(str_detect(gs_name, regex(immune_rx, ignore_case = TRUE))) |>
      distinct()
    
    # Weighted immune score
    term_sizes <- m_df |>
      count(gs_name, name = "term_size")
    
    m_w <- m_df |>
      left_join(term_sizes, by = "gs_name") |>
      mutate(w = 1 / log1p(term_size))
    
    immune_weight <- m_w |>
      group_by(gene_symbol) |>
      summarise(immune_score = sum(w), .groups = "drop") |>
      rename(gene = gene_symbol)
    
    # Join with DE table and compute scaled scores + relevance
    df2 <- df |>
      left_join(immune_weight, by = "gene") |>
      mutate(
        immune_score = ifelse(is.na(immune_score), 0, immune_score)
      )
    
    max_imm <- max(df2$immune_score, na.rm = TRUE)
    
    if (is.finite(max_imm) && max_imm > 0) {
      df2 <- df2 |>
        mutate(
          immune_scaled = 100 * immune_score / max_imm
        )
    } else {
      df2 <- df2 |>
        mutate(
          immune_scaled = 0
        )
    }
    
    df2 |>
      mutate(
        relevance = abs(log2fc) * immune_score
      )
  })
  
  ## --------- 3) Relevance table output ----------
  output$table <- renderTable({
    req(scored())
    scored() |>
      arrange(desc(relevance)) |>
      head(150)
  })
  
  ## --------- 4) Volcano plot (as a reactive object) ----------
  volcano_plot <- reactive({
    req(scored())
    df <- scored()
    topn <- input$topn
    
    up_labels <- df |>
      filter(log2fc > 0) |>
      arrange(desc(relevance)) |>
      slice_head(n = topn) |>
      pull(gene)
    
    down_labels <- df |>
      filter(log2fc < 0) |>
      arrange(desc(relevance)) |>
      slice_head(n = topn) |>
      pull(gene)
    
    label_genes <- unique(c(up_labels, down_labels))
    
    ggplot(df, aes(x = log2fc, y = immune_scaled)) +
      geom_hline(yintercept = 0, linewidth = 0.25, color = "grey80") +
      geom_vline(xintercept = 0, linewidth = 0.25, color = "grey80") +
      geom_point(alpha = 0.7, size = 1.8, color = "steelblue3") +
      ggrepel::geom_text_repel(
        data = dplyr::filter(df, gene %in% label_genes),
        aes(label = gene),
        size = 3.1,
        max.overlaps = 100,
        min.segment.length = 0
      ) +
      labs(
        title = "Immune Relevance Volcano Plot",
        x = "Log2 Fold Change",
        y = "Immune Pathway Score (scaled 0–100)"
      ) +
      coord_cartesian(ylim = c(0, 100)) +
      theme_classic(base_size = 12)
  })
  
  output$volcano <- renderPlot({
    volcano_plot()
  })
  
  ## --------- 5) Download handlers ----------
  output$download_table <- downloadHandler(
    filename = function() {
      "relevance_scores.csv"
    },
    content = function(file) {
      write.csv(scored(), file, row.names = FALSE)
    }
  )
  
  output$download_plot <- downloadHandler(
    filename = function() {
      "immune_volcano.png"
    },
    content = function(file) {
      ggsave(
        filename = file,
        plot = volcano_plot(),
        width = 7,
        height = 5,
        dpi = 300
      )
    }
  )
}

shinyApp(ui, server)
