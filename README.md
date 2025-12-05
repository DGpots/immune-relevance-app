# Immune Relevance App

This is a Shiny app for ranking differentially expressed genes (DEGs) by immune involvement and visualizing them in a volcano plot.

## How to run

1. Clone this repository:
   git clone https://github.com/DGpots/immune-relevance-app.git

2. Open R or RStudio.

3. Run:
   shiny::runApp("immune-relevance-app")

4. Upload any DEG CSV with columns including:
   - gene name
   - log2 fold change
   - padj/FDR/qvalue

5. The app will compute:
   - immune_score (weighted GO:BP involvement)
   - immune_scaled
   - relevance
   - an immune volcano plot (downloadable)
   - a ranked relevance table (downloadable)
