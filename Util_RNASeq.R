

# Deifne metadata
define_metadata <- function(data, selected_columns = "", split = "_", patterns, metadata_labels, 
                            paterns_num = 1){
  # data : expression data
  # selected_columns : The list of character to use to extract the the patterns
  # patterns : The pattern to extract
  # metadata_labels : The metadata labels to give for each pattern
  # paterns_num : The number of character to extract 
  # split : The character after which to perform the splitting. 
  
  if(length(patterns) != length(metadata_labels)) {
    stop("The number of patterns must match the number of metadata labels.")
  }
  
  # Extract the first two letters after the underscore for each column
  metadata <- as.data.frame(rownames(data))
  metadata$condition <- sapply(strsplit(selected_columns, split), function(x) {
    group <- substr(x[2], 1, 2)
    if (grepl("W", group)) {
      return("WT")
    } else if (grepl("K", group)) {
      return("KO")
    } else {
      return(group)
    }
  })
  # Define conditions based on the patterns
  metadata$Days <- sapply(strsplit(selected_columns, split), function(x) {
    group <- substr(x[2], 1, paterns_num)
    matched <- FALSE
    for (i in seq_along(patterns)) {
      if (grepl(patterns[[i]], group)) {
        matched <- TRUE
        return(metadata_labels[[i]])
      }
    }
    if (!matched) {
      return(group)
    }
  })
  metadata <- metadata[,-1]
  rownames(metadata) <- rownames(data)
  return(metadata)
}

# Function to create boxplots for a given gene
Plot_genes <- function(data, metadata, selected_genes, title = "", xlab = "", 
                       ylab = "", filllab = "") {
  # data : expression data
  # metadata : metadata file with condition information
  # selected_genes : List of genes to plot
  # title : main title of the plot
  # xlab : x-axis label
  # ylab : y-axis label
  # filllab : label for the fill legend
  
  # Initialize a list to store the plots
  plots <- vector("list", length = length(selected_genes))
  
  # Loop to create a plot for each selected gene
  for (i in 1:length(selected_genes)) {
    gene <- selected_genes[i]
    # Check if the gene exists in the 'data' dataframe
    if (!gene %in% colnames(data)) {
      stop(paste("Le gène", gene, "n'existe pas dans le dataframe 'data'."))
    }
    gene_values <- data[[gene]]
    
    # Initialize a list to store the plots
    if (length(unique(metadata)) <= 2) {
      data$Condition <- factor(metadata, levels = c("WT", "KO"))
      plots[[i]] <- data %>%
        #mutate(Condition = fct_reorder(Condition, gene_values, .desc = TRUE)) %>%
        ggplot(aes_string(x = 'Condition', y = gene, fill = 'Condition')) + 
        geom_boxplot() + 
        labs(title = title, x = xlab, y = paste0(ylab, " ", gene), fill = filllab) +
        scale_fill_manual(values = c("#00AFBB", "#E7B800"))
    }else{
      data$Condition <- metadata
      plots[[i]] <- data %>%
        #mutate(Condition = fct_reorder(Condition, gene_values, .desc = TRUE)) %>%
        ggplot(aes_string(x = 'Condition', y = gene, fill = 'Condition')) + 
        geom_boxplot() + 
        labs(title = title, x = xlab, y = paste0(ylab, " ", gene), fill = filllab) 
    }
    
  }
  # Arrange plots in a grid
  grid_plot <- do.call(grid.arrange, c(plots, ncol=2))
  
  return(grid_plot)
}

# Function to plot a heatmap using the genes of a specific module 
heatmap_modules_genes <- function(data, metadata, merge, modnumb = 1, title= " ", cluster_col=T, 
                                  cluster_row=T, show_rownames=T, show_colnames=F, scale = "column", 
                                  cutree_cols = 2, cutree_rows = 2, fontsize=6, fontsize_col = 10, 
                                  fontsize_row = 10){
  # data : expression data
  # metadata : metadata file with condition information
  # merge : List containing modules number, color and wgcna results information
  # modnumb : number of selected module 
  
  # Select as data frame only genes from selected modules
  dataMod = data[,names(data)[merge$colors==modnumb]]
  # Pheatmap plot
  pheatmap <- pheatmap(dataMod, cluster_col=cluster_col, cluster_row=cluster_row, show_rownames=show_rownames, 
                       show_colnames=show_colnames, fontsize=fontsize, annotation_row = metadata, 
                       color = rev(RColorBrewer::brewer.pal(n = 10, name = "RdBu")), scale = scale, 
                       cutree_cols = cutree_cols, cutree_rows = cutree_rows, 
                       fontsize_col = fontsize_col, fontsize_row = fontsize_row, main =paste0(title, " ", modnumb))
  
}

# Define cluster based on the heatmap using cutree
df_cluster_genes <- function(data, metadata, merge, modnumb = 1, title= " ",
                             cutree_cols = 2, cutree_rows = 2, wtcluster = 1, kocluster = 2){
  # data : expression data
  # metadata : metadata file with condition information
  # merge : List containing modules number, color and wgcna results information
  # modnumb : number of selected module 
  # wtcluster : The cluster number correlated with WT condition
  # kocluster : The cluster number correlated with KO condition
  
  # Pheatmap plot
  pheatmap <- heatmap_modules_genes(data = datExpr, merge = merge, metadata = metadata_file, 
                                    modnumb = modnumb,  title = "WT vs KO genes in module")
  dev.off()
  
  df <- data.frame(cluster = cutree(pheatmap$tree_col, k = 2))
  df[df$cluster==kocluster,] <- "KO UP"
  df[df$cluster==wtcluster,] <- "WT UP"
  df$genes <- rownames(df)
  
  return(df)
}

# Plot genes from modules 
Plot_modules_genes <- function(data, metadata, df_cluster, gene_to_remove = NULL, cluster_value,  title= "", xlab="",
                               ylab="", filllab="", ntop=6, order=TRUE) {
  # data : expression data
  # metadata : metadata file with condition information
  # df_cluster: dataframe with 2 columns (genes and the cluster). Cluster obtained from df_cluster_genes function (KO UP and WT UP)
  # gene_to_remove : list of genes to remove. For example, genes beginning with numbers considered as numbers by ggplot
  # cluster_value : specific cluster of interest to plot genes from 
  # title : main title of the plot
  # xlab : x-axis label
  # ylab : y-axis label
  # filllab : label for the fill legend
  # ntop : number of top genes to plot
  # order : whether to order the conditions based on expression levels
  # color_palette : color palette for the boxplots
  # ncol : number of columns in the grid arrangement
  
  # Initialisation d'une liste pour stocker les plots
  plots <- vector("list", length = ntop)
  
  # Filter out the specific genes
  if(is.null(gene_to_remove)){
    df_cluster_filtered <- df_cluster 
  } else {
    df_cluster_filtered <- df_cluster %>% 
      filter(!genes %in% gene_to_remove)
  }
  
  # Extract the list of genes where the cluster is the specified value
  list_genes_cluster <- df_cluster_filtered %>%
    filter(cluster == cluster_value)
  
  # Ensure ntop does not exceed the number of available genes
  ntop <- min(ntop, nrow(list_genes_cluster))
  
  # Ensure metadata and data have the same sample order
  data$Condition <- metadata$condition
  
  # Selection of the top ntop genes
  selected_genes <- rownames(list_genes_cluster)[1:ntop]
  
  # Create plots for each selected gene
  for (i in 1:length(selected_genes)) {
    gene <- selected_genes[i]
    gene_values <- data[[gene]]
    
    # Generate dynamic title if required
    plot_title <- if (dynamic_title) paste("Expression of", gene, "in", cluster_value) else title
    
    # Create the plot for the current gene
    plots[[i]] <- data %>%
      mutate(Condition = fct_reorder(Condition, gene_values, .desc = order)) %>%
      ggplot(aes_string(x = 'Condition', y = gene, fill = 'Condition')) + 
      geom_boxplot() + 
      labs(title = plot_title, x = xlab, y = paste0(ylab, " ", gene), fill = filllab) +
      scale_fill_manual(values = color_palette)
  }
  
  # Arrange plots in a grid with specified number of columns
  grid_plot <- do.call(grid.arrange, c(plots, ncol=ncol))
  
  return(grid_plot)
}


# Functions
## Genelist generation for GSEA analysis
Save_Genelist <- function(data, logFC, Module=NULL){
  # data : expression data 
  # LogFC : the log2fc column
  geneList <- logFC
  names(geneList) <- rownames(data)
  if(is.null(Module)){
    geneList <- na.omit(geneList)
    geneList <- sort(geneList, decreasing = TRUE)
  }else{
    geneList <- geneList[names(geneList) %in% names(Module)]
    geneList <- na.omit(geneList)
    geneList <- sort(geneList, decreasing = TRUE)
  }
}

## Genelist generation for GSEA analysis
Save_GSEA_Rnk <- function(data, logF, path = " ", file.name = " "){
  # data : expression data 
  # LogFC : the log2fc column
  # path : path for output
  # file.name : The output file name
  geneList <- logFC
  names(geneList) <- rownames(data)
  geneList <- na.omit(geneList)
  geneList <- sort(geneList, decreasing = TRUE)
  
  data = data.frame(names(geneList),geneList)
  names(data) = c("Genes", "Rank")
  
  write.table(data, paste0(path, file.name), sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
}

## Save gene expression in txt format for GSEA
Save__GSEA_Expr <- function(data, path = " ", file.name = " ", samples.position = c("row", "column")){
  # data : is a data frame with genes samples. Specify the position of each one in the function.
  # file.name : The name of the output data
  # samples.position : if samples are in row or in column in the data frame
  if (samples.position=="row"){
    data = data.frame(t(data))
  }else if (samples.position=="column"){
    data = data.frame(data)
  }else{
    stop(paste("La position", samples.position, "n'est pas correcte. Choisissez entre 'row' et 'column'"))
  }
  columns = names(data)
  data$desc = NA
  data$NAME = rownames(data)
  data = data %>%
    dplyr::select(NAME, desc, all_of(columns))
  genes = nrow(data)
  samples = length(columns)
  # Les nouvelles lignes à ajouter
  new_lines <- c("#1.2", "", paste0(genes,"\t", samples), "")
  # Créer un nouveau fichier avec les nouvelles lignes en haut
  file.name = paste0(file.name, ".gct.txt")
  writeLines(new_lines, paste0(path, file.name))
  # Ajouter le contenu existant en dessous avec des tabulations
  write.table(data, paste0(path, file.name), sep = "\t", row.names = FALSE, col.names = TRUE, append = TRUE, quote = FALSE)
}

# Selection of diff
Select_diff <- function(data, patterns = "CW.vs.CK", filter = TRUE) {
  deseq2_diff <- data %>%
    dplyr::select(gene_symbol,
                  paste0("diffexp_log2fc_", patterns),
                  paste0("diffexp_deseq2_pvalue_", patterns),
                  paste0("diffexp_deseq2_qvalue_", patterns)
    )
  # Renommer les colonnes
  colnames(deseq2_diff) <- c("gene_symbol", "Log2FC", "Pvalue", "Qvalue")
  # Filtrer les lignes où Log2FC n'est pas NA
  deseq2_diff <- deseq2_diff[!is.na(deseq2_diff$Log2FC), ]
  if(filter==TRUE){
    deseq2_diff <- deseq2_diff[deseq2_diff$Pvalue < 0.05,]
  }else{
    deseq2_diff <- deseq2_diff
  }
  # Définir les noms de lignes
  rownames(deseq2_diff) <- deseq2_diff$gene_symbol
  
  return(deseq2_diff)
}

# Plot the heatmap of form GSEA results term.
heatmap_pathway_genes <- function(data, gsea_summary, pathway = "nucleoplasm"){
  # data : expression data
  # gsea_summary : is the summary of gsea results. Ex : gseaDataGO <- summary(gseaData)
  # pathway : The pathway of interest.
  pathway.selected = gsea_summary[gsea_summary$Description==pathway,]
  pathway.selected.genes <- pathway.selected$core_enrichment
  # Séparer les gènes en utilisant le délimiteur "/"
  gene_list <- unlist(strsplit(pathway.selected.genes, "/"))
  gene_list <- sort(gene_list)
  pathway.selected.data <- data[, colnames(data) %in% gene_list]
  pheatmap(pathway.selected.data, scale = "column", main = paste0("Genes in ", pathway, " pathway"))
}

# Heatmap function for plotting the regulatory network from GENIE3
make_heatmap_ggplot <- function(data, y_name, x_name, ntop=30, x_axis_position = "top", legend_position = "top", 
                                color = "blue", legend_title = "score", ...){
  # data : the results of getLinkList function of GENIE3
  if((x_axis_position %in% c("top","bottom")) == FALSE)
    stop("x_axis_position should be top or bottom")
  if((legend_position %in% c("top","bottom","left","right","none")) == FALSE)
    stop("legend_position should be top, bottom, left, right or none")
  if(!is.character(color) | length(color) != 1)
    stop("color should be character vector of length 1")
  
  requireNamespace("dplyr")
  requireNamespace("ggplot2")
  
  # Ensure the data is ordered by the score
  data <- data %>% arrange(desc(weight))
  
  # Select the top n rows
  data_top <- data[1:ntop, ]
  
  # Plot the data
  plot_object <- ggplot(data_top, aes(x = targetGene, y = regulatoryGene, fill = weight)) +
    geom_tile(color = "white", size = 0.5) +
    scale_fill_gradient(low = "whitesmoke", high = color) +
    theme_minimal()
  
  plot_object <- plot_object + theme(
    panel.grid.minor = element_line(color = "transparent"),
    panel.grid.major = element_line(color = "transparent"),
    legend.position = legend_position,
    axis.ticks = element_line(size = 0),
    axis.text.x.top = element_text(angle = 90, hjust = 0, ...),
    axis.text.x = element_text(angle = 90, hjust = 1, ...),
    axis.title = element_text(...),
    axis.text.y = element_text(...)
  )
  
  plot_object <- plot_object +
    scale_x_discrete(position = x_axis_position) +
    xlab(x_name) +
    ylab(y_name) +
    labs(fill = legend_title)
  
  return(plot_object)
}

# Bar Plot of genes regulatory network results from GENIE3
Plot_regulatory_network <- function(data, n = 50, xlab = "GENIE3 Weight rank", 
                                    ylab = "Regulatory -> Target Genes Interaction", color = "purple"){
  
  # data : the results of getLinkList function of GENIE3
  # n : The number of relation to plot
  
  # filter results to top N interactions
  setected_data <- data %>%
    #arrange(weight) %>%
    slice_head(n = n) %>%
    mutate(id = fct_inorder(paste0(regulatoryGene, " -> ", targetGene)))
  
  # visualize median rank
  setected_data %>%
    mutate(id = fct_reorder(id, weight, .desc = FALSE)) %>%
    ggplot(aes(y = id, x = weight, fill = weight)) +
    geom_bar(stat = "identity") + 
    scale_fill_gradient(low = "whitesmoke", high = color) +
    xlab(xlab) + ylab(ylab) +
    theme_cowplot() +
    theme(axis.text.x = element_text(size = 8, angle = 60, hjust = 1, vjust = 1))
}





