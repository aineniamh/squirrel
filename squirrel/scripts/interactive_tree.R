library(ape)
library(stringr)
library(scales)
library(htmlwidgets)

# Silently to not clutter user-visible output
suppressMessages(library(treeio))
suppressMessages(library(plotly))
suppressMessages(library(ggtree))
suppressMessages(library(dplyr))

#' Annotate a treeio object with mutation quality and path data.
#' @param t_data    A tree object from treeio
#' @param mut_file  The *.aln.tree.amino_acid.reconstruction.csv summarizing mutations
#' @return Inputted tree object from treeio with 'data' attribute overwritten
annotate_tree <- function(t_data, mut_df) {

  # Read in tree and assign node labels
  t <- treeio::get.tree(t_data)
  t$node.label <- treeio::get.data(t_data) %>%
    dplyr::arrange(node) %>% 
    dplyr::pull(label)

  # Initialize table for summary data
  labels <- c(t$tip.label, t$node.label) %>% unname()
  root_node <- length(t$tip.label)+1
  tdf <- dplyr::tibble(
    label = labels,
    node = seq_along(labels)
  )
  tdf$edge_index <- sapply(tdf$node, function(i){
    edge_index <- which(tdf$node[i] == t$edge[,2])
    if (length(edge_index) == 0) {edge_index <- NA}
    return(edge_index)
  })%>% unname()
  

  # Establish path of labels visited from this node to root
  tdf$parent_path <- sapply(tdf$node, function(x){
    p <- nodepath(t, root_node, x)
    fixed_labels <- gsub('\\|', '\\\\|', labels[p])
    grep_res <- paste(fixed_labels, collapse = '$|^')
    grep_res <- paste0('^', grep_res, '$')
    return(grep_res)
  }) %>% unname()
  

  # Count cumulative apobec mutations along path
  # Also obtain a ratio to cumulative non-apobec mutations
  tdf$apobec_count <- sapply(tdf$parent_path, function(x){
    dplyr::filter(mut_df, grepl(x, child) & (apobec == 'True')) %>%
      nrow()
  }) %>% unname()

  tdf$total_muts <- sapply(tdf$parent_path, function(x){
    filter(mut_df, grepl(x, child)) %>%
      nrow()
  }) %>% unname()

  tdf <- tdf %>% 
    dplyr::mutate(apobec_ratio = round((apobec_count/total_muts)*100)) %>%
    dplyr::mutate(apobec_ratio = paste0(apobec_ratio, '%'))

  # Quality checking with grantham score
  tdf$mean_grantham <- sapply(tdf$parent_path, function(x){
    filter(mut_df, grepl(x, child)) %>%
      filter(!is.na(score)) %>%
      pull(score) %>%
      mean() %>%
      round()
  }) %>% unname()

  tdf$na_grantham <- sapply(tdf$parent_path, function(x){
    mut_sdf <- filter(mut_df, grepl(x, child))
    ratio <- nrow(filter(mut_sdf, is.na(score)))/nrow(mut_sdf)
    ratio <- paste0(round(ratio*100), '%')
    return(ratio)
  }) %>% unname()


  # Need some manual label fixing for ggplotly tiplab workaround
  tdf$print_label <- sapply(tdf$node, function(x){
    ifelse(x > length(t$tip.label), '', t$tip.label[x])
  })

  # Overwrite tree metadata and get range of apobec vals
  attr(t_data, 'data') <- tdf

  # Also store tree height as a simple attribute
  attr(t_data, 't_height') <- max(node.depth.edgelength(t))

  return(t_data)
}

#' Generate interactive, heated, tree images from squirrel output
#' @param tree_file The *aln.tree file from iqtree
#' @param mut_file  The *.suggested_mask.csv summarizing mutations
#' @return A .html Widget
plot_interactive_tree <- function(tree_file, mut_file) {
  
  # Read in mutation data and extend table to map to tree
  mut_df <- read.csv(mut_file)
  t_data <- treeio::read.beast(tree_file) %>%
    annotate_tree(mut_df)
  tdf <- t_data %>% get.data()

  apobec3_range <- t_data %>%
    treeio::get.data() %>% 
    dplyr::pull(apobec_count) %>%
    range()

  # Dummy Aesthetic set 
  # Mostly to easily define tooltip inclusion
  set_aes <- aes(
    color = apobec_count,
    total_muts = total_muts,
    apobec_ratio = apobec_ratio,
    mean_grantham = mean_grantham,
    na_grantham = na_grantham,
    label = label.y
  )

  # Initial ggtree object
  p <- ggtree::ggtree(t_data, set_aes) +
    ggplot2::geom_text(aes(label=print_label)) +
    ggplot2::theme(
      panel.background = element_rect(fill = "grey90"),
      panel.border = element_rect(linewidth = 1, fill=NA),
      axis.text.x = element_text(size = 12),
      axis.ticks.x = element_line(),
      axis.line.x = element_line(),
      panel.grid.major.x = element_line(colour = 'white', linewidth = 1),
      panel.grid.minor.x = element_line(colour = 'white', linewidth = 0.5, linetype = 1)
    ) +
    ggplot2::scale_color_gradient(
      name='APOBEC3 Mutations', limits=apobec3_range,
      oob=scales::squish,
      low='blue', high='orange2'
    ) 

  # Select tooltip from set_aes
  tt_cols <- c('label', 'apobec_count', 'total_muts', 'apobec_ratio', 'mean_grantham', 'na_grantham')
  
  # Automatically shift x limits in case of long names
  t_height <- attr(t_data, 't_height')
  buffer   <- max(nchar(tdf$label))*(t_height/180)
  x_min <- -(t_height/100)
  x_max <- t_height + buffer

  # ggplotly html widget creation 
  # Some params are modified at this level
  widget <- p %>% 
    plotly::ggplotly(tooltip = tt_cols) %>%
    plotly::layout(xaxis = list(range = c(x_min, x_max))) %>%
    plotly::style(textposition = "right")

  return(widget)
}

# Run and output command. To avoid too many R packages, sticking with base commandArgs()
# Arg 1 - Pruned tree_file from IQTREE
# Arg 2 - Mutation data from reconstruction phase
# Arg 3 - Output file path
message('Building R-Based interactive HTML outputs')
args  <- commandArgs(trailingOnly = TRUE)
fname <- args[3]
plot_interactive_tree(args[1], args[2]) %>% 
  htmlwidgets::saveWidget(fname)