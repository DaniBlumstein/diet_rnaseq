# this first function runs a Walds test between specified groups, then returns a dataframe from it.
getDEdf <- function(dds, treatment, group1, group2){
  res <- results(dds, contrast = c(treatment, group1, group2))
  res.df <- setDT(as.data.frame(res), keep.rownames = "transcript")
  # add gene annotation information
  res.df <- dplyr::left_join(res.df, annos, by = "transcript")
  return(res.df)
}

# this function pulls out particular a priori genes of interest into a new data frame
aprioris <- function(df, annos){
  df.aprioris <- data.frame()
  for (i in 1:length(colors$gene_name)){
    searchterm <- paste0("\\b", colors$gene_name[i], "\\b")
    searchterm <- tolower(searchterm)
    tmp <- df %>% filter(str_detect(Gene, searchterm))
    df.aprioris <- rbind(df.aprioris, tmp)
  }
  return(df.aprioris)
}


# this function extracts all genes below a specified alpha value, prints the number of DE genes to screen, and saves a spreadsheet to file
siggies <- function(df, alpha, csv){
  sigs <- df %>% filter(padj < alpha)
  print(paste0("Number of significant genes: ", length(sigs$padj)))
  write.csv(sigs, csv, row.names = FALSE)
  return(sigs)
}


# significant gene wrapper function
SigGeneWrapper <- function(model_output, alpha, comparison){
  print(paste0("Running ", comparison, " with alpha = ", alpha))
  df <- setDT(as.data.frame(model_output), keep.rownames = "Gene")
  # add annotation data
  #df1 <- dplyr::left_join(df, annos, by = "transcript") 
  # get significant genes
  print("Overall significant genes")
  sigs <- siggies(df, alpha, paste0("results/", comparison, "_genes.csv"))
  # add annotation data
  colordf <- aprioris(df, annos)
  #get significant color genes
  print("Significant color genes")
  color.sigs <- siggies(colordf, alpha, paste0("results/", comparison, "_colorgenes.csv"))
}



### Functions for plotting

# Define ggplot aesthetic.
mytheme <- theme_classic() + theme(axis.text=element_text(size=20), 
                                   axis.title=element_text(size=22,face="bold"), 
                                   plot.title = element_text(size = 22, hjust = 0.5), strip.text.x = element_text(size = 22),
                                   panel.border = element_rect(color = "black", fill = NA, size = 1), legend.position = "none")



# pull in DE genes that were saved as .csvs
DEgenes <- function(comparison, aretheycolors){
  genes <- read.csv(paste0("results/", comparison, "_", aretheycolors, ".csv"))
  #genes <- unique(tmp$Gene)
  return(genes)
}

target <- 7
#### extract the data from the DESeq models
#  requires a trycatch for the plotcounts function
trycatch_plotCounts = function(dds, target) {
  tryCatch(plotCounts(dds, gene = target, intgroup = c("tissue", "trt"), returnData = TRUE, xlab="treatment"),
           warning = function(w) {},
           error = function(e) {print("transcript missing")}) 
}

# genes = list of genes
# dds = deseq2 object
# annos = transcript 2 gene name df
# sampledata = sample ... data ...
extractTadCounts <- function(genes, dds, annos, sampledata){
  # make gene a character
  genes$Gene <- as.character(genes$Gene)
  #### for loop to extract counts
  # make empty data frame
  countdf <- data.frame()
  for (i in 1:nrow(genes)){
    #target <- as.character(genes$transcript[i])
    gene <- genes$Gene[i]
    #print(target)
    print(gene)
    # pull out count data
    tmp <- trycatch_plotCounts(dds, gene)
    if (tmp == "transcript missing") {
      rm(tmp) 
    } else {
      # make rownames a column
      tmp <- setDT(tmp, keep.rownames = TRUE)
      # change name of new first column
      colnames(tmp)[1] <- "sample"
      # add gene name to a column
      tmp$Gene <- gene
      # add transcript id to a column
      #tmp$transcript<- target
      # add all to a data frame
      countdf <- rbind(countdf, tmp)     
    } # end of if statement
  } # end of loop
  # add a column for species
  #colnames(countdf)[1] <- "Name"
  info <- sampledata[,c("sample","tissue")]
  countdf <- dplyr::left_join(countdf, info, by = "sample")
  return(countdf)
}


extractCounts <- function(genes, dds, annos, sampledata){
  # gene2tx <- annos[annos$Gene %in% genes, ]
  #### for loop to extract counts
  # make empty data frame
  countdf <- data.frame()
  for (i in 1:length(genes)){
    genes <- as.character(genes) ######## the bug lives somewhere in here ###################
    #gene <- as.numeric(genes[i]) ######## the bug lives somewhere in here ###################
    #gene_symbol <- genes$gene_symbol[i]
    #print(target)
    print(genes)
    # pull out count data
    tmp <- trycatch_plotCounts(dds, genes)
    if (tmp == "transcript missing") {
      rm(tmp) 
    } else {
      # make rownames a column
      tmp <- setDT(tmp, keep.rownames = TRUE)
      # change name of new first column
      colnames(tmp)[1] <- "individal"
      # add gene name to a column
      tmp$Gene <- genes
      # add transcript id to a column
      #tmp$transcript<- target
      # add gene symbol to a column
      #tmp$gene_symbol <- gene_symbol
      # add all to a data frame
      countdf <- rbind(countdf, tmp)     
    } # end of if statement
  } # end of loop
  # add a column for tissue
  #colnames(countdf)[1] <- "Name"
  info <- sampledata[,c("individal","sex")]
  countdf <- dplyr::left_join(countdf, info, by = "individal")
  return(countdf)
}



# plot expression. This is for each species individually, and will plot all the genes in a dataframe
plotGenestime <- function(df, dir){
  dir.create(dir)
  genes <- unique(df$Gene)
  for (i in 1:length(genes)){
    ith_gene <- genes[i]
    #gene <- df$Gene[i]
    #gene_df <- df[,i]
    gene_df <- df %>% filter(Gene == ith_gene)
    gene <- gene_df$Gene[1]
    # remove problematic characters that break the script
    gene <- gsub(pattern = "/", replacement = " ", gene)
    gene <- gsub(pattern = ":", replacement = " ", gene)
    # this is a test...
    gene <- gsub(pattern = "\\s\\(.*\\)$", replacement = "", gene)
    print(gene)
    
    ggplot(gene_df,  aes(x = age_weeks, y = count)) + 
      ggtitle(gene,  subtitle = "")  + 
      ylab("Normalized Counts") + 
      xlab("Age (weeks)") +
      geom_boxplot(outlier.shape = NA, lwd=.75, aes(fill = age_weeks)) + 
      facet_grid(. ~ morph) +
      mytheme
    
    ggsave(paste0(dir, "/", gene, ".png"), width = 11.38, height = 4.37)
  }
}



# extract the transcript ids from a list of genes. Supply gene names + annotation document
targetedgenes <- function(genes, annos){
  df.targetedgenes <- data.frame()
  for (i in 1:length(genes)){
    searchterm <- paste0("\\b", genes[i], "\\b")
    searchterm <- tolower(searchterm)
    tmp <- annos %>% filter(str_detect(Gene, searchterm))
    tmp$gene_symbol <- genes[i]
    df.targetedgenes <- rbind(df.targetedgenes, tmp)
  }
  return(df.targetedgenes)
}


#Needs a new plotting function:


#  ```{R, gene specific plots}
GeneFigures <- function(gene_df, dir, title){
  dir.create(dir)
  #figs <- unique(df$Gene_name)
  # for (i in 1:length(figs)){
  #   gene <- figs[i]
  #   gene_df <- df[df$Gene_name %in% gene, ]
  
  # order facets by morph
  #gene_df$morph = factor(gene_df$morph, levels=c("banded", "redheaded", "striped", "spotted"))
  # For better plots/visualizations remove all rows with expression < 1 count
  gene_df <- gene_df %>% filter(count > 1)
  
  # plot
  ggplot(gene_df,  aes(x = trt, y = count)) + 
    ggtitle(title)  + 
    mytheme +
    ylab("Normalized Counts") + 
    xlab("trt") +
    geom_boxplot(outlier.shape = NA, lwd=.75, aes(fill = tissue)) + 
    geom_point(position = "jitter") +
    facet_grid(Gene ~ tissue, scales = "free_y") +
    #scale_fill_manual(values = c("orange1", "red1", "yellow2", "chartreuse4")) +
    #theme(strip.text.x = element_blank(), strip.background = element_blank()) 
    theme_classic() + theme(axis.text=element_text(size=20), 
                            axis.title=element_text(size=22,face="bold"), plot.title = element_text(size = 22, hjust = 0.5), 
                            strip.text.y = element_text(size = 22), strip.text.x = element_text(size = 22),
                            panel.border = element_rect(color = "black", fill = NA, size = 1), legend.position = "none")
  
  # }
  ggsave(paste0(dir, "/", title, ".png"), width = 11.38, height = (3 * length(unique(gene_df$Gene))))
}
# quick wrapper to pull in data from multiple genes for both species. Output from this will go into figure making code.
# i think this doesnt work, scrap it.
extractAllCounts <- function(genes, dds1, dds2, annos, sampledata){
  df1 <- extractTadCounts(genes, dds1, annos, sampledata)
  df2 <- extractTadCounts(genes, dds1, annos, sampledata)
  combined.df <- rbind(df1, df2)
  
  return(combined.df)
}



#dani functions:
#PCA graphing function, note there is hard coding in here with the groups.
PCA_graph <- function(dds, tissue)
{
  vsd <- varianceStabilizingTransformation(dds, blind = FALSE)
  # get PC1 and PC2 data
  pcaData <- plotPCA(vsd, intgroup = c("sex", "trt_combo", "lane"), returnData = TRUE)
  # get percent variation
  percentVar <- round(100 * attr(pcaData, "percentVar"))
  # pca code
  ggplot(pcaData, aes(x = PC1, y = PC2, color = trt_combo, shape = trt_combo, name=name)) +
    stat_ellipse(aes(group = trt_combo, linetype = trt_combo), type = "t", level = 0.95, size = 1.25, show.legend = FALSE) +
    #scale_linetype_manual(values=c("twodash", "longdash", "solid"), guide = FALSE) +
    geom_point(size = 3, show.legend = TRUE) + 
    #scale_color_manual(values = c("orange1", "red1", "yellow2", "chartreuse4")) +
    #scale_shape_manual(values = c(8,9,15,16,10,17,18,11,13)) + 
    xlab(paste0("PC1: ", percentVar[1], "% variance")) +
    ylab(paste0("PC2: ", percentVar[2], "% variance")) +
    coord_fixed() + theme_bw() +
    ggtitle(tissue)+
    guides(shape = guide_legend(order = 1),color = guide_legend(order = 2)) 
}
