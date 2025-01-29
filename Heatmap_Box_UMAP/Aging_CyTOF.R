rm(list = ls())
library(reshape2)
library(randomcoloR)
library(pals)
library(ggplot2)
library(Hmisc)
library(stringr)
library(ComplexHeatmap)
library(circlize)
library(corrplot)
library(readxl)
library(ggridges)
library(ggpubr)
library(raster)
library(matrixStats)
library(limma)

####READ and CLUSTER FUNCTIONS####
returnfcs <- function(FDR_cutoff=.05,
                      metaDataFile='~/',
                      panelDataFile='~/',
                      dataDirectory='~/',
                      shape_timepoint=NULL,
                      color_timepoint=NULL){
  #This function generates an fcs file, subtype_markers, colors and shapes for clustering 
  require(scales);require(readxl);require(dplyr);require(flowCore)
  ##directory and metadatafile checking
  if(!dir.exists(dataDirectory)) {stop('ERR: cannot find data directory')}
  if(!file.exists(metaDataFile)) {stop('ERR: cannot find metadata.xlsx or .csv file')}
  ##readin metadata and clean
  ifelse(grepl(metaDataFile,pattern='.xls'),md <- read_excel(metaDataFile),md <- read.csv(metaDataFile,header = TRUE))#must be in xl format or csv
  md$batch <- factor(md$Batch)
  md$timepoint <- factor(md$timepoint)
  md$run <- factor(md$Run)
  md$cancer_type <- factor(md$cancer_type)
  md$match_id <- factor(md$match_id)
  md$response <- factor(md$response)
    
  rownames(md) = md$sample_id;md$sample_id <- md$sample_id
  #Make sure all files in metadata present in datadirectory
  if(!all(md$file_name %in% list.files(dataDirectory)[grep(list.files(dataDirectory),pattern = '.fcs')])){
    print(paste('ERR: not all filenames in metadata present in data folder - missing',md$file_name[!which(md$file_name %in% list.files(dataDirectory)[grep(list.files(dataDirectory),pattern = '.fcs')])],'Subsetting...'))
    md <- md[-c(!which(md$file_name %in% list.files(dataDirectory)[grep(list.files(dataDirectory),pattern = '.fcs')])),]
  }
  ##Define shapes for timepoint
  if(is.null(shape_timepoint)){shape_timepoint <- c(0:25)[1:length(levels(md$timepoint))]}#can specify as long as number is same
  if(length(shape_timepoint)!=length(levels(md$timepoint))){stop(paste0('ERR no. shapes specified is less than no. of timepoint (',length(levels(md$timepoint)),')'))}
  names(shape_timepoint) <- levels(md$timepoint)
  ## Define colors for the timepoint
  if(is.null(color_timepoint)){color_timepoint <- hue_pal()(length(levels(md$timepoint)))}#can specify as long as number is same
  if(length(color_timepoint)!=length(levels(md$timepoint))){stop(paste0('ERR no. shapes specified is less than no. of timepoint (',length(levels(md$timepoint)),')'))}
  ## read fcs
  fcs_raw <- read.flowSet(md$file_name, path = dataDirectory, transformation = FALSE, truncate_max_range = FALSE)
  #sample_ids <- rep(md$sample_id, fsApply(fcs_raw, nrow))
  panel <- read_excel(panelDataFile)
  head(data.frame(panel))
  ## Replace problematic characters
  panel$Metal <- gsub('-', '_', panel$Metal)
  panel_fcs <- pData(parameters(fcs_raw[[1]]))
  panel_fcs$desc <- gsub('-', '_', panel_fcs$desc)
  panel_fcs$desc[is.na(panel_fcs$desc)] <- paste0('NA_',which(is.na(panel_fcs$desc))) #was labelled 'todo'(mapping based on isotope for now just getting rid of NA keeping rownum) 
  # use panel$Antigen to fix description in panel_fcs
  # use metal+isotope as mapping between panel from xlsx and panel from the fcs files
  rownames(panel_fcs) = panel_fcs$name
  panel_fcs[paste0(panel$Metal,panel$Isotope,'Di'),2] <- panel$Antigen
  ## Replace paramater data in flowSet
  pData(parameters(fcs_raw[[1]])) <- panel_fcs
  ## Define variables indicating marker types
  subtype_markers <- panel$Antigen[panel$Subtype == 1]
  functional_markers <- panel$Antigen[panel$Functional == 1]
  if(!all(subtype_markers %in% panel_fcs$desc)){stop('ERR: Not all subtype_markers in panel_fcs$desc (isotopes)')}
  if(!all(functional_markers %in% panel_fcs$desc)){stop('ERR: Not all functional_markers in panel_fcs$desc (isotopes)')}
  ## arcsinh transformation and column subsetting
  fcs <- fsApply(fcs_raw, function(x, cofactor = 5){
    colnames(x) <- panel_fcs$desc
    expr <- exprs(x)
    expr <- asinh(expr[, union(subtype_markers,functional_markers)] / cofactor)
    expr[!is.finite(expr)] <- NA #convert inf to NA
    expr<-na.omit(expr) #remove NA
    exprs(x) <- expr
    x
  })
  sample_ids <- rep(md$sample_id, fsApply(fcs, nrow))
  return(list('fcs'=fcs,
              'subtype_markers'=subtype_markers,
              'functional_markers'=functional_markers,
              'shape_timepoint'=shape_timepoint,
              'color_timepoint'=color_timepoint,
              'sample_ids'=sample_ids,
              'meta_data'=md))
}

clusterfcs <- function(fcs=output$fcs,
                       subtype_markers = output$subtype_markers,
                       seed=1234,plottitle='consensus_plots',
                       numclusters=40){
  ## Cell population identification with FlowSOM and ConsensusClusterPlus
  require(dplyr);require(FlowSOM);require(ConsensusClusterPlus)
  set.seed(seed)
  som <- ReadInput(fcs, transform = FALSE, scale = FALSE) %>% BuildSOM(colsToUse = subtype_markers)
  ## Get the cell clustering into 100 SOM codes
  cell_clustering_som <- som$map$mapping[,1]
  ## Metaclustering into numclusters with ConsensusClusterPlus
  codes <- som$map$codes
  mc <- ConsensusClusterPlus(t(codes), maxK = numclusters, reps = 100,
                             pItem = 0.9, pFeature = 1, title = plottitle, 
                             plot = "png", clusterAlg = "hc", 
                             innerLinkage = "average", finalLinkage = "average",
                             distance = "euclidean", seed = 1234)
  
  ## Get cluster ids for each cell
  code_clustering <- mc[[numclusters]]$consensusClass#metaclusters consensus
  cell_clustering <- code_clustering[cell_clustering_som]#cell clustering from som
  return(list('code_clustering'=code_clustering,'cell_clustering'=cell_clustering,'metaclusters'=mc))
}

####CLUSTER HEATMAP FUNCTIONS ####
plot_clustering_heatmap_wrapper <- function(fcs, cell_clustering, nclusters=40,
                                            color_clusters=clustercolors, cluster_merging = NULL, 
                                            subtype_markers,
                                            clusterMergeFile=NULL,
                                            fileName = 'clusteringheatmap.pdf'){
  require(matrixStats);require(dplyr);require(RColorBrewer);require(pheatmap);require(readxl);require(flowCore);require(scales)
  ## Will output the heatmap object and print it 
  #if((color_clusters)=='auto'){color_clusters <- hue_pal()(nclusters)}
  #get expression
  expr <- fsApply(fcs, exprs);expr <-expr[,subtype_markers]
  ## Scale expression of all markers to values between 0 and 1
  rng <- colQuantiles(expr, probs = c(0.01, 0.99))
  expr01 <- t((t(expr) - rng[, 1]) / (rng[, 2] - rng[, 1]))
  expr01[expr01 < 0] <- 0; expr01[expr01 > 1] <- 1;expr01 <-expr01[,subtype_markers]
  ## Calculate the mean expression##################################################
  pdf(fileName, width=8, height=11) 
  expr_mean <- data.frame(expr, cell_clustering = cell_clustering, check.names = FALSE) %>%
    group_by(cell_clustering) %>% summarize_all(funs(mean))
  expr01_mean <- data.frame(expr01, cell_clustering = cell_clustering, check.names = FALSE) %>%
    group_by(cell_clustering) %>% summarize_all(funs(mean))
  
  ## Calculate cluster frequencies
  
  clustering_table <- as.numeric(table(cell_clustering))
  clustering_prop <- round(clustering_table / sum(clustering_table) * 100, 2)
  
  ## Sort the cell clusters with hierarchical clustering
  
  d <- dist(expr_mean[, colnames(expr)], method = "euclidean")
  cluster_rows <- hclust(d, method = "average")
  
  expr_heat <- as.matrix(expr01_mean[, colnames(expr01)])
  rownames(expr_heat) <- expr01_mean$cell_clustering
  
  ## Colors for the heatmap
  
  color_heat <- colorRampPalette(brewer.pal(n = 9, name = "YlOrBr"))(100)
  #legend_breaks = seq(from = 0, to = 1, by = 0.2)
  labels_row <- paste0(expr01_mean$cell_clustering, " (", clustering_prop ,
                       "%)")
  
  ## Annotation for the original clusters
  
  annotation_row <- data.frame(Cluster = factor(expr01_mean$cell_clustering))
  rownames(annotation_row) <- rownames(expr_heat)
  color_clusters1 <- color_clusters[1:nlevels(annotation_row$Cluster)]
  names(color_clusters1) <- levels(annotation_row$Cluster)
  annotation_colors <- list(Cluster = color_clusters1)
  
  ## Annotation for the merged clusters
  
  if(!is.null(clusterMergeFile)){
    ifelse(grepl(clusterMergeFile,pattern='.xls'),cluster_merging <- read_excel(clusterMergeFile),cluster_merging <- read.csv(clusterMergeFile,header = TRUE))
    cluster_merging$new_cluster <- factor(cluster_merging$new_cluster)
    annotation_row$Merged <- cluster_merging$new_cluster
    color_clusters2 <- color_clusters[1:nlevels(cluster_merging$new_cluster)]
    names(color_clusters2) <- levels(cluster_merging$new_cluster)
    annotation_colors$Merged <- color_clusters2
  }
  
  p <- pheatmap(expr_heat, color = rev(colorRampPalette(brewer.pal(n = 11, name = "RdYlBu"))(100)), 
                cluster_cols = T,
                cluster_rows = T, 
                labels_row = labels_row,
                #scale="column",
                display_numbers = FALSE, number_color = "black",
                fontsize = 9, fontsize_number = 6,  
                #legend_breaks = legend_breaks,
                annotation_row = annotation_row, 
                annotation_colors = annotation_colors,
                cellwidth = 8,
                cellheight = 8
  )
  dev.off() 
  print('Colors:')
  print(color_clusters)
  print(p);return(p)
}

plot_clustering_heatmap_wrapper2 <- function(fcs, cell_clustering, nclusters=40,
                                             color_clusters=clustercolors,
                                             subtype_markers,
                                             fileName = 'clusteringheatmap.pdf'){
  require(matrixStats);require(dplyr);require(RColorBrewer);require(pheatmap);require(readxl);require(flowCore);require(scales)
  ## Will output the heatmap object and print it 
  #if((color_clusters)=='auto'){color_clusters <- hue_pal()(nclusters)}
  #get expression
  expr <- fsApply(fcs, exprs);expr <-expr[,subtype_markers]
  ## Scale expression of all markers to values between 0 and 1
  rng <- colQuantiles(expr, probs = c(0.01, 0.99))
  expr01 <- t((t(expr) - rng[, 1]) / (rng[, 2] - rng[, 1]))
  expr01[expr01 < 0] <- 0; expr01[expr01 > 1] <- 1;expr01 <-expr01[,subtype_markers]
  ## Calculate the mean expression##################################################
  pdf(fileName, width=8, height=11) 
  expr_mean <- data.frame(expr, cell_clustering = cell_clustering, check.names = FALSE) %>%
    group_by(cell_clustering) %>% summarize_all(funs(mean))
  expr01_mean <- data.frame(expr01, cell_clustering = cell_clustering, check.names = FALSE) %>%
    group_by(cell_clustering) %>% summarize_all(funs(mean))
  
  ## Calculate cluster frequencies
  
  clustering_table <- as.numeric(table(cell_clustering))
  clustering_prop <- round(clustering_table / sum(clustering_table) * 100, 2)
  
  ## Sort the cell clusters with hierarchical clustering
  
  d <- dist(expr_mean[, colnames(expr)], method = "euclidean")
  cluster_rows <- hclust(d, method = "average")
  
  expr_heat <- as.matrix(expr01_mean[, colnames(expr01)])
  rownames(expr_heat) <- expr01_mean$cell_clustering
  
  ## Colors for the heatmap
  
  #legend_breaks = seq(from = 0, to = 1, by = 0.2)
  #labels_row <- expr01_mean$cell_clustering
  
  labels_row <- paste0(expr01_mean$cell_clustering, " ")
  
  ## Annotation for the original clusters
  
  annotation_row <- data.frame(Cluster = factor(expr01_mean$cell_clustering))
  rownames(annotation_row) <- rownames(expr_heat)
  color_clusters1 <- color_clusters[1:nlevels(annotation_row$Cluster)]
  names(color_clusters1) <- levels(annotation_row$Cluster)
  annotation_colors <- list(Cluster = color_clusters1)
  
  p <- pheatmap(expr_heat, 
                color = c(rep(magma(100)[1],25),magma(100)[1:100]), 
                cluster_cols = T,
                cluster_rows = F, 
                labels_row = labels_row,
                #scale="column",
                display_numbers = F, 
                number_color = "black",
                fontsize = 9, fontsize_number = 6,  
                #legend_breaks = legend_breaks,
                annotation_row = annotation_row, 
                annotation_colors = annotation_colors,
                cellwidth = 8,
                cellheight = 8,
                border_color = "black",
                annotation_legend = F
  )
  dev.off() 
  print('Colors:')
  print(color_clusters)
  print(p);return(p)
}

####DIAGNOSTICS####
makeDiagnosticPlots = function(exprData, 
                               md = output$meta_data,
                               sample_ids = output$sample_ids,
                               fcs = output$fcs,
                               subtype_markers = output$subtype_markers,
                               color_conditions = clustercolors,
                               shape_conditions = c(1:13),
                               fileName = 'diagnostics.pdf', 
                               tit = '', 
                               fun = mean)
{
  pdf(file = fileName)
  
  # plot 1
  ggdf <- data.frame(sample_id = sample_ids, exprData)
  ggdf <- melt(ggdf, id.var = 'sample_id', value.name = 'expression', 
               variable.name = 'antigen')
  mm <- match(ggdf$sample_id, md$sample_id)
  ggdf$timepoint <- md$timepoint[mm]
  print(ggplot(ggdf, aes(x = expression, color = timepoint, group = sample_id)) + 
          geom_density() +
          facet_wrap(~ antigen, nrow = 4, scales = 'free') + theme_bw() +
          theme(axis.text.x = element_text(angle = 90, hjust = 1), 
                strip.text = element_text(size = 7),
                axis.text = element_text(size = 5)) + 
          scale_color_manual(values = color_timepoint) )
  # plot 2
  
  ## Spot check - number of cells per sample
  cell_table <- table(sample_ids)
  ggdf <- data.frame(sample_id = names(cell_table), 
                     cell_counts = as.numeric(cell_table))
  mm <- match(ggdf$sample_id, md$sample_id)
  ggdf$timepoint <- md$timepoint[mm]
  print(ggplot(ggdf, aes(x = sample_id, y = cell_counts, fill = timepoint)) + 
          geom_bar(stat = 'identity') + 
          geom_text(aes(label = cell_counts), hjust = 0.5, vjust = -0.5, size = 2.5) + 
          theme_bw() +
          theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +  
          scale_fill_manual(values = color_timepoint, drop = FALSE) + 
          scale_x_discrete(drop = FALSE))
  
  # plot 3
  
  ## Define a function that calculates the Non-Redundancy Score per sample
  NRS <- function(x, ncomp = 3){
    pr <- prcomp(x, center = TRUE, scale. = FALSE)
    score <- rowSums(outer(rep(1, ncol(x)), pr$sdev[1:ncomp]^2) * 
                       abs(pr$rotation[,1:ncomp]))
    return(score)
  }
  
  ## Calculate the score
  ## May want to do the same with other markers
  nrs_sample <- fsApply(fcs[, subtype_markers], NRS, use.exprs = TRUE)
  rownames(nrs_sample) <- md$sample_id
  nrs <- colMeans(nrs_sample, na.rm = TRUE)
  
  ## Plot the NRS for ordered markers
  ## May be helpful to look at tissue instead of condition
  subtype_markers_ord <- names(sort(nrs, decreasing = TRUE))
  nrs_sample <- data.frame(nrs_sample)
  nrs_sample$sample_id <- rownames(nrs_sample)
  ggdf <- melt(nrs_sample, id.var = "sample_id",
               value.name = "nrs", variable.name = "antigen")
  ggdf$antigen <- factor(ggdf$antigen, levels = subtype_markers_ord)
  mm <- match(ggdf$sample_id, md$sample_id)
  ggdf$timepoint <- md$timepoint[mm]
  print(ggplot(ggdf, aes(x = antigen, y = nrs)) +
          geom_point(aes(color = timepoint), alpha = 0.9,
                     position = position_jitter(width = 0.3, height = 0)) +
          geom_boxplot(outlier.color = NA, fill = NA) +
          stat_summary(fun.y = "mean", geom = "point", shape = 21, fill = "white") +
          theme_bw() + ggtitle(tit)+ 
          theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) )
  # scale_color_manual(values = color_conditions)
  
  dev.off()
}


####CLUSTER HISTO####

plot_clustering_distr_wrapper <- function(expr = expr, 
                                          cell_clustering){
  # Calculate the median expression
  cell_clustering <- factor(cell_clustering)
  expr_median <- data.frame(expr, cell_clustering = cell_clustering) %>%
    group_by(cell_clustering) %>% summarize_all(funs(median))
  # Sort the cell clusters with hierarchical clustering
  d <- dist(expr_median[, colnames(expr)], method = "euclidean")
  cluster_rows <- hclust(d, method = "average")
  # Calculate cluster frequencies
  freq_clust <- table(cell_clustering)
  freq_clust <- round(as.numeric(freq_clust)/sum(freq_clust)*100, 2)
  cell_clustering <- factor(cell_clustering,
                            labels = levels(cell_clustering))
  ### Data organized per cluster
  ggd <- melt(data.frame(cluster = cell_clustering, expr),
              id.vars = "cluster", value.name = "expression",
              variable.name = "antigen")
  ggd$antigen <- factor(ggd$antigen, levels = colnames(expr))
  ggd$reference <- "no"
  ### The reference data
  ggd_bg <- ggd
  ggd_bg$cluster <- "reference"
  ggd_bg$reference <- "yes"
  
  ggd_plot <- rbind(ggd, ggd_bg)
  ggd_plot$cluster <- factor(ggd_plot$cluster,
                             levels = c(levels(cell_clustering)[rev(cluster_rows$order)], "reference"))
  
  ggplot() +
    geom_density_ridges(data = ggd_plot, aes(x = expression, y = cluster,
                                             color = reference, fill = reference), alpha = 0.3) +
    facet_wrap( ~ antigen, scales = "free_x", nrow = 2) +
    theme_ridges() +
    theme(axis.text = element_text(size = 7),  
          strip.text = element_text(size = 7), legend.position = "none")
  
}

####UMAP####
do_umap <- function(fcs,subtype_markers,sample_ids,cell_clustering,metadata,
                    clusterMergeFile='~_merging.xlsx',
                    seed = 1234, ncells=200,sample_subset=NULL){
  require(umap);require(flowCore);require(readxl)
  expr <- fsApply(fcs, exprs);expr <-expr[,subtype_markers]
  ## Create vector to later find and skip duplicates
  dups <- duplicated(expr[, subtype_markers])
  dups <- which(!(dups))## Find and skip duplicates
  ifelse(grepl(clusterMergeFile,pattern='.xls'),cluster_merging <- read_excel(clusterMergeFile),cluster_merging <- read.csv(clusterMergeFile,header = TRUE))
  ## New clustering1m
  mm <- match(cell_clustering, cluster_merging$original_cluster)
  cell_clustering1m <- cluster_merging$new_cluster[mm]
  ## Create a data frame of sample_ids and cell_clustering1m
  dtf<-data.frame(ids=sample_ids,type=cell_clustering1m)
  ## Data subsampling: create indices by sample
  inds <- split(1:length(sample_ids), sample_ids) #to get original indexes belonging to each cluster
  #inds2 <- split(1:length(sample_ids2), sample_ids2) #to fill in the original indexes that do not have B cells
  samplenames <- names(inds) #create a name vector of the files
  #FOR THIS TO WORK MUST BE IN FORMAT PRE/POST Tx
  # for (i in 1:(length(samplenames)/2)){#1:15 was here because ids in dtf was 30 factors and could be divided into 0 and 6wk for each so changed it
  #   templength <- length(inds2[[i]])
  #   inds2[[i]] <- inds[[i]][dtf$B[dtf$ids==samplenames[i]]] #subsets the "B cell or not a B cell" column for each sample by the name
  #   inds2[[i]] <- inds2[[i]][1:templength]
  # }
  custom.settings = umap.defaults
  custom.settings$seed = seed
  custom.settings$n.neighbors = neighbors
  ####umapindex generation####
  #umap ncells = table of sample ids with how many to downsample to by default col = id, row = ncells
  #sample ids = chr [1:2851129] "4927_0wk" "4927_0wk" "4927_0wk" "4927_0wk" ...
  #^ from ## Generate sample IDs corresponding to each cell in the 'expr' matrix sample_ids <- rep(md$sample_id, fsApply(fcs_raw, nrow))
  #can subset sample_ids and rerun umap 
  #can do timepoint or patient number or pre/post if you know corresp sample ids
  #sample_subset = '02_0wk' or c('02_0wk','2973_0wk') for example picking all 0 wk ones or use regex sample_ids[(grepl(sample_ids,pattern = '0wk'))]
  ifelse(is.null(sample_subset),
         umap_ncells <- pmin(table(sample_ids), ncells),
         umap_ncells <- pmin(table(sample_ids), ncells)[sample_subset]
  )
  if(!is.null(sample_subset)){inds <- inds[sample_subset]}
  umap_inds <- lapply(names(inds), function(i){
    s <- sample(inds[[i]], umap_ncells[i], replace = FALSE)
    intersect(s, dups)
  })
  set.seed(seed)
  umap_inds <- unlist(umap_inds)
  umap_out <- umap(expr[umap_inds, subtype_markers], config = custom.settings, method = 'naive')
  umapRes2D = data.frame(umap1 = umap_out$layout[, 1], umap2 = umap_out$layout[, 2], 
                         expr[umap_inds, subtype_markers],
                         sample_id = sample_ids[umap_inds], cell_clustering = factor(cell_clustering1m[umap_inds]), check.names = FALSE)
  
  #exclude any unassigned cluster post umap if needed--this has to be done by looking at the two columns 
  #metadata$sampleid is just a number in this metadatafile so to make unique ones combine samp_id col with timepoint (0wk)
  #to get in format of umapres2d$sample_id which looks like "02_0wk" do:
  return(umapRes2D)
}

plotUmap <- function(umapRes,seed=1234,neighbors=10,midpoint,color_clusters='auto',code_clustering,subtype_markers=NULL)
{require(umap);require(ggplot2);require(viridis);require(ggrepel)
  if((color_clusters)=='auto'){color_clusters <- hue_pal()(length(unique(code_clustering)))}
  custom.settings = umap.defaults
  custom.settings$seed = seed
  custom.settings$n.neighbors = neighbors
  ggp <- ggplot(umapRes,  aes(x = umap1, y = umap2, color = cell_clustering)) +
    geom_point(size = 1) +
    theme_bw() +
    
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()
    ) +
    
    scale_color_manual(values = color_clusters, name="CLUSTERS") +
    guides(color = guide_legend(override.aes = list(size = 3), ncol = 2))
  
  print(ggp)
  #other options
  print(ggp + facet_wrap(~ timepoint, ncol = 3)+ggtitle('TIMEPOINTS'))
  print(ggp + facet_wrap(~ batch, ncol = 2)+ggtitle('BATCH'))
  
  
  ggp2 <- ggplot(umapRes,  aes(x = umap1, y = umap2)) +
    geom_point(size = 1) +
    theme_bw() +
    
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()
    ) +
    guides(color = guide_legend(override.aes = list(size = 3), ncol = 2))
  print(ggp2 + facet_wrap(~ sample_id, ncol = 8)+ggtitle('SAMPLE'))
  
  
  
  
  ggp3 <- ggplot(umapRes,  aes(x = umap1, y = umap2, color = sample_id)) +
    geom_point(size = 1) +
    theme_bw() +
    
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()
    ) +
    guides(color = guide_legend(override.aes = list(size = 3), ncol = 2))
  
  print(ggp3)
  #can specify which markers to display
  if(!is.null(subtype_markers)){
    for(i in subtype_markers)
    {
      ggp <- ggplot(umapRes,  aes(x = umap1, y = umap2, color = umapRes[,i])) +
        geom_point(size = 1) +
        theme_bw() +
        theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank()) +
        scale_color_gradient2(i, low="dark blue",mid="white",high="dark red", midpoint = mean(unlist(umapRes[,i])))
      print(ggp)
    }
  }
}




#======================
#     RUNNING DATA
#======================


####DATA LOADING AND CLUSTERING####
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
workd<-getwd()

#read and cluster

output <- returnfcs(metaDataFile = paste0(workd,"/Config/","metadata_Aging.xlsx"),
                    panelDataFile = paste0(workd,"/Config/","Aging_panel.xlsx"),
                    dataDirectory = paste0(workd,'/Data'))

### Remove NA and INF Values from Normalized data
expr<-fsApply(output$fcs, exprs)

#check number of cells in expr:
nrow(expr)

#check if expr is free of NA and inf
sum(is.na(expr))
sum(is.infinite(expr))

#Cluster
output <- readRDS("backup_output.rds")
names(output)[8:10] <- c('code_clustering','cell_clustering','metaclusters')

clusterMergeFile = paste0(workd,"/Config/",'Aging_merged.xlsx')
cluster_merging <- read_xlsx(clusterMergeFile)

#set up factor levels

clusterlevels = c("TcEFF",
                  "TcEM",
                  "TcCM",
                  "TcN",
                  "ThCTL",
                  "Th2",
                  "Th2EM",
                  "Th2CM",
                  "Th17",
                  "Treg",
                  "ThN",
                  "DNT",
                  "NK",
                  "B",
                  "Myeloid",
                  "UA")

samplevels <- c("468_0",
                "468_1",
                "963_0",
                "963_1",
                "615_0",
                "615_1",
                "645_0",
                "645_1",
                "732_0",
                "732_1",
                "438_0",
                "438_1",
                "133_0",
                "133_1",
                "682_0",
                "682_1",
                "134_0",
                "134_1",
                "163_0",
                "163_1",
                "582_0",
                "582_1",
                "937_0",
                "937_1",
                "523_0",
                "523_1",
                "523_2",
                "360_0",
                "360_1",
                "239_0",
                "239_1",
                "404_0",
                "404_1",
                "404_2",
                "941_0",
                "941_1",
                "789_0",
                "789_1",
                "600_0",
                "600_1",
                "290_0",
                "290_1",
                "504_0",
                "504_1",
                "064_0",
                "064_1",
                "485_0",
                "485_1",
                "718_0",
                "718_1",
                "718_2",
                "258_0",
                "258_1",
                "890_0",
                "890_1",
                "890_2",
                "065_0",
                "065_1",
                "065_2",
                "617_0",
                "617_1",
                "617_2",
                "680_0",
                "680_1",
                "630_0",
                "630_1",
                "619_0",
                "619_1",
                "154_0",
                "154_1",
                "748_0",
                "748_1",
                "800_0",
                "800_1",
                "882_0",
                "882_1",
                "591_0",
                "591_1",
                "821_0",
                "821_1",
                "570_0",
                "570_1",
                "570_2",
                "108_0",
                "108_1",
                "188_0",
                "188_1",
                "409_0",
                "409_1",
                "409_2",
                "382_0",
                "382_1",
                "111_0",
                "111_1",
                "869_0",
                "869_1",
                "634_0",
                "634_1",
                "242_0",
                "242_1",
                "533_0",
                "533_1",
                "494_0",
                "494_1",
                "997_0",
                "997_1",
                "388_0",
                "388_1",
                "388_2",
                "636_0",
                "636_1",
                "636_2",
                "439_0",
                "439_1",
                "439_2",
                "446_0",
                "446_1",
                "450_0",
                "450_1",
                "527_0",
                "527_1",
                "758_0",
                "758_1",
                "758_2",
                "275_0",
                "275_1",
                "463_0",
                "463_1",
                "463_2",
                "701_0",
                "701_1",
                "959_0",
                "959_1",
                "412_0",
                "412_1",
                "020_0",
                "020_1",
                "453_0",
                "453_1",
                "961_0",
                "961_1",
                "009_0",
                "009_1",
                "543_0",
                "543_1",
                "543_2",
                "406_0",
                "406_1",
                "432_0",
                "432_1",
                "106_0",
                "106_1",
                "573_0",
                "573_1",
                "266_0",
                "266_1",
                "825_0",
                "825_1",
                "505_0",
                "505_1",
                "505_2",
                "750_0",
                "750_1",
                "530_0",
                "530_1",
                "605_0",
                "605_1",
                "875_0",
                "875_1",
                "589_0",
                "589_1",
                "541_0",
                "541_1",
                "488_0",
                "488_1",
                "387_0",
                "387_1",
                "464_0",
                "464_1",
                "464_2",
                "572_0",
                "572_1",
                "374_0",
                "374_1",
                "460_0",
                "460_1",
                "950_0",
                "950_1",
                "806_0",
                "806_1",
                "706_0",
                "706_1",
                "542_0",
                "542_1",
                "542_2",
                "499_0",
                "499_1",
                "499_2",
                "063_0",
                "063_1")

timelevels=c("0","1","2")

clustercolors <- as.character(c(cols25(n=25),alphabet(n=19)))

mm1 <- match(output$cell_clustering, cluster_merging$original_cluster)

cell_clustering1m <- cluster_merging$new_cluster[mm1]

output$cell_clustering1m <- cell_clustering1m

#metacluster heatmap
plot_clustering_heatmap_wrapper2(fcs=output$fcs,
                                 color_clusters = clustercolors,
                                 cell_clustering = factor(output$cell_clustering1m, levels=clusterlevels), 
                                 subtype_markers=output$subtype_markers,
                                 fileName = 'Extended_Aging_clusteringheatmap_final.pdf');dev.off()

####DIFFERENTIAL PLOTS####

#set up count and prop matrices
counts_table <- table(output$cell_clustering1m, output$sample_ids)
props_table <- t(t(counts_table) / colSums(counts_table)) * 100
counts <- as.data.frame.matrix(counts_table)
props <- as.data.frame.matrix(props_table)

write.csv(counts, file='Age_counts.csv')
write.csv(props, file='Age_props.csv')

#set up the data frame for proportional plotting
ggdf <- melt(data.frame(cluster = rownames(props),props, check.names = FALSE),
             id.vars = "cluster", value.name = "proportion", 
             variable.name = "sample_id")
ggdf$sample_id <- factor(ggdf$sample_id, levels=samplevels)
ggdf$cluster <- factor(ggdf$cluster, levels=clusterlevels)
ggdf$timepoint <- factor(output$meta_data$timepoint[match(ggdf$sample_id,output$meta_data$sample_id)], levels=timelevels)
ggdf$patient_id <- factor(output$meta_data$patient_id[match(ggdf$sample_id,output$meta_data$sample_id)])
ggdf$batch <- output$meta_data$Batch[match(ggdf$sample_id,output$meta_data$sample_id)]
ggdf$cancer_type <- factor(output$meta_data$cancer_type[match(ggdf$sample_id,output$meta_data$sample_id)])
ggdf$response <- factor(output$meta_data$response[match(ggdf$sample_id,output$meta_data$sample_id)])
ggdf$Aged <- factor(output$meta_data$Aged[match(ggdf$sample_id,output$meta_data$sample_id)])

#### Fold Change Dataframe #### 
ggdf2 <- ggdf[ggdf$sample_id %nin% c("523_2","404_2","718_2","890_2","065_2","617_2","570_2","409_2","388_2","636_2","439_2","758_2","463_2","543_2","505_2","464_2","542_2","499_2"),]
ggdf_T0 <- ggdf2[ggdf2$timepoint %nin% c("1","2"),]
ggdf_T1 <- ggdf2[ggdf2$timepoint %nin% c("0"),]
ggdf_T0$post_proportion = ggdf_T1$proportion
ggdf_T0$fc = ggdf_T0$post_proportion/ggdf_T0$proportion
ggdf_T0

#### Other Dataframes ####
ggdf_Baseline <- ggdf[ggdf$timepoint %nin% c("1", "2"),]
ggdf_OnTreatment <-ggdf[ggdf$timepoint %nin% c("0", "2"),]

ggdf_Baseline_response <- ggdf_Baseline[ggdf_Baseline$response %nin% c("NA"),]
ggdf_OnTreatment_response <- ggdf_OnTreatment[ggdf_OnTreatment$response %nin% c("NA"),]
ggdf_FC_response <- ggdf_T0[ggdf_T0$response %nin% c("NA"),]

ggdf_Young<- ggdf[ggdf$Aged %nin% c("Yes"),]
ggdf_Young<- ggdf_Young[ggdf_Young$response %nin% c("NA"),]
ggdf_Baseline_Young <- ggdf_Young[ggdf_Young$timepoint %nin% c("1", "2"),]
ggdf_OnTreatment_Young <- ggdf_Young[ggdf_Young$timepoint %nin% c("0", "2"),]
ggdf_FC_Young <-ggdf_T0[ggdf_T0$Aged %nin% c("Yes"),]
ggdf_FC_Young <-ggdf_FC_Young[ggdf_FC_Young$response %nin% c("NA"),]

ggdf_Aged<- ggdf[ggdf$Aged %nin% c("No"),]
ggdf_Aged<- ggdf_Aged[ggdf_Aged$response %nin% c("NA"),]
ggdf_Baseline_Aged <- ggdf_Aged[ggdf_Aged$timepoint %nin% c("1", "2"),]
ggdf_OnTreatment_Aged <- ggdf_Aged[ggdf_Aged$timepoint %nin% c("0", "2"),]
ggdf_FC_Aged <-ggdf_T0[ggdf_T0$Aged %nin% c("No"),]
ggdf_FC_Aged <-ggdf_FC_Aged[ggdf_FC_Aged$response %nin% c("NA"),]

#### Cancer Type Dataframes####
ggdf_Baseline_H <- ggdf_Baseline[ggdf_Baseline$cancer_type %in% c("HCC"),]
ggdf_Baseline_R <- ggdf_Baseline[ggdf_Baseline$cancer_type %in% c("RCC"),]
ggdf_Baseline_A <- ggdf_Baseline[ggdf_Baseline$cancer_type %nin% c("HCC", "RCC"),]
ggdf_OnTreatment_H <- ggdf_OnTreatment[ggdf_OnTreatment$cancer_type %in% c("HCC"),]
ggdf_OnTreatment_R <- ggdf_OnTreatment[ggdf_OnTreatment$cancer_type %in% c("RCC"),]
ggdf_OnTreatment_A <- ggdf_OnTreatment[ggdf_OnTreatment$cancer_type %nin% c("RCC", "HCC"),]
ggdf_T0_H <-ggdf_T0[ggdf_T0$cancer_type %in% c("HCC"),]
ggdf_T0_R <-ggdf_T0[ggdf_T0$cancer_type %in% c("RCC"),]
ggdf_T0_A <-ggdf_T0[ggdf_T0$cancer_type %nin% c("RCC", "HCC"),]
####Comparisons Lists for Stats####
irae_comps <- list(c("No", "Yes"))
age_comps <- list (c("No", "Yes"))
cancer_comps <- list (c("HCC", "RCC"))
response_comps <- list (c("no", "yes"))
####Response Box Plots ####
boxplot <- ggplot(ggdf_Baseline_response, aes(x=response, y=proportion, color=response))+
  geom_boxplot(outlier.size=0, lwd=0.75)+
  scale_color_manual(values=c("#4DBBD5FF", "#440154FF"))+
  geom_jitter(width=0.15)+
  facet_wrap(~cluster,ncol=10,scales="free")+
  ylab("% of Cells")+
  stat_compare_means(method="wilcox.test", 
                     show.legend=FALSE, 
                     method.args = list(var.equal = TRUE), 
                     label="p.value",
                     comparisons = response_comps,
                     label.x.npc=0.37, 
                     size=2.5, aes(group=Aged))+
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(size=8, color="black"),
        axis.title.x = element_blank(),
        axis.line.x = element_line(size=0.25, color="black"),
        axis.line.y = element_line(size=0.25, color="black"),
        axis.ticks.x = element_blank(),
        axis.title.y = element_text(size=8, color="black"),
        strip.background = element_rect(fill=NA),
        strip.text = element_text(size=8, color="black"),
        panel.background = element_rect(fill="white"),
        legend.title = element_blank(),
        legend.key.size = unit(3,'lines'),
        legend.text = element_text(size=12),
        legend.key = element_rect(fill="white"))
pdf("Response_Baseline.pdf", width=25, height=7);boxplot;dev.off()
write.csv(ggdf_Baseline_response, "Baseline_Response.csv")


boxplot <- ggplot(ggdf_OnTreatment_response, aes(x=response, y=proportion, color=response))+
  geom_boxplot(outlier.size=0, lwd=0.75)+
  scale_color_manual(values=c("#4DBBD5FF", "#440154FF"))+
  geom_jitter(width=0.15)+
  facet_wrap(~cluster,ncol=10,scales="free")+
  ylab("% of Cells")+
  stat_compare_means(method="wilcox.test", 
                     show.legend=FALSE, 
                     method.args = list(var.equal = TRUE), 
                     label="p.value",
                     comparisons = response_comps,
                     label.x.npc=0.37, 
                     size=2.5, aes(group=Aged))+
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(size=8, color="black"),
        axis.title.x = element_blank(),
        axis.line.x = element_line(size=0.25, color="black"),
        axis.line.y = element_line(size=0.25, color="black"),
        axis.ticks.x = element_blank(),
        axis.title.y = element_text(size=8, color="black"),
        strip.background = element_rect(fill=NA),
        strip.text = element_text(size=8, color="black"),
        panel.background = element_rect(fill="white"),
        legend.title = element_blank(),
        legend.key.size = unit(3,'lines'),
        legend.text = element_text(size=12),
        legend.key = element_rect(fill="white"))
pdf("Response_OnTreatment.pdf", width=25, height=7);boxplot;dev.off()
write.csv(ggdf_OnTreatment_response, "OnTreatment_Response.csv")

boxplot <- ggplot(ggdf_FC_response, aes(x=response, y=fc, color=response))+
  geom_boxplot(outlier.size=0, lwd=0.75)+
  scale_color_manual(values=c("#4DBBD5FF", "#440154FF"))+
  geom_jitter(width=0.15)+
  facet_wrap(~cluster,ncol=10,scales="free")+
  ylab("Fold Change")+
  stat_compare_means(method="wilcox.test", 
                     show.legend=FALSE, 
                     method.args = list(var.equal = TRUE), 
                     label="p.value",
                     comparisons = response_comps,
                     label.x.npc=0.37, 
                     size=2.5, aes(group=Aged))+
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(size=8, color="black"),
        axis.title.x = element_blank(),
        axis.line.x = element_line(size=0.25, color="black"),
        axis.line.y = element_line(size=0.25, color="black"),
        axis.ticks.x = element_blank(),
        axis.title.y = element_text(size=8, color="black"),
        strip.background = element_rect(fill=NA),
        strip.text = element_text(size=8, color="black"),
        panel.background = element_rect(fill="white"),
        legend.title = element_blank(),
        legend.key.size = unit(3,'lines'),
        legend.text = element_text(size=12),
        legend.key = element_rect(fill="white"))
pdf("Response_FC.pdf", width=25, height=7);boxplot;dev.off()
write.csv(ggdf_FC_response, "FC_Response.csv")

####Aged and Young Response Box Plots ####
boxplot <- ggplot(ggdf_Baseline_Young, aes(x=response, y=proportion, color=response))+
  geom_boxplot(outlier.size=0, lwd=0.75)+
  scale_color_manual(values=c("#4DBBD5FF", "#440154FF"))+
  geom_jitter(width=0.15)+
  facet_wrap(~cluster,ncol=10,scales="free")+
  ylab("% of Cells")+
  stat_compare_means(method="wilcox.test", 
                     show.legend=FALSE, 
                     method.args = list(var.equal = TRUE), 
                     label="p.value",
                     comparisons = response_comps,
                     label.x.npc=0.37, 
                     size=2.5, aes(group=Aged))+
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(size=8, color="black"),
        axis.title.x = element_blank(),
        axis.line.x = element_line(size=0.25, color="black"),
        axis.line.y = element_line(size=0.25, color="black"),
        axis.ticks.x = element_blank(),
        axis.title.y = element_text(size=8, color="black"),
        strip.background = element_rect(fill=NA),
        strip.text = element_text(size=8, color="black"),
        panel.background = element_rect(fill="white"),
        legend.title = element_blank(),
        legend.key.size = unit(3,'lines'),
        legend.text = element_text(size=12),
        legend.key = element_rect(fill="white"))
pdf("Young_Response_Baseline.pdf", width=25, height=7);boxplot;dev.off()
write.csv(ggdf_Baseline_Young, "Young_Baseline_Response.csv")

boxplot <- ggplot(ggdf_Baseline_Aged, aes(x=response, y=proportion, color=response))+
  geom_boxplot(outlier.size=0, lwd=0.75)+
  scale_color_manual(values=c("#4DBBD5FF", "#440154FF"))+
  geom_jitter(width=0.15)+
  facet_wrap(~cluster,ncol=10,scales="free")+
  ylab("% of Cells")+
  stat_compare_means(method="wilcox.test", 
                     show.legend=FALSE, 
                     method.args = list(var.equal = TRUE), 
                     label="p.value",
                     comparisons = response_comps,
                     label.x.npc=0.37, 
                     size=2.5, aes(group=Aged))+
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(size=8, color="black"),
        axis.title.x = element_blank(),
        axis.line.x = element_line(size=0.25, color="black"),
        axis.line.y = element_line(size=0.25, color="black"),
        axis.ticks.x = element_blank(),
        axis.title.y = element_text(size=8, color="black"),
        strip.background = element_rect(fill=NA),
        strip.text = element_text(size=8, color="black"),
        panel.background = element_rect(fill="white"),
        legend.title = element_blank(),
        legend.key.size = unit(3,'lines'),
        legend.text = element_text(size=12),
        legend.key = element_rect(fill="white"))
pdf("Aged_Response_Baseline.pdf", width=25, height=7);boxplot;dev.off()
write.csv(ggdf_Baseline_Aged, "Aged_Baseline_Response.csv")

boxplot <- ggplot(ggdf_OnTreatment_Young, aes(x=response, y=proportion, color=response))+
  geom_boxplot(outlier.size=0, lwd=0.75)+
  scale_color_manual(values=c("#4DBBD5FF", "#440154FF"))+
  geom_jitter(width=0.15)+
  facet_wrap(~cluster,ncol=10,scales="free")+
  ylab("% of Cells")+
  stat_compare_means(method="wilcox.test", 
                     show.legend=FALSE, 
                     method.args = list(var.equal = TRUE), 
                     label="p.value",
                     comparisons = response_comps,
                     label.x.npc=0.37, 
                     size=2.5, aes(group=Aged))+
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(size=8, color="black"),
        axis.title.x = element_blank(),
        axis.line.x = element_line(size=0.25, color="black"),
        axis.line.y = element_line(size=0.25, color="black"),
        axis.ticks.x = element_blank(),
        axis.title.y = element_text(size=8, color="black"),
        strip.background = element_rect(fill=NA),
        strip.text = element_text(size=8, color="black"),
        panel.background = element_rect(fill="white"),
        legend.title = element_blank(),
        legend.key.size = unit(3,'lines'),
        legend.text = element_text(size=12),
        legend.key = element_rect(fill="white"))
pdf("Young_Response_OnTreatment.pdf", width=25, height=7);boxplot;dev.off()
write.csv(ggdf_OnTreatment_Young, "Young_OnTreatment_Response.csv")

boxplot <- ggplot(ggdf_OnTreatment_Aged, aes(x=response, y=proportion, color=response))+
  geom_boxplot(outlier.size=0, lwd=0.75)+
  scale_color_manual(values=c("#4DBBD5FF", "#440154FF"))+
  geom_jitter(width=0.15)+
  facet_wrap(~cluster,ncol=10,scales="free")+
  ylab("% of Cells")+
  stat_compare_means(method="wilcox.test", 
                     show.legend=FALSE, 
                     method.args = list(var.equal = TRUE), 
                     label="p.value",
                     comparisons = response_comps,
                     label.x.npc=0.37, 
                     size=2.5, aes(group=Aged))+
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(size=8, color="black"),
        axis.title.x = element_blank(),
        axis.line.x = element_line(size=0.25, color="black"),
        axis.line.y = element_line(size=0.25, color="black"),
        axis.ticks.x = element_blank(),
        axis.title.y = element_text(size=8, color="black"),
        strip.background = element_rect(fill=NA),
        strip.text = element_text(size=8, color="black"),
        panel.background = element_rect(fill="white"),
        legend.title = element_blank(),
        legend.key.size = unit(3,'lines'),
        legend.text = element_text(size=12),
        legend.key = element_rect(fill="white"))
pdf("Aged_Response_OnTreatment.pdf", width=25, height=7);boxplot;dev.off()
write.csv(ggdf_OnTreatment_Aged, "Aged_OnTreatment_Response.csv")

boxplot <- ggplot(ggdf_FC_Young, aes(x=response, y=fc, color=response))+
  geom_boxplot(outlier.size=0, lwd=0.75)+
  scale_color_manual(values=c("#4DBBD5FF", "#440154FF"))+
  geom_jitter(width=0.15)+
  facet_wrap(~cluster,ncol=10,scales="free")+
  ylab("Fold Change")+
  stat_compare_means(method="p.value", 
                     show.legend=FALSE, 
                     method.args = list(var.equal = TRUE), 
                     label="p.format",
                     comparisons = response_comps,
                     label.x.npc=0.37, 
                     size=2.5, aes(group=Aged))+
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(size=8, color="black"),
        axis.title.x = element_blank(),
        axis.line.x = element_line(size=0.25, color="black"),
        axis.line.y = element_line(size=0.25, color="black"),
        axis.ticks.x = element_blank(),
        axis.title.y = element_text(size=8, color="black"),
        strip.background = element_rect(fill=NA),
        strip.text = element_text(size=8, color="black"),
        panel.background = element_rect(fill="white"),
        legend.title = element_blank(),
        legend.key.size = unit(3,'lines'),
        legend.text = element_text(size=12),
        legend.key = element_rect(fill="white"))
pdf("Young_Response_FC.pdf", width=25, height=7);boxplot;dev.off()
write.csv(ggdf_FC_Young, "Young_FC_Response.csv")

boxplot <- ggplot(ggdf_FC_Aged, aes(x=response, y=fc, color=response))+
  geom_boxplot(outlier.size=0, lwd=0.75)+
  scale_color_manual(values=c("#4DBBD5FF", "#440154FF"))+
  geom_jitter(width=0.15)+
  facet_wrap(~cluster,ncol=10,scales="free")+
  ylab("Fold Change")+
  stat_compare_means(method="wilcox.test", 
                     show.legend=FALSE, 
                     method.args = list(var.equal = TRUE), 
                     label="p.value",
                     comparisons = response_comps,
                     label.x.npc=0.37, 
                     size=2.5, aes(group=Aged))+
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(size=8, color="black"),
        axis.title.x = element_blank(),
        axis.line.x = element_line(size=0.25, color="black"),
        axis.line.y = element_line(size=0.25, color="black"),
        axis.ticks.x = element_blank(),
        axis.title.y = element_text(size=8, color="black"),
        strip.background = element_rect(fill=NA),
        strip.text = element_text(size=8, color="black"),
        panel.background = element_rect(fill="white"),
        legend.title = element_blank(),
        legend.key.size = unit(3,'lines'),
        legend.text = element_text(size=12),
        legend.key = element_rect(fill="white"))
pdf("Aged_Response_FC.pdf", width=25, height=7);boxplot;dev.off()
write.csv(ggdf_FC_Aged, "Aged_FC_Response.csv")



####Box Plots for abundance at baseline HCC RCC separate####
boxplot <- ggplot(ggdf_Baseline_H, aes(x=Aged, y=proportion, color=Aged))+
  geom_boxplot(outlier.size=0, lwd=0.75)+
  scale_color_manual(values=c("#4DBBD5FF", "#440154FF"))+
  geom_jitter(width=0.15)+
  facet_wrap(~cluster,ncol=10,scales="free")+
  ylab("% of Cells")+
  stat_compare_means(method="wilcox.test", 
                     show.legend=FALSE, 
                     method.args = list(var.equal = TRUE), 
                     label="p.value",
                     comparisons = age_comps,
                     label.x.npc=0.37, 
                     size=2.5, aes(group=Aged))+
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(size=8, color="black"),
        axis.title.x = element_blank(),
        axis.line.x = element_line(size=0.25, color="black"),
        axis.line.y = element_line(size=0.25, color="black"),
        axis.ticks.x = element_blank(),
        axis.title.y = element_text(size=8, color="black"),
        strip.background = element_rect(fill=NA),
        strip.text = element_text(size=8, color="black"),
        panel.background = element_rect(fill="white"),
        legend.title = element_blank(),
        legend.key.size = unit(3,'lines'),
        legend.text = element_text(size=12),
        legend.key = element_rect(fill="white"))
pdf("HCC_Baseline.pdf", width=25, height=7);boxplot;dev.off()
write.csv(ggdf_Baseline_H, "HCC_Baseline_Aged.csv")


boxplot <- ggplot(ggdf_Baseline_R, aes(x=Aged, y=proportion, color=Aged))+
  geom_boxplot(outlier.size=0, lwd=0.75)+
  scale_color_manual(values=c("#4DBBD5FF", "#440154FF"))+
  geom_jitter(width=0.15)+
  facet_wrap(~cluster,ncol=10,scales="free")+
  ylab("% of Cells")+
  stat_compare_means(method="wilcox.test", 
                     show.legend=FALSE, 
                     method.args = list(var.equal = TRUE), 
                     label="p.value",
                     comparisons = age_comps,
                     label.x.npc=0.37, 
                     size=2.5, aes(group=Aged))+
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(size=8, color="black"),
        axis.title.x = element_blank(),
        axis.line.x = element_line(size=0.25, color="black"),
        axis.line.y = element_line(size=0.25, color="black"),
        axis.ticks.x = element_blank(),
        axis.title.y = element_text(size=8, color="black"),
        strip.background = element_rect(fill=NA),
        strip.text = element_text(size=8, color="black"),
        panel.background = element_rect(fill="white"),
        legend.title = element_blank(),
        legend.key.size = unit(3,'lines'),
        legend.text = element_text(size=12),
        legend.key = element_rect(fill="white"))
pdf("RCC_Baseline.pdf", width=25, height=7);boxplot;dev.off()
write.csv(ggdf_Baseline_R, "RCC_Baseline_Aged.csv")


####Box Plot for abundance on treatment RCC HCC separate####
boxplot <- ggplot(ggdf_OnTreatment_H, aes(x=Aged, y=proportion, color=Aged))+
  geom_boxplot(outlier.size=0, lwd=0.75)+
  scale_color_manual(values=c("#4DBBD5FF", "#440154FF"))+
  geom_jitter(width=0.15)+
  facet_wrap(~cluster,ncol=10,scales="free")+
  ylab("% of Cells")+
  stat_compare_means(method="wilcox.test", 
                     show.legend=FALSE, 
                     method.args = list(var.equal = TRUE), 
                     label="p.value",
                     comparisons = age_comps,
                     label.x.npc=0.37, 
                     size=2.5, aes(group=Aged))+
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(size=8, color="black"),
        axis.title.x = element_blank(),
        axis.line.x = element_line(size=0.25, color="black"),
        axis.line.y = element_line(size=0.25, color="black"),
        axis.ticks.x = element_blank(),
        axis.title.y = element_text(size=8, color="black"),
        strip.background = element_rect(fill=NA),
        strip.text = element_text(size=8, color="black"),
        panel.background = element_rect(fill="white"),
        legend.title = element_blank(),
        legend.key.size = unit(3,'lines'),
        legend.text = element_text(size=12),
        legend.key = element_rect(fill="white"))
pdf("HCC_OnTreatment.pdf", width=25, height=7);boxplot;dev.off()
write.csv(ggdf_OnTreatment_H, "HCC_OnTreatment_Aged.csv")

boxplot <- ggplot(ggdf_OnTreatment_R, aes(x=Aged, y=proportion, color=Aged))+
  geom_boxplot(outlier.size=0, lwd=0.75)+
  scale_color_manual(values=c("#4DBBD5FF", "#440154FF"))+
  geom_jitter(width=0.15)+
  facet_wrap(~cluster,ncol=10,scales="free")+
  ylab("% of Cells")+
  stat_compare_means(method="wilcox.test", 
                     show.legend=FALSE, 
                     method.args = list(var.equal = TRUE), 
                     label="p.value",
                     comparisons = age_comps,
                     label.x.npc=0.37, 
                     size=2.5, aes(group=Aged))+
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(size=8, color="black"),
        axis.title.x = element_blank(),
        axis.line.x = element_line(size=0.25, color="black"),
        axis.line.y = element_line(size=0.25, color="black"),
        axis.ticks.x = element_blank(),
        axis.title.y = element_text(size=8, color="black"),
        strip.background = element_rect(fill=NA),
        strip.text = element_text(size=8, color="black"),
        panel.background = element_rect(fill="white"),
        legend.title = element_blank(),
        legend.key.size = unit(3,'lines'),
        legend.text = element_text(size=12),
        legend.key = element_rect(fill="white"))
pdf("RCC_OnTreatment.pdf", width=25, height=7);boxplot;dev.off()
write.csv(ggdf_OnTreatment_R, "RCC_OnTreatment_Aged.csv")


####Box Plot for Fold Change RCC HCC separate####
boxplot <- ggplot(ggdf_T0_H, aes(x=Aged, y=fc, color=Aged))+
  geom_boxplot(outlier.size=0, lwd=0.75)+
  scale_color_manual(values=c("#4DBBD5FF", "#440154FF"))+
  geom_jitter(width=0.15)+
  facet_wrap(~cluster,ncol=10,scales="free")+
  ylab("Fold CHange")+
  stat_compare_means(method="wilcox.test", 
                     show.legend=FALSE, 
                     method.args = list(var.equal = TRUE), 
                     label="p.value",
                     comparisons = age_comps,
                     label.x.npc=0.37, 
                     size=2.5, aes(group=Aged))+
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(size=8, color="black"),
        axis.title.x = element_blank(),
        axis.line.x = element_line(size=0.25, color="black"),
        axis.line.y = element_line(size=0.25, color="black"),
        axis.ticks.x = element_blank(),
        axis.title.y = element_text(size=8, color="black"),
        strip.background = element_rect(fill=NA),
        strip.text = element_text(size=8, color="black"),
        panel.background = element_rect(fill="white"),
        legend.title = element_blank(),
        legend.key.size = unit(3,'lines'),
        legend.text = element_text(size=12),
        legend.key = element_rect(fill="white"))
pdf("HCC_FC.pdf", width=25, height=11);boxplot;dev.off()
write.csv(ggdf_T0_H, "HCC_FC_Aged.csv")


boxplot <- ggplot(ggdf_T0_R, aes(x=Aged, y=fc, color=Aged))+
  geom_boxplot(outlier.size=0, lwd=0.75)+
  scale_color_manual(values=c("#4DBBD5FF", "#440154FF"))+
  geom_jitter(width=0.15)+
  facet_wrap(~cluster,ncol=10,scales="free")+
  ylab("Fold CHange")+
  stat_compare_means(method="wilcox.test", 
                     show.legend=FALSE, 
                     method.args = list(var.equal = TRUE), 
                     label="p.value",
                     comparisons = age_comps,
                     label.x.npc=0.37, 
                     size=2.5, aes(group=Aged))+
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(size=8, color="black"),
        axis.title.x = element_blank(),
        axis.line.x = element_line(size=0.25, color="black"),
        axis.line.y = element_line(size=0.25, color="black"),
        axis.ticks.x = element_blank(),
        axis.title.y = element_text(size=8, color="black"),
        strip.background = element_rect(fill=NA),
        strip.text = element_text(size=8, color="black"),
        panel.background = element_rect(fill="white"),
        legend.title = element_blank(),
        legend.key.size = unit(3,'lines'),
        legend.text = element_text(size=12),
        legend.key = element_rect(fill="white"))
pdf("RCC_FC.pdf", width=25, height=11);boxplot;dev.off()
write.csv(ggdf_T0_R, "RCC_FC_Aged.csv")


####Box Plot for All Other Cancer Types####
boxplot <- ggplot(ggdf_Baseline_A, aes(x=Aged, y=proportion, color=Aged))+
  geom_boxplot(outlier.size=0, lwd=0.75)+
  scale_color_manual(values=c("#4DBBD5FF", "#440154FF"))+
  geom_jitter(width=0.15)+
  facet_wrap(~cluster,ncol=10,scales="free")+
  ylab("% of Cells")+
  stat_compare_means(method="wilcox.test", 
                     show.legend=FALSE, 
                     method.args = list(var.equal = TRUE), 
                     label="p.value",
                     comparisons = age_comps,
                     label.x.npc=0.37, 
                     size=2.5, aes(group=Aged))+
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(size=8, color="black"),
        axis.title.x = element_blank(),
        axis.line.x = element_line(size=0.25, color="black"),
        axis.line.y = element_line(size=0.25, color="black"),
        axis.ticks.x = element_blank(),
        axis.title.y = element_text(size=8, color="black"),
        strip.background = element_rect(fill=NA),
        strip.text = element_text(size=8, color="black"),
        panel.background = element_rect(fill="white"),
        legend.title = element_blank(),
        legend.key.size = unit(3,'lines'),
        legend.text = element_text(size=12),
        legend.key = element_rect(fill="white"))
pdf("AllOther_Baseline.pdf", width=25, height=7);boxplot;dev.off()
write.csv(ggdf_Baseline_A, "AllOther_Baseline_Aged.csv")

boxplot <- ggplot(ggdf_OnTreatment_A, aes(x=Aged, y=proportion, color=Aged))+
  geom_boxplot(outlier.size=0, lwd=0.75)+
  scale_color_manual(values=c("#4DBBD5FF", "#440154FF"))+
  geom_jitter(width=0.15)+
  facet_wrap(~cluster,ncol=10,scales="free")+
  ylab("% of Cells")+
  stat_compare_means(method="wilcox.test", 
                     show.legend=FALSE, 
                     method.args = list(var.equal = TRUE), 
                     label="p.value",
                     comparisons = age_comps,
                     label.x.npc=0.37, 
                     size=2.5, aes(group=Aged))+
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(size=8, color="black"),
        axis.title.x = element_blank(),
        axis.line.x = element_line(size=0.25, color="black"),
        axis.line.y = element_line(size=0.25, color="black"),
        axis.ticks.x = element_blank(),
        axis.title.y = element_text(size=8, color="black"),
        strip.background = element_rect(fill=NA),
        strip.text = element_text(size=8, color="black"),
        panel.background = element_rect(fill="white"),
        legend.title = element_blank(),
        legend.key.size = unit(3,'lines'),
        legend.text = element_text(size=12),
        legend.key = element_rect(fill="white"))
pdf("AllOther_OnTreatment.pdf", width=25, height=7);boxplot;dev.off()
write.csv(ggdf_OnTreatment_A, "AllOther_OnTreatment_Aged.csv")

boxplot <- ggplot(ggdf_T0_A, aes(x=Aged, y=proportion, color=Aged))+
  geom_boxplot(outlier.size=0, lwd=0.75)+
  scale_color_manual(values=c("#4DBBD5FF", "#440154FF"))+
  geom_jitter(width=0.15)+
  facet_wrap(~cluster,ncol=10,scales="free")+
  ylab("FC")+
  stat_compare_means(method="wilcox.test", 
                     show.legend=FALSE, 
                     method.args = list(var.equal = TRUE), 
                     label="p.value",
                     comparisons = age_comps,
                     label.x.npc=0.37, 
                     size=2.5, aes(group=Aged))+
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(size=8, color="black"),
        axis.title.x = element_blank(),
        axis.line.x = element_line(size=0.25, color="black"),
        axis.line.y = element_line(size=0.25, color="black"),
        axis.ticks.x = element_blank(),
        axis.title.y = element_text(size=8, color="black"),
        strip.background = element_rect(fill=NA),
        strip.text = element_text(size=8, color="black"),
        panel.background = element_rect(fill="white"),
        legend.title = element_blank(),
        legend.key.size = unit(3,'lines'),
        legend.text = element_text(size=12),
        legend.key = element_rect(fill="white"))
pdf("AllOther_FC.pdf", width=25, height=7);boxplot;dev.off()
write.csv(ggdf_T0_A, "AllOther_FC_Aged.csv")

### Ratios ###
####Plot ratio of TcN over all Tc ####
#make a new table containing counts only for TcN & Tc clusters
counts_table_Tc <- counts_table[rownames(counts_table)[rownames(counts_table) %in% c("TcEFF",
                                                                                     "TcEM",
                                                                                     "TcCM",
                                                                                     "TcN")],]

counts_table_TcN <- counts_table[rownames(counts_table)[rownames(counts_table) %in% c("TcN")],]

#make new object that divides the sum of TcN columns / sum of Tc columns & convert to df
ratio_table <- t(t(counts_table_TcN) / colSums(counts_table_Tc)) 
ratio <- as.data.frame.matrix(ratio_table)
#ratio<-ratio %>% 
#filter_if(~is.numeric(.), all_vars(!is.infinite(.))) #remove INF numbers resulting from division by 0 (samples with no TAM2 cells)
#ratio<-na.omit(ratio) #remove NAs in data from division by 0

ggdf_ratio <- melt(data.frame(sample_id = rownames(ratio),ratio, check.names = FALSE),
                   id.vars = "sample_id", value.name = "ratio")
ggdf_ratio$sample_id <- factor(ggdf_ratio$sample_id, levels=samplevels)
ggdf_ratio$Aged <- factor(output$meta_data$Aged[match(ggdf_ratio$sample_id,output$meta_data$sample_id)])
ggdf_ratio$timepoint <- factor(output$meta_data$timepoint[match(ggdf_ratio$sample_id,output$meta_data$sample_id)])
ggdf_ratio_baseline <- ggdf_ratio[ggdf_ratio$timepoint %nin% c("1", "2"),]
ggdf_ratio_OnTreatment <- ggdf_ratio[ggdf_ratio$timepoint %nin% c("0", "2"),]
ggdf3 <- ggdf_ratio[ggdf_ratio$sample_id %nin% c("523_2","404_2","718_2","890_2","065_2","617_2","570_2","409_2","388_2","636_2","439_2","758_2","463_2","543_2","505_2","464_2","542_2","499_2"),]
ggdf_T0 <- ggdf3[ggdf3$timepoint %nin% c("1","2"),]
ggdf_T1 <- ggdf3[ggdf3$timepoint %nin% c("0"),]
ggdf_T0$post_ratio = ggdf_T1$ratio
ggdf_T0$fc = ggdf_T0$post_ratio/ggdf_T0$ratio
ggdf_T0

boxplot <- ggplot(ggdf_ratio_baseline, aes(x=Aged, y=ratio, color=Aged))+
  geom_boxplot(outlier.size=0, lwd=0.75)+
  scale_color_manual(values=c("#4DBBD5FF", "#440154FF"))+
  geom_jitter(width=0.15)+
  #facet_wrap(~cluster,ncol=10,scales="free")+
  ylab("TcN / Total Tc")+
  stat_compare_means(method="wilcox.test", 
                     show.legend=FALSE, 
                     method.args = list(var.equal = TRUE), 
                     label="p.value",
                     comparisons = age_comps,
                     label.x.npc=0.37, 
                     size=2.5, aes(group=Aged))+
  #scale_shape_manual(values=c(1:8,1:8))+
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(size=12, color="black"),
        axis.title.x = element_blank(),
        axis.line.x = element_line(size=0.25, color="black"),
        axis.line.y = element_line(size=0.25, color="black"),
        axis.ticks.x = element_blank(),
        axis.title.y = element_text(size=15, color="black"),
        strip.background = element_rect(fill=NA),
        strip.text = element_text(size=15, color="black"),
        panel.background = element_rect(fill="white"),
        legend.title = element_blank(),
        legend.key.size = unit(2,'lines'),
        legend.text = element_text(size=12),
        legend.key = element_rect(fill="white"))
pdf("TcN_Tc_Ratio_Baseline.pdf", width=4, height=4);boxplot;dev.off()
write.csv(ggdf_ratio_baseline, "TcN_Tc_Ratio_Baseline.csv")

boxplot <- ggplot(ggdf_ratio_OnTreatment, aes(x=Aged, y=ratio, color=Aged))+
  geom_boxplot(outlier.size=0, lwd=0.75)+
  scale_color_manual(values=c("#4DBBD5FF", "#440154FF"))+
  geom_jitter(width=0.15)+
  #facet_wrap(~cluster,ncol=10,scales="free")+
  ylab("TcN / Total Tc")+
  stat_compare_means(method="wilcox.test", 
                     show.legend=FALSE, 
                     method.args = list(var.equal = TRUE), 
                     label="p.value",
                     comparisons = age_comps,
                     label.x.npc=0.37, 
                     size=2.5, aes(group=Aged))+
  #scale_shape_manual(values=c(1:8,1:8))+
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(size=12, color="black"),
        axis.title.x = element_blank(),
        axis.line.x = element_line(size=0.25, color="black"),
        axis.line.y = element_line(size=0.25, color="black"),
        axis.ticks.x = element_blank(),
        axis.title.y = element_text(size=15, color="black"),
        strip.background = element_rect(fill=NA),
        strip.text = element_text(size=15, color="black"),
        panel.background = element_rect(fill="white"),
        legend.title = element_blank(),
        legend.key.size = unit(2,'lines'),
        legend.text = element_text(size=12),
        legend.key = element_rect(fill="white"))
pdf("TcN_Tc_Ratio_OnTreatment.pdf", width=4, height=4);boxplot;dev.off()
write.csv(ggdf_ratio_OnTreatment, "TcN_Tc_Ratio_OnTreatment.csv")


boxplot <- ggplot(ggdf_T0, aes(x=Aged, y=fc, color=Aged))+
  geom_boxplot(outlier.size=0, lwd=0.75)+
  scale_color_manual(values=c("#4DBBD5FF", "#440154FF"))+
  geom_jitter(width=0.15)+
  #facet_wrap(~cluster,ncol=10,scales="free")+
  ylab("Fold Change of TcN / Total Tc")+
  stat_compare_means(method="wilcox.test", 
                     show.legend=FALSE, 
                     method.args = list(var.equal = TRUE), 
                     label="p.value",
                     comparisons = age_comps,
                     label.x.npc=0.37, 
                     size=2.5, aes(group=Aged))+
  #scale_shape_manual(values=c(1:8,1:8))+
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(size=12, color="black"),
        axis.title.x = element_blank(),
        axis.line.x = element_line(size=0.25, color="black"),
        axis.line.y = element_line(size=0.25, color="black"),
        axis.ticks.x = element_blank(),
        axis.title.y = element_text(size=15, color="black"),
        strip.background = element_rect(fill=NA),
        strip.text = element_text(size=15, color="black"),
        panel.background = element_rect(fill="white"),
        legend.title = element_blank(),
        legend.key.size = unit(2,'lines'),
        legend.text = element_text(size=12),
        legend.key = element_rect(fill="white"))
pdf("TcN_Tc_Ratio_FC.pdf", width=4, height=4);boxplot;dev.off()
write.csv(ggdf_T0, "TcN_Tc_Ratio_FC.csv")


####Plot ratio of TcEFF over all Tc ####
#make a new table containing counts only for TcN & Tc clusters
counts_table_Tc <- counts_table[rownames(counts_table)[rownames(counts_table) %in% c("TcEFF",
                                                                                     "TcEM",
                                                                                     "TcCM",
                                                                                     "TcN")],]

counts_table_TcEFF <- counts_table[rownames(counts_table)[rownames(counts_table) %in% c("TcEFF")],]

#make new object that divides the sum of TcN columns / sum of Tc columns & convert to df
ratio_table <- t(t(counts_table_TcEFF) / colSums(counts_table_Tc)) 
ratio <- as.data.frame.matrix(ratio_table)
#ratio<-ratio %>% 
#filter_if(~is.numeric(.), all_vars(!is.infinite(.))) #remove INF numbers resulting from division by 0 (samples with no TAM2 cells)
#ratio<-na.omit(ratio) #remove NAs in data from division by 0

ggdf_ratio <- melt(data.frame(sample_id = rownames(ratio),ratio, check.names = FALSE),
                   id.vars = "sample_id", value.name = "ratio")
ggdf_ratio$sample_id <- factor(ggdf_ratio$sample_id, levels=samplevels)
ggdf_ratio$Aged <- factor(output$meta_data$Aged[match(ggdf_ratio$sample_id,output$meta_data$sample_id)])
ggdf_ratio$timepoint <- factor(output$meta_data$timepoint[match(ggdf_ratio$sample_id,output$meta_data$sample_id)])
ggdf_ratio_baseline <- ggdf_ratio[ggdf_ratio$timepoint %nin% c("1", "2"),]
ggdf_ratio_OnTreatment <- ggdf_ratio[ggdf_ratio$timepoint %nin% c("0", "2"),]
ggdf3 <- ggdf_ratio[ggdf_ratio$sample_id %nin% c("523_2","404_2","718_2","890_2","065_2","617_2","570_2","409_2","388_2","636_2","439_2","758_2","463_2","543_2","505_2","464_2","542_2","499_2"),]
ggdf_T0 <- ggdf3[ggdf3$timepoint %nin% c("1","2"),]
ggdf_T1 <- ggdf3[ggdf3$timepoint %nin% c("0"),]
ggdf_T0$post_ratio = ggdf_T1$ratio
ggdf_T0$fc = ggdf_T0$post_ratio/ggdf_T0$ratio

boxplot <- ggplot(ggdf_ratio_baseline, aes(x=Aged, y=ratio, color=Aged))+
  geom_boxplot(outlier.size=0, lwd=0.75)+
  scale_color_manual(values=c("#4DBBD5FF", "#440154FF"))+
  geom_jitter(width=0.15)+
  #facet_wrap(~cluster,ncol=10,scales="free")+
  ylab("TcEFF / Total Tc")+
  stat_compare_means(method="wilcox.test", 
                     show.legend=FALSE, 
                     method.args = list(var.equal = TRUE), 
                     label="p.value",
                     comparisons = age_comps,
                     label.x.npc=0.37, 
                     size=2.5, aes(group=Aged))+
  #scale_shape_manual(values=c(1:8,1:8))+
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(size=12, color="black"),
        axis.title.x = element_blank(),
        axis.line.x = element_line(size=0.25, color="black"),
        axis.line.y = element_line(size=0.25, color="black"),
        axis.ticks.x = element_blank(),
        axis.title.y = element_text(size=15, color="black"),
        strip.background = element_rect(fill=NA),
        strip.text = element_text(size=15, color="black"),
        panel.background = element_rect(fill="white"),
        legend.title = element_blank(),
        legend.key.size = unit(2,'lines'),
        legend.text = element_text(size=12),
        legend.key = element_rect(fill="white"))
pdf("TcEFF_Tc_Ratio_Baseline.pdf", width=4, height=4);boxplot;dev.off()
write.csv(ggdf_ratio_baseline, "TcEFF_Tc_Ratio_Baseline.csv")

boxplot <- ggplot(ggdf_ratio_OnTreatment, aes(x=Aged, y=ratio, color=Aged))+
  geom_boxplot(outlier.size=0, lwd=0.75)+
  scale_color_manual(values=c("#4DBBD5FF", "#440154FF"))+
  geom_jitter(width=0.15)+
  #facet_wrap(~cluster,ncol=10,scales="free")+
  ylab("TcEFF / Total Tc")+
  stat_compare_means(method="wilcox.test", 
                     show.legend=FALSE, 
                     method.args = list(var.equal = TRUE), 
                     label="p.value",
                     comparisons = age_comps,
                     label.x.npc=0.37, 
                     size=2.5, aes(group=Aged))+
  #scale_shape_manual(values=c(1:8,1:8))+
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(size=12, color="black"),
        axis.title.x = element_blank(),
        axis.line.x = element_line(size=0.25, color="black"),
        axis.line.y = element_line(size=0.25, color="black"),
        axis.ticks.x = element_blank(),
        axis.title.y = element_text(size=15, color="black"),
        strip.background = element_rect(fill=NA),
        strip.text = element_text(size=15, color="black"),
        panel.background = element_rect(fill="white"),
        legend.title = element_blank(),
        legend.key.size = unit(2,'lines'),
        legend.text = element_text(size=12),
        legend.key = element_rect(fill="white"))
pdf("TcEFF_Tc_Ratio_OnTreatment.pdf", width=4, height=4);boxplot;dev.off()
write.csv(ggdf_ratio_OnTreatment, "TcEFF_Tc_Ratio_OnTreatment.csv")

boxplot <- ggplot(ggdf_T0, aes(x=Aged, y=fc, color=Aged))+
  geom_boxplot(outlier.size=0, lwd=0.75)+
  scale_color_manual(values=c("#4DBBD5FF", "#440154FF"))+
  geom_jitter(width=0.15)+
  #facet_wrap(~cluster,ncol=10,scales="free")+
  ylab("Fold Change of TcEFF / Total Tc")+
  stat_compare_means(method="wilcox.test", 
                     show.legend=FALSE, 
                     method.args = list(var.equal = TRUE), 
                     label="p.value",
                     comparisons = age_comps,
                     label.x.npc=0.37, 
                     size=2.5, aes(group=Aged))+
  #scale_shape_manual(values=c(1:8,1:8))+
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(size=12, color="black"),
        axis.title.x = element_blank(),
        axis.line.x = element_line(size=0.25, color="black"),
        axis.line.y = element_line(size=0.25, color="black"),
        axis.ticks.x = element_blank(),
        axis.title.y = element_text(size=15, color="black"),
        strip.background = element_rect(fill=NA),
        strip.text = element_text(size=15, color="black"),
        panel.background = element_rect(fill="white"),
        legend.title = element_blank(),
        legend.key.size = unit(2,'lines'),
        legend.text = element_text(size=12),
        legend.key = element_rect(fill="white"))
pdf("TcEFF_Tc_Ratio_FC.pdf", width=4, height=4);boxplot;dev.off()
write.csv(ggdf_T0, "TcEFF_Tc_Ratio_FC.csv")

####Plot ratio of TcEM over all Tc ####
#make a new table containing counts only for TcN & Tc clusters
counts_table_Tc <- counts_table[rownames(counts_table)[rownames(counts_table) %in% c("TcEFF",
                                                                                     "TcEM",
                                                                                     "TcCM",
                                                                                     "TcN")],]

counts_table_TcEM <- counts_table[rownames(counts_table)[rownames(counts_table) %in% c("TcEM")],]

#make new object that divides the sum of TcN columns / sum of Tc columns & convert to df
ratio_table <- t(t(counts_table_TcEM) / colSums(counts_table_Tc)) 
ratio <- as.data.frame.matrix(ratio_table)
#ratio<-ratio %>% 
#filter_if(~is.numeric(.), all_vars(!is.infinite(.))) #remove INF numbers resulting from division by 0 (samples with no TAM2 cells)
#ratio<-na.omit(ratio) #remove NAs in data from division by 0

ggdf_ratio <- melt(data.frame(sample_id = rownames(ratio),ratio, check.names = FALSE),
                   id.vars = "sample_id", value.name = "ratio")
ggdf_ratio$sample_id <- factor(ggdf_ratio$sample_id, levels=samplevels)
ggdf_ratio$Aged <- factor(output$meta_data$Aged[match(ggdf_ratio$sample_id,output$meta_data$sample_id)])
ggdf_ratio$timepoint <- factor(output$meta_data$timepoint[match(ggdf_ratio$sample_id,output$meta_data$sample_id)])
ggdf_ratio_baseline <- ggdf_ratio[ggdf_ratio$timepoint %nin% c("1", "2"),]
ggdf_ratio_OnTreatment <- ggdf_ratio[ggdf_ratio$timepoint %nin% c("0", "2"),]
ggdf3 <- ggdf_ratio[ggdf_ratio$sample_id %nin% c("523_2","404_2","718_2","890_2","065_2","617_2","570_2","409_2","388_2","636_2","439_2","758_2","463_2","543_2","505_2","464_2","542_2","499_2"),]
ggdf_T0 <- ggdf3[ggdf3$timepoint %nin% c("1","2"),]
ggdf_T1 <- ggdf3[ggdf3$timepoint %nin% c("0"),]
ggdf_T0$post_ratio = ggdf_T1$ratio
ggdf_T0$fc = ggdf_T0$post_ratio/ggdf_T0$ratio


boxplot <- ggplot(ggdf_ratio_baseline, aes(x=Aged, y=ratio, color=Aged))+
  geom_boxplot(outlier.size=0, lwd=0.75)+
  scale_color_manual(values=c("#4DBBD5FF", "#440154FF"))+
  geom_jitter(width=0.15)+
  #facet_wrap(~cluster,ncol=10,scales="free")+
  ylab("TcEM / Total Tc")+
  stat_compare_means(method="wilcox.test", 
                     show.legend=FALSE, 
                     method.args = list(var.equal = TRUE), 
                     label="p.value",
                     comparisons = age_comps,
                     label.x.npc=0.37, 
                     size=2.5, aes(group=Aged))+
  #scale_shape_manual(values=c(1:8,1:8))+
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(size=12, color="black"),
        axis.title.x = element_blank(),
        axis.line.x = element_line(size=0.25, color="black"),
        axis.line.y = element_line(size=0.25, color="black"),
        axis.ticks.x = element_blank(),
        axis.title.y = element_text(size=15, color="black"),
        strip.background = element_rect(fill=NA),
        strip.text = element_text(size=15, color="black"),
        panel.background = element_rect(fill="white"),
        legend.title = element_blank(),
        legend.key.size = unit(2,'lines'),
        legend.text = element_text(size=12),
        legend.key = element_rect(fill="white"))
pdf("TcEM_Tc_Ratio_Baseline.pdf", width=4, height=4);boxplot;dev.off()
write.csv(ggdf_ratio_baseline, "TcEM_Tc_Ratio_Baseline.csv")

boxplot <- ggplot(ggdf_ratio_OnTreatment, aes(x=Aged, y=ratio, color=Aged))+
  geom_boxplot(outlier.size=0, lwd=0.75)+
  scale_color_manual(values=c("#4DBBD5FF", "#440154FF"))+
  geom_jitter(width=0.15)+
  #facet_wrap(~cluster,ncol=10,scales="free")+
  ylab("TcEM / Total Tc")+
  stat_compare_means(method="wilcox.test", 
                     show.legend=FALSE, 
                     method.args = list(var.equal = TRUE), 
                     label="p.value",
                     comparisons = age_comps,
                     label.x.npc=0.37, 
                     size=2.5, aes(group=Aged))+
  #scale_shape_manual(values=c(1:8,1:8))+
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(size=12, color="black"),
        axis.title.x = element_blank(),
        axis.line.x = element_line(size=0.25, color="black"),
        axis.line.y = element_line(size=0.25, color="black"),
        axis.ticks.x = element_blank(),
        axis.title.y = element_text(size=15, color="black"),
        strip.background = element_rect(fill=NA),
        strip.text = element_text(size=15, color="black"),
        panel.background = element_rect(fill="white"),
        legend.title = element_blank(),
        legend.key.size = unit(2,'lines'),
        legend.text = element_text(size=12),
        legend.key = element_rect(fill="white"))
pdf("TcEM_Tc_Ratio_OnTreatment.pdf", width=4, height=4);boxplot;dev.off()
write.csv(ggdf_ratio_OnTreatment, "TcEM_Tc_Ratio_OnTreatment.csv")

boxplot <- ggplot(ggdf_T0, aes(x=Aged, y=fc, color=Aged))+
  geom_boxplot(outlier.size=0, lwd=0.75)+
  scale_color_manual(values=c("#4DBBD5FF", "#440154FF"))+
  geom_jitter(width=0.15)+
  #facet_wrap(~cluster,ncol=10,scales="free")+
  ylab("Fold Change of TcEM / Total Tc")+
  stat_compare_means(method="wilcox.test", 
                     show.legend=FALSE, 
                     method.args = list(var.equal = TRUE), 
                     label="p.value",
                     comparisons = age_comps,
                     label.x.npc=0.37, 
                     size=2.5, aes(group=Aged))+
  #scale_shape_manual(values=c(1:8,1:8))+
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(size=12, color="black"),
        axis.title.x = element_blank(),
        axis.line.x = element_line(size=0.25, color="black"),
        axis.line.y = element_line(size=0.25, color="black"),
        axis.ticks.x = element_blank(),
        axis.title.y = element_text(size=15, color="black"),
        strip.background = element_rect(fill=NA),
        strip.text = element_text(size=15, color="black"),
        panel.background = element_rect(fill="white"),
        legend.title = element_blank(),
        legend.key.size = unit(2,'lines'),
        legend.text = element_text(size=12),
        legend.key = element_rect(fill="white"))
pdf("TcEM_Tc_Ratio_FC.pdf", width=4, height=4);boxplot;dev.off()
write.csv(ggdf_T0, "TcEM_Tc_Ratio_FC.csv")

####Plot ratio of TcCM over all Tc ####
#make a new table containing counts only for TcN & Tc clusters
counts_table_Tc <- counts_table[rownames(counts_table)[rownames(counts_table) %in% c("TcEFF",
                                                                                     "TcEM",
                                                                                     "TcCM",
                                                                                     "TcN")],]

counts_table_TcCM <- counts_table[rownames(counts_table)[rownames(counts_table) %in% c("TcCM")],]

#make new object that divides the sum of TcN columns / sum of Tc columns & convert to df
ratio_table <- t(t(counts_table_TcCM) / colSums(counts_table_Tc)) 
ratio <- as.data.frame.matrix(ratio_table)
#ratio<-ratio %>% 
#filter_if(~is.numeric(.), all_vars(!is.infinite(.))) #remove INF numbers resulting from division by 0 (samples with no TAM2 cells)
#ratio<-na.omit(ratio) #remove NAs in data from division by 0

ggdf_ratio <- melt(data.frame(sample_id = rownames(ratio),ratio, check.names = FALSE),
                   id.vars = "sample_id", value.name = "ratio")
ggdf_ratio$sample_id <- factor(ggdf_ratio$sample_id, levels=samplevels)
ggdf_ratio$Aged <- factor(output$meta_data$Aged[match(ggdf_ratio$sample_id,output$meta_data$sample_id)])
ggdf_ratio$timepoint <- factor(output$meta_data$timepoint[match(ggdf_ratio$sample_id,output$meta_data$sample_id)])
ggdf_ratio_baseline <- ggdf_ratio[ggdf_ratio$timepoint %nin% c("1", "2"),]
ggdf_ratio_OnTreatment <- ggdf_ratio[ggdf_ratio$timepoint %nin% c("0", "2"),]
ggdf3 <- ggdf_ratio[ggdf_ratio$sample_id %nin% c("523_2","404_2","718_2","890_2","065_2","617_2","570_2","409_2","388_2","636_2","439_2","758_2","463_2","543_2","505_2","464_2","542_2","499_2"),]
ggdf_T0 <- ggdf3[ggdf3$timepoint %nin% c("1","2"),]
ggdf_T1 <- ggdf3[ggdf3$timepoint %nin% c("0"),]
ggdf_T0$post_ratio = ggdf_T1$ratio
ggdf_T0$fc = ggdf_T0$post_ratio/ggdf_T0$ratio

boxplot <- ggplot(ggdf_ratio_baseline, aes(x=Aged, y=ratio, color=Aged))+
  geom_boxplot(outlier.size=0, lwd=0.75)+
  scale_color_manual(values=c("#4DBBD5FF", "#440154FF"))+
  geom_jitter(width=0.15)+
  #facet_wrap(~cluster,ncol=10,scales="free")+
  ylab("TcCM / Total Tc")+
  stat_compare_means(method="wilcox.test", 
                     show.legend=FALSE, 
                     method.args = list(var.equal = TRUE), 
                     label="p.value",
                     comparisons = age_comps,
                     label.x.npc=0.37, 
                     size=2.5, aes(group=Aged))+
  #scale_shape_manual(values=c(1:8,1:8))+
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(size=12, color="black"),
        axis.title.x = element_blank(),
        axis.line.x = element_line(size=0.25, color="black"),
        axis.line.y = element_line(size=0.25, color="black"),
        axis.ticks.x = element_blank(),
        axis.title.y = element_text(size=15, color="black"),
        strip.background = element_rect(fill=NA),
        strip.text = element_text(size=15, color="black"),
        panel.background = element_rect(fill="white"),
        legend.title = element_blank(),
        legend.key.size = unit(2,'lines'),
        legend.text = element_text(size=12),
        legend.key = element_rect(fill="white"))
pdf("TcCM_Tc_Ratio_Baseline.pdf", width=4, height=4);boxplot;dev.off()
write.csv(ggdf_ratio_baseline, "TcCM_Tc_Ratio_Baseline.csv")

boxplot <- ggplot(ggdf_ratio_OnTreatment, aes(x=Aged, y=ratio, color=Aged))+
  geom_boxplot(outlier.size=0, lwd=0.75)+
  scale_color_manual(values=c("#4DBBD5FF", "#440154FF"))+
  geom_jitter(width=0.15)+
  #facet_wrap(~cluster,ncol=10,scales="free")+
  ylab("TcCM / Total Tc")+
  stat_compare_means(method="wilcox.test", 
                     show.legend=FALSE, 
                     method.args = list(var.equal = TRUE), 
                     label="p.value",
                     comparisons = age_comps,
                     label.x.npc=0.37, 
                     size=2.5, aes(group=Aged))+
  #scale_shape_manual(values=c(1:8,1:8))+
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(size=12, color="black"),
        axis.title.x = element_blank(),
        axis.line.x = element_line(size=0.25, color="black"),
        axis.line.y = element_line(size=0.25, color="black"),
        axis.ticks.x = element_blank(),
        axis.title.y = element_text(size=15, color="black"),
        strip.background = element_rect(fill=NA),
        strip.text = element_text(size=15, color="black"),
        panel.background = element_rect(fill="white"),
        legend.title = element_blank(),
        legend.key.size = unit(2,'lines'),
        legend.text = element_text(size=12),
        legend.key = element_rect(fill="white"))
pdf("TcCM_Tc_Ratio_OnTreatment.pdf", width=4, height=4);boxplot;dev.off()
write.csv(ggdf_ratio_OnTreatment, "TcCM_Tc_Ratio_OnTreatment.csv")

boxplot <- ggplot(ggdf_T0, aes(x=Aged, y=fc, color=Aged))+
  geom_boxplot(outlier.size=0, lwd=0.75)+
  scale_color_manual(values=c("#4DBBD5FF", "#440154FF"))+
  geom_jitter(width=0.15)+
  #facet_wrap(~cluster,ncol=10,scales="free")+
  ylab("Fold Change of TcCM / Total Tc")+
  stat_compare_means(method="wilcox.test", 
                     show.legend=FALSE, 
                     method.args = list(var.equal = TRUE), 
                     label="p.value",
                     comparisons = age_comps,
                     label.x.npc=0.37, 
                     size=2.5, aes(group=Aged))+
  #scale_shape_manual(values=c(1:8,1:8))+
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(size=12, color="black"),
        axis.title.x = element_blank(),
        axis.line.x = element_line(size=0.25, color="black"),
        axis.line.y = element_line(size=0.25, color="black"),
        axis.ticks.x = element_blank(),
        axis.title.y = element_text(size=15, color="black"),
        strip.background = element_rect(fill=NA),
        strip.text = element_text(size=15, color="black"),
        panel.background = element_rect(fill="white"),
        legend.title = element_blank(),
        legend.key.size = unit(2,'lines'),
        legend.text = element_text(size=12),
        legend.key = element_rect(fill="white"))
pdf("TcCM_Tc_Ratio_FC.pdf", width=4, height=4);boxplot;dev.off()
write.csv(ggdf_T0, "TcCM_Tc_Ratio_FC.csv")

####Plot ratio of ThN over all Th ####
#make a new table containing counts only for ThN & Th clusters
counts_table_Th <- counts_table[rownames(counts_table)[rownames(counts_table) %in% c("ThCTL",
                                                                                     "Th2",
                                                                                     "Th2EM",
                                                                                     "Th2CM",
                                                                                     "Th17",
                                                                                     "Treg",
                                                                                     "ThN")],]

counts_table_ThN <- counts_table[rownames(counts_table)[rownames(counts_table) %in% c("ThN")],]

#make new object that divides the sum of TcN columns / sum of Tc columns & convert to df
ratio_table <- t(t(counts_table_ThN) / colSums(counts_table_Th)) 
ratio <- as.data.frame.matrix(ratio_table)
#ratio<-ratio %>% 
#filter_if(~is.numeric(.), all_vars(!is.infinite(.))) #remove INF numbers resulting from division by 0 (samples with no TAM2 cells)
#ratio<-na.omit(ratio) #remove NAs in data from division by 0

ggdf_ratio <- melt(data.frame(sample_id = rownames(ratio),ratio, check.names = FALSE),
                   id.vars = "sample_id", value.name = "ratio")
ggdf_ratio$sample_id <- factor(ggdf_ratio$sample_id, levels=samplevels)
ggdf_ratio$Aged <- factor(output$meta_data$Aged[match(ggdf_ratio$sample_id,output$meta_data$sample_id)])
ggdf_ratio$timepoint <- factor(output$meta_data$timepoint[match(ggdf_ratio$sample_id,output$meta_data$sample_id)])
ggdf_ratio_baseline <- ggdf_ratio[ggdf_ratio$timepoint %nin% c("1", "2"),]
ggdf_ratio_OnTreatment <- ggdf_ratio[ggdf_ratio$timepoint %nin% c("0", "2"),]
ggdf3 <- ggdf_ratio[ggdf_ratio$sample_id %nin% c("523_2","404_2","718_2","890_2","065_2","617_2","570_2","409_2","388_2","636_2","439_2","758_2","463_2","543_2","505_2","464_2","542_2","499_2"),]
ggdf_T0 <- ggdf3[ggdf3$timepoint %nin% c("1","2"),]
ggdf_T1 <- ggdf3[ggdf3$timepoint %nin% c("0"),]
ggdf_T0$post_ratio = ggdf_T1$ratio
ggdf_T0$fc = ggdf_T0$post_ratio/ggdf_T0$ratio


boxplot <- ggplot(ggdf_ratio_baseline, aes(x=Aged, y=ratio, color=Aged))+
  geom_boxplot(outlier.size=0, lwd=0.75)+
  scale_color_manual(values=c("#4DBBD5FF", "#440154FF"))+
  geom_jitter(width=0.15)+
  #facet_wrap(~cluster,ncol=10,scales="free")+
  ylab("ThN / Total Th")+
  stat_compare_means(method="wilcox.test", 
                     show.legend=FALSE, 
                     method.args = list(var.equal = TRUE), 
                     label="p.value",
                     comparisons = age_comps,
                     label.x.npc=0.37, 
                     size=2.5, aes(group=Aged))+
  #scale_shape_manual(values=c(1:8,1:8))+
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(size=12, color="black"),
        axis.title.x = element_blank(),
        axis.line.x = element_line(size=0.25, color="black"),
        axis.line.y = element_line(size=0.25, color="black"),
        axis.ticks.x = element_blank(),
        axis.title.y = element_text(size=15, color="black"),
        strip.background = element_rect(fill=NA),
        strip.text = element_text(size=15, color="black"),
        panel.background = element_rect(fill="white"),
        legend.title = element_blank(),
        legend.key.size = unit(2,'lines'),
        legend.text = element_text(size=12),
        legend.key = element_rect(fill="white"))
pdf("ThN_Th_Ratio_Baseline.pdf", width=4, height=4);boxplot;dev.off()
write.csv(ggdf_ratio_baseline, "ThN_Th_Ratio_Baseline.csv")

boxplot <- ggplot(ggdf_ratio_OnTreatment, aes(x=Aged, y=ratio, color=Aged))+
  geom_boxplot(outlier.size=0, lwd=0.75)+
  scale_color_manual(values=c("#4DBBD5FF", "#440154FF"))+
  geom_jitter(width=0.15)+
  #facet_wrap(~cluster,ncol=10,scales="free")+
  ylab("ThN / Total Th")+
  stat_compare_means(method="wilcox.test", 
                     show.legend=FALSE, 
                     method.args = list(var.equal = TRUE), 
                     label="p.value",
                     comparisons = age_comps,
                     label.x.npc=0.37, 
                     size=2.5, aes(group=Aged))+
  #scale_shape_manual(values=c(1:8,1:8))+
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(size=12, color="black"),
        axis.title.x = element_blank(),
        axis.line.x = element_line(size=0.25, color="black"),
        axis.line.y = element_line(size=0.25, color="black"),
        axis.ticks.x = element_blank(),
        axis.title.y = element_text(size=15, color="black"),
        strip.background = element_rect(fill=NA),
        strip.text = element_text(size=15, color="black"),
        panel.background = element_rect(fill="white"),
        legend.title = element_blank(),
        legend.key.size = unit(2,'lines'),
        legend.text = element_text(size=12),
        legend.key = element_rect(fill="white"))
pdf("ThN_Th_Ratio_OnTreatment.pdf", width=4, height=4);boxplot;dev.off()
write.csv(ggdf_ratio_OnTreatment, "ThN_Th_Ratio_OnTreatment.csv")

boxplot <- ggplot(ggdf_T0, aes(x=Aged, y=fc, color=Aged))+
  geom_boxplot(outlier.size=0, lwd=0.75)+
  scale_color_manual(values=c("#4DBBD5FF", "#440154FF"))+
  geom_jitter(width=0.15)+
  #facet_wrap(~cluster,ncol=10,scales="free")+
  ylab("Fold Change of ThN / Total Th")+
  stat_compare_means(method="wilcox.test", 
                     show.legend=FALSE, 
                     method.args = list(var.equal = TRUE), 
                     label="p.value",
                     comparisons = age_comps,
                     label.x.npc=0.37, 
                     size=2.5, aes(group=Aged))+
  #scale_shape_manual(values=c(1:8,1:8))+
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(size=12, color="black"),
        axis.title.x = element_blank(),
        axis.line.x = element_line(size=0.25, color="black"),
        axis.line.y = element_line(size=0.25, color="black"),
        axis.ticks.x = element_blank(),
        axis.title.y = element_text(size=15, color="black"),
        strip.background = element_rect(fill=NA),
        strip.text = element_text(size=15, color="black"),
        panel.background = element_rect(fill="white"),
        legend.title = element_blank(),
        legend.key.size = unit(2,'lines'),
        legend.text = element_text(size=12),
        legend.key = element_rect(fill="white"))
pdf("ThN_Th_Ratio_FC.pdf", width=4, height=4);boxplot;dev.off()
write.csv(ggdf_T0, "ThN_Th_Ratio_FC.csv")


####Plot ratio of ThCTL over all Th####
#make a new table containing counts only for ThN & Th clusters
counts_table_Th <- counts_table[rownames(counts_table)[rownames(counts_table) %in% c("ThCTL",
                                                                                     "Th2",
                                                                                     "Th2EM",
                                                                                     "Th2CM",
                                                                                     "Th17",
                                                                                     "Treg",
                                                                                     "ThN")],]

counts_table_ThCTL <- counts_table[rownames(counts_table)[rownames(counts_table) %in% c("ThCTL")],]

#make new object that divides the sum of TcN columns / sum of Tc columns & convert to df
ratio_table <- t(t(counts_table_ThCTL) / colSums(counts_table_Th)) 
ratio <- as.data.frame.matrix(ratio_table)
#ratio<-ratio %>% 
#filter_if(~is.numeric(.), all_vars(!is.infinite(.))) #remove INF numbers resulting from division by 0 (samples with no TAM2 cells)
#ratio<-na.omit(ratio) #remove NAs in data from division by 0

ggdf_ratio <- melt(data.frame(sample_id = rownames(ratio),ratio, check.names = FALSE),
                   id.vars = "sample_id", value.name = "ratio")
ggdf_ratio$sample_id <- factor(ggdf_ratio$sample_id, levels=samplevels)
ggdf_ratio$Aged <- factor(output$meta_data$Aged[match(ggdf_ratio$sample_id,output$meta_data$sample_id)])
ggdf_ratio$timepoint <- factor(output$meta_data$timepoint[match(ggdf_ratio$sample_id,output$meta_data$sample_id)])
ggdf_ratio_baseline <- ggdf_ratio[ggdf_ratio$timepoint %nin% c("1", "2"),]
ggdf_ratio_OnTreatment <- ggdf_ratio[ggdf_ratio$timepoint %nin% c("0", "2"),]
ggdf3 <- ggdf_ratio[ggdf_ratio$sample_id %nin% c("523_2","404_2","718_2","890_2","065_2","617_2","570_2","409_2","388_2","636_2","439_2","758_2","463_2","543_2","505_2","464_2","542_2","499_2"),]
ggdf_T0 <- ggdf3[ggdf3$timepoint %nin% c("1","2"),]
ggdf_T1 <- ggdf3[ggdf3$timepoint %nin% c("0"),]
ggdf_T0$post_ratio = ggdf_T1$ratio
ggdf_T0$fc = ggdf_T0$post_ratio/ggdf_T0$ratio

boxplot <- ggplot(ggdf_ratio_baseline, aes(x=Aged, y=ratio, color=Aged))+
  geom_boxplot(outlier.size=0, lwd=0.75)+
  scale_color_manual(values=c("#4DBBD5FF", "#440154FF"))+
  geom_jitter(width=0.15)+
  #facet_wrap(~cluster,ncol=10,scales="free")+
  ylab("ThCTL / Total Th")+
  stat_compare_means(method="wilcox.test", 
                     show.legend=FALSE, 
                     method.args = list(var.equal = TRUE), 
                     label="p.value",
                     comparisons = age_comps,
                     label.x.npc=0.37, 
                     size=2.5, aes(group=Aged))+
  #scale_shape_manual(values=c(1:8,1:8))+
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(size=12, color="black"),
        axis.title.x = element_blank(),
        axis.line.x = element_line(size=0.25, color="black"),
        axis.line.y = element_line(size=0.25, color="black"),
        axis.ticks.x = element_blank(),
        axis.title.y = element_text(size=15, color="black"),
        strip.background = element_rect(fill=NA),
        strip.text = element_text(size=15, color="black"),
        panel.background = element_rect(fill="white"),
        legend.title = element_blank(),
        legend.key.size = unit(2,'lines'),
        legend.text = element_text(size=12),
        legend.key = element_rect(fill="white"))
pdf("ThCTL_Th_Ratio_Baseline.pdf", width=4, height=4);boxplot;dev.off()
write.csv(ggdf_ratio_baseline, "ThCTL_Th_Ratio_Baseline.csv")


boxplot <- ggplot(ggdf_ratio_OnTreatment, aes(x=Aged, y=ratio, color=Aged))+
  geom_boxplot(outlier.size=0, lwd=0.75)+
  scale_color_manual(values=c("#4DBBD5FF", "#440154FF"))+
  geom_jitter(width=0.15)+
  #facet_wrap(~cluster,ncol=10,scales="free")+
  ylab("ThCTL / Total Th")+
  stat_compare_means(method="wilcox.test", 
                     show.legend=FALSE, 
                     method.args = list(var.equal = TRUE), 
                     label="p.value",
                     comparisons = age_comps,
                     label.x.npc=0.37, 
                     size=2.5, aes(group=Aged))+
  #scale_shape_manual(values=c(1:8,1:8))+
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(size=12, color="black"),
        axis.title.x = element_blank(),
        axis.line.x = element_line(size=0.25, color="black"),
        axis.line.y = element_line(size=0.25, color="black"),
        axis.ticks.x = element_blank(),
        axis.title.y = element_text(size=15, color="black"),
        strip.background = element_rect(fill=NA),
        strip.text = element_text(size=15, color="black"),
        panel.background = element_rect(fill="white"),
        legend.title = element_blank(),
        legend.key.size = unit(2,'lines'),
        legend.text = element_text(size=12),
        legend.key = element_rect(fill="white"))
pdf("ThCTL_Th_Ratio_OnTreatment.pdf", width=4, height=4);boxplot;dev.off()
write.csv(ggdf_ratio_OnTreatment, "ThCTL_Th_Ratio_OnTreatment.csv")

boxplot <- ggplot(ggdf_T0, aes(x=Aged, y=fc, color=Aged))+
  geom_boxplot(outlier.size=0, lwd=0.75)+
  scale_color_manual(values=c("#4DBBD5FF", "#440154FF"))+
  geom_jitter(width=0.15)+
  #facet_wrap(~cluster,ncol=10,scales="free")+
  ylab("Fold Change of ThCTL / Total Th")+
  stat_compare_means(method="wilcox.test", 
                     show.legend=FALSE, 
                     method.args = list(var.equal = TRUE), 
                     label="p.value",
                     comparisons = age_comps,
                     label.x.npc=0.37, 
                     size=2.5, aes(group=Aged))+
  #scale_shape_manual(values=c(1:8,1:8))+
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(size=12, color="black"),
        axis.title.x = element_blank(),
        axis.line.x = element_line(size=0.25, color="black"),
        axis.line.y = element_line(size=0.25, color="black"),
        axis.ticks.x = element_blank(),
        axis.title.y = element_text(size=15, color="black"),
        strip.background = element_rect(fill=NA),
        strip.text = element_text(size=15, color="black"),
        panel.background = element_rect(fill="white"),
        legend.title = element_blank(),
        legend.key.size = unit(2,'lines'),
        legend.text = element_text(size=12),
        legend.key = element_rect(fill="white"))
pdf("ThCTL_Th_Ratio_FC.pdf", width=4, height=4);boxplot;dev.off()
write.csv(ggdf_T0, "ThN_Th_Ratio_FC.csv")


####Plot ratio of Th2 over all Th####
#make a new table containing counts only for ThN & Th clusters
counts_table_Th <- counts_table[rownames(counts_table)[rownames(counts_table) %in% c("ThCTL",
                                                                                     "Th2",
                                                                                     "Th2EM",
                                                                                     "Th2CM",
                                                                                     "Th17",
                                                                                     "Treg",
                                                                                     "ThN")],]

counts_table_Th2 <- counts_table[rownames(counts_table)[rownames(counts_table) %in% c("Th2")],]

#make new object that divides the sum of TcN columns / sum of Tc columns & convert to df
ratio_table <- t(t(counts_table_Th2) / colSums(counts_table_Th)) 
ratio <- as.data.frame.matrix(ratio_table)
#ratio<-ratio %>% 
#filter_if(~is.numeric(.), all_vars(!is.infinite(.))) #remove INF numbers resulting from division by 0 (samples with no TAM2 cells)
#ratio<-na.omit(ratio) #remove NAs in data from division by 0

ggdf_ratio <- melt(data.frame(sample_id = rownames(ratio),ratio, check.names = FALSE),
                   id.vars = "sample_id", value.name = "ratio")
ggdf_ratio$sample_id <- factor(ggdf_ratio$sample_id, levels=samplevels)
ggdf_ratio$Aged <- factor(output$meta_data$Aged[match(ggdf_ratio$sample_id,output$meta_data$sample_id)])
ggdf_ratio$timepoint <- factor(output$meta_data$timepoint[match(ggdf_ratio$sample_id,output$meta_data$sample_id)])
ggdf_ratio_baseline <- ggdf_ratio[ggdf_ratio$timepoint %nin% c("1", "2"),]
ggdf_ratio_OnTreatment <- ggdf_ratio[ggdf_ratio$timepoint %nin% c("0", "2"),]
ggdf3 <- ggdf_ratio[ggdf_ratio$sample_id %nin% c("523_2","404_2","718_2","890_2","065_2","617_2","570_2","409_2","388_2","636_2","439_2","758_2","463_2","543_2","505_2","464_2","542_2","499_2"),]
ggdf_T0 <- ggdf3[ggdf3$timepoint %nin% c("1","2"),]
ggdf_T1 <- ggdf3[ggdf3$timepoint %nin% c("0"),]
ggdf_T0$post_ratio = ggdf_T1$ratio
ggdf_T0$fc = ggdf_T0$post_ratio/ggdf_T0$ratio

boxplot <- ggplot(ggdf_ratio_baseline, aes(x=Aged, y=ratio, color=Aged))+
  geom_boxplot(outlier.size=0, lwd=0.75)+
  scale_color_manual(values=c("#4DBBD5FF", "#440154FF"))+
  geom_jitter(width=0.15)+
  #facet_wrap(~cluster,ncol=10,scales="free")+
  ylab("Th2 / Total Th")+
  stat_compare_means(method="wilcox.test", 
                     show.legend=FALSE, 
                     method.args = list(var.equal = TRUE), 
                     label="p.value",
                     comparisons = age_comps,
                     label.x.npc=0.37, 
                     size=2.5, aes(group=Aged))+
  #scale_shape_manual(values=c(1:8,1:8))+
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(size=12, color="black"),
        axis.title.x = element_blank(),
        axis.line.x = element_line(size=0.25, color="black"),
        axis.line.y = element_line(size=0.25, color="black"),
        axis.ticks.x = element_blank(),
        axis.title.y = element_text(size=15, color="black"),
        strip.background = element_rect(fill=NA),
        strip.text = element_text(size=15, color="black"),
        panel.background = element_rect(fill="white"),
        legend.title = element_blank(),
        legend.key.size = unit(2,'lines'),
        legend.text = element_text(size=12),
        legend.key = element_rect(fill="white"))
pdf("Th2_Th_Ratio_Baseline.pdf", width=4, height=4);boxplot;dev.off()
write.csv(ggdf_ratio_baseline, "Th2_Th_Ratio_Baseline.csv")

boxplot <- ggplot(ggdf_ratio_OnTreatment, aes(x=Aged, y=ratio, color=Aged))+
  geom_boxplot(outlier.size=0, lwd=0.75)+
  scale_color_manual(values=c("#4DBBD5FF", "#440154FF"))+
  geom_jitter(width=0.15)+
  #facet_wrap(~cluster,ncol=10,scales="free")+
  ylab("Th2 / Total Th")+
  stat_compare_means(method="wilcox.test", 
                     show.legend=FALSE, 
                     method.args = list(var.equal = TRUE), 
                     label="p.value",
                     comparisons = age_comps,
                     label.x.npc=0.37, 
                     size=2.5, aes(group=Aged))+
  #scale_shape_manual(values=c(1:8,1:8))+
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(size=12, color="black"),
        axis.title.x = element_blank(),
        axis.line.x = element_line(size=0.25, color="black"),
        axis.line.y = element_line(size=0.25, color="black"),
        axis.ticks.x = element_blank(),
        axis.title.y = element_text(size=15, color="black"),
        strip.background = element_rect(fill=NA),
        strip.text = element_text(size=15, color="black"),
        panel.background = element_rect(fill="white"),
        legend.title = element_blank(),
        legend.key.size = unit(2,'lines'),
        legend.text = element_text(size=12),
        legend.key = element_rect(fill="white"))
pdf("Th2_Th_Ratio_OnTreatment.pdf", width=4, height=4);boxplot;dev.off()
write.csv(ggdf_ratio_OnTreatment, "Th2_Th_Ratio_OnTreatment.csv")

boxplot <- ggplot(ggdf_T0, aes(x=Aged, y=fc, color=Aged))+
  geom_boxplot(outlier.size=0, lwd=0.75)+
  scale_color_manual(values=c("#4DBBD5FF", "#440154FF"))+
  geom_jitter(width=0.15)+
  #facet_wrap(~cluster,ncol=10,scales="free")+
  ylab("Fold Change of Th2 / Total Th")+
  stat_compare_means(method="wilcox.test", 
                     show.legend=FALSE, 
                     method.args = list(var.equal = TRUE), 
                     label="p.value",
                     comparisons = age_comps,
                     label.x.npc=0.37, 
                     size=2.5, aes(group=Aged))+
  #scale_shape_manual(values=c(1:8,1:8))+
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(size=12, color="black"),
        axis.title.x = element_blank(),
        axis.line.x = element_line(size=0.25, color="black"),
        axis.line.y = element_line(size=0.25, color="black"),
        axis.ticks.x = element_blank(),
        axis.title.y = element_text(size=15, color="black"),
        strip.background = element_rect(fill=NA),
        strip.text = element_text(size=15, color="black"),
        panel.background = element_rect(fill="white"),
        legend.title = element_blank(),
        legend.key.size = unit(2,'lines'),
        legend.text = element_text(size=12),
        legend.key = element_rect(fill="white"))
pdf("Th2_Th_Ratio_FC.pdf", width=4, height=4);boxplot;dev.off()
write.csv(ggdf_T0, "Th2_Th_Ratio_FC.csv")

####Plot ratio of Th2EM over all Th####
#make a new table containing counts only for ThN & Th clusters
counts_table_Th <- counts_table[rownames(counts_table)[rownames(counts_table) %in% c("ThCTL",
                                                                                     "Th2",
                                                                                     "Th2EM",
                                                                                     "Th2CM",
                                                                                     "Th17",
                                                                                     "Treg",
                                                                                     "ThN")],]

counts_table_Th2EM <- counts_table[rownames(counts_table)[rownames(counts_table) %in% c("Th2EM")],]

#make new object that divides the sum of TcN columns / sum of Tc columns & convert to df
ratio_table <- t(t(counts_table_Th2EM) / colSums(counts_table_Th)) 
ratio <- as.data.frame.matrix(ratio_table)
#ratio<-ratio %>% 
#filter_if(~is.numeric(.), all_vars(!is.infinite(.))) #remove INF numbers resulting from division by 0 (samples with no TAM2 cells)
#ratio<-na.omit(ratio) #remove NAs in data from division by 0

ggdf_ratio <- melt(data.frame(sample_id = rownames(ratio),ratio, check.names = FALSE),
                   id.vars = "sample_id", value.name = "ratio")
ggdf_ratio$sample_id <- factor(ggdf_ratio$sample_id, levels=samplevels)
ggdf_ratio$Aged <- factor(output$meta_data$Aged[match(ggdf_ratio$sample_id,output$meta_data$sample_id)])
ggdf_ratio$timepoint <- factor(output$meta_data$timepoint[match(ggdf_ratio$sample_id,output$meta_data$sample_id)])
ggdf_ratio_baseline <- ggdf_ratio[ggdf_ratio$timepoint %nin% c("1", "2"),]
ggdf_ratio_OnTreatment <- ggdf_ratio[ggdf_ratio$timepoint %nin% c("0", "2"),]
ggdf3 <- ggdf_ratio[ggdf_ratio$sample_id %nin% c("523_2","404_2","718_2","890_2","065_2","617_2","570_2","409_2","388_2","636_2","439_2","758_2","463_2","543_2","505_2","464_2","542_2","499_2"),]
ggdf_T0 <- ggdf3[ggdf3$timepoint %nin% c("1","2"),]
ggdf_T1 <- ggdf3[ggdf3$timepoint %nin% c("0"),]
ggdf_T0$post_ratio = ggdf_T1$ratio
ggdf_T0$fc = ggdf_T0$post_ratio/ggdf_T0$ratio

boxplot <- ggplot(ggdf_ratio_baseline, aes(x=Aged, y=ratio, color=Aged))+
  geom_boxplot(outlier.size=0, lwd=0.75)+
  scale_color_manual(values=c("#4DBBD5FF", "#440154FF"))+
  geom_jitter(width=0.15)+
  #facet_wrap(~cluster,ncol=10,scales="free")+
  ylab("Th2EM / Total Th")+
  stat_compare_means(method="wilcox.test", 
                     show.legend=FALSE, 
                     method.args = list(var.equal = TRUE), 
                     label="p.value",
                     comparisons = age_comps,
                     label.x.npc=0.37, 
                     size=2.5, aes(group=Aged))+
  #scale_shape_manual(values=c(1:8,1:8))+
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(size=12, color="black"),
        axis.title.x = element_blank(),
        axis.line.x = element_line(size=0.25, color="black"),
        axis.line.y = element_line(size=0.25, color="black"),
        axis.ticks.x = element_blank(),
        axis.title.y = element_text(size=15, color="black"),
        strip.background = element_rect(fill=NA),
        strip.text = element_text(size=15, color="black"),
        panel.background = element_rect(fill="white"),
        legend.title = element_blank(),
        legend.key.size = unit(2,'lines'),
        legend.text = element_text(size=12),
        legend.key = element_rect(fill="white"))
pdf("Th2EM_Th_Ratio_baseline.pdf", width=4, height=4);boxplot;dev.off()
write.csv(ggdf_ratio_baseline, "Th2EM_Th_Ratio_Baseline.csv")

boxplot <- ggplot(ggdf_ratio_OnTreatment, aes(x=Aged, y=ratio, color=Aged))+
  geom_boxplot(outlier.size=0, lwd=0.75)+
  scale_color_manual(values=c("#4DBBD5FF", "#440154FF"))+
  geom_jitter(width=0.15)+
  #facet_wrap(~cluster,ncol=10,scales="free")+
  ylab("Th2EM / Total Th")+
  stat_compare_means(method="wilcox.test", 
                     show.legend=FALSE, 
                     method.args = list(var.equal = TRUE), 
                     label="p.value",
                     comparisons = age_comps,
                     label.x.npc=0.37, 
                     size=2.5, aes(group=Aged))+
  #scale_shape_manual(values=c(1:8,1:8))+
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(size=12, color="black"),
        axis.title.x = element_blank(),
        axis.line.x = element_line(size=0.25, color="black"),
        axis.line.y = element_line(size=0.25, color="black"),
        axis.ticks.x = element_blank(),
        axis.title.y = element_text(size=15, color="black"),
        strip.background = element_rect(fill=NA),
        strip.text = element_text(size=15, color="black"),
        panel.background = element_rect(fill="white"),
        legend.title = element_blank(),
        legend.key.size = unit(2,'lines'),
        legend.text = element_text(size=12),
        legend.key = element_rect(fill="white"))
pdf("Th2EM_Th_Ratio_OnTreatment.pdf", width=4, height=4);boxplot;dev.off()
write.csv(ggdf_ratio_OnTreatment, "Th2EM_Th_Ratio_OnTreatment.csv")

boxplot <- ggplot(ggdf_T0, aes(x=Aged, y=fc, color=Aged))+
  geom_boxplot(outlier.size=0, lwd=0.75)+
  scale_color_manual(values=c("#4DBBD5FF", "#440154FF"))+
  geom_jitter(width=0.15)+
  #facet_wrap(~cluster,ncol=10,scales="free")+
  ylab("Fold Change of Th2EM / Total Th")+
  stat_compare_means(method="wilcox.test", 
                     show.legend=FALSE, 
                     method.args = list(var.equal = TRUE), 
                     label="p.value",
                     comparisons = age_comps,
                     label.x.npc=0.37, 
                     size=2.5, aes(group=Aged))+
  #scale_shape_manual(values=c(1:8,1:8))+
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(size=12, color="black"),
        axis.title.x = element_blank(),
        axis.line.x = element_line(size=0.25, color="black"),
        axis.line.y = element_line(size=0.25, color="black"),
        axis.ticks.x = element_blank(),
        axis.title.y = element_text(size=15, color="black"),
        strip.background = element_rect(fill=NA),
        strip.text = element_text(size=15, color="black"),
        panel.background = element_rect(fill="white"),
        legend.title = element_blank(),
        legend.key.size = unit(2,'lines'),
        legend.text = element_text(size=12),
        legend.key = element_rect(fill="white"))
pdf("Th2EM_Th_Ratio_FC.pdf", width=4, height=4);boxplot;dev.off()
write.csv(ggdf_T0, "Th2EM_Th_Ratio_FC.csv")

####Plot ratio of Th2CM over all Th####
#make a new table containing counts only for ThN & Th clusters
counts_table_Th <- counts_table[rownames(counts_table)[rownames(counts_table) %in% c("ThCTL",
                                                                                     "Th2",
                                                                                     "Th2EM",
                                                                                     "Th2CM",
                                                                                     "Th17",
                                                                                     "Treg",
                                                                                     "ThN")],]

counts_table_Th2CM <- counts_table[rownames(counts_table)[rownames(counts_table) %in% c("Th2CM")],]

#make new object that divides the sum of TcN columns / sum of Tc columns & convert to df
ratio_table <- t(t(counts_table_Th2CM) / colSums(counts_table_Th)) 
ratio <- as.data.frame.matrix(ratio_table)
#ratio<-ratio %>% 
#filter_if(~is.numeric(.), all_vars(!is.infinite(.))) #remove INF numbers resulting from division by 0 (samples with no TAM2 cells)
#ratio<-na.omit(ratio) #remove NAs in data from division by 0

ggdf_ratio <- melt(data.frame(sample_id = rownames(ratio),ratio, check.names = FALSE),
                   id.vars = "sample_id", value.name = "ratio")
ggdf_ratio$sample_id <- factor(ggdf_ratio$sample_id, levels=samplevels)
ggdf_ratio$Aged <- factor(output$meta_data$Aged[match(ggdf_ratio$sample_id,output$meta_data$sample_id)])
ggdf_ratio$timepoint <- factor(output$meta_data$timepoint[match(ggdf_ratio$sample_id,output$meta_data$sample_id)])
ggdf_ratio_baseline <- ggdf_ratio[ggdf_ratio$timepoint %nin% c("1", "2"),]
ggdf_ratio_OnTreatment <- ggdf_ratio[ggdf_ratio$timepoint %nin% c("0", "2"),]
ggdf3 <- ggdf_ratio[ggdf_ratio$sample_id %nin% c("523_2","404_2","718_2","890_2","065_2","617_2","570_2","409_2","388_2","636_2","439_2","758_2","463_2","543_2","505_2","464_2","542_2","499_2"),]
ggdf_T0 <- ggdf3[ggdf3$timepoint %nin% c("1","2"),]
ggdf_T1 <- ggdf3[ggdf3$timepoint %nin% c("0"),]
ggdf_T0$post_ratio = ggdf_T1$ratio
ggdf_T0$fc = ggdf_T0$post_ratio/ggdf_T0$ratio

boxplot <- ggplot(ggdf_ratio_baseline, aes(x=Aged, y=ratio, color=Aged))+
  geom_boxplot(outlier.size=0, lwd=0.75)+
  scale_color_manual(values=c("#4DBBD5FF", "#440154FF"))+
  geom_jitter(width=0.15)+
  #facet_wrap(~cluster,ncol=10,scales="free")+
  ylab("Th2CM / Total Th")+
  stat_compare_means(method="wilcox.test", 
                     show.legend=FALSE, 
                     method.args = list(var.equal = TRUE), 
                     label="p.value",
                     comparisons = age_comps,
                     label.x.npc=0.37, 
                     size=2.5, aes(group=Aged))+
  #scale_shape_manual(values=c(1:8,1:8))+
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(size=12, color="black"),
        axis.title.x = element_blank(),
        axis.line.x = element_line(size=0.25, color="black"),
        axis.line.y = element_line(size=0.25, color="black"),
        axis.ticks.x = element_blank(),
        axis.title.y = element_text(size=15, color="black"),
        strip.background = element_rect(fill=NA),
        strip.text = element_text(size=15, color="black"),
        panel.background = element_rect(fill="white"),
        legend.title = element_blank(),
        legend.key.size = unit(2,'lines'),
        legend.text = element_text(size=12),
        legend.key = element_rect(fill="white"))
pdf("Th2CM_Th_Ratio_Baseline.pdf", width=4, height=4);boxplot;dev.off()
write.csv(ggdf_ratio_baseline, "Th2CM_Th_Ratio_Baseline.csv")

boxplot <- ggplot(ggdf_ratio_OnTreatment, aes(x=Aged, y=ratio, color=Aged))+
  geom_boxplot(outlier.size=0, lwd=0.75)+
  scale_color_manual(values=c("#4DBBD5FF", "#440154FF"))+
  geom_jitter(width=0.15)+
  #facet_wrap(~cluster,ncol=10,scales="free")+
  ylab("Th2CM / Total Th")+
  stat_compare_means(method="wilcox.test", 
                     show.legend=FALSE, 
                     method.args = list(var.equal = TRUE), 
                     label="p.value",
                     comparisons = age_comps,
                     label.x.npc=0.37, 
                     size=2.5, aes(group=Aged))+
  #scale_shape_manual(values=c(1:8,1:8))+
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(size=12, color="black"),
        axis.title.x = element_blank(),
        axis.line.x = element_line(size=0.25, color="black"),
        axis.line.y = element_line(size=0.25, color="black"),
        axis.ticks.x = element_blank(),
        axis.title.y = element_text(size=15, color="black"),
        strip.background = element_rect(fill=NA),
        strip.text = element_text(size=15, color="black"),
        panel.background = element_rect(fill="white"),
        legend.title = element_blank(),
        legend.key.size = unit(2,'lines'),
        legend.text = element_text(size=12),
        legend.key = element_rect(fill="white"))
pdf("Th2CM_Th_Ratio_OnTreatment.pdf", width=4, height=4);boxplot;dev.off()
write.csv(ggdf_ratio_OnTreatment, "Th2CM_Th_Ratio_OnTreatment.csv")

boxplot <- ggplot(ggdf_T0, aes(x=Aged, y=fc, color=Aged))+
  geom_boxplot(outlier.size=0, lwd=0.75)+
  scale_color_manual(values=c("#4DBBD5FF", "#440154FF"))+
  geom_jitter(width=0.15)+
  #facet_wrap(~cluster,ncol=10,scales="free")+
  ylab("Fold Change of Th2CM / Total Th")+
  stat_compare_means(method="wilcox.test", 
                     show.legend=FALSE, 
                     method.args = list(var.equal = TRUE), 
                     label="p.value",
                     comparisons = age_comps,
                     label.x.npc=0.37, 
                     size=2.5, aes(group=Aged))+
  #scale_shape_manual(values=c(1:8,1:8))+
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(size=12, color="black"),
        axis.title.x = element_blank(),
        axis.line.x = element_line(size=0.25, color="black"),
        axis.line.y = element_line(size=0.25, color="black"),
        axis.ticks.x = element_blank(),
        axis.title.y = element_text(size=15, color="black"),
        strip.background = element_rect(fill=NA),
        strip.text = element_text(size=15, color="black"),
        panel.background = element_rect(fill="white"),
        legend.title = element_blank(),
        legend.key.size = unit(2,'lines'),
        legend.text = element_text(size=12),
        legend.key = element_rect(fill="white"))
pdf("Th2CM_Th_Ratio_FC.pdf", width=4, height=4);boxplot;dev.off()
write.csv(ggdf_T0, "Th2CM_Th_Ratio_FC.csv")

####Plot ratio of Th17 over all Th####
#make a new table containing counts only for ThN & Th clusters
counts_table_Th <- counts_table[rownames(counts_table)[rownames(counts_table) %in% c("ThCTL",
                                                                                     "Th2",
                                                                                     "Th2EM",
                                                                                     "Th2CM",
                                                                                     "Th17",
                                                                                     "Treg",
                                                                                     "ThN")],]

counts_table_Th17 <- counts_table[rownames(counts_table)[rownames(counts_table) %in% c("Th17")],]

#make new object that divides the sum of TcN columns / sum of Tc columns & convert to df
ratio_table <- t(t(counts_table_Th17) / colSums(counts_table_Th)) 
ratio <- as.data.frame.matrix(ratio_table)
#ratio<-ratio %>% 
#filter_if(~is.numeric(.), all_vars(!is.infinite(.))) #remove INF numbers resulting from division by 0 (samples with no TAM2 cells)
#ratio<-na.omit(ratio) #remove NAs in data from division by 0

ggdf_ratio <- melt(data.frame(sample_id = rownames(ratio),ratio, check.names = FALSE),
                   id.vars = "sample_id", value.name = "ratio")
ggdf_ratio$sample_id <- factor(ggdf_ratio$sample_id, levels=samplevels)
ggdf_ratio$Aged <- factor(output$meta_data$Aged[match(ggdf_ratio$sample_id,output$meta_data$sample_id)])
ggdf_ratio$timepoint <- factor(output$meta_data$timepoint[match(ggdf_ratio$sample_id,output$meta_data$sample_id)])
ggdf_ratio_baseline <- ggdf_ratio[ggdf_ratio$timepoint %nin% c("1", "2"),]
ggdf_ratio_OnTreatment <- ggdf_ratio[ggdf_ratio$timepoint %nin% c("0", "2"),]
ggdf3 <- ggdf_ratio[ggdf_ratio$sample_id %nin% c("523_2","404_2","718_2","890_2","065_2","617_2","570_2","409_2","388_2","636_2","439_2","758_2","463_2","543_2","505_2","464_2","542_2","499_2"),]
ggdf_T0 <- ggdf3[ggdf3$timepoint %nin% c("1","2"),]
ggdf_T1 <- ggdf3[ggdf3$timepoint %nin% c("0"),]
ggdf_T0$post_ratio = ggdf_T1$ratio
ggdf_T0$fc = ggdf_T0$post_ratio/ggdf_T0$ratio

boxplot <- ggplot(ggdf_ratio_baseline, aes(x=Aged, y=ratio, color=Aged))+
  geom_boxplot(outlier.size=0, lwd=0.75)+
  scale_color_manual(values=c("#4DBBD5FF", "#440154FF"))+
  geom_jitter(width=0.15)+
  #facet_wrap(~cluster,ncol=10,scales="free")+
  ylab("Th17 / Total Th")+
  stat_compare_means(method="wilcox.test", 
                     show.legend=FALSE, 
                     method.args = list(var.equal = TRUE), 
                     label="p.value",
                     comparisons = age_comps,
                     label.x.npc=0.37, 
                     size=2.5, aes(group=Aged))+
  #scale_shape_manual(values=c(1:8,1:8))+
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(size=12, color="black"),
        axis.title.x = element_blank(),
        axis.line.x = element_line(size=0.25, color="black"),
        axis.line.y = element_line(size=0.25, color="black"),
        axis.ticks.x = element_blank(),
        axis.title.y = element_text(size=15, color="black"),
        strip.background = element_rect(fill=NA),
        strip.text = element_text(size=15, color="black"),
        panel.background = element_rect(fill="white"),
        legend.title = element_blank(),
        legend.key.size = unit(2,'lines'),
        legend.text = element_text(size=12),
        legend.key = element_rect(fill="white"))
pdf("Th17_Th_Ratio_Baseline.pdf", width=4, height=4);boxplot;dev.off()
write.csv(ggdf_ratio_baseline, "Th17_Th_Ratio_Baseline.csv")

boxplot <- ggplot(ggdf_ratio_OnTreatment, aes(x=Aged, y=ratio, color=Aged))+
  geom_boxplot(outlier.size=0, lwd=0.75)+
  scale_color_manual(values=c("#4DBBD5FF", "#440154FF"))+
  geom_jitter(width=0.15)+
  #facet_wrap(~cluster,ncol=10,scales="free")+
  ylab("Th17 / Total Th")+
  stat_compare_means(method="wilcox.test", 
                     show.legend=FALSE, 
                     method.args = list(var.equal = TRUE), 
                     label="p.value",
                     comparisons = age_comps,
                     label.x.npc=0.37, 
                     size=2.5, aes(group=Aged))+
  #scale_shape_manual(values=c(1:8,1:8))+
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(size=12, color="black"),
        axis.title.x = element_blank(),
        axis.line.x = element_line(size=0.25, color="black"),
        axis.line.y = element_line(size=0.25, color="black"),
        axis.ticks.x = element_blank(),
        axis.title.y = element_text(size=15, color="black"),
        strip.background = element_rect(fill=NA),
        strip.text = element_text(size=15, color="black"),
        panel.background = element_rect(fill="white"),
        legend.title = element_blank(),
        legend.key.size = unit(2,'lines'),
        legend.text = element_text(size=12),
        legend.key = element_rect(fill="white"))
pdf("Th17_Th_Ratio_OnTreatment.pdf", width=4, height=4);boxplot;dev.off()
write.csv(ggdf_ratio_OnTreatment, "Th17_Th_Ratio_OnTreatment.csv")

boxplot <- ggplot(ggdf_T0, aes(x=Aged, y=fc, color=Aged))+
  geom_boxplot(outlier.size=0, lwd=0.75)+
  scale_color_manual(values=c("#4DBBD5FF", "#440154FF"))+
  geom_jitter(width=0.15)+
  #facet_wrap(~cluster,ncol=10,scales="free")+
  ylab("Fold Change of Th17 / Total Th")+
  stat_compare_means(method="wilcox.test", 
                     show.legend=FALSE, 
                     method.args = list(var.equal = TRUE), 
                     label="p.value",
                     comparisons = age_comps,
                     label.x.npc=0.37, 
                     size=2.5, aes(group=Aged))+
  #scale_shape_manual(values=c(1:8,1:8))+
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(size=12, color="black"),
        axis.title.x = element_blank(),
        axis.line.x = element_line(size=0.25, color="black"),
        axis.line.y = element_line(size=0.25, color="black"),
        axis.ticks.x = element_blank(),
        axis.title.y = element_text(size=15, color="black"),
        strip.background = element_rect(fill=NA),
        strip.text = element_text(size=15, color="black"),
        panel.background = element_rect(fill="white"),
        legend.title = element_blank(),
        legend.key.size = unit(2,'lines'),
        legend.text = element_text(size=12),
        legend.key = element_rect(fill="white"))
pdf("Th17_Th_Ratio_FC.pdf", width=4, height=4);boxplot;dev.off()
write.csv(ggdf_T0, "Th17_Th_Ratio_FC.csv")

####Plot ratio of Treg over all Th####
#make a new table containing counts only for ThN & Th clusters
counts_table_Th <- counts_table[rownames(counts_table)[rownames(counts_table) %in% c("ThCTL",
                                                                                     "Th2",
                                                                                     "Th2EM",
                                                                                     "Th2CM",
                                                                                     "Th17",
                                                                                     "Treg",
                                                                                     "ThN")],]

counts_table_Treg <- counts_table[rownames(counts_table)[rownames(counts_table) %in% c("Treg")],]

#make new object that divides the sum of TcN columns / sum of Tc columns & convert to df
ratio_table <- t(t(counts_table_Treg) / colSums(counts_table_Th)) 
ratio <- as.data.frame.matrix(ratio_table)
#ratio<-ratio %>% 
#filter_if(~is.numeric(.), all_vars(!is.infinite(.))) #remove INF numbers resulting from division by 0 (samples with no TAM2 cells)
#ratio<-na.omit(ratio) #remove NAs in data from division by 0

ggdf_ratio <- melt(data.frame(sample_id = rownames(ratio),ratio, check.names = FALSE),
                   id.vars = "sample_id", value.name = "ratio")
ggdf_ratio$sample_id <- factor(ggdf_ratio$sample_id, levels=samplevels)
ggdf_ratio$Aged <- factor(output$meta_data$Aged[match(ggdf_ratio$sample_id,output$meta_data$sample_id)])
ggdf_ratio$timepoint <- factor(output$meta_data$timepoint[match(ggdf_ratio$sample_id,output$meta_data$sample_id)])
ggdf_ratio_baseline <- ggdf_ratio[ggdf_ratio$timepoint %nin% c("1", "2"),]
ggdf_ratio_OnTreatment <- ggdf_ratio[ggdf_ratio$timepoint %nin% c("0", "2"),]
ggdf3 <- ggdf_ratio[ggdf_ratio$sample_id %nin% c("523_2","404_2","718_2","890_2","065_2","617_2","570_2","409_2","388_2","636_2","439_2","758_2","463_2","543_2","505_2","464_2","542_2","499_2"),]
ggdf_T0 <- ggdf3[ggdf3$timepoint %nin% c("1","2"),]
ggdf_T1 <- ggdf3[ggdf3$timepoint %nin% c("0"),]
ggdf_T0$post_ratio = ggdf_T1$ratio
ggdf_T0$fc = ggdf_T0$post_ratio/ggdf_T0$ratio

boxplot <- ggplot(ggdf_ratio_baseline, aes(x=Aged, y=ratio, color=Aged))+
  geom_boxplot(outlier.size=0, lwd=0.75)+
  scale_color_manual(values=c("#4DBBD5FF", "#440154FF"))+
  geom_jitter(width=0.15)+
  #facet_wrap(~cluster,ncol=10,scales="free")+
  ylab("Treg / Total Th")+
  stat_compare_means(method="wilcox.test", 
                     show.legend=FALSE, 
                     method.args = list(var.equal = TRUE), 
                     label="p.value",
                     comparisons = age_comps,
                     label.x.npc=0.37, 
                     size=2.5, aes(group=Aged))+
  #scale_shape_manual(values=c(1:8,1:8))+
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(size=12, color="black"),
        axis.title.x = element_blank(),
        axis.line.x = element_line(size=0.25, color="black"),
        axis.line.y = element_line(size=0.25, color="black"),
        axis.ticks.x = element_blank(),
        axis.title.y = element_text(size=15, color="black"),
        strip.background = element_rect(fill=NA),
        strip.text = element_text(size=15, color="black"),
        panel.background = element_rect(fill="white"),
        legend.title = element_blank(),
        legend.key.size = unit(2,'lines'),
        legend.text = element_text(size=12),
        legend.key = element_rect(fill="white"))
pdf("Treg_Th_Ratio_Baseline.pdf", width=4, height=4);boxplot;dev.off()
write.csv(ggdf_ratio_baseline, "Treg_Th_Ratio_Baseline.csv")

boxplot <- ggplot(ggdf_ratio_OnTreatment, aes(x=Aged, y=ratio, color=Aged))+
  geom_boxplot(outlier.size=0, lwd=0.75)+
  scale_color_manual(values=c("#4DBBD5FF", "#440154FF"))+
  geom_jitter(width=0.15)+
  #facet_wrap(~cluster,ncol=10,scales="free")+
  ylab("Treg / Total Th")+
  stat_compare_means(method="wilcox.test", 
                     show.legend=FALSE, 
                     method.args = list(var.equal = TRUE), 
                     label="p.value",
                     comparisons = age_comps,
                     label.x.npc=0.37, 
                     size=2.5, aes(group=Aged))+
  #scale_shape_manual(values=c(1:8,1:8))+
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(size=12, color="black"),
        axis.title.x = element_blank(),
        axis.line.x = element_line(size=0.25, color="black"),
        axis.line.y = element_line(size=0.25, color="black"),
        axis.ticks.x = element_blank(),
        axis.title.y = element_text(size=15, color="black"),
        strip.background = element_rect(fill=NA),
        strip.text = element_text(size=15, color="black"),
        panel.background = element_rect(fill="white"),
        legend.title = element_blank(),
        legend.key.size = unit(2,'lines'),
        legend.text = element_text(size=12),
        legend.key = element_rect(fill="white"))
pdf("Treg_Th_Ratio_OnTreatment.pdf", width=4, height=4);boxplot;dev.off()
write.csv(ggdf_ratio_OnTreatment, "Treg_Th_Ratio_OnTreatment.csv")

boxplot <- ggplot(ggdf_T0, aes(x=Aged, y=fc, color=Aged))+
  geom_boxplot(outlier.size=0, lwd=0.75)+
  scale_color_manual(values=c("#4DBBD5FF", "#440154FF"))+
  geom_jitter(width=0.15)+
  #facet_wrap(~cluster,ncol=10,scales="free")+
  ylab("Fold Change of Treg / Total Th")+
  stat_compare_means(method="wilcox.test", 
                     show.legend=FALSE, 
                     method.args = list(var.equal = TRUE), 
                     label="p.value",
                     comparisons = age_comps,
                     label.x.npc=0.37, 
                     size=2.5, aes(group=Aged))+
  #scale_shape_manual(values=c(1:8,1:8))+
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(size=12, color="black"),
        axis.title.x = element_blank(),
        axis.line.x = element_line(size=0.25, color="black"),
        axis.line.y = element_line(size=0.25, color="black"),
        axis.ticks.x = element_blank(),
        axis.title.y = element_text(size=15, color="black"),
        strip.background = element_rect(fill=NA),
        strip.text = element_text(size=15, color="black"),
        panel.background = element_rect(fill="white"),
        legend.title = element_blank(),
        legend.key.size = unit(2,'lines'),
        legend.text = element_text(size=12),
        legend.key = element_rect(fill="white"))
pdf("Treg_Th_Ratio_FC.pdf", width=4, height=4);boxplot;dev.off()
write.csv(ggdf_T0, "Treg_Th_Ratio_FC.csv")

####Func Markers Age####
fmlistplot <- c("CCR6",
                "CCR5",
                "TBET",
                "CD25",
                "OX40",
                "CCR3",
                "TIM3",
                "CCR10",
                "PDL1",
                "GATA3",
                "CCR7",
                "CD28",
                "CTLA4",
                "FOXP3",
                "CXCR3",
                "CCR4",
                "CD27",
                "RORY",
                "GZMB",
                "KI67",
                "PD1",
                "LAG3",
                "TIGIT",
                "CD137",
                "HLADR")


exprtbl <- 
  data.frame(fsApply(output$fcs,exprs)[, union(output$subtype_markers,output$functional_markers)],
             sample_id = output$sample_ids, cluster = output$cell_clustering1m) %>%
  group_by(sample_id, cluster) %>%
  summarize_all(funs(mean))

Age_levels=c("No", "Yes")
Time_levels=c("0", "1", "2")
response_levels=c("no", "yes", "NA")

ggdf_Func <-melt(exprtbl, id.var=c("cluster","sample_id"))
ggdf_Func$cluster <- factor(ggdf_Func$cluster, levels=clusterlevels)
ggdf_Func$Aged<- factor(output$meta_data$Aged[match(ggdf_Func$sample_id,output$meta_data$sample_id)],levels=Age_levels)
ggdf_Func$timepoint <- factor (output$meta_data$timepoint[match(ggdf_Func$sample_id,output$meta_data$sample_id)],levels=Time_levels)
ggdf_Func$sample_id <- factor(ggdf_Func$sample_id, levels = samplevels)
ggdf_Func$response <- factor(output$meta_data$response[match(ggdf_Func$sample_id,output$meta_data$sample_id)],levels=response_levels)
ggdf_Func<- ggdf_Func[ggdf_Func$cluster %nin% c("UA", "ThCTL"),]

##substitute CD137 for 41BB##
#Replace 41BB for CD137 in DF and Make New Column
ggdf_Func1 <- gsub('X41BB', "CD137", ggdf_Func$variable)
#Replace old olumn with new colummn
ggdf_Func$variable <- ggdf_Func1

##substitute HLADR for HLA-DR##
ggdf_Func1 <- gsub('HLA.DR', "HLADR", ggdf_Func$variable)
ggdf_Func$variable <- ggdf_Func1

##substitute RORY for RORYT##
ggdf_Func1 <- gsub('RORYT', "RORY", ggdf_Func$variable)
ggdf_Func$variable <- ggdf_Func1


ggdf_Func_Baseline<-ggdf_Func[ggdf_Func$timepoint %nin% c("1", "2"),]
ggdf_Func_OnTreatment<-ggdf_Func[ggdf_Func$timepoint %nin% c("0", "2"),]

ggdf3 <- ggdf_Func[ggdf_Func$sample_id %nin% c("523_2","404_2","718_2","890_2","065_2","617_2","570_2","409_2","388_2","636_2","439_2","758_2","463_2","543_2","505_2","464_2","542_2","499_2"),]
ggdf_Func_T0 <- ggdf3[ggdf3$timepoint %nin% c("1","2"),]
ggdf_Func_T1 <- ggdf3[ggdf3$timepoint %nin% c("0"),]
ggdf_Func_T0$post_value = ggdf_Func_T1$value
ggdf_Func_T0$fc = ggdf_Func_T0$post_value/ggdf_Func_T0$value
ggdf_Func_FC <- ggdf_Func_T0

pdf("Baseline_funcmarkers_Age.pdf",width=20,height=12)
for(i in 1:length(fmlistplot)){
  ggp <- ggplot(ggdf_Func_Baseline[ggdf_Func_Baseline$variable==fmlistplot[i],], aes(x=Aged, y=value, fill=Aged))+
    geom_boxplot(outlier.size=0, lwd=0.25, color="black")+
    geom_jitter(width=0.15, size=0.5)+
    facet_wrap(~cluster,ncol=8,scales="free")+
    ylab("MMI")+
    stat_compare_means(method="wilcox.test", 
                       show.legend=FALSE, 
                       method.args = list(var.equal = TRUE), 
                       #label="p.signif",
                       comparisons = age_comps,
                       label.x.npc=0.37, 
                       size=2.5)+
    ggtitle(fmlistplot[i])+
    theme(axis.text.x = element_blank(),
          axis.text.y = element_text(size=10, color="black"),
          axis.title.x = element_blank(),
          axis.line.x = element_line(size=0.25, color="black"),
          axis.line.y = element_line(size=0.25, color="black"),
          axis.ticks.x = element_blank(),
          axis.title.y = element_text(size=12, color="black"),
          strip.background = element_rect(fill=NA),
          strip.text = element_text(size=7, color="black"),
          panel.background = element_rect(fill="white"),
          legend.title = element_blank(),
          legend.key.size = unit(2,'lines'),
          legend.text = element_text(size=8),
          legend.key = element_rect(fill="white")  
    )
  print(ggp)
}
dev.off()
write.csv(ggdf_Func_Baseline, "Func_Markers_Baseline_Aged.csv")

pdf("OnTreatment_funcmarkers_Age.pdf",width=20,height=12)
for(i in 1:length(fmlistplot)){
  ggp <- ggplot(ggdf_Func_OnTreatment[ggdf_Func_OnTreatment$variable==fmlistplot[i],], aes(x=Aged, y=value, fill=Aged))+
    geom_boxplot(outlier.size=0, lwd=0.25, color="black")+
    geom_jitter(width=0.15, size=0.5)+
    facet_wrap(~cluster,ncol=8,scales="free")+
    ylab("MMI")+
    stat_compare_means(method="wilcox.test", 
                       show.legend=FALSE, 
                       method.args = list(var.equal = TRUE), 
                       #label="p.signif",
                       comparisons = age_comps,
                       label.x.npc=0.37, 
                       size=2.5)+
    ggtitle(fmlistplot[i])+
    theme(axis.text.x = element_blank(),
          axis.text.y = element_text(size=10, color="black"),
          axis.title.x = element_blank(),
          axis.line.x = element_line(size=0.25, color="black"),
          axis.line.y = element_line(size=0.25, color="black"),
          axis.ticks.x = element_blank(),
          axis.title.y = element_text(size=12, color="black"),
          strip.background = element_rect(fill=NA),
          strip.text = element_text(size=7, color="black"),
          panel.background = element_rect(fill="white"),
          legend.title = element_blank(),
          legend.key.size = unit(2,'lines'),
          legend.text = element_text(size=8),
          legend.key = element_rect(fill="white")  
    )
  print(ggp)
}
dev.off()
write.csv(ggdf_Func_OnTreatment, "Func_Markers_OnTreatment_Aged.csv")


pdf("FC_funcmarkers_Age.pdf",width=20,height=12)
for(i in 1:length(fmlistplot)){
  ggp <- ggplot(ggdf_Func_FC[ggdf_Func_FC$variable==fmlistplot[i],], aes(x=Aged, y=fc, fill=Aged))+
    geom_boxplot(outlier.size=0, lwd=0.25, color="black")+
    geom_jitter(width=0.15, size=0.5)+
    facet_wrap(~cluster,ncol=8,scales="free")+
    ylab("Fold Change MMI")+
    stat_compare_means(method="wilcox.test", 
                       show.legend=FALSE, 
                       method.args = list(var.equal = TRUE), 
                       #label="p.signif",
                       comparisons = age_comps,
                       label.x.npc=0.37, 
                       size=2.5)+
    ggtitle(fmlistplot[i])+
    theme(axis.text.x = element_blank(),
          axis.text.y = element_text(size=10, color="black"),
          axis.title.x = element_blank(),
          axis.line.x = element_line(size=0.25, color="black"),
          axis.line.y = element_line(size=0.25, color="black"),
          axis.ticks.x = element_blank(),
          axis.title.y = element_text(size=12, color="black"),
          strip.background = element_rect(fill=NA),
          strip.text = element_text(size=7, color="black"),
          panel.background = element_rect(fill="white"),
          legend.title = element_blank(),
          legend.key.size = unit(2,'lines'),
          legend.text = element_text(size=8),
          legend.key = element_rect(fill="white")  
    )
  print(ggp)
}
dev.off()
write.csv(ggdf_Func_FC, "Func_Markers_Baseline_Aged.csv")

#### Functional Markers Aged and Young responder vs non-responders####

ggdf_Func_Baseline<-ggdf_Func[ggdf_Func$timepoint %nin% c("1", "2"),]

ggdf_Func_Baseline_Aged <- ggdf_Func_Baseline[ggdf_Func_Baseline$Aged %nin% c("No"),]
ggdf_Func_Baseline_Aged <- ggdf_Func_Baseline_Aged[ggdf_Func_Baseline_Aged$response %nin% c("NA"),]

ggdf_Func_Baseline_Young <- ggdf_Func_Baseline[ggdf_Func_Baseline$Aged %nin% c("Yes"),]
ggdf_Func_Baseline_Young <- ggdf_Func_Baseline_Young[ggdf_Func_Baseline_Young$response %nin% c("NA"),]

ggdf_Func_OnTreatment<-ggdf_Func[ggdf_Func$timepoint %nin% c("0", "2"),]

ggdf_Func_OnTreatment_Aged <- ggdf_Func_OnTreatment[ggdf_Func_OnTreatment$Aged %nin% c("No"),]
ggdf_Func_OnTreatment_Aged <- ggdf_Func_OnTreatment_Aged[ggdf_Func_OnTreatment_Aged$response %nin% c("NA"),]

ggdf_Func_OnTreatment_Young <- ggdf_Func_OnTreatment[ggdf_Func_OnTreatment$Aged %nin% c("Yes"),]
ggdf_Func_OnTreatment_Young <- ggdf_Func_OnTreatment_Young[ggdf_Func_OnTreatment_Young$response %nin% c("NA"),]

ggdf_Func_FC <- ggdf_Func_T0

ggdf_Func_FC_Aged <- ggdf_Func_FC[ggdf_Func_FC$Aged %nin% c("No"),]
ggdf_Func_FC_Aged <- ggdf_Func_FC_Aged[ggdf_Func_FC_Aged$response %nin% c("NA"),]

ggdf_Func_FC_Young <- ggdf_Func_FC[ggdf_Func_FC$Aged %nin% c("Yes"),]
ggdf_Func_FC_Young <- ggdf_Func_FC_Young[ggdf_Func_FC_Young$response %nin% c("NA"),]

pdf("Baseline_funcmarkers_Aged_Response.pdf",width=20,height=12)
for(i in 1:length(fmlistplot)){
  ggp <- ggplot(ggdf_Func_Baseline_Aged[ggdf_Func_Baseline_Aged$variable==fmlistplot[i],], aes(x=response, y=value, fill=response))+
    geom_boxplot(outlier.size=0, lwd=0.25, color="black")+
    geom_jitter(width=0.15, size=0.5)+
    facet_wrap(~cluster,ncol=8,scales="free")+
    ylab("MMI")+
    stat_compare_means(method="wilcox.test", 
                       show.legend=FALSE, 
                       method.args = list(var.equal = TRUE), 
                       label="p.format",
                       comparisons = response_comps,
                       label.x.npc=0.37, 
                       size=2.5)+
    ggtitle(fmlistplot[i])+
    theme(axis.text.x = element_blank(),
          axis.text.y = element_text(size=10, color="black"),
          axis.title.x = element_blank(),
          axis.line.x = element_line(size=0.25, color="black"),
          axis.line.y = element_line(size=0.25, color="black"),
          axis.ticks.x = element_blank(),
          axis.title.y = element_text(size=12, color="black"),
          strip.background = element_rect(fill=NA),
          strip.text = element_text(size=7, color="black"),
          panel.background = element_rect(fill="white"),
          legend.title = element_blank(),
          legend.key.size = unit(2,'lines'),
          legend.text = element_text(size=8),
          legend.key = element_rect(fill="white")  
    )
  print(ggp)
}
dev.off()
write.csv(ggdf_Func_Baseline_Aged, "Aged_Func_Markers_Baseline_Response.csv")

pdf("Baseline_funcmarkers_Young_Response.pdf",width=20,height=12)
for(i in 1:length(fmlistplot)){
  ggp <- ggplot(ggdf_Func_Baseline_Young[ggdf_Func_Baseline_Young$variable==fmlistplot[i],], aes(x=response, y=value, fill=response))+
    geom_boxplot(outlier.size=0, lwd=0.25, color="black")+
    geom_jitter(width=0.15, size=0.5)+
    facet_wrap(~cluster,ncol=8,scales="free")+
    ylab("MMI")+
    stat_compare_means(method="wilcox.test", 
                       show.legend=FALSE, 
                       method.args = list(var.equal = TRUE), 
                       label="p.format",
                       comparisons = response_comps,
                       label.x.npc=0.37, 
                       size=2.5)+
    ggtitle(fmlistplot[i])+
    theme(axis.text.x = element_blank(),
          axis.text.y = element_text(size=10, color="black"),
          axis.title.x = element_blank(),
          axis.line.x = element_line(size=0.25, color="black"),
          axis.line.y = element_line(size=0.25, color="black"),
          axis.ticks.x = element_blank(),
          axis.title.y = element_text(size=12, color="black"),
          strip.background = element_rect(fill=NA),
          strip.text = element_text(size=7, color="black"),
          panel.background = element_rect(fill="white"),
          legend.title = element_blank(),
          legend.key.size = unit(2,'lines'),
          legend.text = element_text(size=8),
          legend.key = element_rect(fill="white")  
    )
  print(ggp)
}
dev.off()
write.csv(ggdf_Func_Baseline_Young, "Young_Func_Markers_Baseline_Response.csv")


pdf("OnTreatment_funcmarkers_Aged_Response.pdf",width=20,height=12)
for(i in 1:length(fmlistplot)){
  ggp <- ggplot(ggdf_Func_OnTreatment_Aged[ggdf_Func_OnTreatment_Aged$variable==fmlistplot[i],], aes(x=response, y=value, fill=response))+
    geom_boxplot(outlier.size=0, lwd=0.25, color="black")+
    geom_jitter(width=0.15, size=0.5)+
    facet_wrap(~cluster,ncol=8,scales="free")+
    ylab("MMI")+
    stat_compare_means(method="wilcox.test", 
                       show.legend=FALSE, 
                       method.args = list(var.equal = TRUE), 
                       label="p.format",
                       comparisons = response_comps,
                       label.x.npc=0.37, 
                       size=2.5)+
    ggtitle(fmlistplot[i])+
    theme(axis.text.x = element_blank(),
          axis.text.y = element_text(size=10, color="black"),
          axis.title.x = element_blank(),
          axis.line.x = element_line(size=0.25, color="black"),
          axis.line.y = element_line(size=0.25, color="black"),
          axis.ticks.x = element_blank(),
          axis.title.y = element_text(size=12, color="black"),
          strip.background = element_rect(fill=NA),
          strip.text = element_text(size=7, color="black"),
          panel.background = element_rect(fill="white"),
          legend.title = element_blank(),
          legend.key.size = unit(2,'lines'),
          legend.text = element_text(size=8),
          legend.key = element_rect(fill="white")  
    )
  print(ggp)
}
dev.off()
write.csv(ggdf_Func_OnTreatment_Aged, "Aged_Func_Markers_OnTreatment_Response.csv")

pdf("OnTreatment_funcmarkers_Young_Response.pdf",width=20,height=12)
for(i in 1:length(fmlistplot)){
  ggp <- ggplot(ggdf_Func_OnTreatment_Young[ggdf_Func_OnTreatment_Young$variable==fmlistplot[i],], aes(x=response, y=value, fill=response))+
    geom_boxplot(outlier.size=0, lwd=0.25, color="black")+
    geom_jitter(width=0.15, size=0.5)+
    facet_wrap(~cluster,ncol=8,scales="free")+
    ylab("MMI")+
    stat_compare_means(method="wilcox.test", 
                       show.legend=FALSE, 
                       method.args = list(var.equal = TRUE), 
                       label="p.format",
                       comparisons = response_comps,
                       label.x.npc=0.37, 
                       size=2.5)+
    ggtitle(fmlistplot[i])+
    theme(axis.text.x = element_blank(),
          axis.text.y = element_text(size=10, color="black"),
          axis.title.x = element_blank(),
          axis.line.x = element_line(size=0.25, color="black"),
          axis.line.y = element_line(size=0.25, color="black"),
          axis.ticks.x = element_blank(),
          axis.title.y = element_text(size=12, color="black"),
          strip.background = element_rect(fill=NA),
          strip.text = element_text(size=7, color="black"),
          panel.background = element_rect(fill="white"),
          legend.title = element_blank(),
          legend.key.size = unit(2,'lines'),
          legend.text = element_text(size=8),
          legend.key = element_rect(fill="white")  
    )
  print(ggp)
}
dev.off()
write.csv(ggdf_Func_OnTreatment_Young, "Young_Func_Markers_OnTreatment_Response.csv")


pdf("FC_funcmarkers_Aged_Response.pdf",width=20,height=12)
for(i in 1:length(fmlistplot)){
  ggp <- ggplot(ggdf_Func_FC_Aged[ggdf_Func_FC_Aged$variable==fmlistplot[i],], aes(x=response, y=fc, fill=response))+
    geom_boxplot(outlier.size=0, lwd=0.25, color="black")+
    geom_jitter(width=0.15, size=0.5)+
    facet_wrap(~cluster,ncol=8,scales="free")+
    ylab("Fold Change MMI")+
    stat_compare_means(method="wilcox.test", 
                       show.legend=FALSE, 
                       method.args = list(var.equal = TRUE), 
                       label="p.format",
                       comparisons = response_comps,
                       label.x.npc=0.37, 
                       size=2.5)+
    ggtitle(fmlistplot[i])+
    theme(axis.text.x = element_blank(),
          axis.text.y = element_text(size=10, color="black"),
          axis.title.x = element_blank(),
          axis.line.x = element_line(size=0.25, color="black"),
          axis.line.y = element_line(size=0.25, color="black"),
          axis.ticks.x = element_blank(),
          axis.title.y = element_text(size=12, color="black"),
          strip.background = element_rect(fill=NA),
          strip.text = element_text(size=7, color="black"),
          panel.background = element_rect(fill="white"),
          legend.title = element_blank(),
          legend.key.size = unit(2,'lines'),
          legend.text = element_text(size=8),
          legend.key = element_rect(fill="white")  
    )
  print(ggp)
}
dev.off()
write.csv(ggdf_Func_FC_Aged, "Aged_Func_Markers_FC_Response.csv")


pdf("FC_funcmarkers_Young_Response.pdf",width=20,height=12)
for(i in 1:length(fmlistplot)){
  ggp <- ggplot(ggdf_Func_FC_Young[ggdf_Func_FC_Young$variable==fmlistplot[i],], aes(x=response, y=fc, fill=response))+
    geom_boxplot(outlier.size=0, lwd=0.25, color="black")+
    geom_jitter(width=0.15, size=0.5)+
    facet_wrap(~cluster,ncol=8,scales="free")+
    ylab("Fold Change MMI")+
    stat_compare_means(method="wilcox.test", 
                       show.legend=FALSE, 
                       method.args = list(var.equal = TRUE), 
                       label="p.format",
                       comparisons = response_comps,
                       label.x.npc=0.37, 
                       size=2.5)+
    ggtitle(fmlistplot[i])+
    theme(axis.text.x = element_blank(),
          axis.text.y = element_text(size=10, color="black"),
          axis.title.x = element_blank(),
          axis.line.x = element_line(size=0.25, color="black"),
          axis.line.y = element_line(size=0.25, color="black"),
          axis.ticks.x = element_blank(),
          axis.title.y = element_text(size=12, color="black"),
          strip.background = element_rect(fill=NA),
          strip.text = element_text(size=7, color="black"),
          panel.background = element_rect(fill="white"),
          legend.title = element_blank(),
          legend.key.size = unit(2,'lines'),
          legend.text = element_text(size=8),
          legend.key = element_rect(fill="white")  
    )
  print(ggp)
}
dev.off()
write.csv(ggdf_Func_FC_Young, "Young_Func_Markers_FC_Response.csv")

####Umap####
umapRes <- readRDS("backup_umap.rds")

mm <- match(as.character(umapRes$sample_id), as.character(output[["meta_data"]]$sample_id))
umapRes$timepoint <- factor(output[["meta_data"]]$timepoint[mm], levels=timelevels)
umapRes$batch <- output[["meta_data"]]$batch[mm]
umapRes$sample_id <- factor(output[["meta_data"]]$sample_id[mm], levels=samplevels)
umapRes$cell_clustering = factor(umapRes$cell_clustering, levels=clusterlevels)
pdf('Aging_umaps.pdf',width=10,height=10)
plotUmap(umapRes = umapRes,
         code_clustering=cell_clustering1m,
         color_clusters = "auto",
         subtype_markers = output$subtype_markers)
dev.off()
