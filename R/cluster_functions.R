firstKneeValue <- function(data=spliced, lower_num=2000){
  br.out <- DropletUtils::barcodeRanks(data)
  plot(br.out$rank, br.out$total, log="xy", xlab="Rank", ylab="Total")
  o <- order(br.out$rank)
  first_knee <- DropletUtils::barcodeRanks(data, lower=lower_num, exclude.from = 50)
  knee <- S4Vectors::metadata(first_knee)$knee
  plot(br.out$rank, br.out$total, log="xy", xlab="Rank", ylab="Total")
  o <- order(br.out$rank)
  abline(h=knee, col="dodgerblue", lty=2)
  lines(br.out$rank[o], br.out$fitted[o], col="red")
  return(first_knee)
}

findKnees <- function(data2=spliced, first_knee){

  all_knees <- NULL
  knee <- S4Vectors::metadata(first_knee)$knee
  for (i in 1:3)
  {
    all_knees[[length(all_knees) + 1]] <- knee
    br <- DropletUtils::barcodeRanks(data2, lower=0, exclude.from = knee)
    knee <- S4Vectors::metadata(br)$knee
    if (knee < 10) {
      all_knees[[length(all_knees) + 1]] <- knee
      break
    }
  }

  all_knee_df <- do.call(rbind.data.frame, all_knees)
  colnames(all_knee_df)[1] = "knee_val"

  sorted_knees <- as.data.frame(sort(all_knee_df$knee_val, decreasing = TRUE))

  if (nrow(sorted_knees) == 3) {
    if (sorted_knees[1,] - sorted_knees[2,] <= 50) {
      knee1 = sorted_knees[1,]
      sorted_knees[2,] <- S4Vectors::metadata(first_knee)$inflection - (S4Vectors::metadata(first_knee)$inflection * 0.5)
      knee2 = sorted_knees[2,]
      knee3 = NA
      warning('Second knee unsuccessfully found, first knee inflection point used instead.\nInflection value: ', knee2)
    }
    else if (sorted_knees[2,] - sorted_knees[3,] >= 30) {
      knee1 = sorted_knees[1,]
      knee2 = sorted_knees[2,]
      knee3 = sorted_knees[3,]
    }
    else{
      knee1 = sorted_knees[1,]
      knee2 = sorted_knees[2,]
      knee3 = NA
    }
  }
  if (nrow(sorted_knees) == 2) {
    knee1 = sorted_knees[1,]
    knee2 = sorted_knees[2,]
    knee3 = NA
  }

  if (nrow(sorted_knees) == 1) {
    knee1 = sorted_knees[1,]
    knee2 = NA
    knee3 = NA
  }

  if (nrow(sorted_knees > 1)){
    if (knee1 < knee2 ) {
    stop("knee 1 and knee 2 not correct")
  }
  if (!is.na(knee3)) {
    if (knee2 < knee3) {
      stop("knee 2 and knee 3 not correct")
    }
    if (knee1 < knee3 ) {
      stop("knee 1 and knee 3 not correct")
    }
  }
    }
  knee_list <- I(list(knee1, knee2, knee3))
  knee_df <- as.data.frame(knee_list)

  rownames(knee_df) <- c("knee1","knee2","knee3")
  return(knee_df)
}

plotKnees <- function(first_knee, knee_data){
  plot(first_knee$rank, first_knee$total, log="xy", xlab="Rank", ylab="Total")
  abline(h=knee_data[1,], col="dodgerblue", lty=2)
  if (!is.na(knee_data[2,])) {
    abline(h=knee_data[2,], col="dodgerblue", lty=2)
  }
  if (!is.na(knee_data[3,])) {
    abline(h=knee_data[3,], col="dodgerblue", lty=2)
  }
}

filterKnees <- function(data3=spliced, knee_data){
  #filter if between k1k2 or k2k3
  br1 <- DropletUtils::barcodeRanks(
    data3,
    lower = 100,
    fit.bounds = NULL,
    exclude.from = 50,
    df = 10,
    BPPARAM = BiocParallel::SerialParam()
  )

  #all barcodes between first and last knee
  if (is.na(knee_data[2,])) {
    tot <- DropletUtils::emptyDrops(data3, retain = knee_data[1,])
  }else if (is.na(knee_data[3,])) {
    tot <- DropletUtils::emptyDrops(data3, lower = knee_data[2,], retain = knee_data[1,])
  }else {
    tot <- DropletUtils::emptyDrops(data3, lower = knee_data[3,], retain = knee_data[1,])
  }

  #adding col name to row names
  e.out <- cbind(cellID = rownames(tot), tot)
  br.out <- cbind(cellID = rownames(br1), br1)

  is.cell <- e.out$FDR <= 0.001

  barcodes <- as.data.frame(colnames(data3))
  colnames(barcodes) <- c('barcodes')
  is.higherknee1 <- (br.out$total > knee_data[1,])

  if (!is.na(knee_data[3,])) {
    is.betweenknee1andknee2 <- ((br.out$total > knee_data[2,]) & (br.out$total <= knee_data[1,]))
    is.betweenknee2andknee3 <- ((br.out$total >= knee_data[3,]) & (br.out$total <= knee_data[2,]))
    knee_groups_df <- cbind(barcodes, is.betweenknee1andknee2, is.betweenknee2andknee3, is.higherknee1)
  }
  else if (!is.na(knee_data[2,])) {
    is.betweenknee1andknee2 <- ((br.out$total > knee_data[2,]) & (br.out$total <= knee_data[1,]))
    knee_groups_df <- cbind(barcodes, is.betweenknee1andknee2, is.higherknee1)
  }
   else{
    knee_groups_df <- cbind(barcodes,is.higherknee1)
  }
  return(knee_groups_df)
}

combineSplicedUnsplicedData <- function(spliced_data=spliced, unspliced_data=unspliced, log_transform = TRUE){
  total_spliced <- Matrix::colSums(spliced_data)
  total_unspliced <- Matrix::colSums(unspliced_data)

  if (log_transform == TRUE) {
    log_total_spliced <- log10(total_spliced+1)
    log_total_unspliced <- log10(total_unspliced+1)
    comb <-cbind(log_total_spliced, log_total_unspliced)
    ds_comb <- as.data.frame(comb)
  }
  else{
    comb <-cbind(total_spliced, total_unspliced)
    ds_comb <- as.data.frame(comb)
  }
  return(ds_comb)
}

filterUMIData <- function(ds_comb, filter_num = 1){
  filt_comb_data_1 <- dplyr::filter(ds_comb, ds_comb[,1] > filter_num & ds_comb[,2] > filter_num)
  plot(filt_comb_data_1[,1], filt_comb_data_1[,2], xlab = "Total Spliced", ylab = 'Total Unspliced')
  filt_comb_data_1[is.na(filt_comb_data_1)] <- 0
  return(filt_comb_data_1)
}

clusterGroups <- function(filter_data, sort_knee, clusters = NULL ){

  if (is.null(clusters)) {
    if (!is.na(sort_knee[3,])) {
      clusters = 3
  }else{clusters = 2}
  }
  else{
    km <- kmeans(filter_data, centers=clusters, nstart=45)
  }

  km <- kmeans(filter_data, centers=clusters, nstart=45)
  print(factoextra::fviz_cluster(km, data = filter_data, stand = FALSE, geom = "point", show.clust.cent=TRUE, xlab = "Total Spliced", ylab = "Total Unspliced"))
  filt_data <- as.data.frame(filter_data)
  filt_data$cluster <- km$cluster

  cluster_data <- list(km, filt_data)
  return(cluster_data)
}

kneeGroupsPlot <- function(ds_comb, knee_groups){
  big_data <- cbind(ds_comb, knee_groups)
  big_data$col <- "grey"
  big_data[big_data$is.betweenknee1andknee2, "col"] <- "green"
  big_data[big_data$is.betweenknee2andknee3, "col"] <- "red"
  big_data[big_data$is.higherknee1, "col"] <- "orange"

  knee_plot <- ggplot2::ggplot(big_data, ggplot2::aes(x = log_total_spliced, y = log_total_unspliced, color = col)) +
    ggplot2::geom_point() +
    ggplot2::ggtitle("Knee group plot") + ggplot2::xlab("Total Spliced") + ggplot2::ylab("Total Unspliced") +
    ggplot2::scale_color_identity(name = "Legend",
                         labels = c(orange = "Above Knee 1", green = "Between Knee 1\nand Knee 2", red = "Between Knee 2\nand Knee 3",
                                    grey = "Other"), guide = "legend")
  print(knee_plot)
}

compareKneeClusterPlot <- function(filt_comb_data_2, cluster, knee_plot) {
  km <- cluster[[1]]
  cluster_plot <- factoextra::fviz_cluster(km, data = filt_comb_data_2, stand = FALSE, geom = "point", pointsize = 1, show.clust.cent=TRUE, xlab = "Total Spliced", ylab = "Total Unspliced")

  print(gridExtra::grid.arrange(cluster_plot, knee_plot, ncol=2, nrow = 1))
}

barcodeInFirstCluster <- function(kmeans_data, cluster_data){
  centroid <- kmeans_data$centers
  if (is.null(centroid)) {
    centroid <- kmeans_data[[1]]$centers
  }
  sorted_clust <- centroid[order(centroid[,1],decreasing = TRUE),]
  sorted_clust.out <- as.data.frame(cbind(cluster = rownames(sorted_clust), sorted_clust))
  order_cluster <- as.data.frame(sorted_clust.out[,1])
  firstcluster <- as.numeric(order_cluster[1,])
  cluster_data_out <- as.data.frame(cluster_data)
  first_clust <- dplyr::filter(cluster_data_out, cluster == firstcluster)
  bars_in_first_clust <- as.data.frame(rownames(first_clust))
  colnames(bars_in_first_clust) <- c('barcodes')
  return(bars_in_first_clust)
}

meanExpressionSecondClust <- function(kmeans_data1, cluster_data1){

  centroid <- kmeans_data1$centers
  if (is.null(centroid)) {
    centroid <- kmeans_data1[[1]]$centers
  }
  sorted_clust <- centroid[order(centroid[,1],decreasing = TRUE),]
  sorted_clust.out <- as.data.frame(cbind(cluster = rownames(sorted_clust), sorted_clust))
  order_cluster <- as.data.frame(sorted_clust.out[,1])
  secondcluster <- as.numeric(order_cluster[2,])
  cluster_data_out <- as.data.frame(cluster_data1)
  bars_in_second_clust <-  dplyr::filter(cluster_data_out, cluster == secondcluster)
  named_bars_second_clust <- as.data.frame(cbind(barcode = rownames(bars_in_second_clust), bars_in_second_clust))

  spliced_sec_clust <- spliced[,named_bars_second_clust$barcode]
  unspliced_sec_clust <- unspliced[,named_bars_second_clust$barcode]

  spliced_gene_count <- Matrix::rowSums(spliced_sec_clust)
  unspliced_gene_count <- Matrix::rowSums(unspliced_sec_clust)
  plot(log10(spliced_gene_count+1),log10(unspliced_gene_count)+1, main = 'Second Cluster Gene Expression', xlab = 'Spliced', ylab = 'Unspliced')
  sec_clust_gene_exp <- as.data.frame(cbind(spliced_gene_count,unspliced_gene_count))
  return(sec_clust_gene_exp)
}

















