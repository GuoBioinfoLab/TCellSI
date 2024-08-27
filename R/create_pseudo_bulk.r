#' @title TCSS_Calculate
#' @description You can use this function to calculate scores for T cell states.
#'
#' @param object Sample expression
#' @param reference Reference inputs
#' @param nbin How many gene bins
#' @param ctrl How many genes to select from each bins
#' @param seed Random number seed set
#' @return features.scores.df
#' @import ggplot2
#' @import dplyr
#' @importFrom stats na.omit rnorm
#' @export
#' @examples
#'
#' # Example usage:
#' \donttest{
#' sample_expression <- TCellSI::exampleSample
#' TCSS_Calculate(sample_expression)
#' }

create_pseudo_bulk <- function(annotation_data, expression_data, cluster_col, cell_id_col, n_clusters = 18, factor = 5, sampling_rate = 0.6) {
  # 计算伪散装数据列表
  pseudo_bulk_list <- lapply(1:n_clusters, function(i) {
    pseudo_bulk_df <- as.data.frame(lapply(1:round(table(annotation_data[[cluster_col]])[i] / factor), function(j) {
      sp <- sample(annotation_data[which(annotation_data[[cluster_col]] == names(table(annotation_data[[cluster_col]]))[i]), ][[cell_id_col]], round(table(annotation_data[[cluster_col]])[i] / factor) * sampling_rate)
      pseudo_bulk <- as.data.frame(apply(expression_data[, which(colnames(expression_data) %in% sp)], 1, mean))
      return(pseudo_bulk)
    }))
    return(pseudo_bulk_df)
  })
  pseudo_bulk <- as.data.frame(matrix(nrow = nrow(expression_data), ncol = 0))
  for (i in 1:n_clusters) {
    colnames(pseudo_bulk_list[[i]]) <- rep(paste0(names(table(annotation_data[[cluster_col]]))[i], "_bulk"), round(table(annotation_data[[cluster_col]])[i] / factor))
    pseudo_bulk <- cbind(pseudo_bulk, pseudo_bulk_list[[i]])
  }
  return(pseudo_bulk)
}
