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

TCSS_Calculate <- function(object, reference = ref_data, nbin = 50, ctrl = 100, seed = 1, ref = TRUE) {
  set.seed(seed = seed)
  all_sample_HK <- mean(rowMeans(object[which(rownames(object) %in% HKgenes), ]))
  all_sample_TPM1 <- apply(object, 2, function(i) {
    return(i * all_sample_HK / mean(na.omit(i[HKgenes])))
  })
  helper_message_shown <- FALSE
  # 计算特征得分列表
  scores_list <- lapply(seq_len(length(markers)), function(k) {
    if (k %in% c(1, 2, 3, 7, 8, 9, 10)) {
      message(paste0("Calculating of ", names(markers)[k]), " state")
    } else if (k == 4 && !helper_message_shown) {
      message("Calculating of Helper state")
      helper_message_shown <- TRUE
    }
    features <- markers[[k]]
    # 检查缺失的特征
    missing_features_input_sample <- setdiff(features, rownames(all_sample_TPM1))
    if (length(missing_features_input_sample) > 0) {
      warning("The following features are not present in the object: ",
              paste(missing_features_input_sample, collapse = ", "),
              call. = FALSE, immediate. = TRUE)
    }
    features <- intersect(features, rownames(all_sample_TPM1))
    if (ref) {
      missing_features_ref_data <- setdiff(features, rownames(reference))
      if (length(missing_features_ref_data) > 0) {
        warning("The following features are not present in the reference: ",
                paste(missing_features_ref_data, collapse = ", "),
                call. = FALSE, immediate. = TRUE)
      }
      features <- intersect(features, rownames(reference))
      if (length(features) < 3) {
        stop("The number of features is less than 3")
      }
    }
    all_sample_TPM2 <- all_sample_TPM1[which(!rownames(all_sample_TPM1) %in% features), ]
    features.scores.vec <- sapply(seq_len(ncol(all_sample_TPM2)), function(i) {
      features.exp <- all_sample_TPM1[features, i]
      data.avg <- all_sample_TPM2[, i]
      data.avg <- sort(data.avg)
      data.cut <- ggplot2::cut_number(data.avg + rnorm(n = length(data.avg)) / 1e+6, n = nbin, labels = FALSE)
      names(data.cut) <- data.avg
      ctrl.use <- lapply(seq_len(nbin), function(j) {
        as.numeric(names(sample(data.cut[which(data.cut == j)], ctrl, replace = FALSE)))
      })
      sample_feature <- c(features.exp, unlist(ctrl.use))
      features.use <- rank(sample_feature,ties.method = "first") / length(sample_feature)
      feature.score <- features.use[features]
      ctrl.score <- sort(features.use)[1:length(features.exp)]
      return(mean(feature.score) - mean(ctrl.score))
    })
    # 如果 ref = TRUE，计算与参考数据的比例
    if (ref) {
      features_sample <- all_sample_TPM1[features, ]
      features_ref <- reference[features, names(markers)[k]]
      div_percent <- features_sample / features_ref
      flag_lower <- (features_sample - features_ref) < 0
      number_lower <- colSums(flag_lower)
      percentage <- sapply(seq_len(ncol(div_percent)), function(i) {
      ifelse(number_lower[i] == 0, 1,
             (sum(div_percent[, i][flag_lower[, i]]) + sum(flag_lower[, i] == FALSE)) / length(features))
      })
      features.scores.vec <- features.scores.vec * percentage
    }
    if (length(features.scores.vec) > 0) {
      return(as.data.frame(t(features.scores.vec)))
    } else {
      return(NULL)
    }
  })
  # 过滤掉 NULL 结果并绑定数据框
  features.scores.df <- dplyr::bind_rows(scores_list)
  colnames(features.scores.df) <- colnames(object)
  rownames(features.scores.df) <- names(markers)

  # 如果 ref = TRUE，调整特定 marker 的得分
  if (ref) {
    for (i in c(2, 4:6)) {
      for (j in 1:ncol(all_sample_TPM1)) {
        pq <- as.numeric(reference["CD4", names(markers)[i]])
        percentage_T <- ifelse(all_sample_TPM1["CD4", j] < pq,
                               all_sample_TPM1["CD4", j] / pq, 1)
        features.scores.df[i, j] <- features.scores.df[i, j] * percentage_T
      }
    }
  }
  vec <- sapply(seq_len(ncol(features.scores.df)), function(i) {
    max(features.scores.df[4, i], features.scores.df[5, i], features.scores.df[6, i])
  })
  features.scores.df <- rbind(features.scores.df, vec)
  rownames(features.scores.df)[11] <- "Helper"
  features.scores.df <- features.scores.df[-4:-6, ]
  features.scores.df <- features.scores.df[c("Quiescence", "Regulating", "Proliferation",
                                             "Helper", "Cytotoxicity",
                                             "Progenitor_exhaustion",
                                             "Terminal_exhaustion", "Senescence"), ]

  return(features.scores.df)
}
