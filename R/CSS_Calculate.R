#' @title CSS_Calculate
#' @description You can use this function to calculate scores for other cell states that interest you.
#'
#' @param object Sample expression
#' @param reference Reference inputs
#' @param markers A list of multiple cell states containing specific gene sets
#' @param nbin How many gene bins
#' @param ctrl How many genes to select from each bins
#' @param seed Random number seed set
#' @return features.scores.df
#' @import ggplot2
#' @import dplyr
#' @importFrom stats na.omit rnorm
#' @export

CSS_Calculate <- function(object, ref = TRUE, reference = NULL, markers, nbin = 50, ctrl = 100, seed = 1) {
  set.seed(seed = seed)
  all_sample_HK <- mean(rowMeans(object[which(rownames(object) %in% HKgenes), ]))
  all_sample_TPM1 <- apply(object, 2, function(i) {
    return(i * all_sample_HK / mean(na.omit(i[HKgenes])))
  })

  features.scores.df <- as.data.frame(dplyr::bind_rows(lapply(seq_len(length(markers)), function(k) {
    message(paste0("Calculating ", names(markers)[k]), " state : ( ", k, "/", length(markers), " )")
    features <- markers[[k]]
    missing_features_input_sample <- setdiff(features, rownames(all_sample_TPM1))
    if (length(missing_features_input_sample) > 0) {
      warning("The following features are not present in the object: ",
        paste(missing_features_input_sample, collapse = ", "), ifelse(test = FALSE,
          yes = ", attempting to find updated synonyms",
          no = ", not searching for symbol synonyms"
        ),
        call. = FALSE, immediate. = TRUE
      )
    }
    features <- intersect(features, rownames(all_sample_TPM1))
    
    if (ref) {
      # Perform reference-related calculations only if ref = TRUE
      missing_features_ref_data <- setdiff(features, rownames(reference))
      if (length(x = missing_features_ref_data) > 0) {
        warning("The following features are not present in the reference: ",
          paste(missing_features_ref_data, collapse = ", "), ifelse(test = FALSE,
            yes = ", attempting to find updated synonyms",
            no = ", not searching for symbol synonyms"
          ),
          call. = FALSE, immediate. = TRUE
        )
      }
      features <- intersect(features, rownames(reference))
    }
    
    if (length(features) < 3) {
      stop("The number of features is less than 3")
    }
    
    all_sample_TPM2 <- all_sample_TPM1[which(!rownames(all_sample_TPM1) %in% features), ]
    features.scores.vec <- pbapply::pbsapply(seq_len(ncol(all_sample_TPM2)), function(i) {
      features.exp <- all_sample_TPM1[features, i]
      data.avg <- all_sample_TPM2[, i]
      data.avg <- sort(data.avg)
      data.cut <- ggplot2::cut_number(data.avg + rnorm(n = length(data.avg)) / 1e+6, n = nbin, labels = FALSE, right = FALSE)
      names(data.cut) <- data.avg
      ctrl.use <- lapply(seq_len(nbin), function(j) {
        as.numeric(names(sample(data.cut[which(data.cut == j)], ctrl, FALSE)))
      })
      sample_feature <- c(features.exp, unlist(ctrl.use))
      features.use <- rank(sample_feature) / length(sample_feature)
      feature.score <- features.use[features]
      ctrl.score <- sort(features.use)[1:length(features.exp)]
      features.scores.use <- mean(feature.score) - mean(ctrl.score)
      return(features.scores.use)
    })
    
    if (ref) {
      features_sample <- all_sample_TPM1[features, ]
      features_ref <- reference[features, names(markers)[k]]
      div_percent <- features_sample / features_ref
      flag_lower <- (features_sample - features_ref) < 0
      number_lower <- colSums(flag_lower)
      percentage <- sapply(seq_len(ncol(div_percent)), function(i) {
        ifelse(number_lower[i] == 0, 1, sum(div_percent[, i][flag_lower[, i]], (length(flag_lower[, i]) - sum(flag_lower[, i]))) / length(features))
      })
      features.scores.vec <- features.scores.vec * percentage
    }
    
    return(features.scores.vec)
  })))
  
  rownames(features.scores.df) <- names(markers)
  return(features.scores.df)
}

