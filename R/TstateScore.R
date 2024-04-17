#' @title TstateScore
#' @description xxxxxx.
#' @details xxxxxx.
#' @param object xxxxxxx.
#' @param reference xxxxxxx.
#' @param nbin xxxxxxx.
#' @param ctrl xxxxxxx.
#' @param seed xxxxxxx.
#' @return xxxxxxxx.
#' @import ggplot2
#' @import dplyr
#' @importFrom stats na.omit rnorm
NULL
#' @export

Tstate_calcScore <- function (object, reference = ref_data, nbin=50, ctrl = 100, seed = 1) {
    set.seed(seed = seed)
    all_sample_HK <- mean(rowMeans(object[which(rownames(object) %in% HKgenes), ]))
    all_sample_TPM1 <- apply(object,2,function(i) { return(i * all_sample_HK/mean(na.omit(i[HKgenes]))) })

    features.scores.df <- as.data.frame(dplyr::bind_rows(lapply(seq_len(length(markers)),function(k) {
        message(paste0("Calculating of ",names(markers)[k])," state : ( ",k,"/",length(markers)," )")
        features <- markers[[k]]
        missing_features_input_sample <- setdiff(features, rownames(all_sample_TPM1))
        if (length(missing_features_input_sample) > 0) {
            warning("The following features are not present in the object: ",
                    paste(missing_features_input_sample, collapse = ", "), ifelse(test = FALSE,
                                                                     yes = ", attempting to find updated synonyms",
                                                                     no = ", not searching for symbol synonyms"),
                    call. = FALSE, immediate. = TRUE)
        }
        features <- intersect(features,rownames(all_sample_TPM1))
        missing_features_ref_data <- setdiff(features, rownames(x = reference))
        if (length(x = missing_features_ref_data) > 0) {
            warning("The following features are not present in the reference: ",
                    paste(missing_features_ref_data, collapse = ", "), ifelse(test = FALSE,
                                                                     yes = ", attempting to find updated synonyms",
                                                                     no = ", not searching for symbol synonyms"),
                    call. = FALSE, immediate. = TRUE)
        }
        features <- intersect(features,rownames(reference))
        if (length(features) < 3) { stop("The number of features less than 3") }
        all_sample_TPM2 <- all_sample_TPM1[which(!rownames(all_sample_TPM1) %in% features),]
        features.scores.vec <- pbapply::pbsapply(seq_len(ncol(all_sample_TPM2)),function(i) {
            features.exp <- all_sample_TPM1[features,i]
            data.avg <- all_sample_TPM2[,i]
            data.avg <- sort(data.avg)
            data.cut <- ggplot2::cut_number(data.avg + rnorm(n = length(data.avg))/1e+6,n = nbin, labels = FALSE, right = FALSE)
            names(data.cut) <- data.avg
            ctrl.use <- lapply(seq_len(nbin), function(j) {
                as.numeric(names(sample(data.cut[which(data.cut == j)], ctrl, FALSE)))
            })
            sample_feature <- c(features.exp,unlist(ctrl.use))
            features.use <- rank(sample_feature)/length(sample_feature)
            feature.score <- features.use[features]
            ctrl.score <- sort(features.use)[1:length(features.exp)]
            features.scores.use <- mean(feature.score) - mean(ctrl.score)
            return(features.scores.use)
        })
        features_sample = all_sample_TPM1[features,]
        features_ref = reference[features,names(markers)[k]]
        div_percent <- features_sample/features_ref
        flag_lower <- (features_sample - features_ref) < 0
        #number_lower <- colSums(flag_lower)
        percentage <- sapply(seq_len(ncol(div_percent)),function(i){
            ifelse(number_lower[i] == 0,1,sum( div_percent[,i][flag_lower[,i]],(length(flag_lower[, i]) - sum(flag_lower[, i]))) / length(features))
        })
        features.scores.vec <- features.scores.vec * percentage
        return(features.scores.vec)
    })))
    rownames(features.scores.df) <- names(markers)
    for (i in c(6,8:10)) {
        for (j in 1:ncol(all_sample_TPM1)) {
            pq <- as.numeric(reference["CD4", names(markers)[i]] )
            if (all_sample_TPM1["CD4",j] < pq) {
                percentage_T <- all_sample_TPM1["CD4", j]/pq
            }else {
                percentage_T <- 1
            }
            features.scores.df[i,j] <- features.scores.df[i,j]*percentage_T
        }
    }
    vec <- sapply(seq_len(ncol(features.scores.df)),function(i) {
        max(features.scores.df[8,i],features.scores.df[9,i],features.scores.df[10,i])
    })
    features.scores.df <- rbind(features.scores.df,vec)
    rownames(features.scores.df)[11] <- "Helper"
    features.scores.df <- features.scores.df[-8:-10,]
    return(features.scores.df)
}
