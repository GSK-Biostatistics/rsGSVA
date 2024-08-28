#' Compare the observed expression data to the training data distribution
#'
#' @param X matrix of expression data. Rows are genes and columns are subjects
#' @param train_expr matrix of reference expression data. Rows are genes and columns are subjects
#' @param gs Either a list of gene sets or indices of genes to compare the distributions. If NULL, then will default to using all genes
#' @param tail Either 'both', 'left', or 'right'. Specifies which side of the distribution to estimate probabilities.
#'
#' @author Lorin Towle-Miller
#'
#' @description
#' This function estimates the probability of observing the expression data as compared to the training data's empirical
#' distribution. The distribution is estimated for each gene separately, and then each of the observed values is compared
#' against the estimated distribution. To evaluate a directional bias, the tail of distribution may be specified.
#'
#' If tail="both", the calculation closely mimics a two-sided p-value calculation where the probability equals
#' 2*Pr_training(-|observed|), where Pr_training() denotes the empirical CDF from the training data.
#'
#' If tail="right", the probability equals 1-Pr_training(observed). I.e., this is the probability that the training data
#' is greater than the observed value.
#'
#' If tail="left", the probability equals Pr_training(observed). I.e., this is the probability that the training data
#' is less than the observed value.
#'
#' @return
#' returns object of S3 class 'ecdf_comparison' with the following:
#' - probabilities: Matrix with rows for subjects and columns for features containing the estimated probabilities that the input data falls on the emperical distribution of the reference data
#' - tail: User input for which tail to use
#' - gs: The set of genes (if any) that was input by the user
#'
#' @examples
#' \dontrun{
#' # Simulate reference data
#' reference_data <- simulate_data()
#' # Simulate input data (similar distribution)
#' input_data <- simulate_data()
#' # Simulate input data (convert to count data so different than reference)
#' input_data_counts <- simulate_data()
#' input_data_counts$x <- ceiling(10^input_data_counts$x)
#'
#' # Estimate probabilities for the values from comparable distributions
#' gene_probs_similar <- ecdf_comparison(X = input_data$x,
#'                                       train_expr = reference_data$x,
#'                                       tail = "both")
#' # Plot the gene probability summaries
#' plot(gene_probs_similar)
#'
#' # Estimate probabilities for the values from different distributions
#' # i.e., counts versus normally distributed
#' gene_probs_different <- ecdf_comparison(X = input_data_counts$x,
#'                                         train_expr = reference_data$x,
#'                                         tail = "both")
#' # Plot the gene probability summaries
#' plot(gene_probs_different)
#' }
#'
#' @export
ecdf_comparison <- function(X, train_expr, gs = NULL, tail = "both"){

  if(!is.null(gs))
    gs <- unique(unlist(gs))
  else
    gs <- 1:nrow(X)

  train_expr <- t(train_expr[gs,])
  X <- t(X[gs,])

  # Check the location on the distribution for each gene
  gene.density <- sapply(1:ncol(train_expr), function(i) {
    f <- ecdf(train_expr[,i])
    if(tail=="both"){
      sapply(f(X[,i]), FUN = function(j){2*min(c(1-j, j))})
    } else if(tail=="left"){
      sapply(f(X[,i]), FUN = function(j){j})
    } else if(tail=="right"){
      sapply(f(X[,i]), FUN = function(j){1-j})
    }
  })

  # If the genes that was passed was a character, then update the column names from the returned dataset
  if(!is.null(gs) && is.numeric(gs)){
    colnames(gene.density) <- gs
  }

  # Return the object
  gene.density <- list(probabilities = gene.density,
                       tail = tail,
                       gs = gs)

  class(gene.density) <- "ecdf_comparison"

  return(gene.density)
}


#' Plot the gene probability summaries that the observed expression data belongs
#' to the training data distributions
#'
#' @param x ecdf_comparison object
#' @param summary.func Function to summarize across subjects for each gene
#'
#' @author Lorin Towle-Miller
#'
#' @description
#' This function summarizes the probabilities across subjects for each gene and
#' plots the summary probabilities to a histogram. The red reference line at
#' 50% refers to a negligible difference between the observed and training data.
#'
#' @export
plot.ecdf_comparison <- function(x, summary.func = mean){
  summary_values <- apply(x$probabilities, 2, summary.func)

  summary_function_name <- deparse(summary.func)[2]
  summary_function_name <- gsub("UseMethod(\"", "", summary_function_name, fixed = T)
  summary_function_name <- gsub("\")", "", summary_function_name, fixed = T)

  # Check if there is a significant percentage of values below or above 0.5
  percent_below <- sum(summary_values < 0.5)/length(summary_values)
  percent_above <- sum(summary_values > 0.5)/length(summary_values)

  hist(summary_values, main = paste0("Histogram of ", summary_function_name, " probabilities"),
       xlim = c(0,1),
       xlab = ifelse(x$tail=="both", "2*Pr_training(-|observed|)",
                     ifelse(x$tail=="left", "Pr_training(observed)",
                            "1-Pr_training(observed)")))
  abline(v = 0.5, col = "red")
  # Add the user messages
  if(percent_below > 0.95){
    text(x=0.8, y = length(summary_values)*0.1, paste0("Caution!\n", round(percent_below, 3)*100, "% genes below 0.5"),
         col = "red")
  } else if(percent_above > 0.95){
    text(x=0.2, y = length(summary_values)*0.1, paste0("Caution!\n", round(percent_below, 3)*100, "% genes above 0.5"),
         col = "red")
  } else if(percent_below > percent_above){
    text(x=0.8, y = length(summary_values)*0.1, paste0(round(percent_below, 3)*100, "% genes below 0.5"),
         col = "darkgreen")
  } else {
    text(x=0.2, y = length(summary_values)*0.1, paste0(round(percent_above, 3)*100, "% genes above 0.5"),
         col = "darkgreen")
  }
}
