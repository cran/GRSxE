#' Testing gene-environment interactions
#'
#' Fitting and evaluating GRS (genetic risk scores) for testing the
#' presence of GxE (gene-environment) interactions.
#'
#' The GRS is usually constructed through random forests for taking gene-gene
#' interactions into account and using its OOB (out-of-bag) prediction
#' mechanism. Alternatively, a classical GRS construction approach can be
#' employed by fitting an elastic net.
#' Bagging can also be applied to fit multiple elastic net models to also be
#' able to perform OOB predictions.
#'
#' The advantage of OOB predictions is that they allow the GRS model to be
#' constructed on the full available data, while performing unbiased
#' predictions also on the full available data.
#' Thus, both the GRS construction and the GxE interaction testing can utilize
#' all observations.
#'
#' If desired, sampling can be performed without replacement in contrast to
#' the classical bagging approach that utilizes bootstrap sampling.
#'
#' Potentially confounding variables can also be supplied that will then be
#' adjusted for in the GxE interaction testing.
#'
#' This function uses a GLM (generalized linear model) for modelling the
#' marginal genetic effect, marginal environmental effect, the GRSxE interaction
#' effect, and potential confounding effects.
#' The fitted GLM is returned, which can be, e.g., inspected via
#' \code{summary(...)} to retrieve the Wald test p-values for the individual
#' terms. The p-value corresponding to the \code{G:E} term is the p-value
#' for testing the presence of a GRSxE interaction.
#'
#' @param X Matrix or data frame of genetic variables such as SNPs usually
#'   coded as 0-1-2.
#' @param y Numeric vector of the outcome/phenotype. Binary outcomes such
#'   as a disease status should be coded as 0-1 (control-case).
#' @param E Numeric vector of the environmental exposure.
#' @param C Optional data frame containing potentially confounding
#'   variables to be adjusted for.
#' @param test.type Testing type. The standard setting is \code{"bagging"},
#'   which employs its OOB (out-of-bag) prediction mechanism such that the full
#'   data can be used for both training the GRS and testing the GxE interaction.
#'   Alternatively, this can be set to \code{"holdout"}, which requires
#'   splitting the available data into a training data set and test data set.
#'   For that, \code{test.ind} needs to be set to the data indices used for
#'   testing.
#' @param B The number of bagging iterations if \code{test.type = "bagging"} is
#'   used. Also used as the number of trees grown in the random forest if
#'   \code{grs.type = "rf"} is set.
#' @param replace Should sampling with or without replacement be performed?
#'   Only used if \code{test.type = "bagging"} is set.
#' @param subsample Subsample fraction if \code{test.type = "bagging"} is used.
#' @param test.ind Vector of indices in the supplied data for testing the GxE
#'   interaction. Only used if \code{test.type = "holdout"} is set.
#'   The standard setting corresponds to a random 50:50 training-test split.
#' @param grs.type Type of GRS to be constructed. Either \code{"rf"} for a
#'   random forest or \code{"elnet"} for an elastic net.
#' @param grs.args Optional list of arguments passed to the GRS fitting
#'   procedure.
#'
#' @return An object of class \code{glm} is returned, in which \code{G:E}
#'   describes the GRSxE term.
#' @example examples/example.R
#' @references
#' \itemize{
#'   \item Lau, M., Schikowski, T. & Schwender, H. (2022).
#'   Efficient gene-environment interaction testing through
#'   bootstrap aggregating. To be submitted.
#'   \item Lau, M., Wigmann C., Kress S., Schikowski, T. & Schwender, H. (2022).
#'   Evaluation of tree-based statistical learning methods for constructing
#'   genetic risk scores. BMC Bioinformatics 23:97.
#'   \doi{https://doi.org/10.1186/s12859-022-04634-w}
#'   \item Breiman, L. (1996).
#'   Bagging predictors. Machine Learning 24:123–140.
#'   \doi{https://doi.org/10.1007/BF00058655}
#'   \item Breiman, L. (2001).
#'   Random Forests. Machine Learning 45:5–32.
#'   \doi{https://doi.org/10.1023/A:1010933404324}
#'   \item Friedman J., Hastie T. & Tibshirani R. (2010).
#'   Regularization Paths for Generalized Linear Models via Coordinate Descent.
#'   Journal of Statistical Software 33(1):1–22.
#'   \doi{https://doi.org/10.18637/jss.v033.i01}
#' }
#'
#' @importFrom ranger ranger
#' @importFrom glmnet cv.glmnet predict.glmnet
#' @importFrom stats as.formula predict
#' @importFrom utils modifyList
#' @export
GRSxE <- function(X, y, E, C = NULL,
                  test.type = "bagging", B = 500,
                  replace = TRUE, subsample = ifelse(replace, 1, 0.632),
                  test.ind = sample(nrow(X), floor(nrow(X)/2)),
                  grs.type = "rf", grs.args = list()) {
  N <- nrow(X)
  train.ind <- -test.ind

  if(any(!(y %in% c(0, 1)))) {
    mode <- "reg"
    rf.probability <- FALSE
    elnet.family <- "gaussian"
  } else {
    mode <- "bin"
    rf.probability <- TRUE
    elnet.family <- "binomial"
  }

  if(test.type == "holdout") {
    train <- data.frame(y = y[train.ind], X[train.ind,,drop=FALSE])
    test <- data.frame(y = y[test.ind], X[test.ind,,drop=FALSE])
  } else { # Bagging
    train <- data.frame(y = y, X) -> test
    train.ind <- 1:N -> test.ind
  }

  if(grs.type == "rf") {
    if(mode == "bin") train$y <- as.factor(train$y)
    grs.args <- modifyList(list(formula = y ~ ., data = train, min.node.size = floor(0.05 * nrow(train)),
                                num.trees = B, probability = rf.probability, replace = replace,
                                sample.fraction = subsample), grs.args)
    grs.model <- do.call(ranger, grs.args)
    if(test.type == "holdout") {
      grs <- predict(grs.model, test[,-1])$predictions
    } else {
      grs <- grs.model$predictions
    }
    if(mode == "bin") {
      grs <- grs[,2]
      grs <- logit(grs)
    }
  } else if (grs.type == "elnet") {
    if(test.type == "bagging") {
      N.sample <- round(subsample * N)
      grs <- rep(0, N)
      n.trees <- rep(0, N)
      for(i in 1:B) {
        bag <- sample(N, N.sample, replace = replace)
        oob <- setdiff(1:N, bag)
        grs.args <- modifyList(list(x = as.matrix(train[bag, -1]), y = train[bag, "y"], family = elnet.family, alpha = 0.5),
                               grs.args)
        grs.model <- do.call(cv.glmnet, grs.args)
        grs[oob] <- grs[oob] + predict(grs.model, as.matrix(train[oob, -1]), type = "link", s = grs.model$lambda.min)
        n.trees[oob] <- n.trees[oob] + 1
      }
      grs <- as.numeric(grs/n.trees)
    } else {
      grs.args <- modifyList(list(x = as.matrix(train[,-1]), y = train$y, family = elnet.family, alpha = 0.5),
                             grs.args)
      grs.model <- do.call(cv.glmnet, grs.args)
      grs <- predict(grs.model, as.matrix(test[,-1]), type="link", s=grs.model$lambda.min)
      grs <- as.numeric(grs)
    }
  } else {
    stop("GRS type not supported!")
  }

  if(!is.null(C)) C.test <- C[test.ind,,drop=FALSE] else C.test <- NULL
  gxe.mod <- GxE(grs, y[test.ind], E[test.ind], C = C.test)
  return(gxe.mod)
}

#' Testing individual gene-environment interactions
#'
#' Function for testing univariate GxE interactions, e.g., using single SNPs
#' or a GRS.
#'
#' This function uses a GLM (generalized linear model) for modelling the
#' marginal genetic effect, marginal environmental effect, the GxE interaction
#' effect, and potential confounding effects.
#' The fitted GLM is returned, which can be, e.g., inspected via
#' \code{summary(...)} to retrieve the Wald test p-values for the individual
#' terms. The p-value corresponding to the \code{G:E} term is the p-value
#' for testing the presence of a GxE interaction.
#'
#' @param G Numeric vector of a genetic variable such as a GRS (genetic risk
#'   score) or a SNP coded as 0-1-2.
#' @param y Numeric vector of the outcome/phenotype. Binary outcomes such
#'   as a disease status should be coded as 0-1 (control-case).
#' @param E Numeric vector of the environmental exposure.
#' @param C Optional data frame containing potentially confounding
#'   variables to be adjusted for.
#'
#' @return An object of class \code{glm} is returned, in which \code{G:E}
#'   describes the GxE term.
#' @references
#' \itemize{
#'   \item Lau, M., Schikowski, T. & Schwender, H. (2022).
#'   Efficient gene-environment interaction testing through
#'   bootstrap aggregating. To be submitted.
#' }
#'
#' @importFrom stats gaussian binomial glm
#' @export
GxE <- function(G, y, E, C = NULL) {
  if(any(!(y %in% c(0, 1)))) {
    glm.family <- gaussian(link = "identity")
  } else {
    glm.family <- binomial(link = "logit")
  }
  N <- length(y)
  if(is.null(C)) C <- matrix(nrow = N, ncol = 0)
  dat <- data.frame(G = as.numeric(G), y = as.numeric(y), E = as.numeric(E), C)
  form <- paste("y ~ G * E", paste("+", colnames(C), collapse = " ", recycle0 = TRUE))
  gxe.mod <- glm(as.formula(form), data = dat, family = glm.family)
  return(gxe.mod)
}

logit <- function(p) {
  lower_tol <- 1e-3
  upper_tol <- 1-1e-3
  p[p < lower_tol] <- lower_tol
  p[p > upper_tol] <- upper_tol
  log(p/(1-p))
}

