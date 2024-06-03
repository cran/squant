#######################################################################
########################## SQUANT: YAN SUN ############################
#######################################################################

#' The SQUANT method
#'
#' \code{squant} conducts subgroup identification based on quantitative criteria.
#'
#' This is the main function of SQUANT to train subgroup signatures.
#' This method can handle continuous, binary and survival endpoint for both
#' prognostic and predictive case. For the predictive case, the method aims at
#' identifying a subgroup for which treatment is better than control by at
#' least a pre-specified or auto-selected constant. For the prognostic case,
#' the method aims at identifying a subgroup that is at least better than a
#' pre-specified/auto-selected constant. The derived signature is a linear
#' combination of predictors, and the selected subgroup are subjects with
#' the signature > 0. The false discover rate when no true subgroup exists
#' is strictly controlled at a user-specified level.
#'
#' @param yvar A character. The response variable name in the \code{data}.
#'   The corresponding column in the \code{data} should be numeric.
#' @param censorvar A character or NULL. The event indicator variable name in
#'   the \code{data}. The corresponding column in the \code{data} should
#'   be 0(censor) or 1(event). Use NULL when it is not a time-to-event case.
#' @param xvars A vector of characters. The covariates (predictors) variable
#'   names in the \code{data}. The corresponding columns in the \code{data}
#'   should be numeric.
#' @param trtvar A character or NULL. The trt variable name in the \code{data}
#'   for the predictive case. The corresponding column in the \code{data}
#'   should contain the treatment assignments, and can be either numeric
#'   or character. Use NULL for the prognostics case.
#' @param trtcd The code for the treatment arm for the predictive case,
#'   e.g. trtcd="treatment" or trtcd=1, etc.
#' @param data The data frame for training.
#' @param type The response type. Use "s" for survival, "b" for binary, and
#'   "c" for continuous.
#' @param weight The weight of every observation, has to be a numeric vector>0
#'   or NULL (equivalent to all 1).
#' @param dir A character, "larger" or "smaller".
#'   When dir == "larger", larger response is preferred for the target subgroup.
#'   In the predictive case, it means selecting patients satisfying
#'                  E(Y|X,TRT)-E(Y|X,CTRL)>=quant.
#'   In the prognostic case, it means selecting patients satisfying
#'                  E(Y|X)>=quant.
#'   When dir == "smaller", smaller response is preferred for the target subgroup.
#'   In the predictive case, it means selecting patients satisfying
#'                  E(Y|X,CTRL)-E(Y|X,TRT)>=quant.
#'   In the prognostic case, it means selecting patients satisfying
#'                  E(Y|X)<=quant.
#' @param quant A numeric value or NULL. The quantitative subgroup selection criterion.
#'   Please see \code{dir}. When NULL, the program will automatically select the best
#'   quant based on cross validation.
#' @param xvars.keep A character vector. The names of variables that we want to keep
#'   in the final model.
#' @param alpha The same alpha as in \code{glmnet}. alpha=1 is the lasso penalty.
#' @param fold A numeric value. The number of folds for internal cross validation
#'   for variable selection.
#' @param n.cv A numeric value. The number of different values of \code{quant} used
#'   for cross validation. It's also the number of CV to conduct variable selection.
#' @param FDR A numeric value. The level of FDR control for variable selection and
#'   the entire training process.
#' @param progress a logical value (TRUE/FALSE), whether to display the program progress.
#' @return An object of "squant". A list containing the following elements.
#'   \item{squant.fit}{The fitted signature from training, which is the coefficients of
#'     the linear combination of predictors plus an intercept.}
#'   \item{data.pred}{The training data with the predicted subgroup in the last column.}
#'   \item{performance}{The output of eval_squant (excluding the data.pred).
#'     The performance of subgroup identification. In the predictive
#'     case, the performance includes the interaction p value, the p value of the
#'     trt difference in the selected positive group, the p value of the trt difference
#'     in the unselected negative group (all adjusted for prognostic markers if any) and
#'     the stats for each arm in each group. In the prognostic case, the performance
#'     includes p value of group comparison and the stats of each group.}
#'   \item{d.sel}{Closely related to quant.Please see element: \code{interpretation}.}
#'   \item{min.diff, threshold}{Please see \code{interpretation}.}
#'   \item{xvars.top}{The ordered variable importance list.}
#'   \item{FDR.min}{The minimum achievable FDR threshold so that a signature
#'     can be derived. This is useful when a pre-specified \code{FDR} does not lead to
#'     a signature, in which case the \code{FDR.min} can be used instead.}
#'   \item{prog.adj}{Prognostic effect contributed by xvars.adj for each subject (predictive case only).}
#'   \item{xvars.adj}{Important prognostic markers to adjust in the model (predictive case only).}
#'   \item{interpretation1}{Interpretation of the result.}
#'   \item{interpretation2}{Interpretation of the result.}
#' @example example_squant/example.R
#' @references Yan Sun, Samad Hedayat. Subgroup Identification based on Quantitative Objectives. (submitted)
#' @import stats graphics utils survival glmnet ggplot2
#' @importFrom methods is
#' @export
squant = function(yvar, censorvar=NULL, xvars, trtvar=NULL, trtcd=1, data, type="c", weight=NULL, dir="larger",
                  quant=NULL, xvars.keep=NULL, alpha=1, fold=5, n.cv = 50, FDR = 0.15, progress=TRUE){


  if(is.null(trtvar)){
    result = squant.prog(yvar=yvar, censorvar=censorvar, xvars=xvars, data=data, type=type, weight=weight, dir=dir,
                         threshold=quant, xvars.keep=xvars.keep, alpha=alpha, fold=fold, n.cv = n.cv,
                         FDR = FDR, progress=progress)
  }else if(is.character(trtvar) && length(trtvar)==1){
    result = squant.pred(yvar=yvar, censorvar=censorvar, xvars=xvars, trtvar=trtvar, trtcd=trtcd, data=data, type=type,
                         weight=weight, dir=dir, min.diff=quant, xvars.keep=xvars.keep, alpha=alpha, fold=fold, n.cv = n.cv,
                         FDR = FDR, progress=progress)
  }else{
    stop("trtvar should be either NULL or a character of length 1.")
  }

  class(result) = c("squant", "list")

  result

}




#' @export
print.squant = function(x, ...){
  #x: an object of class "squant" (a list)
  cat(x$interpretation1, fill=TRUE)
  cat(x$interpretation2, fill=TRUE)
}



#' SQUANT prediction
#'
#' \code{predict} assigns subgroup for each individual in a new dataset.
#'
#' This function assigns subgroup for each individual in a new dataset
#' based on the derived signature contained within the squant object.
#'
#' @param object The squant object, the signature of which will be applied
#'   to the specified data. The output of \code{squant} function.
#' @param data The data frame for prediction.
#' @param ... Ignored.
#' @return A data frame with the predicted subgroup in the last column.
#' @export
predict.squant = function(object, data, ...){
  if(class(object)[1]!="squant" || class(object)[2]!="list") stop("object should be of class squant and list.")
  data.pred = predict_squant(squant.fit=object$squant.fit, data=data)$data.pred
  data.pred
}




########### performance evaluation #################

#' SQUANT performance evaluation
#'
#' \code{eval_squant} evaluates the subgroup identification performance.
#'
#' This function evaluates the subgroup identification performance through
#' applying the derived signature (the squant object) to a specified dataset.
#' Note that when the specified dataset is the same as the training set,
#' the performance is usually over-optimistic and is subject to over-fitting.
#' Ideally, use an independent testing set to have an honest evaluation of
#' the performance.
#'
#' @param yvar A character. The response variable name in the \code{data}.
#'   The corresponding column in the \code{data} should be numeric.
#' @param censorvar A character or NULL. The event indicator variable name in
#'   the \code{data}. The corresponding column in the \code{data} should
#'   be 0(censor) or 1(event). Use NULL when it is not a time-to-event case.
#' @param trtvar A character or NULL. The trt variable name in the \code{data}
#'   for the predictive case. The corresponding column in the \code{data}
#'   should contain the treatment assignments, and can be either numeric
#'   or character. Use NULL for the prognostics case.
#' @param trtcd The code for the treatment arm for the predictive case,
#'   e.g. trtcd="treatment" or trtcd=1, etc.
#' @param dir A character, "larger" or "smaller".
#'   When dir == "larger", larger response is preferred for the target subgroup.
#'   In the predictive case, it means the derived signature from \code{squant}
#'   selects patients satisfying
#'                  E(Y|X,TRT)-E(Y|X,CTRL)>=quant.
#'   In the prognostic case, it means the derived signature from \code{squant}
#'   selects patients satisfying
#'                  E(Y|X)>=quant.
#'   When dir == "smaller", smaller response is preferred for the target subgroup.
#'   In the predictive case, it means the derived signature from \code{squant}
#'   selects patients satisfying
#'                  E(Y|X,CTRL)-E(Y|X,TRT)>=quant.
#'   In the prognostic case, it means the derived signature from \code{squant}
#'   selects patients satisfying
#'                  E(Y|X)<=quant.
#' @param type The response type. Use "s" for survival, "b" for binary, and
#'   "c" for continuous.
#' @param data The data frame for performance evaluation of the derived signature.
#' @param squant.out The squant object, the signature of which will be applied
#'   to the specified data. The output of \code{squant} function.
#' @param brief A logical value, TRUE or FALSE. When TRUE, only the most important p value
#'   will be reported.
#' @return An object of "eval_squant". A list containing the following elements.
#'   \item{inter.pval}{Treatment*subgroup Interaction p value (predictive case only,
#'     adjusted for prognostic markers if any).}
#'   \item{pos.group.pval}{The p value of the trt difference in the selected positive
#'     group (predictive case only, adjusted for prognostic markers if any).}
#'   \item{neg.group.pval}{The p value of the trt difference in the negative group
#'     (predictive case only, adjusted for prognostic markers if any).}
#'   \item{pval}{The p value of group comparison (prognostic case only).}
#'   \item{group.stats}{The performance of each arm by group (predictive case) or
#'     the performance of each group (prognostic case).}
#'   \item{data.pred}{The data with the predicted subgroup in the last column.}
#' @example example_squant/example.R
#' @export
eval_squant = function(yvar, censorvar, trtvar, trtcd=1, dir, type, data, squant.out, brief=FALSE){

  if(class(squant.out)[1]!="squant" || class(squant.out)[2]!="list") stop("squant.out should be of class squant and list.")
  data.pred = predict_squant(squant.fit=squant.out$squant.fit, data=data)$data.pred

  if(is.null(trtvar)){
    performance = eval.squant.prog(yvar=yvar, censorvar=censorvar, dir=dir, type=type, data.pred=data.pred, brief=brief)

  }else if(is.character(trtvar) && length(trtvar)==1){
    performance = eval.squant.pred(yvar=yvar, censorvar=censorvar, trtvar=trtvar, trtcd=trtcd, dir=dir, type=type,
                                   data.pred=data.pred, xvars.adj = squant.out$xvars.adj, brief=brief)

  }else{
    stop("trtvar should be either NULL or a character of length 1.")
  }

  performance$data.pred = data.pred

  class(performance) = c("eval_squant", "list")
  performance

}


#' @export
print.eval_squant = function(x, ...){
  #x: an object of class "eval_squant" (a list)
  cat("Apply the derived signature to the specified data.", fill=TRUE)
  if("inter.pval" %in% names(x)){
    cat(paste("Interaction p value (adjusted for prognostic markers if any):", x$inter.pval), fill=TRUE)
  }else if("pval" %in% names(x)){
    cat(paste("Group comparison p value:", x$pval), fill=TRUE)
  }

}




######### plot the result #########

#' Plot SQUANT result
#'
#' \code{plot} plots the subgroup identification performance.
#'
#' An interaction plot is plotted for the predictive case and a group
#' plot is plotted for the prognostic case.
#'
#' @param x A squant object. The output of \code{squant} function.
#' @param trt.name The name used on plot for the treatment arm.
#' @param ctrl.name The name used on plot for the control arm.
#' @param ... Ignored.
#' @return A ggplot.
#' @example example_squant/example.R
#' @export
plot.squant = function(x, trt.name="Trt", ctrl.name="Ctrl", ...){

  group.stats = x$performance$group.stats
  if(nrow(group.stats)==2){
    fig = plotsquant.prog(group.stats=group.stats)
  }else if(nrow(group.stats)==4){
    fig = plotsquant.pred(group.stats=group.stats, trt.name=trt.name, ctrl.name=ctrl.name)
  }else{
    stop("Wrong format of group.stats.")
  }

  fig

}


#' Plot eval_squant result
#'
#' \code{plot} plots the subgroup identification performance.
#'
#' An interaction plot is plotted for the predictive case and a group
#' plot is plotted for the prognostic case.
#'
#' @param x An eval_squant object. The output of \code{eval_squant} function.
#' @param trt.name The name used on plot for the treatment arm.
#' @param ctrl.name The name used on plot for the control arm.
#' @param ... Ignored.
#' @return A ggplot.
#' @example example_squant/example.R
#' @export
plot.eval_squant = function(x, trt.name="Trt", ctrl.name="Ctrl", ...){

  group.stats = x$group.stats
  if(nrow(group.stats)==2){
    fig = plotsquant.prog(group.stats=group.stats)
  }else if(nrow(group.stats)==4){
    fig = plotsquant.pred(group.stats=group.stats, trt.name=trt.name, ctrl.name=ctrl.name)
  }else{
    stop("Wrong format of group.stats.")
  }

  fig

}

















