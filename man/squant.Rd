% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SQUANT.R
\name{squant}
\alias{squant}
\title{The SQUANT method}
\usage{
squant(
  yvar,
  censorvar = NULL,
  xvars,
  trtvar = NULL,
  trtcd = 1,
  data,
  type = "c",
  weight = NULL,
  dir = "larger",
  quant = NULL,
  xvars.keep = NULL,
  alpha = 1,
  fold = 5,
  n.cv = 50,
  FDR = 0.15,
  progress = TRUE
)
}
\arguments{
\item{yvar}{A character. The response variable name in the \code{data}.
The corresponding column in the \code{data} should be numeric.}

\item{censorvar}{A character or NULL. The event indicator variable name in
the \code{data}. The corresponding column in the \code{data} should
be 0(censor) or 1(event). Use NULL when it is not a time-to-event case.}

\item{xvars}{A vector of characters. The covariates (predictors) variable
names in the \code{data}. The corresponding columns in the \code{data}
should be numeric.}

\item{trtvar}{A character or NULL. The trt variable name in the \code{data}
for the predictive case. The corresponding column in the \code{data}
should contain the treatment assignments, and can be either numeric
or character. Use NULL for the prognostics case.}

\item{trtcd}{The code for the treatment arm for the predictive case,
e.g. trtcd="treatment" or trtcd=1, etc.}

\item{data}{The data frame for training.}

\item{type}{The response type. Use "s" for survival, "b" for binary, and
"c" for continuous.}

\item{weight}{The weight of every observation, has to be a numeric vector>0
or NULL (equivalent to all 1).}

\item{dir}{A character, "larger" or "smaller".
When dir == "larger", larger response is preferred for the target subgroup.
In the predictive case, it means selecting patients satisfying
               E(Y|X,TRT)-E(Y|X,CTRL)>=quant.
In the prognostic case, it means selecting patients satisfying
               E(Y|X)>=quant.
When dir == "smaller", smaller response is preferred for the target subgroup.
In the predictive case, it means selecting patients satisfying
               E(Y|X,CTRL)-E(Y|X,TRT)>=quant.
In the prognostic case, it means selecting patients satisfying
               E(Y|X)<=quant.}

\item{quant}{A numeric value or NULL. The quantitative subgroup selection criterion.
Please see \code{dir}. When NULL, the program will automatically select the best
quant based on cross validation.}

\item{xvars.keep}{A character vector. The names of variables that we want to keep
in the final model.}

\item{alpha}{The same alpha as in \code{glmnet}. alpha=1 is the lasso penalty.}

\item{fold}{A numeric value. The number of folds for internal cross validation
for variable selection.}

\item{n.cv}{A numeric value. The number of different values of \code{quant} used
for cross validation. It's also the number of CV to conduct variable selection.}

\item{FDR}{A numeric value. The level of FDR control for variable selection and
the entire training process.}

\item{progress}{a logical value (TRUE/FALSE), whether to display the program progress.}
}
\value{
An object of "squant". A list containing the following elements.
  \item{squant.fit}{The fitted signature from training, which is the coefficients of
    the linear combination of predictors plus an intercept.}
  \item{data.pred}{The training data with the predicted subgroup in the last column.}
  \item{performance}{The output of eval_squant (excluding the data.pred).
    The performance of subgroup identification. In the predictive
    case, the performance includes the interaction p value, the p value of the
    trt difference in the selected positive group, the p value of the trt difference
    in the unselected negative group (all adjusted for prognostic markers if any) and
    the stats for each arm in each group. In the prognostic case, the performance
    includes p value of group comparison and the stats of each group.}
  \item{d.sel}{Closely related to quant.Please see element: \code{interpretation}.}
  \item{min.diff, threshold}{Please see \code{interpretation}.}
  \item{xvars.top}{The ordered variable importance list.}
  \item{FDR.min}{The minimum achievable FDR threshold so that a signature
    can be derived. This is useful when a pre-specified \code{FDR} does not lead to
    a signature, in which case the \code{FDR.min} can be used instead.}
  \item{prog.adj}{Prognostic effect contributed by xvars.adj for each subject (predictive case only).}
  \item{xvars.adj}{Important prognostic markers to adjust in the model (predictive case only).}
  \item{interpretation1}{Interpretation of the result.}
  \item{interpretation2}{Interpretation of the result.}
}
\description{
\code{squant} conducts subgroup identification based on quantitative criteria.
}
\details{
This is the main function of SQUANT to train subgroup signatures.
This method can handle continuous, binary and survival endpoint for both
prognostic and predictive case. For the predictive case, the method aims at
identifying a subgroup for which treatment is better than control by at
least a pre-specified or auto-selected constant. For the prognostic case,
the method aims at identifying a subgroup that is at least better than a
pre-specified/auto-selected constant. The derived signature is a linear
combination of predictors, and the selected subgroup are subjects with
the signature > 0. The false discover rate when no true subgroup exists
is strictly controlled at a user-specified level.
}
\examples{
#toy example#
set.seed(888)
x=as.data.frame(matrix(rnorm(200),100,2))
names(x) = c("x1", "x2")
trt = sample(0:1, size=100, replace=TRUE)
y= 2*x[,2]*trt+rnorm(100)
data = cbind(y=y, trt=trt, x)
res = squant(yvar="y", censorvar=NULL, xvars=c("x1", "x2"),
             trtvar="trt", trtcd=1, data=data, type="c", weight=NULL,
             dir="larger", quant=NULL, xvars.keep=NULL, alpha=1,
             fold=2, n.cv = 10, FDR = 0.1, progress=FALSE)


\donttest{
#predictive case with continuous response#
set.seed(888)
x=as.data.frame(matrix(rnorm(20000),200,100))
names(x) = paste("x", 1:100,sep="")
trt = sample(0:1, size=200, replace=TRUE)
y=x[,1]+x[,2]*trt+rnorm(200)
data = cbind(y=y, trt=trt, x)
res = squant(yvar="y", censorvar=NULL, xvars=paste("x", 1:100,sep=""),
             trtvar="trt", trtcd=1, data=data, type="c", weight=NULL,
             dir="larger", quant=NULL, xvars.keep=NULL, alpha=1,
             fold=5, n.cv = 50, FDR = 0.1)
res
#fitted signature#
res$squant.fit
#performance of the identified subgroup#
#including:
#  interaction p value,
#  p valve of trt difference in positive group,
#  p value of trt difference in negative group,
#  and stats for each arm in each group.
res$performance
#interpretation#
res$interpretation1
res$interpretation2

#evaluation of prediction performance#
eval.res = eval_squant(yvar="y", censorvar=NULL, trtvar="trt", trtcd=1, dir="larger",
                       type="c", data=data, squant.out=res, brief=FALSE)
#plot the subgroups#
plot(res, trt.name="Trt", ctrl.name="Ctrl")
plot(eval.res, trt.name="Trt", ctrl.name="Ctrl")



#prognostic case with survival response#
set.seed(888)
x=as.data.frame(matrix(rnorm(20000),200,100))
names(x) = paste("x", 1:100,sep="")
y=10*(10+x[,1]+rnorm(200))
data = cbind(y=y, x)
data$event = sample(c(rep(1,150),rep(0,50)))
res = squant(yvar="y", censorvar="event", xvars=paste("x", 1:100,sep=""),
             trtvar=NULL, trtcd=NULL, data=data, type="s", weight=NULL,
             dir="larger", quant=NULL, xvars.keep=NULL, alpha=1,
             fold=5, n.cv = 50, FDR = 0.1)
res

#fitted signature#
res$squant.fit
#performance of the identified subgroup#
res$performance
#evaluation of prediction performance#
eval.res = eval_squant(yvar="y", censorvar="event", trtvar=NULL, trtcd=NULL, dir="larger",
                       type="s", data=data, squant.out=res, brief=FALSE)

#plot the subgroups#
plot(res, trt.name=NULL, ctrl.name=NULL)
plot(eval.res, trt.name=NULL, ctrl.name=NULL)

}
}
\references{
Yan Sun, Samad Hedayat. Subgroup Identification based on Quantitative Objectives. (submitted)
}
