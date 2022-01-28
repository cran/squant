#######################################################################
################## SQUANT General: YAN SUN ############################
#######################################################################


############## beta fitting using working data ###################
fit.squant = function (data.wk, xvars.sel){
  ##data.wk: working data
  ##xvars.sel: selected predictors

  fmla = as.formula(paste("squant.y.psdo ~",paste(c("1",xvars.sel),collapse="+")))
  fit=glm(fmla, family="quasibinomial", data=data.wk,
          weights=data.wk[,"squant.logit.weight"], control=glm.control(maxit=100))
  beta=coef(fit)
  squant.fit=data.frame(variables=names(beta), coef=beta, stringsAsFactors=FALSE)

  squant.fit
}




##### order selected xvars #######
order.xvars = function(beta.sel.all, xvars, fold){
  #beta.sel.all: selected betas that meet FDR
  #xvars: covariates variable names
  #fold: # of CV folds

  if(dim(beta.sel.all)[1]!=fold*length(xvars)) stop("dim(beta.sel.all)[1]!=fold*length(xvars).")
  beta.sel.all = abs(beta.sel.all)
  for(i in 1:fold){
    beta.sel.fold.i = beta.sel.all[((i-1)*length(xvars)+1):(i*length(xvars)),,drop=FALSE]
    beta.sel.fold.i = t(t(beta.sel.fold.i)/colSums(beta.sel.fold.i))
    beta.sel.all[((i-1)*length(xvars)+1):(i*length(xvars)),] = beta.sel.fold.i
  }
  beta.sel.all[is.nan(beta.sel.all)]=0

  xvars.weight = rep(0, length(xvars))
  names(xvars.weight) = xvars
  for(xvar.i in xvars){
    beta.sel.i = beta.sel.all[rownames(beta.sel.all)%in%xvar.i,]
    xvars.weight[xvar.i] = mean(beta.sel.i)
  }

  if(length(xvars.weight)>=2){
    xvars.weight = sort(xvars.weight)
    xvars.weight.inc = xvars.weight[2:length(xvars.weight)]-xvars.weight[1:(length(xvars.weight)-1)]
    max.inc = max(xvars.weight.inc)
    idx.inc.cut = which(xvars.weight.inc > max.inc/4)[1]+1
    if(idx.inc.cut==2) idx.inc.cut=1
    xvars.weight = xvars.weight[idx.inc.cut:length(xvars.weight)]
  }


  xvars.weight.cutoff = max(xvars.weight)/10
  xvars.weight = xvars.weight[xvars.weight > xvars.weight.cutoff]
  weight.ordered = sort(xvars.weight, decreasing = TRUE)
  xvars.ordered = data.frame(xvars.ordered=names(weight.ordered), Importance=weight.ordered,
                             stringsAsFactors = FALSE)
  xvars.ordered

}




############## predictition based on tranining result #################
predict_squant=function(squant.fit, data){
  #squant.fit: effect coefficients and their names, output of squant function.
  #data: a dataset that contains required variables to make predictions

  if(!is.data.frame(data)) stop("data needs to be a data frame.")
  if(sum(names(data)%in%c("(Intercept)", "squant.subgroup"))>0) {
    stop("variables in data cannot have any of the following names: squant.subgroup, (Intercept)")
  }

  data=cbind("(Intercept)"=1,data)
  xvars=squant.fit[,"variables"]
  if(!all(xvars%in%names(data))) stop("data does not contain all required variables to make predictions")

  x=data.matrix(data[,xvars,drop=FALSE])
  sgn=((x%*%squant.fit[,"coef"])>0)+0
  data.pred=cbind(data,squant.subgroup=sgn)

  data.pred = data.pred[!names(data.pred)%in%"(Intercept)"]
  list(squant.subgroup=as.vector(sgn), data.pred=data.pred)
}




################ progress status ########################
pg_status = function(x, end=FALSE){
  #x: a numeric value between 0 and 1 to relfect the progress
  #end: whether or not to end the progress stats line (a new line)

  if(end){
    cat("\n")
    flush.console()
  }else{
    progress = paste("Progress: ",round(x*100), "%", sep="")
    cat("\r", progress)
    flush.console()
  }
}





