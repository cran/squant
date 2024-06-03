########################################################################################
################## SQUANT Prognostic Case: YAN SUN #####################################
########################################################################################


######## data pre-processing ######
prep.prog = function(yvar, censorvar, xvars, data, weight, dir){
  #yvar: response variable name
  #censorvar: event indicator variable name 0-censor 1-event
  #xvars: covariates variable names
  #data: the data frame
  #weight: w(x), weight of every observation, has to be >0.
  #dir:
  #   "larger", prefer larger response  E(Y|X) - 0 >= 2d
  #   "smaller", prefer smaller response  0 -E(Y|X) >= 2d

  if(sum(c(yvar,censorvar,xvars)%in%c("squant.trt","squant.weight", "squant.y.pseudo", "squant.logit.weight",
                                      "(Intercept)", "squant.subgroup", "squant.pseudo.int"))>0) {
    stop("variables cannot have any of the following names: squant.trt, squant.weight,
         squant.y.pseudo, squant.logit.weight, squant.subgroup, (Intercept), squant.pseudo.int.")
  }

  if(!is.null(weight) && !is.numeric(weight)) stop("weight needs to be a numeric vector or NULL")
  if(is.null(weight)) weight=1
  if(any(weight <=0)) stop("weight cannot be <= 0")

  #missing data handling
  data=data[c(yvar,censorvar,xvars)]
  data=cbind(data,squant.weight=weight)
  data=data[complete.cases(data),]

  #change one arm to two arms
  if(dir=="larger"){
    data2 = data
    data[["squant.trt"]] = 1
    data2[[yvar]] = 0
    data2[["squant.trt"]] = -1
  }

  if(dir=="smaller"){
    data2=data
    data[["squant.trt"]] = -1
    data2[[yvar]]=0
    data2[["squant.trt"]] = 1
  }

  data = rbind(data, data2)

  data
}



### working data generation ####
wk.data.macro.prog = function(data, d, yvar){
  #data: the data frame
  #d: threshold to qualify for sig+, E(Y|X,T=1)-E(Y|X,T=-1)>=2d
  #yvar: response variable name

  off.set = min(data[[yvar]]-data[["squant.trt"]]*d)

  logit.weight=(data[[yvar]]-data[["squant.trt"]]*d-off.set)*data[["squant.weight"]]
  logit.weight[logit.weight<=0] = min(logit.weight[logit.weight>0])/10000

  data[["squant.logit.weight"]]=logit.weight
  data[["squant.y.psdo"]]=(data[["squant.trt"]]+1)/2

  data

}



wk.data.prog = function(data, d, yvar){
  #data: the data frame
  #d: threshold to qualify for sig+, E(Y|X,T=1)-E(Y|X,T=-1)>=2d
  #yvar: response variable name

  data2 = data
  data2[["squant.trt"]]=-1*data2[["squant.trt"]]
  data2[[yvar]] = -1*data2[[yvar]]
  data.wk1 = wk.data.macro.prog(data=data, d=d, yvar=yvar)
  data.wk2 = wk.data.macro.prog(data=data2, d=d, yvar=yvar)
  data.wk = rbind(data.wk1, data.wk2)
  data.wk

}



############## performance evaluation ##############
eval.squant.prog = function(yvar, censorvar, dir, type, data.pred, brief=FALSE){
  #yvar: response variable name
  #censorvar: event indicator variable name 0-censor 1-event
  #dir:
  #   "larger", prefer larger response, i.e., E(Y|X)-0>=2d
  #   "smaller", prefer smaller response,  i.e.,  0-E(Y|X)>=2d
  #type: "c" = continuous, "b" = binary, "s" = survival
  #data.pred: the data frame of prediction, output of predict_squant
  #brief: whether to only calculate interaction p values or also a whole bunch of other statistics

  if(!is.data.frame(data.pred)) stop("data.pred needs to be a data frame.")
  if(length(yvar)>1) stop("yvar can only contain one variable.")
  if(!is.null(censorvar) && is.na(censorvar)) censorvar=NULL

  pval=NULL
  group.stats=NULL

  #missing data handling
  data.pred=data.pred[c(yvar,censorvar,"squant.subgroup")]
  data.pred=data.pred[complete.cases(data.pred),]

  ###### one-sided p value comparing sig+ and sig- group ######
  if(type=="c"){
    fml = as.formula(paste(yvar, "~ squant.subgroup", sep=""))
    res.subgrp = try(summary(lm(fml, data=data.pred))$coefficients["squant.subgroup",c("Estimate","Pr(>|t|)")],silent=TRUE)

  }else if(type=="b"){
    fml = as.formula(paste(yvar, "~ squant.subgroup", sep=""))
    res.subgrp = try(summary(glm(fml, family="binomial",data=data.pred))$coefficients["squant.subgroup",c("Estimate","Pr(>|z|)")],silent=TRUE)

  }else if(type=="s"){
    fml = as.formula(paste("Surv(",yvar,",",censorvar, ")~", "squant.subgroup" , sep=""))
    res.subgrp = try(summary(coxph(fml, data=data.pred))$coefficients["squant.subgroup",c("coef","Pr(>|z|)")], silent=TRUE)
    if(!is(res.subgrp, "try-error")) res.subgrp["coef"]=-res.subgrp["coef"]

  }

  if(!is(res.subgrp, "try-error") && !is.na(res.subgrp[1])){
    if(dir=="larger" & res.subgrp[1]>=0){
      pval = res.subgrp[2]/2
    }else if(dir=="larger" & res.subgrp[1]<0){
      pval = 1-res.subgrp[2]/2
    }else if(dir=="smaller" & res.subgrp[1]>=0){
      pval = 1-res.subgrp[2]/2
    }else if(dir=="smaller" & res.subgrp[1]<0){
      pval = res.subgrp[2]/2
    }
  }else{
    pval = NA
  }



  if(!brief & !is.na(pval)){
    #### group stats ####
    sigpos = data.pred[data.pred$squant.subgroup==1,]
    signeg = data.pred[data.pred$squant.subgroup==0,]

    n.all = c(nrow(sigpos), nrow(signeg))

    if(type=="c"){
      mean.all = c(mean(sigpos[[yvar]]), mean(signeg[[yvar]]))
      sd.all = c(sd(sigpos[[yvar]]), sd(signeg[[yvar]]))
      group.stats = data.frame(N=n.all, Mean=mean.all, SD=sd.all, stringsAsFactors=FALSE, row.names = c("Sig+", "Sig-"))

    }else if(type=="b"){
      event.rate.all = c(mean(sigpos[[yvar]]==1), mean(signeg[[yvar]]==1))
      group.stats = data.frame(N=n.all, Event.Rate=event.rate.all, stringsAsFactors=FALSE, row.names = c("Sig+", "Sig-"))

    }else if(type=="s"){
      fml = as.formula(paste("Surv(",yvar,",",censorvar, ")~squant.subgroup", sep=""))
      table.all = summary(survfit(fml, data=data.pred), rmean="common")$table
      name.pos = "squant.subgroup=1"
      name.neg = "squant.subgroup=0"

      table.all = table.all[c(name.pos, name.neg),]

      if("*rmean"%in%colnames(table.all)){
        mean.all = table.all[,"*rmean"]
      }else if("rmean"%in%colnames(table.all)){
        mean.all = table.all[,"rmean"]
      }

      if("*se(rmean)"%in%colnames(table.all)){
        se.all = table.all[,"*se(rmean)"]
      }else if("se(rmean)"%in%colnames(table.all)){
        se.all = table.all[,"se(rmean)"]
      }

      median.all = table.all[,"median"]

      fml = as.formula(paste("Surv(",yvar,",",censorvar, ")~ squant.subgroup" , sep=""))
      fit = coxph(fml, data=data.pred)
      hr = exp(coef(fit))
      hr.all = c(hr, 1)
      group.stats = data.frame(N=n.all, Mean.Surv = mean.all, SE.Mean.Surv = se.all, Median.Surv = median.all, HR = hr.all,
                               stringsAsFactors=FALSE, row.names = c("Sig+", "Sig-"))
    }
  }


  list(pval=pval, group.stats=group.stats)

}




###### cv for variable selection: different lambda and a specific d #########
cv.squant.prog = function(yvar, censorvar, xvars, d, data, type, dir, xvars.keep=NULL, alpha=1, fold=5){
  #yvar: response variable name
  #censorvar: event indicator variable name 0-censor 1-event
  #xvars: covariates variable names
  #d: threshold to qualify for sig+, E(Y|X,T=1)-E(Y|X,T=-1)>=2d
  #data: the data frame after prep function
  #type: "c" = continuous, "b" = binary, "s" = survival
  #dir:
  #   "larger", prefer larger response, i.e., E(Y|X)-0>=2d
  #   "smaller", prefer smaller response,  i.e.,  0-E(Y|X)>=2d
  #xvars.keep: the covariates that need to stay in the model
  #alpha: same alpha as in glmnet (controlling elastic net/lasso)
  #fold: # of CV folds

  if(dir=="larger") {
    cd.org = 1
  }else{
    cd.org = -1
  }

  data.wk = wk.data.prog(data=data, d=d, yvar=yvar)

  x.sd = apply(data.wk[xvars], 2, sd, na.rm=T)
  x.sd[x.sd==0] = 1
  x.std = data.matrix(scale(data.wk[xvars], center = T, scale = x.sd))
  x.std = cbind(squant.pseudo.int=1, x.std)
  penalty.factor = rep(1, dim(x.std)[2])
  penalty.factor[colnames(x.std)%in%xvars.keep] = 0

  lasso.fit = glmnet(x=x.std, y=factor(data.wk[,"squant.y.psdo"]), family="binomial", weights=data.wk[,"squant.logit.weight"],
                     alpha=alpha, nlambda=100, standardize=FALSE, intercept=TRUE, penalty.factor=penalty.factor)
  lambda.seq = lasso.fit$lambda

  if(all(is.na(lambda.seq))) {
    lambda.seq = c(0.2, 0.1, 0)
  }else if (any(is.na(lambda.seq))) {
    lambda.seq[is.na(lambda.seq)]=max(lambda.seq,na.rm = TRUE)+0.1
  }

  lambda.seq = c(max(lambda.seq, na.rm = T)*1.1, lambda.seq, min(lambda.seq, na.rm = T)*0.9)
  lambda.seq = sort(unique(lambda.seq), decreasing = TRUE)
  lambda.seq = lambda.seq[1:min(25, length(lambda.seq))]


  N = dim(data)[1]
  cv.idx = rep(NA, N)

  if(type!="b"){
    n.org = sum(data[["squant.trt"]]==cd.org)
    cv.idx.org = sample(rep(1:fold,ceiling(n.org/fold))[1:n.org])

    cv.idx[data[["squant.trt"]]==cd.org] = cv.idx.org
    cv.idx[data[["squant.trt"]]==-cd.org] = cv.idx.org

  }else{
    n.org.1 = sum(data[["squant.trt"]]==cd.org&data[[yvar]]==1)
    n.org.0 = sum(data[["squant.trt"]]==cd.org&data[[yvar]]==0)
    cv.idx.org.1 = sample(rep(1:fold,ceiling(n.org.1/fold))[1:n.org.1])
    cv.idx.org.0 = sample(rep(1:fold,ceiling(n.org.0/fold))[1:n.org.0])
    cv.idx[data[["squant.trt"]]==cd.org&data[[yvar]]==1] = cv.idx.org.1
    cv.idx[data[["squant.trt"]]==cd.org&data[[yvar]]==0] = cv.idx.org.0
    cv.idx[data[["squant.trt"]]==-cd.org] = cv.idx[data[["squant.trt"]]==cd.org]

  }

  data.pred = data
  data.pred[,paste("squant.subgroup", 1:length(lambda.seq),sep="")] = NA
  beta.all = matrix(NA, nrow=length(xvars)*fold, ncol=length(lambda.seq)+1,dimnames = list(rep(xvars,fold), c(paste("lambda",1:length(lambda.seq),sep=""),"fold")))

  for(i in 1:fold){
    idx.train.i = c(1:N)[cv.idx!=i]
    idx.test.i = c(1:N)[cv.idx==i]

    data.wk.i = wk.data.prog(data=data[idx.train.i,], d=d, yvar=yvar)

    x.sd.i = apply(data.wk.i[xvars], 2, sd, na.rm=T)
    x.sd.i[x.sd.i==0] = 1
    x.std.i = data.matrix(scale(data.wk.i[xvars], center = T, scale = x.sd.i))
    x.std.i = cbind(squant.pseudo.int=1, x.std.i)

    lasso.fit.i = glmnet(x=x.std.i, y=factor(data.wk.i[,"squant.y.psdo"]), family="binomial", weights=data.wk.i[,"squant.logit.weight"],
                         alpha=alpha, lambda=lambda.seq, standardize=FALSE, intercept=TRUE, penalty.factor=penalty.factor)
    beta.i = as.matrix(coef(lasso.fit.i, s=lambda.seq))[xvars,,drop=FALSE]
    beta.all[((i-1)*length(xvars)+1):(i*length(xvars)),1:length(lambda.seq)] = beta.i
    beta.all[((i-1)*length(xvars)+1):(i*length(xvars)),"fold"] = i

    for(j in 1:length(lambda.seq)){
      xvars.sel.ij = xvars[beta.i[,j]!=0]
      squant.fit.ij = fit.squant(data.wk=data.wk.i, xvars.sel=xvars.sel.ij)
      data.pred[idx.test.i, paste("squant.subgroup",j,sep="")] = predict_squant(squant.fit=squant.fit.ij,
                                                                                data=data[idx.test.i,])$squant.subgroup
    }
  }

  pval.all = rep(NA, length(lambda.seq))

  for(i in 1:length(lambda.seq)){
    data.pred.i = data.pred[c(yvar,censorvar,"squant.trt",paste("squant.subgroup",i,sep=""))]
    names(data.pred.i) = c(yvar,censorvar,"squant.trt","squant.subgroup")
    data.pred.i = data.pred.i[data.pred.i[["squant.trt"]]==cd.org,]
    pval.all[i] = eval.squant.prog(yvar=yvar, censorvar=censorvar, dir=dir, type=type,
                                    data.pred=data.pred.i, brief=TRUE)$pval
  }


  list(beta.all = beta.all, lambda.seq=lambda.seq, pval.all = pval.all)

}









##### prognostic plot ######
plotsquant.prog = function(group.stats){
  #group.stats: output of squant function: $performance$group.stats; Or output of eval.squant: $group.stats

  if(!is.data.frame(group.stats)) stop("group.stats needs to be a data frame.")

  group.stats = group.stats[c("Sig+", "Sig-"),]
  group.stats[,"group"] = factor(c("Sig+ Group", "Sig- Group"), levels=c("Sig- Group","Sig+ Group"))

  if("Mean.Surv" %in% names(group.stats)){
    type="s"
  }else if ("Mean" %in% names(group.stats)){
    type="c"
  }else if ("Event.Rate" %in% names(group.stats)){
    type="b"
  }

  if(type=="s"){
    group.stats[,"resp"] = group.stats[,"Mean.Surv"]
    group.stats[,"se"] = group.stats[,"SE.Mean.Surv"]
    ylabel = "Restricted Mean Survival Time (with SE)"
  }else if (type == "c"){
    group.stats[,"resp"] = group.stats[,"Mean"]
    group.stats[,"se"] = group.stats[,"SD"]/sqrt(group.stats[,"N"])
    ylabel = "Mean (with SE)"
  }else if (type == "b"){
    group.stats[, "resp"] = group.stats[,"Event.Rate"]
    group.stats[, "se"] = sqrt(group.stats[,"Event.Rate"]*(1-group.stats[,"Event.Rate"])/group.stats[,"N"])
    ylabel = "Response/Event Rate (with SE)"
  }

  group=NULL
  resp=NULL
  se=NULL
  fig = ggplot(group.stats, aes(x=group, y=resp, group=1, color=group))+
    labs(color="")+
    geom_errorbar(aes(ymin=resp-se, ymax=resp+se),width=0, size=1.2)+
    xlab("")+
    ylab(ylabel)+
    geom_line(size=1.2, colour="black")+
    theme_bw(base_size=15)+
    theme(text=element_text(size=15), axis.text=element_text(size=15), legend.text=element_text(size=15), legend.position="none")

  fig

}









########## main function ################
squant.prog = function(yvar, censorvar=NULL, xvars, data, type="c", weight=NULL, dir="larger",
                       threshold=NULL, xvars.keep=NULL, alpha=1, fold=5, n.cv = 50,
                       FDR = 0.15, progress=FALSE){

  #yvar: response variable name
  #censorvar: event indicator variable name 0-censor 1-event
  #xvars: covariates variable names (predictors)
  #data: the data frame for training
  #type: response type
  #       "s" survival; "b" binary; "c" continuous
  #weight: w(x), weight of every observation, has to be >0 or NULL (all 1).
  #dir:
  #   "larger", prefer larger response E(Y|X)>=2d
  #   "smaller", prefer smaller response E(Y|X)<=-2d
  #threshold: subgroup selection objective: threshold = 2d or -2d
  ### when preferring larger response, we want to identify E(Y|X)>=threshold
  ### when preferring smaller response, we want to identify E(Y|X)<=threshold
  ### when NULL, the program will automatically select the best value
  #xvars.keep: the covariates that need to stay in the final model
  #alpha: same alpha as in glmnet (alpha=1 is the lasso penalty)
  #fold: # of folds for internal CV for variable selection
  #n.cv: # of different "d" used for CV (i.e, the # of CV) to conduct variable selection
  #FDR: FDR control of the training process
  #progress: a logical value (TRUE/FALSE), whether to display the program progress.


  #### input verification #####
  if(!is.data.frame(data)) stop("data needs to a data frame.")
  if(is.null(yvar)||is.null(xvars)) stop("yvar, xvars canot be NULL.")
  if(length(yvar)>1) stop("yvar can only contain one variable.")
  if(!dir %in% c("larger", "smaller")) stop("dir needs to be either larger or smaller.")
  if(fold<2) stop("fold has to be >=2 for cross validation.")
  if(alpha<=0 | alpha>1) stop("alpha needs to be >0 and <=1.")
  if(type=="s" && length(censorvar)==0) stop("event indicator (censorvar) is missing!")
  if(!type %in% c("c","b","s")) stop("type needs to be c, b or s.")
  if(n.cv < 10) stop("n.cv needs to be >= 10")

  if(!is.null(censorvar) && is.na(censorvar)) censorvar=NULL
  if(!all(make.names(c(yvar, censorvar, xvars), unique=TRUE)%in%c(yvar, censorvar, xvars)))
    stop(paste("Some variable (column) names are not syntactically valid or duplicated. \n",
               "Consider using make.names() to change the names automatically.", sep=""))
  if(!all(sapply(data[c(yvar,censorvar, xvars)], function(x) is.numeric(x))))
    stop("There are non-numeric columns.")



  if(length(threshold)>1) {
    threshold = threshold[1]
    warning("Only the 1st element of threshold will be used.")
  }else if(length(threshold)==1&&is.na(threshold)){
    threshold = NULL
  }

  if(type=="b"){
      data[[yvar]] = as.numeric(as.factor(data[[yvar]]))
      data[[yvar]] = (data[[yvar]]>1)*1
      if(!all(c(0,1)%in%data[[yvar]])) stop("Response only has one unique value.")
  }


  #### data pre-processing ######
  res.final = NULL
  data.org = data
  data= prep.prog(yvar=yvar, censorvar=censorvar, xvars=xvars, data=data, weight=weight, dir=dir)


  #### generate sequence of d #######
  if(type=="c" | type == "s"){
    quant.y = quantile(data.org[[yvar]],probs=c(0.25,0.75),na.rm=TRUE)
    threshold.seq = seq(from=quant.y[1], to=quant.y[2], length.out = n.cv)
  }else{
    threshold.seq = seq(from=0.1, to=0.9, length.out = n.cv)
  }
  d.seq = threshold.seq/2
  if(dir=="smaller") d.seq = -1*d.seq

  #### progress bar setup
  if(progress) pg_status(0, end=FALSE)

  #####  cv for variable selection:generate ordered selected xvars #####
  beta.all = NULL
  pval.all = NULL
  d.all = NULL
  for(d.i in d.seq){
    cv.i = cv.squant.prog(yvar=yvar, censorvar=censorvar, xvars=xvars, d=d.i, data=data,
                     type=type, dir=dir, xvars.keep=xvars.keep, alpha=alpha, fold=fold)
    beta.i = cv.i$beta.all
    pval.i = cv.i$pval.all

    beta.all = cbind(beta.all, beta.i[,!colnames(beta.i)%in%"fold", drop=FALSE])
    pval.all = c(pval.all, pval.i)
    d.all = c(d.all, rep(d.i, length(pval.i)))

    if(progress) pg_status(which(d.seq==d.i)/(n.cv+1), end=FALSE)
  }



  pval.adj = p.adjust(pval.all, method = "BH")
  pval.sel.idx = which(pval.adj <= FDR)

  for(FDR.min.i in sort(unique(pval.adj),decreasing = FALSE)){
    idx.i = which(pval.adj <= FDR.min.i)
    if(length(unique(d.all[idx.i])) > n.cv*0.15){
      FDR.min = FDR.min.i
      break
    }
  }



  if(length(pval.sel.idx)>0 && length(unique(d.all[pval.sel.idx])) > n.cv*0.15){
    pval.sel.all = pval.all[pval.sel.idx]
    d.sel.all = d.all[pval.sel.idx]
    beta.sel.all = beta.all[,pval.sel.idx, drop=FALSE]
    xvars.top = order.xvars(beta.sel.all=beta.sel.all, xvars=xvars, fold=fold)
    xvars.ordered = unique(c(xvars.keep, xvars.top$xvars.ordered$xvars.ordered))

    if(length(xvars.ordered)>0){

      xvars.sel.final = xvars.ordered

      ##### select the best d #####
      if(is.null(threshold)){
        res.d.sel = NULL

        for(d.i in unique(d.sel.all)){
          data.wk.i = wk.data.prog(data=data, d=d.i, yvar=yvar)
          squant.fit.i = fit.squant(data.wk=data.wk.i, xvars.sel=xvars.sel.final)
          data.pred.i = predict_squant(squant.fit=squant.fit.i, data=data.org)$data.pred
          performance.i = eval.squant.prog(yvar=yvar, censorvar=censorvar, dir=dir, type=type,
                                         data.pred=data.pred.i, brief=T)
          res.d.sel = rbind(res.d.sel, data.frame(d.sel = d.i, pval = performance.i$pval))

        }

        d.sel = res.d.sel[which.min(res.d.sel$pval), "d.sel"]

        if(length(d.sel)==0) d.sel = d.sel.all[which.min(pval.sel.all)]

      }else{
        d.sel = threshold/2
        if(dir=="smaller") d.sel = -1*d.sel
      }




      #### final fit and performance evaluation ###
      data.wk = wk.data.prog(data=data, d=d.sel, yvar=yvar)
      squant.fit = fit.squant(data.wk=data.wk, xvars.sel=xvars.sel.final)
      data.pred = predict_squant(squant.fit=squant.fit, data=data.org)$data.pred
      performance = eval.squant.prog(yvar=yvar, censorvar=censorvar, dir=dir, type=type,
                                     data.pred=data.pred, brief=FALSE)

      #### put together the results ####
      squant.fit.print = squant.fit
      squant.fit.print$coef = signif(squant.fit.print$coef, 3)

      interpretation1.1 = paste("Selected Positive Subgroup:",squant.fit.print$coef[squant.fit.print$variables=="(Intercept)"],"+")
      interpretation1.2 = paste(squant.fit.print$coef[squant.fit.print$variables!="(Intercept)"],
                                squant.fit.print$variables[squant.fit.print$variables!="(Intercept)"],collapse=" + ", sep="*")
      interpretation1 = paste(interpretation1.1, interpretation1.2,"> 0")
      interpretation1 = gsub(pattern="+ -", replacement="- ", x=interpretation1, fixed = TRUE)

      if(dir=="larger"){
        interpretation2 = "Subgroup Selection Objective: E(Y|X) >= threshold (i.e, 2*d.sel)"
      }else{
        interpretation2 = "Subgroup Selection Objective: E(Y|X) <= threshold (i.e, -2*d.sel)"
      }

      res.final = list(squant.fit = squant.fit, data.pred = data.pred, performance=performance, d.sel=d.sel,
                       threshold=ifelse(dir=="larger",2*d.sel, -2*d.sel), xvars.top = xvars.top$xvars.ordered.all,
                       FDR.min = FDR.min, interpretation1=interpretation1, interpretation2=interpretation2)

      if(is.na(performance$pval)||performance$pval>FDR) {
        res.final=list(squant.fit = NULL, data.pred = NULL, performance = NULL, d.sel=d.sel,
                       threshold=ifelse(dir=="larger",2*d.sel, -2*d.sel), xvars.top = xvars.top$xvars.ordered.all,
                       FDR.min = FDR.min,
                       interpretation1="No significant subgroup can be identified!", interpretation2=NULL)
      }




    }
  }

  if(is.null(res.final)) res.final = list(squant.fit = NULL, data.pred = NULL, performance = NULL, d.sel=NULL,
                                          threshold=NULL, xvars.top = NULL, FDR.min = FDR.min,
                                          interpretation1="No significant subgroup can be identified!", interpretation2=NULL)

  if(progress){
    pg_status(1, end=FALSE)
    pg_status(1, end=TRUE)
  }

  res.final

}












