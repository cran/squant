#######################################################################
################## SQUANT Predictive Case: YAN SUN ####################
#######################################################################


######## data pre-processing ######
prep.pred = function(yvar, censorvar, xvars, trtvar, trtcd, data, weight, dir){
  #yvar: response variable name
  #censorvar: event indicator variable name 0-censor 1-event
  #xvars: covariates variable names
  #trtvar: trt variable names
  #trtcd: code for the treatment arm, e.g. trtcd="treatment"
  #data: the data frame
  #weight: w(x), weight of every observation, has to be a vector >0 or NULL.
  #dir:
  #   "larger", prefer larger response after taking the drug E(Y|X,TRT)-E(Y|X,CTRL)>=2d
  #   "smaller", prefer smaller response after taking the drug E(Y|X,CTRL)-E(Y|X,TRT)>=2d

  if(sum(c(yvar,censorvar,xvars,trtvar)%in%c("squant.weight", "squant.y.pseudo", "squant.logit.weight",
                                             "(Intercept)", "squant.subgroup", "squant.pseudo.int"))>0) {
    stop("Variables cannot have any of the following names: squant.weight,
         squant.y.pseudo, squant.logit.weight, squant.subgroup, (Intercept), squant.pseudo.int.")
  }

  if(!is.null(weight) && !is.numeric(weight)) stop("weight needs to be a numeric vector or NULL")
  if(is.null(weight)) weight=1
  if(any(weight <=0)) stop("weight cannot be <= 0")

  #missing data handling
  data=data[c(yvar,censorvar,trtvar,xvars)]
  data=cbind(data,squant.weight=weight)
  data=data[complete.cases(data),]

  #coding trt
  if(dir=="larger") data[[trtvar]]=c(-1,1)[(data[[trtvar]]==trtcd)+1]
  if(dir=="smaller") data[[trtvar]]=c(1,-1)[(data[[trtvar]]==trtcd)+1]


  data
}



### working data generation ####
wk.data.macro.pred = function(data, d, yvar, trtvar){
  #data: the data frame
  #d: threshold to qualify for sig+, E(Y|X,T=1)-E(Y|X,T=-1)>=2d
  #yvar: response variable name
  #trtvar: trt variable names

  off.set = min(data[[yvar]]-data[[trtvar]]*d)


  logit.weight=(data[[yvar]]-data[[trtvar]]*d-off.set)*data[["squant.weight"]]
  logit.weight[logit.weight<=0] = min(logit.weight[logit.weight>0])/10000

  data[["squant.logit.weight"]]=logit.weight
  data[["squant.y.psdo"]]=(data[[trtvar]]+1)/2

  data

}



wk.data.pred = function(data, d, yvar, trtvar){
  #data: the data frame
  #d: threshold to qualify for sig+, E(Y|X,T=1)-E(Y|X,T=-1)>=2d
  #yvar: response variable name
  #trtvar: trt variable names

  data2 = data
  data2[[trtvar]]=-1*data2[[trtvar]]
  data2[[yvar]] = -1*data2[[yvar]]
  data.wk1 = wk.data.macro.pred(data=data, d=d, yvar=yvar, trtvar=trtvar)
  data.wk2 = wk.data.macro.pred(data=data2, d=d, yvar=yvar, trtvar=trtvar)
  data.wk = rbind(data.wk1, data.wk2)
  data.wk

}



############## performance evaluation ##############
eval.squant.pred = function(yvar, censorvar, trtvar, trtcd, dir, type, data.pred, xvars.adj = NULL, brief=FALSE){
  #yvar: response variable name
  #censorvar: event indicator variable name 0-censor 1-event
  #trtvar: trt variable names
  #trtcd: code for the treatment arm, e.g. trtcd="treatment"
  #dir:
  #   "larger", prefer larger response after taking the drug, i.e., E(Y|X,TRT)-E(Y|X,CTRL)>=2d
  #   "smaller", prefer smaller response after taking the drug, i.e.,  E(Y|X,CTRL)-E(Y|X,TRT)>=2d
  #type: "c" = continuous, "b" = binary, "s" = survival
  #data.pred: the data frame of prediction, output of predict_squant
  #xvars.adj: adjusted prognostic covariates in the model
  #brief: whether to only calculate interaction p values or also a whole bunch of other statistics

  if(!is.data.frame(data.pred)) stop("data.pred needs to be a data frame.")
  if(length(yvar)>1) stop("yvar can only contain one variable.")
  if(!is.null(censorvar) && is.na(censorvar)) censorvar=NULL

  inter.pval=NULL
  pos.pval=NULL
  neg.pval=NULL
  group.stats=NULL

  #missing data handling
  data.pred=data.pred[c(yvar,censorvar,trtvar,"squant.subgroup", xvars.adj)]
  data.pred=data.pred[complete.cases(data.pred),]
  data.pred[[trtvar]] = (data.pred[[trtvar]]==trtcd)*1


  ##### interaction one-sided p value ######
  if(!is.null(xvars.adj)){
    xvars.adj.fml = paste("+",xvars.adj, collapse="", sep="")
  }else{
    xvars.adj.fml = NULL
  }

  if(type=="c"){
    fml = as.formula(paste(yvar, "~", trtvar, "*squant.subgroup", xvars.adj.fml, sep=""))
    res.inter = try(summary(lm(fml, data=data.pred))$coefficients[paste(trtvar,":squant.subgroup",sep=""),c("Estimate","Pr(>|t|)")],silent=TRUE)
  }else if(type=="b"){
    fml = as.formula(paste(yvar, "~", trtvar, "*squant.subgroup", xvars.adj.fml, sep=""))
    res.inter = try(summary(glm(fml, family="binomial", data=data.pred))$coefficients[paste(trtvar,":squant.subgroup",sep=""),c("Estimate","Pr(>|z|)")],silent=TRUE)
  }else if(type=="s"){
    fml = as.formula(paste("Surv(",yvar,",",censorvar, ")~", trtvar, "*squant.subgroup", xvars.adj.fml, sep=""))
    res.inter = try(summary(coxph(fml, data=data.pred))$coefficients[paste(trtvar,":squant.subgroup",sep=""),c("coef","Pr(>|z|)")], silent=TRUE)
    if(!is(res.inter, "try-error")) res.inter["coef"]=-res.inter["coef"]
  }

  if(!is(res.inter, "try-error")&& !is.na(res.inter[1])){
    if(dir=="larger" & res.inter[1]>=0){
      inter.pval = res.inter[2]/2
    }else if(dir=="larger" & res.inter[1]<0){
      inter.pval = 1-res.inter[2]/2
    }else if(dir=="smaller" & res.inter[1]>=0){
      inter.pval = 1-res.inter[2]/2
    }else if(dir=="smaller" & res.inter[1]<0){
      inter.pval = res.inter[2]/2
    }
  }else{
    inter.pval = NA
  }

  if(!brief & !is.na(inter.pval)){
    ###### trt one-sided p value in sig+ and sig- group ######
    sigpos = data.pred[data.pred$squant.subgroup==1,]
    signeg = data.pred[data.pred$squant.subgroup==0,]
    if(type=="c"){
      fml = as.formula(paste(yvar, "~", trtvar, xvars.adj.fml, sep=""))
      res.pos = try(summary(lm(fml, data=sigpos))$coefficients[trtvar,c("Estimate","Pr(>|t|)")],silent=TRUE)
      res.neg = try(summary(lm(fml, data=signeg))$coefficients[trtvar,c("Estimate","Pr(>|t|)")],silent=TRUE)
    }else if(type=="b"){
      fml = as.formula(paste(yvar, "~", trtvar, xvars.adj.fml, sep=""))
      res.pos = try(summary(glm(fml, family="binomial",data=sigpos))$coefficients[trtvar,c("Estimate","Pr(>|z|)")],silent=TRUE)
      res.neg = try(summary(glm(fml, family="binomial",data=signeg))$coefficients[trtvar,c("Estimate","Pr(>|z|)")],silent=TRUE)
    }else if(type=="s"){
      fml = as.formula(paste("Surv(",yvar,",",censorvar, ")~", trtvar, xvars.adj.fml, sep=""))
      res.pos = try(summary(coxph(fml, data=sigpos))$coefficients[trtvar,c("coef","Pr(>|z|)")], silent=TRUE)
      res.neg = try(summary(coxph(fml, data=signeg))$coefficients[trtvar,c("coef","Pr(>|z|)")], silent=TRUE)
      if(!is(res.pos, "try-error")) res.pos["coef"]=-res.pos["coef"]
      if(!is(res.neg, "try-error")) res.neg["coef"]=-res.neg["coef"]
    }

    if(!is(res.pos, "try-error") && !is.na(res.pos[1])){
      if(dir=="larger" & res.pos[1]>=0){
        pos.pval = res.pos[2]/2
      }else if(dir=="larger" & res.pos[1]<0){
        pos.pval = 1-res.pos[2]/2
      }else if(dir=="smaller" & res.pos[1]>=0){
        pos.pval = 1-res.pos[2]/2
      }else if(dir=="smaller" & res.pos[1]<0){
        pos.pval = res.pos[2]/2
      }
    }else{
      pos.pval = NA
    }

    if(!is(res.neg, "try-error") && !is.na(res.neg[1])){
      if(dir=="larger" & res.neg[1]>=0){
        neg.pval = res.neg[2]/2
      }else if(dir=="larger" & res.neg[1]<0){
        neg.pval = 1-res.neg[2]/2
      }else if(dir=="smaller" & res.neg[1]>=0){
        neg.pval = 1-res.neg[2]/2
      }else if(dir=="smaller" & res.neg[1]<0){
        neg.pval = res.neg[2]/2
      }
    }else{
      neg.pval = NA
    }


    #### group stats ####
    sigpos.trt = data.pred[data.pred$squant.subgroup==1&data.pred[[trtvar]]==1,]
    sigpos.ctrl = data.pred[data.pred$squant.subgroup==1&data.pred[[trtvar]]==0,]
    signeg.trt = data.pred[data.pred$squant.subgroup==0&data.pred[[trtvar]]==1,]
    signeg.ctrl = data.pred[data.pred$squant.subgroup==0&data.pred[[trtvar]]==0,]

    n.all = c(nrow(sigpos.trt), nrow(sigpos.ctrl), nrow(signeg.trt), nrow(signeg.ctrl))

    if(type=="c"){
      mean.all = c(mean(sigpos.trt[[yvar]]), mean(sigpos.ctrl[[yvar]]), mean(signeg.trt[[yvar]]), mean(signeg.ctrl[[yvar]]))
      sd.all = c(sd(sigpos.trt[[yvar]]), sd(sigpos.ctrl[[yvar]]), sd(signeg.trt[[yvar]]), sd(signeg.ctrl[[yvar]]))
      group.stats = data.frame(N=n.all, Mean=mean.all, SD=sd.all, stringsAsFactors=FALSE, row.names = c("Sig+.Trt", "Sig+.Ctrl", "Sig-.Trt", "Sig-.Ctrl"))
    }else if(type=="b"){
      event.rate.all = c(mean(sigpos.trt[[yvar]]==1), mean(sigpos.ctrl[[yvar]]==1), mean(signeg.trt[[yvar]]==1), mean(signeg.ctrl[[yvar]]==1))
      group.stats = data.frame(N=n.all, Event.Rate=event.rate.all, stringsAsFactors=FALSE, row.names = c("Sig+.Trt", "Sig+.Ctrl", "Sig-.Trt", "Sig-.Ctrl"))

    }else if(type=="s"){
      fml = as.formula(paste("Surv(",yvar,",",censorvar, ")~squant.subgroup+",trtvar, sep=""))
      table.all = summary(survfit(fml, data=data.pred), rmean="common")$table
      name.pos.trt = paste("squant.subgroup=1, ", trtvar, "=1",sep="")
      name.pos.ctrl = paste("squant.subgroup=1, ", trtvar, "=0",sep="")
      name.neg.trt = paste("squant.subgroup=0, ", trtvar, "=1",sep="")
      name.neg.ctrl = paste("squant.subgroup=0, ", trtvar, "=0",sep="")
      table.all = table.all[c(name.pos.trt, name.pos.ctrl, name.neg.trt, name.neg.ctrl),]

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



      fml = as.formula(paste("Surv(",yvar,",",censorvar, ")~", trtvar, xvars.adj.fml, sep=""))
      hr.pos = exp(summary(coxph(fml, data=sigpos))$coefficients[trtvar,"coef"])
      hr.neg = exp(summary(coxph(fml, data=signeg))$coefficients[trtvar,"coef"])
      hr.all = c(hr.pos, 1, hr.neg, 1)
      group.stats = data.frame(N=n.all, Mean.Surv = mean.all, SE.Mean.Surv = se.all, Median.Surv = median.all, HR = hr.all,
                               stringsAsFactors=FALSE, row.names = c("Sig+.Trt", "Sig+.Ctrl", "Sig-.Trt", "Sig-.Ctrl"))
    }


  }

  list(inter.pval=inter.pval, pos.group.pval=pos.pval, neg.group.pval=neg.pval, group.stats=group.stats)

}




###### cv for variable selection: different lambda and a specific d #########
cv.squant.pred = function(yvar, censorvar, xvars, trtvar, d, data, type, xvars.keep=NULL,
                          prog.adj=0, xvars.adj=NULL, alpha=1, fold=5){
  #yvar: response variable name
  #censorvar: event indicator variable name 0-censor 1-event
  #xvars: covariates variable names
  #trtvar: trt variable names
  #d: threshold to qualify for sig+, E(Y|X,T=1)-E(Y|X,T=-1)>=2d
  #data: the data frame after prep function
  #type: "c" = continuous, "b" = binary, "s" = survival
  #xvars.keep: the covariates that need to stay in the model
  #prog.adj: prognostic effect adjustment (output of prog.rm)
  #xvars.adj: prognostic covariates to adjust in the model (output of prog.rm)
  #alpha: same alpha as in glmnet (controlling elastic net/lasso)
  #fold: # of CV folds

  data.wk = data
  data.wk[[yvar]] = data.wk[[yvar]]-prog.adj
  data.wk = wk.data.pred(data=data.wk, d=d, yvar=yvar, trtvar=trtvar)

  x.sd = apply(data.wk[xvars], 2, sd, na.rm=T)
  x.sd[x.sd==0] = 1
  x.std = data.matrix(scale(data.wk[xvars], center = T, scale=x.sd))
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
    n.trt = sum(data[[trtvar]]==1)
    n.ctrl = sum(data[[trtvar]]==-1)
    cv.idx.trt = sample(rep(1:fold,ceiling(n.trt/fold))[1:n.trt])
    cv.idx.ctrl = sample(rep(1:fold,ceiling(n.ctrl/fold))[1:n.ctrl])
    cv.idx[data[[trtvar]]==1] = cv.idx.trt
    cv.idx[data[[trtvar]]==-1] = cv.idx.ctrl

  }else{
    n.trt.1 = sum(data[[trtvar]]==1&data[[yvar]]==1)
    n.trt.0 = sum(data[[trtvar]]==1&data[[yvar]]==0)
    n.ctrl.1 = sum(data[[trtvar]]==-1&data[[yvar]]==1)
    n.ctrl.0 = sum(data[[trtvar]]==-1&data[[yvar]]==0)

    cv.idx.trt.1 = sample(rep(1:fold,ceiling(n.trt.1/fold))[1:n.trt.1])
    cv.idx.trt.0 = sample(rep(1:fold,ceiling(n.trt.0/fold))[1:n.trt.0])
    cv.idx.ctrl.1 = sample(rep(1:fold,ceiling(n.ctrl.1/fold))[1:n.ctrl.1])
    cv.idx.ctrl.0 = sample(rep(1:fold,ceiling(n.ctrl.0/fold))[1:n.ctrl.0])

    cv.idx[data[[trtvar]]==1&data[[yvar]]==1] = cv.idx.trt.1
    cv.idx[data[[trtvar]]==1&data[[yvar]]==0] = cv.idx.trt.0
    cv.idx[data[[trtvar]]==-1&data[[yvar]]==1] = cv.idx.ctrl.1
    cv.idx[data[[trtvar]]==-1&data[[yvar]]==0] = cv.idx.ctrl.0
  }


  data.wk = data
  data.wk[[yvar]] = data.wk[[yvar]]-prog.adj
  data.pred = data
  data.pred[,paste("squant.subgroup", 1:length(lambda.seq),sep="")] = NA
  beta.all = matrix(NA, nrow=length(xvars)*fold, ncol=length(lambda.seq)+1,dimnames = list(rep(xvars,fold), c(paste("lambda",1:length(lambda.seq),sep=""),"fold")))
  for(i in 1:fold){
    idx.train.i = c(1:N)[cv.idx!=i]
    idx.test.i = c(1:N)[cv.idx==i]

    data.wk.i = wk.data.pred(data=data.wk[idx.train.i,], d=d, yvar=yvar, trtvar=trtvar)


    x.sd.i = apply(data.wk.i[xvars], 2, sd, na.rm=T)
    x.sd.i[x.sd.i==0] = 1
    x.std.i = data.matrix(scale(data.wk.i[xvars], center=T, scale=x.sd.i))
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


  inter.pval.all = rep(NA, length(lambda.seq))
  for(i in 1:length(lambda.seq)){
    data.pred.i = data.pred[c(yvar,censorvar,trtvar,paste("squant.subgroup",i,sep=""), xvars.adj)]
    names(data.pred.i) = c(yvar,censorvar,trtvar,"squant.subgroup", xvars.adj)
    inter.pval.all[i] = eval.squant.pred(yvar=yvar, censorvar=censorvar, trtvar=trtvar, trtcd=1, dir="larger", type=type,
                                    data.pred=data.pred.i, xvars.adj = xvars.adj, brief=TRUE)$inter.pval
  }

  list(beta.all = beta.all, lambda.seq=lambda.seq, inter.pval.all = inter.pval.all)

}


######################################### functions to remove prognostic factors #############################
####### cv for removing prognostic effect ########
cv.prog.rm = function(yvar, censorvar, xvars, trtvar, d, data, type, alpha=1, fold=5){
  #yvar: response variable name
  #censorvar: event indicator variable name 0-censor 1-event
  #xvars: covariates variable names
  #trtvar: trt variable names
  #d: prognostic effect threshold to qualify for sig+, E(Y|X,T=1)-E(Y|X,T=-1)>=2d
  #data: the data frame after prep function
  #type: "c" = continuous, "b" = binary, "s" = survival
  #alpha: same alpha as in glmnet (controlling elastic net/lasso)
  #fold: # of CV folds

  data.wk = data
  data.wk[[yvar]] = data.wk[[yvar]]*data.wk[[trtvar]]
  data.wk = wk.data.pred(data=data.wk, d=d, yvar=yvar, trtvar=trtvar)


  x.sd = apply(data.wk[xvars], 2, sd, na.rm=T)
  x.sd[x.sd==0] = 1
  x.std = data.matrix(scale(data.wk[xvars], center = T, scale = x.sd))
  x.std = cbind(squant.pseudo.int=1, x.std)

  lasso.fit = glmnet(x=x.std, y=factor(data.wk[,"squant.y.psdo"]), family="binomial", weights=data.wk[,"squant.logit.weight"],
                     alpha=alpha, nlambda=100, standardize=FALSE, intercept=TRUE)
  lambda.seq = lasso.fit$lambda

  if(all(is.na(lambda.seq))) {
    lambda.seq = c(0.2, 0.1, 0)
  }else if (any(is.na(lambda.seq))) {
    lambda.seq[is.na(lambda.seq)]=max(lambda.seq,na.rm = TRUE)+0.1
  }

  lambda.seq = c(max(lambda.seq, na.rm = T)*1.1, lambda.seq, min(lambda.seq, na.rm = T)*0.9)
  lambda.seq = sort(unique(lambda.seq), decreasing = TRUE)
  lambda.seq = lambda.seq[1:min(10, length(lambda.seq))]


  N = dim(data)[1]
  cv.idx = rep(NA, N)

  if(type!="b"){
    n.trt = sum(data[[trtvar]]==1)
    n.ctrl = sum(data[[trtvar]]==-1)
    cv.idx.trt = sample(rep(1:fold,ceiling(n.trt/fold))[1:n.trt])
    cv.idx.ctrl = sample(rep(1:fold,ceiling(n.ctrl/fold))[1:n.ctrl])
    cv.idx[data[[trtvar]]==1] = cv.idx.trt
    cv.idx[data[[trtvar]]==-1] = cv.idx.ctrl

  }else{
    n.trt.1 = sum(data[[trtvar]]==1&data[[yvar]]==1)
    n.trt.0 = sum(data[[trtvar]]==1&data[[yvar]]==0)
    n.ctrl.1 = sum(data[[trtvar]]==-1&data[[yvar]]==1)
    n.ctrl.0 = sum(data[[trtvar]]==-1&data[[yvar]]==0)

    cv.idx.trt.1 = sample(rep(1:fold,ceiling(n.trt.1/fold))[1:n.trt.1])
    cv.idx.trt.0 = sample(rep(1:fold,ceiling(n.trt.0/fold))[1:n.trt.0])
    cv.idx.ctrl.1 = sample(rep(1:fold,ceiling(n.ctrl.1/fold))[1:n.ctrl.1])
    cv.idx.ctrl.0 = sample(rep(1:fold,ceiling(n.ctrl.0/fold))[1:n.ctrl.0])

    cv.idx[data[[trtvar]]==1&data[[yvar]]==1] = cv.idx.trt.1
    cv.idx[data[[trtvar]]==1&data[[yvar]]==0] = cv.idx.trt.0
    cv.idx[data[[trtvar]]==-1&data[[yvar]]==1] = cv.idx.ctrl.1
    cv.idx[data[[trtvar]]==-1&data[[yvar]]==0] = cv.idx.ctrl.0
  }


  data.pred = data
  data.pred[,paste("squant.subgroup", 1:length(lambda.seq),sep="")] = NA
  beta.all = matrix(NA, nrow=length(xvars)*fold, ncol=length(lambda.seq)+1,dimnames = list(rep(xvars,fold), c(paste("lambda",1:length(lambda.seq),sep=""),"fold")))
  for(i in 1:fold){
    idx.train.i = c(1:N)[cv.idx!=i]
    idx.test.i = c(1:N)[cv.idx==i]

    data.wk.i = data[idx.train.i,]
    data.wk.i[[yvar]] = data.wk.i[[yvar]]*data.wk.i[[trtvar]]
    data.wk.i = wk.data.pred(data=data.wk.i, d=d, yvar=yvar, trtvar=trtvar)


    x.sd.i = apply(data.wk.i[xvars], 2, sd, na.rm=T)
    x.sd.i[x.sd.i==0] = 1
    x.std.i = data.matrix(scale(data.wk.i[xvars], center=T, scale=x.sd.i))
    x.std.i = cbind(squant.pseudo.int=1, x.std.i)

    lasso.fit.i = glmnet(x=x.std.i, y=factor(data.wk.i[,"squant.y.psdo"]), family="binomial", weights=data.wk.i[,"squant.logit.weight"],
                     alpha=alpha, lambda=lambda.seq, standardize=FALSE, intercept=TRUE)
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


  trt.pval.all = rep(NA, length(lambda.seq))
  ctrl.pval.all = rep(NA, length(lambda.seq))
  for(i in 1:length(lambda.seq)){
    data.pred.i = data.pred[c(yvar,censorvar,trtvar,paste("squant.subgroup",i,sep=""))]
    names(data.pred.i) = c(yvar,censorvar,trtvar,"squant.subgroup")
    eval.prog.rm.res.i = eval.prog.rm(yvar=yvar, censorvar=censorvar, trtvar=trtvar, trtcd=1, dir="larger", type=type,
                                      data.pred=data.pred.i)
    trt.pval.all[i] = eval.prog.rm.res.i$trt.pval
    ctrl.pval.all[i] = eval.prog.rm.res.i$ctrl.pval
  }

  list(beta.all = beta.all, lambda.seq=lambda.seq, trt.pval.all = trt.pval.all, ctrl.pval.all = ctrl.pval.all)

}




eval.prog.rm = function(yvar, censorvar, trtvar, trtcd, dir, type, data.pred){
  #yvar: response variable name
  #censorvar: event indicator variable name 0-censor 1-event
  #trtvar: trt variable names
  #trtcd: code for the treatment arm, e.g. trtcd="treatment"
  #dir:
  #   "larger", prefer larger response after taking the drug, i.e., E(Y|X,TRT)-E(Y|X,CTRL)>=2d
  #   "smaller", prefer smaller response after taking the drug, i.e.,  E(Y|X,CTRL)-E(Y|X,TRT)>=2d
  #type: "c" = continuous, "b" = binary, "s" = survival
  #data.pred: the data frame of prediction, output of predict_squant


  if(!is.data.frame(data.pred)) stop("data.pred needs to be a data frame.")
  if(length(yvar)>1) stop("yvar can only contain one variable.")
  if(!is.null(censorvar) && is.na(censorvar)) censorvar=NULL

  trt.pval=NULL
  ctrl.pval=NULL

  #missing data handling
  data.pred=data.pred[c(yvar,censorvar,trtvar,"squant.subgroup")]
  data.pred=data.pred[complete.cases(data.pred),]
  data.pred[[trtvar]] = (data.pred[[trtvar]]==trtcd)*1

  ###### sig group one-sided p value in trt and ctrl group ######
  trt_grp = data.pred[data.pred[[trtvar]]%in%1,]
  ctrl_grp = data.pred[data.pred[[trtvar]]%in%0,]

  if(type=="c"){
    fml = as.formula(paste(yvar, "~ squant.subgroup", sep=""))
    res.trt = try(summary(lm(fml, data=trt_grp))$coefficients["squant.subgroup",c("Estimate","Pr(>|t|)")],silent=TRUE)
    res.ctrl = try(summary(lm(fml, data=ctrl_grp))$coefficients["squant.subgroup",c("Estimate","Pr(>|t|)")],silent=TRUE)
  }else if(type=="b"){
    fml = as.formula(paste(yvar, "~ squant.subgroup", sep=""))
    res.trt = try(summary(glm(fml, family="binomial",data=trt_grp))$coefficients["squant.subgroup",c("Estimate","Pr(>|z|)")],silent=TRUE)
    res.ctrl = try(summary(glm(fml, family="binomial",data=ctrl_grp))$coefficients["squant.subgroup",c("Estimate","Pr(>|z|)")],silent=TRUE)
  }else if(type=="s"){
    fml = as.formula(paste("Surv(",yvar,",",censorvar, ")~squant.subgroup" , sep=""))
    res.trt = try(summary(coxph(fml, data=trt_grp))$coefficients["squant.subgroup",c("coef","Pr(>|z|)")], silent=TRUE)
    res.ctrl = try(summary(coxph(fml, data=ctrl_grp))$coefficients["squant.subgroup",c("coef","Pr(>|z|)")], silent=TRUE)
    if(!is(res.trt, "try-error")) res.trt["coef"]=-res.trt["coef"]
    if(!is(res.ctrl, "try-error")) res.ctrl["coef"]=-res.ctrl["coef"]
  }

  if(!is(res.trt, "try-error") && !is.na(res.trt[1])){
    if(dir=="larger" & res.trt[1]>=0){
      trt.pval = res.trt[2]/2
    }else if(dir=="larger" & res.trt[1]<0){
      trt.pval = 1-res.trt[2]/2
    }else if(dir=="smaller" & res.trt[1]>=0){
      trt.pval = 1-res.trt[2]/2
    }else if(dir=="smaller" & res.trt[1]<0){
      trt.pval = res.trt[2]/2
    }
  }else{
    trt.pval = NA
  }

  if(!is(res.ctrl, "try-error") && !is.na(res.ctrl[1])){
    if(dir=="larger" & res.ctrl[1]>=0){
      ctrl.pval = res.ctrl[2]/2
    }else if(dir=="larger" & res.ctrl[1]<0){
      ctrl.pval = 1-res.ctrl[2]/2
    }else if(dir=="smaller" & res.ctrl[1]>=0){
      ctrl.pval = 1-res.ctrl[2]/2
    }else if(dir=="smaller" & res.ctrl[1]<0){
      ctrl.pval = res.ctrl[2]/2
    }
  }else{
    ctrl.pval = NA
  }



  list(trt.pval = trt.pval, ctrl.pval = ctrl.pval)

}






prog.rm = function(yvar, censorvar=NULL, xvars, trtvar, data, type="c", alpha=1,
                   fold=5, n.cv = 50, FDR = 0.15){

  #yvar: response variable name
  #censorvar: event indicator variable name 0-censor 1-event
  #xvars: covariates variable names (predictors)
  #trtvar: trt variable name
  #data: the data frame after prep function
  #type: response type, "s" survival; "b" binary; "c" continuous
  #alpha: same alpha as in glmnet (alpha=1 is the lasso penalty)
  #fold: # of folds for internal CV for variable selection
  #n.cv: # of different "d" used for CV (i.e, the # of CV) to conduct variable selection
  #FDR: FDR control of the training process


  res.final = NULL
  #### generate sequence of d #######
  if(type=="c" | type == "s"){
    quant.y = quantile(data[[yvar]],probs=c(0.1,0.9))
    d.seq = seq(from=quant.y[1], to=quant.y[2], length.out = n.cv)
  }else{
    d.seq = seq(from=0.1, to=0.9, length.out = n.cv)
  }


  #####  cv for variable selection:generate ordered selected xvars #####
  beta.all = NULL
  trt.pval.all = NULL
  ctrl.pval.all = NULL
  d.all = NULL
  for(d.i in d.seq){
    cv.i = cv.prog.rm(yvar=yvar, censorvar=censorvar, xvars=xvars, trtvar=trtvar,
                      d=d.i, data=data, type=type, alpha=alpha, fold=fold)

    beta.i = cv.i$beta.all
    trt.pval.i = cv.i$trt.pval.all
    ctrl.pval.i = cv.i$ctrl.pval.all

    beta.all = cbind(beta.all, beta.i[,!colnames(beta.i)%in%"fold", drop=FALSE])
    trt.pval.all = c(trt.pval.all, trt.pval.i)
    ctrl.pval.all = c(ctrl.pval.all, ctrl.pval.i)
    d.all = c(d.all, rep(d.i, length(trt.pval.i)))

  }


  trt.pval.adj = p.adjust(trt.pval.all, method = "BH")
  ctrl.pval.adj = p.adjust(ctrl.pval.all, method = "BH")
  pval.sel.idx = which(trt.pval.adj <= FDR & ctrl.pval.adj <= FDR)


  if(length(pval.sel.idx)>0){
    trt.pval.sel.all = trt.pval.all[pval.sel.idx]
    ctrl.pval.sel.all = ctrl.pval.all[pval.sel.idx]
    d.sel.all = d.all[pval.sel.idx]
    beta.sel.all = beta.all[,pval.sel.idx, drop=FALSE]
    xvars.top = order.xvars(beta.sel.all=beta.sel.all, xvars=xvars, fold=fold)
    xvars.ordered = xvars.top$xvars.ordered$xvars.ordered

    if(length(xvars.ordered)>0){
      xvars.sel.final = xvars.ordered

      data.pred.all = data
      data.pred.all[,paste("squant.subgroup", 1:length(d.seq),sep="")] = NA
      data.wk = data
      data.wk[[yvar]] = data.wk[[yvar]]*data.wk[[trtvar]]
      for(i in 1:length(d.seq)){
        d.i = d.seq[i]
        data.wk.i = wk.data.pred(data=data.wk, d=d.i, yvar=yvar, trtvar=trtvar)

        squant.fit.i = fit.squant(data.wk=data.wk.i, xvars.sel=xvars.sel.final)
        data.pred.all[, paste("squant.subgroup",i,sep="")] = predict_squant(squant.fit=squant.fit.i, data=data)$squant.subgroup

      }


      slide.window1 = diag(1, length(d.seq))
      slide.window2 = rbind(0, slide.window1)
      slide.window2 = slide.window2[-nrow(slide.window2),]
      slide.window3 = rbind(slide.window1,0)
      slide.window3 = slide.window3[-1,]
      slide.window = (slide.window1 + slide.window2 + slide.window3)[,-1]

      d.subgrp.all = data.matrix(data.pred.all[,paste("squant.subgroup", 1:length(d.seq),sep="")]) %*% slide.window
      d.subgrp = apply(d.subgrp.all, 1, function(x) max(c(which(x>=2),0)))
      prog.adj = d.seq[d.subgrp+1]

      res.final = list(prog.adj = prog.adj, xvars.adj = xvars.sel.final)



    }
  }


  if(is.null(res.final)) res.final = list(prog.adj = 0, xvars.adj = NULL)




  res.final

}


##################################################################################################################



##### interaction plot ######
plotsquant.pred = function(group.stats, trt.name="Trt", ctrl.name="Ctrl"){
  #group.stats: output of squant function: $performance$group.stats; Or output of eval.squant: $group.stats
  #trt.name: trt name to be used in the graph
  #ctrl.name: ctrl name to be used in the graph

  if(!is.data.frame(group.stats)) stop("group.stats needs to be a data frame.")

  group.stats = group.stats[c("Sig+.Trt", "Sig+.Ctrl", "Sig-.Trt", "Sig-.Ctrl"),]
  group.stats[,"trt"] = factor(c(trt.name, ctrl.name, trt.name, ctrl.name), levels=c(trt.name, ctrl.name))
  group.stats[,"group"] = factor(c("Sig+ Group", "Sig+ Group", "Sig- Group","Sig- Group"), levels=c("Sig- Group","Sig+ Group"))

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
  trt=NULL
  se=NULL
  fig = ggplot(group.stats, aes(x=group, y=resp, group=trt, color=trt, linetype=trt))+
    labs(color="", linetype="")+
    geom_errorbar(aes(ymin=resp-se, ymax=resp+se),width=0, size=1.2, position=position_dodge(width=0.2))+
    xlab("")+
    ylab(ylabel)+
    geom_line(size=1.2, position=position_dodge(width=0.2))+
    theme_bw(base_size=15)+
    theme(text=element_text(size=15), axis.text=element_text(size=15), legend.text=element_text(size=15))

  fig

}





########## main function ################
squant.pred = function(yvar, censorvar=NULL, xvars, trtvar, trtcd=1, data, type="c", weight=NULL, dir="larger",
                  min.diff=NULL, xvars.keep=NULL, alpha=1, fold=5, n.cv = 50, FDR = 0.15, progress=FALSE){

  #yvar: response variable name
  #censorvar: event indicator variable name 0-censor 1-event
  #xvars: covariates variable names (predictors)
  #trtvar: trt variable name
  #trtcd: code for the treatment arm, e.g. trtcd="treatment"
  #data: the data frame for training
  #type: response type
  #       "s" survival; "b" binary; "c" continuous
  #weight: w(x), weight of every observation, has to be a vector >0 or NULL (all 1).
  #dir:
  #   "larger", prefer larger response after taking the drug E(Y|X,TRT)-E(Y|X,CTRL)>=2d
  #   "smaller", prefer smaller response after taking the drug E(Y|X,CTRL)-E(Y|X,TRT)>=2d
  #min.diff: subgroup selection objective, the smallest benefit after switching to Trt from ctrl: min.diff = 2d
        ### min.diff is usually a value >=0, when preferring larger response, we look at E(Y|X,TRT)-E(Y|X,CTRL)>=min.diff
        ### when preferring smaller response, we look at E(Y|X,CTRL)-E(Y|X,TRT)>=min.diff
        ### min.diff can be a value < 0 if we allow worse improvment than standard of care in exchange for better safety profile
        ### when NULL, the program will automatically select the best value
  #xvars.keep: the covariates that need to stay in the final model
  #alpha: same alpha as in glmnet (alpha=1 is the lasso penalty)
  #fold: # of folds for internal CV for variable selection
  #n.cv: # of different "d" used for CV (i.e, the # of CV) to conduct variable selection
  #FDR: FDR control of the training process
  #progress: a logical value (TRUE/FALSE), whether to display the program progress.


  #### input verification #####
  if(!is.data.frame(data)) stop("data needs to a data frame.")
  if(is.null(yvar)||is.null(xvars)||is.null(trtvar)) stop("yvar, xvars, trtvar canot be NULL.")
  if(length(yvar)>1) stop("yvar can only contain one variable.")
  if(!trtcd %in% unique(data[[trtvar]])) stop("trtcd is not assigned a correct value.")
  if(!dir %in% c("larger", "smaller")) stop("dir needs to be either larger or smaller.")
  if(fold<2) stop("fold has to be >=2 for cross validation.")
  if(alpha<=0 | alpha>1) stop("alpha needs to be >0 and <=1.")
  if(type=="s" && length(censorvar)==0) stop("event indicator (censorvar) is missing!")
  if(!type %in% c("c","b","s")) stop("type needs to be c, b or s.")
  if(n.cv < 10) stop("n.cv needs to be >= 10")

  if(!is.null(censorvar) && is.na(censorvar)) censorvar=NULL
  if(!all(make.names(c(yvar, censorvar, xvars, trtvar), unique=TRUE)%in%c(yvar, censorvar, xvars, trtvar)))
    stop(paste("Some variable (column) names are not syntactically valid or duplicated. \n",
               "Consider using make.names() to change the names automatically.", sep=""))
  if(!all(sapply(data[c(yvar,censorvar,xvars)], function(x) is.numeric(x))))
     stop("There are non-numeric columns.")



  if(length(min.diff)>1) {
    min.diff = min.diff[1]
    warning("Only the 1st element of min.diff will be used.")
  }else if(length(min.diff)==1&&is.na(min.diff)){
    min.diff = NULL
  }

  if(type=="b"){
    data[[yvar]] = as.numeric(as.factor(data[[yvar]]))
    data[[yvar]] = (data[[yvar]]>1)*1
    if(!all(c(0,1)%in%data[[yvar]])) stop("Response only has one unique value.")
  }


  #### data pre-processing ######
  res.final = NULL
  data.org = data
  data= prep.pred(yvar=yvar, censorvar=censorvar, xvars=xvars, trtvar=trtvar, trtcd=trtcd, data=data, weight=weight, dir=dir)

  #### generate sequence of d #######
  if(type=="c" | type == "s"){
    quant.y = quantile(data[[yvar]],probs=c(0.25,0.75))
    quant.diff = abs(quant.y[2]-quant.y[1])
    diff.seq = seq(from=-quant.diff/4, to=quant.diff, length.out = n.cv)
  }else{
    diff.seq = seq(from=-0.2, to=0.8, length.out = n.cv)
  }
  d.seq = diff.seq/2

  #### progress ####'
  if(progress) pg_status(0, end=FALSE)

  ##### estimate prognostic effect and predictors
  prog.rm.res = prog.rm(yvar=yvar, censorvar=censorvar, xvars=xvars, trtvar=trtvar,
                        data=data, type=type, alpha=alpha, fold=fold, n.cv = 50, FDR = 0.15)

  if(progress) pg_status(0.2, end=FALSE)


  #####
  #####  cv for variable selection:generate ordered selected xvars #####
  beta.all = NULL
  inter.pval.all = NULL
  d.all = NULL
  for(d.i in d.seq){
    cv.i = cv.squant.pred(yvar=yvar, censorvar=censorvar, xvars=xvars, trtvar=trtvar, d=d.i,
                     data=data, type=type, xvars.keep=xvars.keep, prog.adj=prog.rm.res$prog.adj,
                     xvars.adj=prog.rm.res$xvars.adj, alpha=alpha, fold=fold)
    beta.i = cv.i$beta.all
    inter.pval.i = cv.i$inter.pval.all

    beta.all = cbind(beta.all, beta.i[,!colnames(beta.i)%in%"fold", drop=FALSE])
    inter.pval.all = c(inter.pval.all, inter.pval.i)
    d.all = c(d.all, rep(d.i, length(inter.pval.i)))

    if(progress) pg_status(0.2 + (which(d.seq==d.i)/(n.cv+1))*0.8, end=FALSE)

  }


  pval.adj = p.adjust(inter.pval.all, method = "BH")
  pval.sel.idx = which(pval.adj <= FDR)

  for(FDR.min.i in sort(unique(pval.adj),decreasing = FALSE)){
    idx.i = which(pval.adj <= FDR.min.i)
    if(length(unique(d.all[idx.i]))>n.cv*0.15){
      FDR.min = FDR.min.i
      break
    }
  }



  if(length(pval.sel.idx)>0 && length(unique(d.all[pval.sel.idx])) > n.cv*0.15){
    pval.sel.all = inter.pval.all[pval.sel.idx]
    d.sel.all = d.all[pval.sel.idx]
    beta.sel.all = beta.all[,pval.sel.idx, drop=FALSE]
    xvars.top = order.xvars(beta.sel.all=beta.sel.all, xvars=xvars, fold=fold)
    xvars.ordered = unique(c(xvars.keep, xvars.top$xvars.ordered$xvars.ordered))

    if(length(xvars.ordered)>0){
      xvars.sel.final = xvars.ordered

      ##### select the best d #####
      if(is.null(min.diff)){
        res.d.sel=NULL
        for(d.i in unique(d.sel.all)){
          data.wk.i = data
          data.wk.i[[yvar]] = data.wk.i[[yvar]]-prog.rm.res$prog.adj
          data.wk.i = wk.data.pred(data=data.wk.i, d=d.i, yvar=yvar, trtvar=trtvar)
          squant.fit.i = fit.squant(data.wk=data.wk.i, xvars.sel=xvars.sel.final)
          data.pred.i = predict_squant(squant.fit=squant.fit.i, data=data.org)$data.pred
          performance.i = eval.squant.pred(yvar=yvar, censorvar=censorvar, trtvar=trtvar, trtcd=trtcd, dir=dir, type=type,
                                         data.pred=data.pred.i, xvars.adj = prog.rm.res$xvars.adj, brief=T)
          res.d.sel = rbind(res.d.sel, data.frame(d.sel = d.i, pval = performance.i$inter.pval))

        }
        d.sel = res.d.sel[which.min(res.d.sel$pval), "d.sel"]
        if(length(d.sel)==0) d.sel = d.sel.all[which.min(pval.sel.all)]

      }else{
        d.sel = min.diff/2
      }

      #### final fit and performance evaluation ###

      data.wk = data
      data.wk[[yvar]] = data.wk[[yvar]]-prog.rm.res$prog.adj
      data.wk = wk.data.pred(data=data.wk, d=d.sel, yvar=yvar, trtvar=trtvar)
      squant.fit = fit.squant(data.wk=data.wk, xvars.sel=xvars.sel.final)
      data.pred = predict_squant(squant.fit=squant.fit, data=data.org)$data.pred
      performance = eval.squant.pred(yvar=yvar, censorvar=censorvar, trtvar=trtvar, trtcd=trtcd, dir=dir, type=type,
                                     data.pred=data.pred, xvars.adj = prog.rm.res$xvars.adj, brief=FALSE)

      #### put together the results ####
      squant.fit.print = squant.fit
      squant.fit.print$coef = signif(squant.fit.print$coef, 3)

      interpretation1.1 = paste("Selected Positive Subgroup:",squant.fit.print$coef[squant.fit.print$variables=="(Intercept)"],"+")
      interpretation1.2 = paste(squant.fit.print$coef[squant.fit.print$variables!="(Intercept)"],
                                squant.fit.print$variables[squant.fit.print$variables!="(Intercept)"],collapse=" + ", sep="*")
      interpretation1 = paste(interpretation1.1, interpretation1.2,"> 0")
      interpretation1 = gsub(pattern="+ -", replacement="- ", x=interpretation1, fixed = TRUE)

      if(dir=="larger"){
        interpretation2 = "Subgroup Selection Objective: E(Y|X,TRT==trtcd)-E(Y|X,TRT!=trtcd) >= min.diff (i.e, 2*d.sel)"
      }else{
        interpretation2 = "Subgroup Selection Objective: E(Y|X,TRT!=trtcd)-E(Y|X,TRT==trtcd) >= min.diff (i.e, 2*d.sel)"
      }

      res.final = list(squant.fit = squant.fit, data.pred = data.pred, performance=performance, d.sel=d.sel, min.diff=d.sel*2,
                       xvars.top = xvars.top$xvars.ordered.all, FDR.min = FDR.min, prog.adj = prog.rm.res$prog.adj,
                       xvars.adj = prog.rm.res$xvars.adj,
                       interpretation1=interpretation1, interpretation2=interpretation2)

      if(is.na(performance$inter.pval)||performance$inter.pval>FDR) {
        res.final=list(squant.fit = NULL, data.pred = NULL, performance=NULL, d.sel=d.sel, min.diff=d.sel*2,
                       xvars.top = xvars.top$xvars.ordered.all, FDR.min = FDR.min, prog.adj = prog.rm.res$prog.adj,
                       xvars.adj = prog.rm.res$xvars.adj,
                       interpretation1="No significant subgroup can be identified!", interpretation2=NULL)
      }


    }
  }


  if(is.null(res.final)) res.final = list(squant.fit = NULL, data.pred = NULL, performance=NULL, d.sel=NULL, min.diff=NULL,
                                          xvars.top = NULL, FDR.min = FDR.min, prog.adj = NULL, xvars.adj = NULL,
                                          interpretation1="No significant subgroup can be identified!", interpretation2=NULL)


  if(progress) {
    pg_status(1, end=FALSE)
    pg_status(1, end=TRUE)
  }

  res.final

}









