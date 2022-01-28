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

  p1=round(mean(data[[trtvar]]==1),digits=3)
  logit.weight.num=(data[[yvar]]-data[[trtvar]]*d-off.set)*data[["squant.weight"]]
  logit.weight.den=data[[trtvar]]*p1+(1-data[[trtvar]])/2
  logit.weight=logit.weight.num/logit.weight.den
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
eval.squant.pred = function(yvar, censorvar, trtvar, trtcd, dir, type, data.pred, brief=FALSE){
  #yvar: response variable name
  #censorvar: event indicator variable name 0-censor 1-event
  #trtvar: trt variable names
  #trtcd: code for the treatment arm, e.g. trtcd="treatment"
  #dir:
  #   "larger", prefer larger response after taking the drug, i.e., E(Y|X,TRT)-E(Y|X,CTRL)>=2d
  #   "smaller", prefer smaller response after taking the drug, i.e.,  E(Y|X,CTRL)-E(Y|X,TRT)>=2d
  #type: "c" = continuous, "b" = binary, "s" = survival
  #data.pred: the data frame of prediction, output of predict_squant
  #brief: whether to only calculate interaction p values or also a whole bunch of other statistics

  if(!is.data.frame(data.pred)) stop("data.pred needs to be a data frame.")
  if(length(yvar)>1) stop("yvar can only contain one variable.")
  if(!is.null(censorvar) && is.na(censorvar)) censorvar=NULL

  inter.pval=NULL
  pos.pval=NULL
  neg.pval=NULL
  group.stats=NULL

  #missing data handling
  data.pred=data.pred[c(yvar,censorvar,trtvar,"squant.subgroup")]
  data.pred=data.pred[complete.cases(data.pred),]
  data.pred[[trtvar]] = (data.pred[[trtvar]]==trtcd)*1


  ##### interaction onse-sided p value ######
  if(type=="c"){
    fml = as.formula(paste(yvar, "~", trtvar, "*squant.subgroup", sep=""))
    res.inter = try(summary(lm(fml, data=data.pred))$coefficients[paste(trtvar,":squant.subgroup",sep=""),c("Estimate","Pr(>|t|)")],silent=TRUE)
  }else if(type=="b"){
    fml = as.formula(paste(yvar, "~", trtvar, "*squant.subgroup", sep=""))
    res.inter = try(summary(glm(fml, family="binomial",data=data.pred))$coefficients[paste(trtvar,":squant.subgroup",sep=""),c("Estimate","Pr(>|z|)")],silent=TRUE)
  }else if(type=="s"){
    fml = as.formula(paste("Surv(",yvar,",",censorvar, ")~", trtvar, "*squant.subgroup", sep=""))
    res.inter = try(summary(coxph(fml, data=data.pred))$coefficients[paste(trtvar,":squant.subgroup",sep=""),c("coef","Pr(>|z|)")], silent=TRUE)
    if(class(res.inter)!="try-error") res.inter["coef"]=-res.inter["coef"]
  }

  if(class(res.inter)!= "try-error" && !is.na(res.inter[1])){
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
      fml = as.formula(paste(yvar, "~", trtvar, sep=""))
      res.pos = try(summary(lm(fml, data=sigpos))$coefficients[trtvar,c("Estimate","Pr(>|t|)")],silent=TRUE)
      res.neg = try(summary(lm(fml, data=signeg))$coefficients[trtvar,c("Estimate","Pr(>|t|)")],silent=TRUE)
    }else if(type=="b"){
      fml = as.formula(paste(yvar, "~", trtvar, sep=""))
      res.pos = try(summary(glm(fml, family="binomial",data=sigpos))$coefficients[trtvar,c("Estimate","Pr(>|z|)")],silent=TRUE)
      res.neg = try(summary(glm(fml, family="binomial",data=signeg))$coefficients[trtvar,c("Estimate","Pr(>|z|)")],silent=TRUE)
    }else if(type=="s"){
      fml = as.formula(paste("Surv(",yvar,",",censorvar, ")~", trtvar , sep=""))
      res.pos = try(summary(coxph(fml, data=sigpos))$coefficients[trtvar,c("coef","Pr(>|z|)")], silent=TRUE)
      res.neg = try(summary(coxph(fml, data=signeg))$coefficients[trtvar,c("coef","Pr(>|z|)")], silent=TRUE)
      if(class(res.pos)!="try-error") res.pos["coef"]=-res.pos["coef"]
      if(class(res.neg)!="try-error") res.neg["coef"]=-res.neg["coef"]
    }

    if(class(res.pos)!= "try-error" && !is.na(res.pos[1])){
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

    if(class(res.neg)!= "try-error" && !is.na(res.neg[1])){
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

      fml = as.formula(paste("Surv(",yvar,",",censorvar, ")~", trtvar , sep=""))
      fit.pos = coxph(fml, data=sigpos)
      fit.neg = coxph(fml, data=signeg)
      hr.pos = exp(coef(fit.pos))
      hr.neg = exp(coef(fit.neg))
      hr.all = c(hr.pos, 1, hr.neg, 1)
      group.stats = data.frame(N=n.all, Mean.Surv = mean.all, SE.Mean.Surv = se.all, Median.Surv = median.all, HR = hr.all,
                               stringsAsFactors=FALSE, row.names = c("Sig+.Trt", "Sig+.Ctrl", "Sig-.Trt", "Sig-.Ctrl"))
    }


  }

  list(inter.pval=inter.pval, pos.group.pval=pos.pval, neg.group.pval=neg.pval, group.stats=group.stats)

}




###### cv for variable selection: different lambda and a specific d #########
cv.squant.pred = function(yvar, censorvar, xvars, trtvar, d, data, type, xvars.keep=NULL, alpha=1, fold=5){
  #yvar: response variable name
  #censorvar: event indicator variable name 0-censor 1-event
  #xvars: covariates variable names
  #trtvar: trt variable names
  #d: threshold to qualify for sig+, E(Y|X,T=1)-E(Y|X,T=-1)>=2d
  #data: the data frame after prep function
  #type: "c" = continuous, "b" = binary, "s" = survival
  #xvars.keep: the covariates that need to stay in the model
  #alpha: same alpha as in glmnet (controlling elastic net/lasso)
  #fold: # of CV folds

  data.wk = wk.data.pred(data=data, d=d, yvar=yvar, trtvar=trtvar)
  x.std = data.matrix(scale(data.wk[xvars]))
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
    lambda.seq = sort(lambda.seq, decreasing = TRUE)
  }

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

    data.wk.i = wk.data.pred(data=data[idx.train.i,], d=d, yvar=yvar, trtvar=trtvar)
    x.std.i = data.matrix(scale(data.wk.i[xvars]))
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
    data.pred.i = data.pred[c(yvar,censorvar,trtvar,paste("squant.subgroup",i,sep=""))]
    names(data.pred.i) = c(yvar,censorvar,trtvar,"squant.subgroup")
    inter.pval.all[i] = eval.squant.pred(yvar=yvar, censorvar=censorvar, trtvar=trtvar, trtcd=1, dir="larger", type=type,
                                    data.pred=data.pred.i, brief=TRUE)$inter.pval
  }

  list(beta.all = beta.all, lambda.seq=lambda.seq, inter.pval.all = inter.pval.all)

}



####### cv for final fitting: an ordered vector of covariates and a specific d ########
cv.fit.pred = function(yvar, censorvar, xvars.ordered, trtvar, d, data, type, fold=5){
  #yvar: response variable name
  #censorvar: event indicator variable name 0-censor 1-event
  #xvars.ordered: ordered covariates variable names based on cv.squant variable selection result
  #trtvar: trt variable names
  #d: threshold to qualify for sig+, E(Y|X,T=1)-E(Y|X,T=-1)>=2d
  #data: the data frame after prep function
  #type: "c" = continuous, "b" = binary, "s" = survival
  #fold: # of CV folds

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
  data.pred[,paste("squant.subgroup", 1:length(xvars.ordered),sep="")] = NA

  for(i in 1:fold){
    idx.train.i = c(1:N)[cv.idx!=i]
    idx.test.i = c(1:N)[cv.idx==i]
    data.wk.i = wk.data.pred(data=data[idx.train.i,], d=d, yvar=yvar, trtvar=trtvar)

    for(j in 1:length(xvars.ordered)){
      xvars.sel.ij = xvars.ordered[1:j]
      squant.fit.ij = fit.squant(data.wk=data.wk.i, xvars.sel=xvars.sel.ij)
      data.pred[idx.test.i, paste("squant.subgroup",j,sep="")] = predict_squant(squant.fit=squant.fit.ij,
                                                                                data=data[idx.test.i,])$squant.subgroup
    }
  }

  inter.pval.all = rep(NA, length(xvars.ordered))
  for(i in 1:length(inter.pval.all)){
    data.pred.i = data.pred[c(yvar,censorvar,trtvar,paste("squant.subgroup",i,sep=""))]
    names(data.pred.i) = c(yvar,censorvar,trtvar,"squant.subgroup")
    inter.pval.all[i] = eval.squant.pred(yvar=yvar, censorvar=censorvar, trtvar=trtvar, trtcd=1, dir="larger", type=type,
                                    data.pred=data.pred.i, brief=TRUE)$inter.pval
  }

  n.xvars.sel = which.min(inter.pval.all)
  list(inter.pval.all = inter.pval.all, n.xvars.sel = n.xvars.sel, data.pred = data.pred)

}





##### interaction plot ######
plot.squant.pred = function(group.stats, trt.name="Trt", ctrl.name="Ctrl"){
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
  if(!all(sapply(data[c(yvar,censorvar,trtvar,xvars)], function(x) is.numeric(x))))
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


  #####  cv for variable selection:generate ordered selected xvars #####
  beta.all = NULL
  inter.pval.all = NULL
  d.all = NULL
  for(d.i in d.seq){
    cv.i = cv.squant.pred(yvar=yvar, censorvar=censorvar, xvars=xvars, trtvar=trtvar, d=d.i,
                     data=data, type=type, xvars.keep=xvars.keep, alpha=alpha, fold=fold)
    beta.i = cv.i$beta.all
    inter.pval.i = cv.i$inter.pval.all

    beta.all = cbind(beta.all, beta.i[,!colnames(beta.i)%in%"fold", drop=FALSE])
    inter.pval.all = c(inter.pval.all, inter.pval.i)
    d.all = c(d.all, rep(d.i, length(inter.pval.i)))

    if(progress) pg_status(which(d.seq==d.i)/(n.cv+1), end=FALSE)

  }


  pval.adj = p.adjust(inter.pval.all, method = "BH")
  pval.sel.idx = which(pval.adj <= FDR)

  for(FDR.min.i in sort(unique(pval.adj))){
    idx.i = which(pval.adj <= FDR.min.i)
    if(length(unique(d.all[idx.i]))>=n.cv*0.15){
      FDR.min = FDR.min.i
      break
    }
  }

  if(length(pval.sel.idx)>0 && length(unique(d.all[pval.sel.idx]))>=n.cv*0.15){
    pval.sel.all = inter.pval.all[pval.sel.idx]
    d.sel.all = d.all[pval.sel.idx]
    beta.sel.all = beta.all[,pval.sel.idx, drop=FALSE]
    xvars.top = order.xvars(beta.sel.all=beta.sel.all, xvars=xvars, fold=fold)
    xvars.ordered = unique(c(xvars.keep, xvars.top$xvars.ordered))

    if(length(xvars.ordered)>0){
      ##### select the best d #####
      if(is.null(min.diff)){
        d.sel = d.sel.all[which.min(pval.sel.all)]
      }else{
        d.sel = min.diff/2
      }

      ##### cv for final fit ##########
      n.cv.final = 25
      inter.pval.all = matrix(NA, nrow=n.cv.final, ncol=length(xvars.ordered))
      names.pred.subgroup = paste("squant.subgroup", 1:length(xvars.ordered),sep="")
      pred.subgroup.all = matrix(NA, nrow=dim(data)[1]*n.cv.final, ncol=length(xvars.ordered))
      colnames(pred.subgroup.all) = names.pred.subgroup

      for(i in 1:n.cv.final){
        cv.i = cv.fit.pred(yvar=yvar, censorvar=censorvar, xvars.ordered=xvars.ordered, trtvar=trtvar,
                      d=d.sel, data=data, type=type, fold=fold)
        inter.pval.all[i,] = cv.i$inter.pval.all
        pred.subgroup.all[((i-1)*dim(data)[1]+1):(i*dim(data)[1]),] = data.matrix(cv.i$data.pred[, names.pred.subgroup, drop=FALSE])

        if(progress) pg_status((n.cv+i/(n.cv.final+1))/(n.cv+1), end=FALSE)

      }


      median.inter.pval.all = apply(inter.pval.all, 2, median)
      min.idx = which.min(median.inter.pval.all)

      if(length(min.idx)==1){
        n.xvars.sel = min.idx
        error.candidate = colMeans(abs(pred.subgroup.all[,1:min.idx,drop=FALSE]-pred.subgroup.all[,min.idx]))
        xvars.sel.final = xvars.ordered[1:which(error.candidate < 0.05)[1]]
        xvars.sel.final = unique(c(xvars.keep, xvars.sel.final))

        #### final fit and performance evaluation ###
        data.wk = wk.data.pred(data=data, d=d.sel, yvar=yvar, trtvar=trtvar)
        squant.fit = fit.squant(data.wk=data.wk, xvars.sel=xvars.sel.final)
        data.pred = predict_squant(squant.fit=squant.fit, data=data.org)$data.pred
        performance = eval.squant.pred(yvar=yvar, censorvar=censorvar, trtvar=trtvar, trtcd=trtcd, dir=dir, type=type,
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
          interpretation2 = "Subgroup Selection Objective: E(Y|X,TRT==trtcd)-E(Y|X,TRT!=trtcd) >= min.diff (i.e, 2*d.sel)"
        }else{
          interpretation2 = "Subgroup Selection Objective: E(Y|X,TRT!=trtcd)-E(Y|X,TRT==trtcd) >= min.diff (i.e, 2*d.sel)"
        }

        res.final = list(squant.fit = squant.fit, data.pred = data.pred, performance=performance, d.sel=d.sel, min.diff=d.sel*2,
                         xvars.top = xvars.top, FDR.min = FDR.min, interpretation1=interpretation1, interpretation2=interpretation2)

        if(is.na(performance$inter.pval)||performance$inter.pval>FDR) res.final=NULL

      }
    }
  }


  if(is.null(res.final)) res.final = list(squant.fit = NULL, data.pred = NULL, performance=NULL, d.sel=NULL, min.diff=NULL,
                                          xvars.top = NULL, FDR.min = FDR.min,
                                          interpretation1="No significant subgroup can be identified!", interpretation2=NULL)


  if(progress) {
    pg_status(1, end=FALSE)
    pg_status(1, end=TRUE)
  }

  res.final

}









