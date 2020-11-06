# BAYESIAN SPECIFICATION CURVE ANALYSIS
#           R FUNCTIONS
# 
# (c) Christoph Semken and David Rossell


#### BSCA #### 

#BMA-based Specification Curve Analysis
single_bsca= function(fit, coefidx=2, omitvars=c(1, coefidx), bbma, 
                 bmodels, maxmodels=100, alpha=0.05, nsamples=2000, col=1:length(coefidx), lty=1:length(coefidx), 
                 pch, plotlegend= ifelse(length(coefidx)>1,TRUE,FALSE), legend.cex=1.5, legend.pos='bottomright', 
                 show.score= TRUE, sort.by= 'score', sort.asc= FALSE, x.labels, 
                 var.labels, y.labels, height.vars=0.4) {
  # Input
  # - msfit: object of class "msfit", as returned by modelSelection
  # - coefidx: index of coefficients for which the plot should be created
  # - omitvars: variables to be ommitted from the variable configurations in the bottom panel, typically the variables in coefidx + the intercept
  # - bbma: optionally, BMA posterior samples on the regression coefficients. If not provided, bbma = rnlp(msfit,nsamples=2000).
  # - bmodels: optionally, coefficient estimates and posterior intervals for individual models. If not provided, bmodels= coefByModel(msfit,alpha=alpha)
  # - maxmodels: maximum number of models for which the coefficient estimates should be plotted
  # - alpha: the plot provides posterior intervals (only used if bbma or bmodels are not provided)
  # - nsamples: number of posterior samples to be obtained when estimating BMA posterior means and intervals. Only used if bbma not provided
  # - col: colours to be used for each coefidx
  # - lty: line type to be used for each coefidx
  # - pch: plotting character to be used for each coefidx
  # - plotlegend: set to TRUE to plot legend indicating the names of the treatment variables in coefidx
  # - legend.cex: size of legend, passed on to legend
  # - legend.pos: legend position (first argument in function legend)
  # - show.score: show score panel between coefficients and variable configurations if TRUE
  # - sort.by: sort coefficients by "score" or "coef" value
  # - sort.asc: sort in ascending order if TRUE
  # - x.labels: labels for plotted variables
  # - var.labels: labels for variables in the bottom panel
  # - height.vars: y position of start of bottom variable panel (as distance from x axis)
  #
  # Output: BMA Specification Curve Analysis plot
  if (missing(coefidx)) stop("You must specify the coefficient index for which you want to create the plot")
  if (!any(class(fit) == "msfit") & !any(class(fit) == "bas")) stop("fit must be of type msfit or bas")
  if (!sort.by %in% c('score', 'coef')) stop("sort.by must be 'score' or 'coef'")
  if (sort.by=="coef" & length(coefidx)>1) stop("Sorting by coefficient size can only be used with one coefficient")
  
  if (any(class(fit) == "bas")) {
    if ((exists('betaprior', where=fit) && fit$call$betaprior != call('bic.prior') ||
         exists('prior', where=fit) && fit$call$prior != 'BIC')
        && fit$call$modelprior != call('beta.binomial')) {
      stop('BAS results are currently only supported with `betaprior=bic.prior()` or `prior="BIC"` and `modelprior=beta.binomial()`')
    }
    # convert BAS fit to mombf fit
    msfit = as.msfit(fit)
    ppvec = msfit$postProb
    if (missing(bmodels)) bmodels = list(
      postmean = msfit$postmean, 
      ci.low = msfit$postmean + qnorm(alpha/2) * sqrt(msfit$postvar),
      ci.up = msfit$postmean - qnorm(alpha/2) * sqrt(msfit$postvar)
    )
    names(ppvec) = rownames(msfit$postmean)
    # limit to maxmodels
    maxmodels = min(maxmodels, length(ppvec))
    topmodels = sort(ppvec, index.return = TRUE, decreasing = TRUE)$ix[1:maxmodels]
    ppvec = ppvec[topmodels]
    bmodels = lapply(bmodels, function(x) x[topmodels,])
  } else {
    msfit = fit
    #Obtain posterior model probabilities
    pp = postProb(msfit)
    ppvec= pp$pp
    names(ppvec)= pp$modelid
  }
  
  if (missing(bbma)) bbma=rnlp(msfit=msfit, niter=nsamples, pp='norm')
  if (missing(bmodels)) bmodels= coefByModel(msfit,maxmodels=maxmodels)
  if (missing(pch)) pch= 16:(16+length(coefidx)-1)
  if (missing(x.labels)) x.labels=colnames(bbma)[coefidx]
  if (missing(var.labels)) var.labels=colnames(msfit$xstd)[-omitvars]
  
  if (length(var.labels) != ncol(msfit$xstd) - length(omitvars)) {
    stop("var.labels must contain one label for each non-omitted variable")
  }
  if (length(x.labels) != length(coefidx)) stop("x.labels must have the same length as coefidx")
  
  #Obtain BMA posterior samples
  bbmaest= cbind(colMeans(bbma[,coefidx,drop=FALSE]), t(apply(bbma[,coefidx,drop=FALSE], 2, quantile, probs=c(alpha/2,1-alpha/2))))
  
  #Obtain parameter estimates under each model
  pm= bmodels$postmean[,coefidx,drop=FALSE]
  ci.low= bmodels$ci.low[,coefidx,drop=FALSE]
  ci.up= bmodels$ci.up[,coefidx,drop=FALSE]
  nmodels= nrow(pm)
  if (any(names(pm) != names(ppvec)[1:nmodels])) { ppvec= ppvec[rownames(pm)] } else { ppvec= ppvec[1:nmodels] }
  
  #Set up plotting areas
  y1.top = ifelse(show.score, height.vars + 0.1, height.vars)
  fig1= c(x1=0.005,x2=0.09,y1=y1.top,y2=1)  #BMA posterior
  fig2= c(x1=0.11,x2=0.995,y1=y1.top,y2=1)  #Model-specific estimates & CI's
  fig3= c(x1=0.11,x2=0.995,y1=height.vars,y2=height.vars+0.11)  #Model scores
  fig4= c(x1=0.11,x2=0.995,y1=0.10,y2=height.vars+0.01) #Variable configurations
  
  #Set up graphical parameters
  if (length(col)==1) col= rep(col,length(coefidx))
  if (length(lty)==1) lty= rep(lty,length(coefidx))
  par(cex=0.7, mai=c(0.1,0.1,0.2,0.01), # make labels and margins smaller
      bg = "transparent") # switch off background to avoid obscuring adjacent plots in cowplot (https://github.com/wilkelab/cowplot/issues/69)
  d= lapply(coefidx, function(z) density(bbma[,z]))
  incld= sapply(1:length(d), function(i) bbmaest[i, 2] != bbmaest[i, 3])
  if (!any(incld)) incld = rep(TRUE, length(d))
  ymax= max(sapply(d[incld], function(z) max(z$y)))
  xlim= sapply(d, function(z) range(z$x)); xlim= c(min(xlim[1,]),max(xlim[2,]))
  xlim[1]= min(xlim[1], ci.low); xlim[2]= max(xlim[2], ci.up)
  
  #Plot BMA posterior
  par(fig=fig1) # define area (full area is x1=0, x2=1, y1=0, y2=1)
  if (incld[1]) {x=-d[[1]]$y; y=d[[1]]$x} else {x=y=0}
  plot(x,y,type='l',xaxt='n',yaxt='n',bty='n',xlim=c(-ymax,0),ylim=xlim,col=col[1])
  points(rep(0,nrow(bbmaest)),bbmaest[,1], col=col, pch=pch, cex=2)
  points(rep(0,nrow(bbmaest)),bbmaest[,2], col=col, pch='-', cex=2)
  points(rep(0,nrow(bbmaest)),bbmaest[,3], col=col, pch='-', cex=2)
  segments(x0=0, y0=bbmaest[,2], y1=bbmaest[,3], lwd=1.5, col=col, lty=lty, pch=pch)
  if (length(coefidx)>1) {
    for (i in 1:length(d)) {
      # plot the density unless CI.low and CI.up are the same
      if (incld[i]) {
        lines(-1-d[[i]]$y,d[[i]]$x,col=col[i],lty=lty[i])
      }
    }
  }
  
  #Sort [ppvec, pm, ci.low, ci.up]
  if (sort.by == "coef") {
    new.order= sort(pm, index.return = TRUE)$ix
    ppvec= ppvec[new.order]
    pm= pm[new.order, 1:dim(pm)[2], drop=FALSE]
    ci.low= ci.low[new.order, 1:dim(pm)[2], drop=FALSE]
    ci.up= ci.up[new.order, 1:dim(pm)[2], drop=FALSE]
  }
  
  #Plot point estimate and 95% intervals under each model
  par(fig=fig2, new=TRUE)
  cex= 2 * (ppvec / max(ppvec))
  ylim= xlim; xlim= c(1,nmodels)
  plot(1:nmodels, pm[,1], xlim=xlim, ylim=ylim, xaxt='n', pch=pch[1], col=col[1], cex=cex, yaxt=ifelse(missing(y.labels), 's', 'n'))
  segments(x0=1:nmodels, y0=ci.low[,1], y1=ci.up[,1], pch=pch[1], col=col[1])
  if (length(coefidx)>1) {
    for (i in 2:length(d)) {
      points(1:nmodels, pm[,i], pch=pch[i], col=col[i], cex=cex)
      segments(x0=1:nmodels, y0=ci.low[,i], y1=ci.up[,i], pch=pch[i], col=col[i], lty=lty[i])
    }
  }
  par(bg = "white")
  if (plotlegend) legend(legend.pos, x.labels, lty=lty, pch=pch, col=col, cex=legend.cex)
  par(bg = "transparent")
  if (!missing(y.labels)) axis(2, names(y.labels), y.labels)
  
  #Plot model scores
  if (show.score) {
    par(fig=fig3, new=TRUE)
    par(mai=c(0.1,0.1,0,.01)) #make upper margin smaller
    plot(1:nmodels, ppvec, xlim=xlim, type='l', xaxt='n', yaxp=c(0, round(max(ppvec),1)-0.1, 1))
    points(1:nmodels, ppvec, pch=16)
    mtext('Score',side=2,line=2)
  }
  
  #Plot variable configurations
  par(fig=fig4, new=TRUE)
  modelvars= strsplit(rownames(pm), split=',')
  ylim= c(0,ncol(msfit$xstd)+1-length(omitvars))
  ypos= (ylim[2]-1):1; names(ypos)= setdiff(1:ncol(msfit$xstd),omitvars)
  plot(NA,NA,ylim=ylim,xlim=xlim,yaxt='n',xlab='')
  for (i in 1:nmodels) {
    selvars= setdiff(modelvars[[i]],as.character(omitvars))
    points(rep(i,length(selvars)), ypos[selvars], pch=15)
  }
  axis(2, at=length(var.labels):1, labels=var.labels, las=1)
  mtext('Specification number',side=1,line=2)
  
  #ylim= c(0,ncol(msfit$xstd)+1)
  #plot(NA,NA,ylim=ylim,xlim=xlim,yaxt='n',xlab='')
  #for (i in 1:nmodels) points(rep(i,length(modelvars[[i]])), ylim[2]-as.numeric(modelvars[[i]]), pch=15)
  #axis(2, at=ncol(msfit$xstd):1, labels= colnames(msfit$xstd), las=1)
  #mtext('Specification number',side=1,line=2)
}


multi_bsca= function(b, ms=NULL, conversion=NULL, treatments=NULL, y.scale='Coefficient', x.name='Treatment', y.name='Outcome', 
                     y.wrap=21, x.wrap=10, legend.position='right', add.ate=FALSE, ate.name='ATE', 
                     add.global.ate=FALSE, global.ate.name='  Global ATE', shapes=c(15:19, 3:10)) {
  # convert b to mombf if necessary
  if (typeof(b[[1]]) == 'list') b = lapply(b, as.mscoef)
  
  # add treatment averages
  if (add.ate) {
    if (is.null(ms)) stop('need to provide model fits ms if adding averages')
    xvars = grep(paste0('^`?',paste(treatments, collapse='|'),'((?!:).)*$'),rownames(b[[1]]), perl=TRUE)
    for (i in seq_along(b)) {
      ate = getATE(ms[[i]], xvars)
      b[[i]] = rbind(b[[i]], ate)
      rownames(b[[i]])[nrow(b[[i]])] = ate.name
    }
  }
  
  # prepare models df
  models = lapply(b, function(r) {
    df = as.data.frame(r)
    names(df)[2:3] = c('conf.low', 'conf.high')
    df$term = rownames(r)
    return(df)
  })
  models = bind_rows(models, .id = 'model')
  models$term = gsub('`', '', models$term)
  
  # limit treatments
  if (!is.null(treatments)) {
    models = models[models$term %in% c(treatments, ate.name), ]
  }
  
  # turn into ordered factors
  for (v in c('model', 'term')) {
    models[[v]] = factor(models[[v]], levels = models[[v]][!duplicated(models[[v]])])
  }
  
  # add outcome average
  if (add.global.ate) {
    global.ate = models[models['term'] == ate.name, !names(models) == 'model'] %>% group_by(term) %>% summarise_all(mean)
    global.ate
  }
  
  # convert coefficients
  if (!is.null(conversion)) {
    for (v in c('estimate', 'conf.low', 'conf.high')) {
      models[[v]] = conversion(models[[v]])
      if (add.global.ate) global.ate[[v]] = conversion(global.ate[[v]])
    }
  }
  
  # wrap text
  models$term = str_wrap(models$term, x.wrap)
  models$model = str_wrap(models$model, y.wrap)
  
  # turn into ordered factors (again)
  for (v in c('model', 'term')) {
    models[[v]] = factor(models[[v]], levels = models[[v]][!duplicated(models[[v]])])
  }
  
  g = ggplot(models, aes(term, estimate, 
                         shape=model, linetype=model, colour=model)) + 
    geom_point(size=3, position = position_dodge(width = 0.5)) +
    geom_errorbar(aes(ymin=conf.low, ymax=conf.high), 
                  lwd=1, width=0, position = position_dodge(width = 0.5)) + 
    theme_bw() + theme(legend.position=legend.position) + 
    ylab(y.scale) + xlab(x.name) + scale_shape_manual(name=y.name, values=shapes) + 
    scale_colour_discrete(name=y.name) + scale_linetype_discrete(name=y.name)
  if (add.global.ate) {
    h = global.ate$estimate
    g = g + geom_hline(yintercept=h, lty=2, colour='grey50', size=0.25) + 
      annotate('text', x=0, y=h, label=global.ate.name, hjust=0, vjust=-1, colour='grey50', size=8/.pt) 
    # geom_text(x=0.2, y=0.5, label='global.ate.name', hjust=0,  colour='black', size=(theme_get()$text$size*theme_get()$legend.text$size/.pt), fontface='plain', family='Times New Roman')
  }
  return(g)
}



#### Convert Results #### 

as.msfit = function(basfit, priors=DEFAULT_PRIORS) {
  x=basfit$X; y=basfit$Y
  center=TRUE; scale=TRUE
  p= ncol(x); n= length(y)
  
  #Standardize (y,x) to mean 0 and variance 1
  if (!is.vector(y)) { y <- as.double(as.vector(y)) } else { y <- as.double(y) }
  if (!is.matrix(x)) x <- as.matrix(x)
  mx= colMeans(x); sx= sqrt(colMeans(x^2) - mx^2) * sqrt(n/(n-1))
  ct= (sx==0)
  if (any(is.na(ct))) stop('x contains NAs, this is currently not supported, please remove the NAs')
  if (sum(ct)>1) stop('There are >1 constant columns in x (e.g. two intercepts)')
  if (!center) { my=0; mx= rep(0,p) } else { my= mean(y) }
  if (!scale) { sy=1; sx= rep(1,p) } else { sy= sd(y) }
  ystd= (y-my)/sy; xstd= x; xstd[,!ct]= t((t(x[,!ct]) - mx[!ct])/sx[!ct])
  stdconstants= rbind(c(my,sy),cbind(mx,sx)); colnames(stdconstants)= c('shift','scale')
  stdconstants[1,] = c(0, 1)
  
  bascoef = coef(basfit)
  family = ifelse(exists("family", where=basfit), paste(basfit$family$family, basfit$family$link), 'normal')
  modelid = sapply(basfit$which, function(x) {x=x+1; paste(x,collapse=',')})
  postmean = bascoef$conditionalmean; rownames(postmean) = modelid
  postvar = bascoef$conditionalsd^2; rownames(postvar) = modelid
  
  fit = new(
    "msfit",
    list(
      postSample=matrix(logical(), nrow=0, ncol=p),
      margpp=basfit$probne0,
      # postMode=postMode,
      # postModeProb=postModeProb,
      postProb=bascoef$postprobs,
      postmean=postmean,
      postvar=postvar,
      family=family,
      p=p,
      enumerate=TRUE,
      priors=list(
        priorCoef=bicprior(), 
        priorGroup=bicprior(), 
        priorDelta=modelbbprior(1,1), 
        priorConstraints=modelbbprior(1,1), 
        priorVar=igprior(0,0), 
        priorSkew=NULL
      ),
      ystd=y,
      xstd=xstd,
      # groups=groups,
      # constraints=constraints,
      stdconstants=stdconstants,
      outcometype='glm',
      # call=call,
      models=data.frame(
        modelid=modelid,
        family=rep(family,length(basfit$postprobs)),
        pp=basfit$postprobs
      )
    )
  )
}


as.mscoef = function(bascoef, alpha=0.05) {
  if (alpha != 0.05) stop('Set alpha=0.05 or change colnames in as.mscoef')
  # Warning: wrong interval; confint throws:
  # Error in .HPDinterval(betas, prob = level) : obj must have nsamp > 1
  mean = bascoef$postmean; sd = bascoef$postsd; c = qnorm(1-alpha/2)
  coef = cbind(mean, mean-c*sd, mean+c*sd, bascoef$probne0)
  rownames(coef) = bascoef[['namesx']]
  colnames(coef) = c('estimate', '2.5%', '97.5%', 'margpp')
  return(coef)
}


#### Display Results #### 

#Extract average treatment effect. Point estimate, 95% interval and posterior probability that the ATE is non-zero
getATE= function(ms, xvars, fun='mean', alpha=0.05) {
  b = rnlp(msfit=ms, niter=10^5)  # sample coefficients
  if (fun=='mean') {
    b.ate = rowMeans(b[,xvars])
  } else if (fun=='median') {
    b.ate = apply(b[,xvars],1,'median')
  } else { stop("fun must be either 'mean' or 'median'") }
  return(c(estimate=mean(b.ate), ci=quantile(b.ate, probs=c(alpha/2, 1-alpha/2)), margpp= mean(b.ate != 0)))
}

#Extract MLE and CI's from glm fit
getci <- function(fit, alpha=.05) {
   b= summary(fit)$coef
   ans= cbind(b[,1], b[,1] + qnorm(alpha/2) * b[,2], b[,1] - qnorm(alpha/2) * b[,2])
   return(ans)
}


#Estimate OR from BMA output associated to the variable treat and values treatvals
#Example of use (for numeric variables): getOR(b, treat='Electronic', treatvals=1:6)
#Alternatively (for factors): getOR(b, treat='TV')
getOR= function(b, treat, treatvals=NULL, digits=NULL) {
  # convert b to mombf if necessary
  if (typeof(b) == 'list') b = as.mscoef(b)
  # match when variable name starts with "[`]treat", but not when followed by a ":" (interaction term)
  btreat= b[grep(paste0('^`?',treat,'((?!:).)*$'),rownames(b), perl=TRUE),1:3]
  if (is.null(treatvals)) {
    treateffect= btreat
  } else {
    treateffect= treatvals * matrix(rep(btreat,length(treatvals)),ncol=3,byrow=TRUE)
  }
  or= exp(treateffect)
  colnames(or)= c('OR','CI.low','CI.up')
  if (is.null(treatvals)) {
    rownames(or)= rownames(btreat)
  } else {
    rownames(or)= paste(treat, treatvals)
  }
  if (!is.null(digits)) {
    or= round(or, digits=digits)
  }
  return(or)
}


#Plot bivariate cross-tabulation along with univariate marginals
plotxtab= function(x, y, xlab='', ylab='') {
    require(MASS)
    require(RColorBrewer)
    h1= table(x)
    h2= table(y)
    top= max(c(max(h1),max(h2)))
    k= table(x,y)
    # margins
    oldpar= par()
    par(mar=c(3,3,1,1))
    layout(matrix(c(2,0,1,3),2,2,byrow=T),c(3,1), c(1,3))
    # colors
    rf= colorRampPalette(rev(brewer.pal(11,'Spectral')))
    r= rf(32)
    # plot
    image(k, col=r, xaxt='n', yaxt='n')
    par(mar=c(0,2,1,0))
    barplot(h1, axes=F, ylim=c(0, top), space=0, col='red', main=xlab)
    par(mar=c(2,0,0.5,1))
    barplot(h2, axes=F, xlim=c(0, top), space=0, col='red', horiz=T, main=ylab)
}


#### Formatting #### 

# Output an objects representation, with every line preceeded by four spaces 
# to create a preformatted block in "as-is" chunks
cat.asis <- function(x, pre='') {
  cat(pre, paste(sprintf('    %s\n', capture.output(show(x))), collapse=''))
}


# Create model matrix from posterior probabilities
as.matrix.pp = function(pp, nummodels, numvars){
  modelids = as.character(pp[1:nummodels,'modelid'])
  # expand modelids into boolean matrix
  models = matrix(FALSE, nrow=nummodels, numvars)
  for (r in 1:nummodels) {
    for (c in strsplit(modelids[r], ',')) {
      models[r, as.numeric(c)] = TRUE
    }
  }
  return(models)
}



#### Simulate data #### 

simmultivdata= function(n, beta, alpha, Sigma, seed) {
  # Simulate multiple outcomes
  # - n: sample size
  # - beta: p x r matrix with regression coefficients for p treatments on r outcomes
  # - alpha: q x r matrix with regression coefficients for q controls on r outcomes
  # - Sigma: r x r error correlation matrix
  # - seed: random number generator seed, passed on to set.seed
  require(mvtnorm)
  if (!is.matrix(beta) | !is.matrix(alpha) | !is.matrix(Sigma)) stop("beta, alpha and Sigma must be matrices")
  if (ncol(beta) != ncol(alpha)) stop("beta and alpha must have the same number of columns")
  if (ncol(beta) != ncol(Sigma)) stop("beta and Sigma must have the same number of columns")
  set.seed(seed)
  p= nrow(beta); q= ifelse(missing(alpha), 0, nrow(alpha)); r= ncol(beta)
  #Simulate the value of treatment(s) and controls
  x = matrix(rnorm(n*p),nrow=n) 
  xnames= paste('x',1:nrow(beta),sep='')
  if (q>0) { #if there are control covariates
    z = rowMeans(x) + matrix(rnorm(n*q),nrow=n) #correlated with treatments
    znames= paste('z',1:nrow(alpha),sep='')
  }
  #Simulate the values of the outcomes
  e= rmvnorm(n, sigma=Sigma)
  y= matrix(NA, nrow=n, ncol=r)
  for (i in 1:r) {
    if (q>0) {
      y[,i] = x %*% matrix(beta[,i],ncol=1) + z %*% matrix(alpha[,i],ncol=1) + e[,i]
    } else {
      y[,i] = x %*% matrix(beta[,i],ncol=1) + e[,i]
    }
  }
  #Format output
  ynames= paste('y',1:r,sep='')
  if (q>0) {
    ans= data.frame(y,1,x,z)
    colnames(ans)= c(ynames,'Intercept',xnames,znames)
  } else {
    ans= data.frame(y,1,x)
    colnames(ans)= c(ynames,'Intercept',xnames)
  }
  return(ans)
}
