## code yield gaps NPC

#################################################
### packages loading and set working directory ##
#################################################

rm(list=ls())
gc()

library(data.table)
library(plyr)
library(reshape2)
library(plm)
library(msm)
library(Matrix)
library(sandwich)
library(tseries)
library(car)
library(lmtest)

wd ="/home/NET1/baldoed/yield_gaps_capri"
setwd(wd)


#############################################
##################### functions definition ##
#############################################


## function to compute variance inflation factor
get_vif_organic = function(ff=ff,dataset=est_yield$data,ORGANICclass=2){
  
  ## inputs
  # ff: function to be estimated 
  # dataset: data to perform estimation
  # ORGANICclass: class to be used as dependent variable
  
  ## output
  # VIF score
  
  part1 = paste('organic',ORGANICclass,sep='')
  dataset[[part1]]=0
  dataset[[part1]][dataset$ORGANIC==ORGANICclass]=1
  part2 = as.character(ff)[[3]]
  part2 = gsub('as.factor\\(ORGANIC\\) \\+ ','',part2)
  
  ff_org = as.formula(paste(part1,part2,sep='~'))
  mod_org = lm(ff_org,dataset)
  sum_mod_org = summary(mod_org)
  vif_organic = 1/(1-sum_mod_org$r.squared)
  names(vif_organic) = paste('vif_',part1,sep='')
  return(vif_organic)
}


## function to perform delta method to obtain standard errors of yield gaps (exp(b)-1)*100
get_sd_coefs = function(coefs,model){
  
  ### inputs
  # coefs: estimated coefficients
  # model: estimated model
  
  ### output
  # diagonal of matrix of standard errors of transformed coefficients
  
  bb = coefs[!is.na(coefs)]
  vv = vcov(model)
  vv = vv[which(rownames(vv) %in% names(bb)), which(rownames(vv) %in% names(bb))]
  
  if(length(bb)==4){vard = deltamethod(list(~exp(x1),~100*(exp(x2)-1),~100*(exp(x3)-1),~100*(exp(x4)-1)),bb,vv,ses=FALSE)}
  if(length(bb)==3){vard = deltamethod(list(~exp(x1),~100*(exp(x2)-1),~100*(exp(x3)-1)),bb,vv,ses=FALSE)}
  if(length(bb)==2){vard = deltamethod(list(~exp(x1),~100*(exp(x2)-1)),bb,vv,ses=FALSE)}
  if(length(bb)==1 & grepl('Int',names(bb)[1])){vard = deltamethod(list(~exp(x1)),bb,vv,ses=FALSE)}
  if(length(bb)==1 & grepl('ORG',names(bb)[1])){vard = deltamethod(list(~100*(exp(x1)-1)),bb,vv,ses=FALSE)}
  
  out = try(sqrt(diag(vard)))
  
  return(out)
}


## function to estimate gaps
estimate_model_delta = function(dataset = temp_data, VAR, fixed_effect=FALSE, ff=NULL,NMIN,planB=TRUE,STEP=FALSE){
  ### inputs
  # dataset: the dataset containing yield and other info
  # VAR: the yield variable
  # fixed-effect (not implemented): estimate the model including individual fixed-effects
  # ff: regression formula
  # NMIN: minimum number of observations to perform estimations
  # planB: if TREU estimate a model with a reduced number of dummy variables. When regions are flat, it removes dummies related to altitude class
  # STEP (not implemented): use the standard stepwise model selection function implemented in base R

  ### outputs
  # estimates: yield gap estimates for all classes of organic farms (exp(b)-1)*100
  # pvals: pvalues of estimates of regression coefficients (b)
  # nobs: number of observations used for estimation
  # model: model estimates
  # data: the final dataset used for estimation
  # comment: comment and details about estimation results
  # nobs_organic: number of organic farms used for estimation
  # sds: standard errors of transformed coefficients (exp(b)-1)*100
  # sds0: standard errors of original coefficients (b)
  # MODEL: full output of estimated model
  
  require(plm)
  
  dataset_prod =  subset(dataset,dataset[[VAR]] > 0 & is.finite(dataset[[VAR]]))
  dataset_prod0 =  subset(dataset,!dataset[[VAR]]>0 )
  
  if(nrow(dataset_prod)< NMIN | length(as.character(levels(as.factor(dataset_prod$ORGANIC))))<2 ){
    
    coefs_yprod_final = rep(NA,4)
    sum_yprod = 'error'
    pvals_yprod = rep(NA,4)
    sds_yprod = rep(NA,4)
    sds_yprod0 = rep(NA,4)
    mats_yprod=rep(NA,4)
    number_obs = nrow(dataset_prod)
    number_org = sum(as.character(dataset_prod$ORGANIC)=='2')
    comment = 'too few observations/no organic variability'
    MODEL = NA
    MODEL_PANEL = NA
    vif_organic=NA

  }else{
    
    dataset_prod$ORGANIC = as.factor(dataset_prod$ORGANIC)
    
    if(fixed_effect){
      
      
      ff_fe1 = NULL
      mod_yprod =try(plm(ff_fe1,
                         data=dataset_prod,
                         model='within',
                         effect='individual',index=c('id','YEAR')))
      
    }else{
      
      if(STEP){dataset_prod = subset(dataset_prod,is.finite(share_irr))}
      
      part1 = as.character(ff)[[2]]
      part1 = gsub('dat.prodc.geo','dataset_prod',part1)
      part1 = gsub('temp_data','dataset_prod',part1)
      part2 = as.character(ff)[[3]]
      if(length(as.character(levels(as.factor(dataset_prod$IRRSYS))))<2 & grepl('IRRSYS',part2)){
        part2 = gsub('\\+ as.factor\\(IRRSYS\\) \\+','+',part2)
      }else{
        part2 = part2
      }
      
      ff = as.formula(paste(part1,part2,sep='~'))
      
      ## model yields (by FADN region)
      mod_yprod =try(lm(ff,data=dataset_prod))
      mod_yprod_panel = try(plm(ff,data=dataset_prod,model='pooling',index=c('YEAR','ID')))
      
      if(class(mod_yprod) != 'try-error'){
        vif_organic = get_vif_organic(ff,dataset_prod,ORGANICclass=2)
      }
      
      if(STEP & !class(mod_yprod)=='try-error'){ mod_yprod = try(step(mod_yprod,scope = list(lower = ~ as.factor(ORGANIC)+as.factor(YEAR))))}
      
    }
    
    
    
    sum_yprod = summary(mod_yprod)
    number_obs = try(length(sum_yprod$residuals));if(class(number_obs)=='try-error'){number_obs=NA}
    number_org = sum(as.character(dataset_prod$ORGANIC)=='2')
    
    if(class(mod_yprod) == 'try-error'){
      
      coefs_yprod_final = rep(NA,4)
      pvals_yprod = rep(NA,4)
      sds_yprod = rep(NA,4)
      sds_yprod0 = rep(NA,4)
      mats_yprod = rep(NA,4)
      number_obs = nrow(dataset_prod)
      number_org = sum(as.character(dataset_prod$ORGANIC)=='2')
      sum_yprod = 'error'
      comment = 'no model could be estimated'
      MODEL = NA
      MODEL_PANEL = NA
      vif_organic  = NA

      if(planB){ 
        if(fixed_effect){
          
          ff_fe1=NULL
          mod_yprod1 =try(plm(ff_fe1,
                              data=dataset_prod,
                              model='within',
                              effect='individual',index=c('id','YEAR')))
          
        }else{
          
          if(STEP){dataset_prod = subset(dataset_prod,is.finite(share_irr))}
          
          part1 = as.character(ff)[[2]]
          part2 = as.character(ff)[[3]]
          if(length(as.character(levels(as.factor(dataset_prod$IRRSYS))))<2 & grepl('IRRSYS',part2)){
            part2 = gsub('\\+ as.factor\\(ALTITUDE\\) \\+ as.factor\\(ANC3\\) \\+ as.factor\\(IRRSYS\\) \\+','+',part2)
          }else{
            part2 = gsub('\\+ as.factor\\(ALTITUDE\\) \\+ as.factor\\(ANC3\\) \\+','+',part2)
          }
          
          
          ff1 = as.formula(paste(part1,part2,sep='~'))
          
          mod_yprod1 =try(lm(ff1,data=dataset_prod))
          mod_yprod_panel1 = try(plm(ff1,data=dataset_prod,model='pooling',index=c('YEAR','ID')))
          if(class(mod_yprod1) != 'try-error'){
            vif_organic = get_vif_organic(ff1,dataset_prod,ORGANICclass=2)
          }
          
          
          if(STEP & !class(mod_yprod1)=='try-error'){ mod_yprod1 = try(step(mod_yprod1,scope = list(lower = ~ as.factor(ORGANIC)+as.factor(YEAR))))}
          
          
        }
        
        
        sum_yprod = summary(mod_yprod1)
        number_obs = try(length(sum_yprod$residuals));if(class(number_obs)=='try-error'){number_obs=NA}
        number_org = sum(as.character(dataset_prod$ORGANIC)=='2')
        
        
        if(class(mod_yprod1) == 'try-error'){
          
          coefs_yprod_final = rep(NA,4)
          pvals_yprod = rep(NA,4)
          sds_yprod = rep(NA,4)
          sds_yprod0 = rep(NA,4)
          mats_yprod = rep(NA,4)
          number_obs = nrow(dataset_prod)
          number_org = sum(as.character(dataset_prod$ORGANIC)=='2')
          sum_yprod = 'error'
          comment = 'no model could be estimated'
          MODEL = NA
          MODEL_PANEL = NA
          vif_organic = NA

        }else{
          
          MODEL = mod_yprod1
          MODEL_PANEL = mod_yprod_panel1
          
          coef_yprod  = coef(mod_yprod1)
          coef_yprod1 = coef_yprod[grep('\\(Intercept\\)',names(coef_yprod))];if(length(coef_yprod1)==0){coef_yprod1 =NA}
          coef_yprod2 = coef_yprod[grep('\\(ORGANIC\\)2',names(coef_yprod))];if(length(coef_yprod2)==0){coef_yprod2 =NA}
          coef_yprod3 = coef_yprod[grep('\\(ORGANIC\\)3',names(coef_yprod))];if(length(coef_yprod3)==0){coef_yprod3 =NA}
          coef_yprod4 = coef_yprod[grep('\\(ORGANIC\\)4',names(coef_yprod))];if(length(coef_yprod4)==0){coef_yprod4 =NA}
          
          
          mat_yprod  = try(summary(mod_yprod1)$coefficients)
          mat_yprod1 = mat_yprod[grep('\\(Intercept\\)',rownames(mat_yprod)),2];if(length(mat_yprod1)==0){mat_yprod1 =NA}
          mat_yprod2 = mat_yprod[grep('\\(ORGANIC\\)2',rownames(mat_yprod)),2];if(length(mat_yprod2)==0){mat_yprod2 =NA}
          mat_yprod3 = mat_yprod[grep('\\(ORGANIC\\)3',rownames(mat_yprod)),2];if(length(mat_yprod3)==0){mat_yprod3 =NA}
          mat_yprod4 = mat_yprod[grep('\\(ORGANIC\\)4',rownames(mat_yprod)),2];if(length(mat_yprod4)==0){mat_yprod4 =NA}
          mats_yprod = c(mat_yprod1,mat_yprod2,mat_yprod3,mat_yprod4)
          
          
          pval_yprod1 = sum_yprod$coefficients[grep('\\(Intercept\\)',rownames(sum_yprod$coefficients)),4];if(length(pval_yprod1 )==0){pval_yprod1 =NA}
          pval_yprod2 = sum_yprod$coefficients[grep('\\(ORGANIC\\)2',rownames(sum_yprod$coefficients)),4];if(length(pval_yprod2)==0){pval_yprod2 =NA}
          pval_yprod3 = sum_yprod$coefficients[grep('\\(ORGANIC\\)3',rownames(sum_yprod$coefficients)),4];if(length(pval_yprod3)==0){pval_yprod3 =NA}
          pval_yprod4 = sum_yprod$coefficients[grep('\\(ORGANIC\\)4',rownames(sum_yprod$coefficients)),4];if(length( pval_yprod4 )==0){pval_yprod4 =NA}
          pvals_yprod  = c(pval_yprod1,pval_yprod2,pval_yprod3,pval_yprod4)
          
          
          coefs_yprod  = c(coef_yprod1, coef_yprod2,coef_yprod3,coef_yprod4)
          coef_yprod_corr = (100*(exp(coefs_yprod[2:length(coefs_yprod)])-1))
          coefs_yprod_final = c(exp(coef_yprod1),coef_yprod_corr)
          
          which_sds = which(!is.na(coefs_yprod))
          sds_yprod = rep(NA,4)
          sds_yprod0 = try(get_sd_coefs(coefs=coefs_yprod,model=mod_yprod1))
          if(class(sds_yprod0)=='try-error'){
            mats_yprod = rep(NA,4)
            NULL
          }else{
            sds_yprod[which_sds] = sds_yprod0
            mats_yprod = mats_yprod
          }
          
          comment = 'excluding ALTITUDE and ANC3'
        }
      }
      
      
      
      
    }else{
      
      MODEL = mod_yprod
      MODEL_PANEL = mod_yprod_panel
      
      coef_yprod  = coef(mod_yprod)
      coef_yprod1 = coef_yprod[grep('\\(Intercept\\)',names(coef_yprod))];if(length(coef_yprod1)==0){coef_yprod1 =NA}
      coef_yprod2 = coef_yprod[grep('\\(ORGANIC\\)2',names(coef_yprod))];if(length(coef_yprod2)==0){coef_yprod2 =NA}
      coef_yprod3 = coef_yprod[grep('\\(ORGANIC\\)3',names(coef_yprod))];if(length(coef_yprod3)==0){coef_yprod3 =NA}
      coef_yprod4 = coef_yprod[grep('\\(ORGANIC\\)4',names(coef_yprod))];if(length(coef_yprod4)==0){coef_yprod4 =NA}
      
      mat_yprod  = try(summary(mod_yprod)$coefficients)
      mat_yprod1 = mat_yprod[grep('\\(Intercept\\)',rownames(mat_yprod)),2];if(length(mat_yprod1)==0){mat_yprod1 =NA}
      mat_yprod2 = mat_yprod[grep('\\(ORGANIC\\)2',rownames(mat_yprod)),2];if(length(mat_yprod2)==0){mat_yprod2 =NA}
      mat_yprod3 = mat_yprod[grep('\\(ORGANIC\\)3',rownames(mat_yprod)),2];if(length(mat_yprod3)==0){mat_yprod3 =NA}
      mat_yprod4 = mat_yprod[grep('\\(ORGANIC\\)4',rownames(mat_yprod)),2];if(length(mat_yprod4)==0){mat_yprod4 =NA}
      mats_yprod = c(mat_yprod1,mat_yprod2,mat_yprod3,mat_yprod4)
      
      
      pval_yprod1 = sum_yprod$coefficients[grep('\\(Intercept\\)',rownames(sum_yprod$coefficients)),4];if(length(pval_yprod1 )==0){pval_yprod1 =NA}
      pval_yprod2 = sum_yprod$coefficients[grep('\\(ORGANIC\\)2',rownames(sum_yprod$coefficients)),4];if(length(pval_yprod2)==0){pval_yprod2 =NA}
      pval_yprod3 = sum_yprod$coefficients[grep('\\(ORGANIC\\)3',rownames(sum_yprod$coefficients)),4];if(length(pval_yprod3)==0){pval_yprod3 =NA}
      pval_yprod4 = sum_yprod$coefficients[grep('\\(ORGANIC\\)4',rownames(sum_yprod$coefficients)),4];if(length( pval_yprod4 )==0){pval_yprod4 =NA}
      pvals_yprod  = c(pval_yprod1,pval_yprod2,pval_yprod3,pval_yprod4)
      
      coefs_yprod  = c(coef_yprod1, coef_yprod2,coef_yprod3,coef_yprod4)
      coef_yprod_corr = (100*(exp(coefs_yprod[2:length(coefs_yprod)])-1))
      coefs_yprod_final = c(exp(coef_yprod1),coef_yprod_corr)
      
      which_sds = which(!is.na(coefs_yprod))
      sds_yprod = rep(NA,4)
      sds_yprod0 = try(get_sd_coefs(coefs=coefs_yprod,model=mod_yprod))
      if(class(sds_yprod0)=='try-error'){
        mats_yprod = rep(NA,4)
        NULL
      }else{
        sds_yprod[which_sds] = sds_yprod0
        mats_yprod = mats_yprod
      }
      
      comment='ok'
    }
  }
  
  out = list(coefs_yprod_final,pvals_yprod,number_obs,sum_yprod,dataset_prod,comment,number_org,sds_yprod,mats_yprod,MODEL,MODEL_PANEL,vif_organic)
  names(out) = c('estimates','pvals','nobs','model','data','comment','nobs_organic','sds','sds0','MODEL','MODEL_PANEL','VIF_ORGANIC')
  return(out)
}




#############################################
############# data loading and cleaning #####
#############################################

tableAB = as.data.frame(readRDS('tableAB.RDS')) ## load general information on the holdings
tableSE = as.data.frame(readRDS('tableSE.RDS')) ## load standard results on the holdings
tableI = as.data.frame(readRDS('tableI.RDS'))   ## load production data of the holdings
productions_sel =  c('CBRL','CFODMZ','CPOT','CWHTC','COAT','CLEG','CMZ','CRYE','CWHTD','CPEA')

### define outliers of production data: define outliers at COUNTRY level
clean_outliers = TRUE
geo_clean = 'COUNTRY'
if(clean_outliers){
  res_cleanI = list()
  idx = 1
  for(j in 1:length(productions_sel)){
    temp = subset(tableI,CROP==productions_sel[j])
    if(nrow(temp)>0){
      set_geos_clean = as.character(unique(temp[[geo_clean]]))
      for(k in 1:length(set_geos_clean)){
        print(paste(productions_sel[j],set_geos_clean[k]))
        temp1 =  subset(temp,temp[[geo_clean]]==set_geos_clean[k]) 
        qy = quantile(temp1[['YIELD']],seq(0,1,by=0.01),na.rm=TRUE)
        temp1$outlier1perc = 0
        temp1$outlier1perc[temp1$YIELD < qy[2] | temp1$YIELD > qy[100]]=1
        temp1$outlier5perc = 0
        temp1$outlier5perc[temp1$YIELD < qy[6] | temp1$YIELD > qy[96]]=1
        temp1$outlier10perc = 0
        temp1$outlier10perc[temp1$YIELD < qy[11] | temp1$YIELD > qy[91]]=1
        res_cleanI[[idx]] = temp1
        idx = idx + 1
      }
    }
  }
  tableIoutliers = do.call('rbind',res_cleanI)
}else{
  tableIoutliers = tableI
}




#############################################
##################### yield gap estimation ##
#############################################
## define intersetion between REGION and TF8 class
tableAB$REGION.TF8 = paste(tableAB$REGION, tableAB$TF8,sep='.')

## recode the ORGANIC variable
tableAB$ORGANICold = tableAB$ORGANIC
tableAB$ORGANIC[tableAB$ORGANIC==4] = 3

## compute Shannon index and join with tableSE
tablea = tableI[,c('ID','YEAR','CROP','A')]
tableac = reshape2::dcast(YEAR+ID~CROP,value.var='A',sum,data=tablea)
tableac[,3:ncol(tableac)] = tableac[,3:ncol(tableac)]/rowSums(tableac[,3:ncol(tableac)])
tableac$SHANNON_CROP_DIVERSITY = apply(tableac[,3:ncol(tableac)],1,function(x) -sum(x[x>0]*log(x[x>0])))
tableSE = join(tableSE,tableac[,c('YEAR','ID','SHANNON_CROP_DIVERSITY')],by=c('YEAR','ID'),type='left',match='first')


## define parameters for the yield gap loop
geo.agg = 'REGION.TF8' # obtain yield gaps for every FADN region and TF8 combination
nmin = 20          # minimum number of observations to perform estimations
outliers_remove = '1perc' # define outliers: bottom 1%, top 1% of yields at national level


## start loop for yield gap estimation: for every product perform estimations for every FADN region
list_coefs = list()
list_crop = list()
IDX = 1
for(i in 1:length(productions_sel)){
  
  
  dat.prod = subset(tableIoutliers[,c('YEAR','ID','CROP','YIELD','PRICE','outlier1perc','outlier5perc','outlier10perc')],CROP==productions_sel[i])
  if(outliers_remove == '1perc'){dat.prod = subset(dat.prod, outlier1perc < 1)} ## remove outliers

  dat.prodm = reshape2::melt(dat.prod,id.vars=c('YEAR','ID','CROP'))
  dat.prodc = reshape2::dcast(dat.prodm,YEAR+ID~CROP+variable,value.var='value',sum,na.rm=TRUE)
  dat.prodc = join(dat.prodc,tableAB[,c('ID','YEAR','NUTS3','COUNTRY','REGION.TF8',
                                        'NUTS1','NUTS2','NUTS0','REGION',
                                        'ORGANIC','ALTITUDE','ANC3','IRRSYS','TF14')],by=c('YEAR','ID'),type='left',match='first')
  dat.prodc = join(dat.prodc,tableSE[,c('ID','YEAR','SE030','SE025','SE015','SE010','SE295',
                                        'SE135','SE131','SE275','SE284','SE206',
                                        'SE605','SE410','SE441','SE490',
                                        'SE080','SE285','SE290','SE455',
                                        'SHANNON_CROP_DIVERSITY')],by=c('YEAR','ID'),type='left',match='first')

  dat.prodc$share_rented = dat.prodc$SE030/dat.prodc$SE025
  dat.prodc$share_unpaid = dat.prodc$SE015/dat.prodc$SE010
  dat.prodc$share_crops = dat.prodc$SE135/dat.prodc$SE131
  dat.prodc$share_crops[is.na(dat.prodc$share_crops)]=0
  
  dat.prodc$ratio_KL = (dat.prodc$SE455/dat.prodc$SE010)/1000
  dat.prodc$subsidies_grossincome = dat.prodc$SE605/dat.prodc$SE410
  dat.prodc$subsidies_grossincome[is.na(dat.prodc$subsidies_grossincome)]= 0
  dat.prodc$assets_ha = dat.prodc$SE441/dat.prodc$SE025
  dat.prodc$assets_to_liabilities = dat.prodc$SE441/(dat.prodc$SE441+dat.prodc$SE490)
  dat.prodc$ORGANIC = as.character(dat.prodc$ORGANIC)
  dat.prodc$REGION = as.character(dat.prodc$REGION)
  dat.prodc$ALTITUDE = as.character(dat.prodc$ALTITUDE)
  dat.prodc$ANC3 = as.character(dat.prodc$ANC3)
  dat.prodc$IRRSYS = as.character(dat.prodc$IRRSYS)
  dat.prodc$TF14 = as.character(dat.prodc$TF14)
  dat.prodc$fertilizers_ha = dat.prodc$SE295/dat.prodc$SE025
  dat.prodc = subset(dat.prodc,is.finite(subsidies_grossincome))
  dat.prodc$seeds_ha = (dat.prodc$SE285+dat.prodc$SE290)/dat.prodc$SE025
  dat.prodc$SE080[is.na(dat.prodc$SE080)] = 0
  
  varyield = paste(productions_sel[i],'_YIELD',sep='')
  
  # Get the list of FADN regions-TF8 combinations. 
  # Add also the list of REGION whose gaps will be
  # used in the case REGION-TF8 combinations do not
  # have enough observations
  set_geoagg0 = unique(dat.prodc[[geo.agg]])
  set_regions = unique(dat.prodc[['REGION']])
  set_geoagg = c(set_geoagg0,set_regions)
    
  print('########################')
  print(productions_sel[i])
  
  ## nested loop for yield gap estimation: split the sample by region and perform estimations
  for(j in 1:length(set_geoagg)){
    
    setname = paste(productions_sel[i],set_geoagg[j],sep='|')  
    
    print(set_geoagg[j])
    
    if(grepl('\\.',set_geoagg[j])){  
      dat.prodc.geo = subset(dat.prodc,dat.prodc[[geo.agg]]==set_geoagg[j])
    }else{
      dat.prodc.geo = subset(dat.prodc,dat.prodc$REGION==set_geoagg[j])
    }
    
    
    if(nrow(dat.prodc.geo)<nmin | length(as.character(levels(as.factor(dat.prodc.geo$ORGANIC))))<2){
      
      list_crop[[IDX]] = list()
      list_crop[[IDX]][[1]] = data.frame(var = rep('yield',4),
                                         geo = rep(set_geoagg[j],4),
                                         crop = rep(productions_sel[i],4),
                                         variable = rep(NA,4),
                                         estimate = rep(NA,4),
                                         sds = rep(NA,4),
                                         pvals = rep(NA,4),
                                         pvals.robust = rep(NA,4),
                                         nobs = rep(NA,4),
                                         nobs_organic = rep(NA,4),
                                         comment = rep(NA,4),
                                         normality.pval.sh = rep(NA,4),
                                         heterosk.pval = rep(NA,4),
                                         share.leverage.larger2mean = rep(NA,4),
                                         share.leverage.larger3mean =  rep(NA,4),
                                         sercorr.pval = rep(NA,4),
                                         r2 = rep(NA,4),
                                         adjr2 = rep(NA,4),
                                         mae = rep(NA,4),
                                         rmse =rep(NA,4),
                                         vif_organic = rep(NA,4))
      
      res_stat  = data.frame(geo = rep(set_geoagg[j],1),
                             crop = rep(productions_sel[i],1),
                             y = NA,
                             yhat = NA,
                             err = NA,
                             ORGANIC = NA)
      
      list_crop[[IDX]][[2]] = res_stat
      list_coefs[[IDX]] = list(setname,NA)
      names(list_coefs)[IDX] = setname
      
      IDX = IDX + 1
      
    }else{
      
      ### MODEL ESTIMATION
      ######################
      
      # define regression equation
      ff_yield =  log(dat.prodc.geo[[varyield]]) ~
        as.factor(ORGANIC)+
        as.factor(YEAR)+
        as.factor(NUTS3)+
        as.factor(TF14)+
        as.factor(ALTITUDE)+
        as.factor(ANC3)+
        as.factor(IRRSYS)+
        I(SE080/SE025)+               
        share_unpaid+
        subsidies_grossincome+
        share_rented+
        assets_to_liabilities+
        I(SE025^2)+
        SE025+
        ratio_KL+
        assets_ha+
        seeds_ha+
        fertilizers_ha+
        share_crops+
        SHANNON_CROP_DIVERSITY
      
      # perform estimation of yield gaps
      ##################################
      est_yield  = estimate_model_delta(dataset=dat.prodc.geo,
                                        VAR=varyield,
                                        fixed_effect = FALSE,
                                        ff=ff_yield,
                                        NMIN=nmin,
                                        planB=TRUE,
                                        STEP=FALSE)
      
 
      # robust estimation (autocorrelation & heteroskedasticity/ heteroskedasticity only)
      ###################################################################################
      
      # adjust standard errors: heteroskedasticity and autocorrelation
      robust=TRUE
      if(!is.null(robust) & sum(is.na(est_yield$estimates)) < 4){
      
        # heteroskedasticity and autocorrelation 
        mod2 = try(coeftest(est_yield$MODEL_PANEL, vcov = vcovHC(est_yield$MODEL_PANEL, method= 'arellano',cluster = c("group"))))
        # heteroskedasticity only
        if(class(mod2)=='try-error'){mod2 = try(coeftest(est_yield$MODEL_PANEL, vcov = vcovHC(est_yield$MODEL_PANEL, method= 'white1')))}
        if(class(mod2)=='try-error'){mod2 = try(coeftest(est_yield$MODEL, vcov = vcovHC(est_yield$MODEL, type= 'HC0')))}
        
        if(class(mod2)=='try-error'){
          pvals_robust = rep(NA,4)
        }else{
          mod3 = as.data.frame(mod2[rownames(mod2) %in% c('(Intercept)','as.factor(ORGANIC)2','as.factor(ORGANIC)3','as.factor(ORGANIC)4'),])
          mod3$variable = rownames(mod3)
          org_coefs = data.frame(variable = c('(Intercept)','as.factor(ORGANIC)2','as.factor(ORGANIC)3','as.factor(ORGANIC)4'))
          org_coefs = join(org_coefs,mod3,by='variable',type='left',match='first')
          pvals_robust = org_coefs[,grep('Pr',names(org_coefs))]
        }
        
       }
      
      ## residuals diagnostics
      ########################
      
      if(sum(!is.na(est_yield$estimate))>0){
        if(abs(sum(residuals(est_yield$MODEL),na.rm=TRUE))> 0){ 
          
          ## sumary, R2 and R2 adjusted
          sum_mod = summary(est_yield$MODEL)     
          r2 = sum_mod$r.squared
          r2c = sum_mod$adj.r.squared
          
          ## Variance inflation factor
          if(is.finite(est_yield$VIF_ORGANIC)){
            vif_organic = est_yield$VIF_ORGANIC
          }else{
            vif_model = try(vif(est_yield$MODEL))
            if(class(vif_model)[1] != 'try-error'){
              vif_organic = vif_model[rownames(vif_model)=='as.factor(ORGANIC)','GVIF']
            }else{
              vif_organic = NA
            }
          }
          
          
          ## get residuals
          residui = residuals(est_yield$MODEL)
      
          ## normality with Shapiro test
          norm_errors = try(shapiro.test(residui))
          if(class(norm_errors)[1] == 'try-error'){
            norm.pval.sh = NA
          }else{
            norm.pval.sh =  c(norm_errors$p.value)
          }
      
          ## leverage (2/3 times higher than average leverage value)
          lev = hatvalues(est_yield$MODEL)
          share_leverage2mean = mean(lev>2*mean(lev,na.rm=TRUE))
          share_leverage3mean = mean(lev>3*mean(lev,na.rm=TRUE))
      
          ## heteroskedasticity with Breusch-Pagan test
          heterosk_errors = bptest(est_yield$MODEL)
          heterosk.pval = heterosk_errors$p.value
      
          ## serial correlation with Breusch-Godfrey test
          sercorr = try(pbgtest(est_yield$MODEL_PANEL, order = 2))
          if(class(sercorr)[1] == 'try-error'){
           sercorr.pval = NA
          }else{
           sercorr.pval = sercorr$p.value
          }
      
          
          # when estimation fails, set parameters of diagnostics to NA
          }else{
            norm.pval.sh = NA
            heterosk.pval = NA
            share_leverage2mean = NA
            share_leverage3mean = NA
            sercorr.pval = NA
            r2 = NA
            adjr2 = NA
            vif_organic = NA
            pvals_robust = NA
          }
        
        # when estimation fails, set parameters of diagnostics to NA
        }else{
          norm.pval.sh = NA
          heterosk.pval = NA
          share_leverage2mean = NA
          share_leverage3mean = NA
          sercorr.pval = NA
          r2 = NA
          adjr2 = NA
          vif_organic = NA
          pvals_robust = NA
        }
      
  
      ### COLLECT MODEL RESULTS
      #########################
      
     ## collect prediction errors
     y = try(est_yield$MODEL$model[[1]])
     if(class(y)[1] != 'try-error'){
       yhat = try(predict(est_yield$MODEL))
       diff_yield = try((yhat-y)/y)
       ORG = est_yield$MODEL$model[['as.factor(ORGANIC)']]
       mae = mean(abs(y-yhat))
       rmse = sqrt(mean((y-yhat)^2))
     }else{
       y = NA
       yhat=NA
       diff_yield=NA
       ORG=NA
       mae = NA
       rmse = NA
     }
     
     ## collect estimates and model statistics
     res_yield = data.frame(var = rep('yield',4),
                            geo = rep(set_geoagg[j],4),
                            crop = rep(productions_sel[i],4),
                            variable = c('(Intercept)','as.factor(ORGANIC)2','as.factor(ORGANIC)3','as.factor(ORGANIC)4'),
                            estimate = est_yield$estimates,
                            sds = est_yield$sds,
                            pvals = est_yield$pval,
                            pvals.robust = pvals_robust,
                            nobs = est_yield$nobs,
                            nobs_organic = est_yield$nobs_organic,
                            comment = est_yield$comment,
                            normality.pval.sh = norm.pval.sh,
                            heterosk.pval = heterosk.pval,
                            share.leverage.larger2mean = share_leverage2mean,
                            share.leverage.larger3mean = share_leverage3mean,
                            sercorr.pval = sercorr.pval,
                            r2 = r2,
                            adjr2 = r2c,
                            mae= mae,
                            rmse = rmse,
                            vif_organic = vif_organic)
     
     rownames(res_yield) = 1:nrow(res_yield)
     list_crop[[IDX]] = list()
     list_crop[[IDX]][[1]] = res_yield
     
     
     res_stat  = data.frame(geo = rep(set_geoagg[j],length(y)),
                            crop = rep(productions_sel[i],length(y)),
                            y = y,
                            yhat = yhat,
                            err = diff_yield,
                            ORGANIC = ORG)
     list_crop[[IDX]][[2]] = res_stat
      
     ## collect full model estimates
     list_coefs[[IDX]] = list(setname,
                               est_yield$model)
     names(list_coefs)[IDX] = setname
      
     IDX = IDX + 1
    }
  }
}



##########################################################
## Aggregate REGION-TF8 yield gaps at REGION level #######
##########################################################

tab_yields = do.call('rbind',lapply(list_crop,function(x) x[[1]]))
## select only relevant organic gaps and ensure compliance with confidentiality agreement with DG AGRI (N organic > 15).
tab_yields_comp = subset(tab_yields,!is.na(estimate) & variable == 'as.factor(ORGANIC)2' & nobs_organic>15) 
tab_yields_comp$geo[!grepl('\\.',tab_yields_comp$geo)] = paste(tab_yields_comp$geo[!grepl('\\.',tab_yields_comp$geo)],'.ALL',sep='')
tab_yields_comp$pvals.robust[is.na(tab_yields_comp$pvals.robust)] = tab_yields_comp$pvals[is.na(tab_yields_comp$pvals.robust)]


### aggregate at regional level
tab_yields_comp$SET = paste(tab_yields_comp$crop,tab_yields_comp$geo,sep='|')
tab_yields_comp$REGION = unlist(lapply(strsplit(tab_yields_comp$geo,split='\\.'),function(x) x[[1]]))
tab_yields_comp$TF8 = unlist(lapply(strsplit(tab_yields_comp$geo,split='\\.'),function(x) x[[2]]))
tab_yields_comp$REGION[tab_yields_comp$TF8 == 'ALL'] = paste(tab_yields_comp$REGION[tab_yields_comp$TF8 == 'ALL'],'.ALL',sep='')

med_reg_gaps = aggregate(tab_yields_comp$estimate,list(tab_yields_comp$REGION,tab_yields_comp$crop),median,na.rm=TRUE)
names(med_reg_gaps) = c('REGION','crop','median_estimate')
avg_reg_numest = aggregate(tab_yields_comp$estimate,list(tab_yields_comp$REGION,tab_yields_comp$crop),length)
names(avg_reg_numest) = c('REGION','crop','num_gap_estimates')
avg_reg_numb = aggregate(tab_yields_comp[,c('nobs','nobs_organic')],list(tab_yields_comp$REGION,tab_yields_comp$crop),sum,na.rm=TRUE)
names(avg_reg_numb) = c('REGION','crop','sum_nobs','sum_nobs_organic')
avg_reg = join(med_reg_gaps,avg_reg_numb,by=c('REGION','crop'),type='left',match='first')
avg_reg = join(avg_reg,avg_reg_numest,by=c('REGION','crop'),type='left',match='first')
avg_reg_comp = subset(avg_reg,sum_nobs_organic>15)
avg_reg_comp$SPEC = 'TF8'
avg_reg_comp$SPEC[grep('.ALL',avg_reg_comp$REGION)] = 'ALL'
avg_reg_comp$REGION = unlist(lapply(strsplit(avg_reg_comp$REGION ,split='\\.'),function(x) x[[1]]))
avg_reg_comp = avg_reg_comp[,c('REGION','SPEC','crop','median_estimate','sum_nobs','sum_nobs_organic','num_gap_estimates')]



##########################################################
###################### Export yield gap estimates ########
##########################################################

write.csv(avg_reg_comp,file='yield_gap_estimates.csv',row.names=FALSE)
