########################################################################################

#TITLE:     Study 1 - GP Model 2 Imputation Evaluation 

#AUTHOR:    Matthew Barstead, M.S.

#CONTACT:   mbarstead@gmail.com

########################################################################################
#Ensuring random starts for two different threads on same computer... 
set.seed(as.numeric(Sys.time()))

#Original Processing: 

dat.time.orig<-read.csv('/media/mbarsted/HDD1/Sim_data_2/Raw_coded.csv',
                        stringsAsFactors = F)

#cleaning up extra rows
dat.time.orig<-dat.time.orig[-91:-nrow(dat.time.orig),]

#Taking only cases with at least 180 continuous seconds worth of data: 
dat.time.orig<-dat.time.orig[dat.time.orig$Window_Length>=180,]

table(dat.time.orig$Population) #not bad - about 50/50 
#40 sections of adult data 
#39 sections of child data

#Converting "Adult" to "adult" and "Child" to "Child" 
dat.time.orig$Population<-ifelse(dat.time.orig$Population=='Adult', 'adult', 'child')

wd<-'/media/mbarsted/HDD1/Sim_data_2'

#Bringing in file information worksheet
file_info<-read.csv(paste0(wd, '/Physio_Tracking - Sheet1.csv'), stringsAsFactors = F)

#Converting ID to common structure:
file_info$ID<-substr(file_info$Video_1, start=1, stop=6)

#Task timing file
task.time<-read.table(paste0(wd, '/Timing_File_forVizEdit.txt'), 
                      header = T, sep = '\t')
colnames(task.time)[1]<-'ID'

#Creating a test case folder - will run five imputations under each condition on 5 different files: 
#All test files will be with adult files. 
data.folder<-paste0(wd, '/Good_3min')

file.names<-list.files(data.folder)

#Getting appropriate timing values for sections of "good" data
DF.segment<-merge(dat.time.orig, task.time[,1:2], by='ID')

#Will not be able to use case 051_T2 as there appears to be a task/timing issue there (loses 1 child case)
DF.segment<-DF.segment[DF.segment$ID!="051_T2",]
DF.segment$Start.adj<-DF.segment$Start+DF.segment$Video
DF.segment$Stop.adj<-DF.segment$Stop+DF.segment$Video
DF.segment$time.min<-DF.segment$Start.adj+31
DF.segment$time.max<-DF.segment$Stop.adj-31

DF.segment<-merge(DF.segment, file_info[,c(1,14:17)], by='ID')

####################################################################################################################
####################################################################################################################
#Functions for identifying heart beats
#===========================================================================
#Function 1 - Finding Peakings Using Specified bandwidth: 
findpeaks <- function (x, m = 3){
  shape <- diff(sign(diff(x, na.pad = FALSE)))
  pks <- sapply(which(shape < 0), FUN = function(i){
    z <- i - m + 1
    z <- ifelse(z > 0, z, 1)
    w <- i + m + 1
    w <- ifelse(w < length(x), w, length(x))
    if(all(x[c(z : i, (i + 2) : w)] <= x[i + 1])) return(i + 1) else return(numeric(0))
  })
  pks <- unlist(pks)
  pks
}
#===========================================================================
#Function 2 - Summing IBIs from Raw PPG file: 
time.sum<-function(x){
  Z<-rep(NA, length(x))
  for(i in 1:length(x)){
    Z[i]<-ifelse(i==1, x[i], x[i]-x[i-1])
  }
  return(Z)  
}
#===========================================================================
#Function 2b - Summing Time from IBIs
IBI.sum<-function(x){
  Z<-rep(NA, length(x))
  for(i in 1:length(x)){
    Z[i]<-sum(x[1:i])
  }
  return(Z)
}
#===========================================================================
#Function 3 - Iterative function for getting IBIs
iter.IBI<-function(x, ds=2000){
  #browser()
  require(psych)
  s<-round(seq(round(ds/50), round(ds/2), length.out = 200))
  Z<-data.frame(rep(NA, length(s)), 
                rep(NA, length(s)), 
                rep(NA, length(s)), 
                rep(NA, length(s)), 
                rep(NA, length(s)), 
                rep(NA, length(s)))
  for(i in 1:length(s)){
    IBI<-findpeaks(x, s[i])
    time<-time.sum(IBI)/ds
    Z[i,1]<-s[i]
    Z[i,2]<-sd(time)
    Z[i,3]<-max(time)-min(time)
    Z[i,4]<-rmssd(time)
    Z[i,5]<-mean(acf(time, lag.max = length(time)/20, plot = F)$acf)
    Z[i,6]<-s[i]/ds
  }
  colnames(Z)<-c('BW', 'SD', 'Range', 'RMSSD', 'AC', 'BW(s)')
  Z<-Z[order(Z$RMSSD, decreasing = F),]
  IBI.fin<-findpeaks(x, m=Z[1,1])-1
  IBI.fin<-IBI.fin/ds
  IBI.done<-time.sum(IBI.fin)
  IBI.comp<-list(IBI.done, Z)
  names(IBI.comp)<-c('IBI.done', 'Z')
  return(IBI.comp)
}
#===========================================================================
#Function 4 - Obtaining Time Values for IBI
sum.rev<-function(x){
  Z<-rep(NA, length(x))
  for(i in 1:length(x)){
    Z[i]<-ifelse(i==1, x[i], sum(x[1:(i-1)])+x[i])
  }
  return(Z)
}

###+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#setting up models:
#See dissertation for original specifications (slightly modifying based on JH comments during dissertation defense)
chains<-2
impute.min<-2
impute.max<-8
warmup<-c(2000)
n.iter<-c(2500)
Data.sample.Hz<-c(2000)
impute.sample.Hz<-c(125, 250)
GP.sample.Hz<-c(4, 8, 12)
impute.fac<-c(1,2,3)

#Going to need to bring in header information
#For now all test cases have 15 rows of header information: 

######################################################################################################################
library(ggplot2)
library(signal)
library(psych)
library(bayesplot)
library(MCMCvis)
library(rstan)
library(rstanarm)
library(astsa)
library(benchmarkme)
options(mc.cores=parallel::detectCores())
rstan_options(auto_write = TRUE)
#Run the complete simulation from here
#Basing this on 5 simulations each to start: 
n.sims<-20

#Basic simulation vectors
ID.vec<-vector()
impute.ID<-vector()
CPU.vec<-vector()
Cores.vec<-vector()
RAM.vec<-vector()
pop.vec<-vector()
start.vec<-vector()
stop.vec<-vector()
IBI.rmssd.15<-vector()
IBI.rmssd.30<-vector()
IBI.rmssd.45<-vector()
IBI.rmssd.60<-vector()
IBI.MR.rmssd.15<-vector()
IBI.MR.rmssd.30<-vector()
IBI.MR.rmssd.45<-vector()
IBI.MR.rmssd.60<-vector()
IBI.HD.rmssd.15<-vector()
IBI.HD.rmssd.30<-vector()
IBI.HD.rmssd.45<-vector()
IBI.HD.rmssd.60<-vector()
IBI.GP.rmssd.15<-vector()
IBI.GP.rmssd.30<-vector()
IBI.GP.rmssd.45<-vector()
IBI.GP.rmssd.60<-vector()
N.IBIs.vec<-vector()
runtime<-vector()
warmup.vec<-vector()
iter.vec<-vector()
Hz.vec<-vector()
impute.Hz.vec<-vector()
GP.Hz.vec<-vector()
impute.fac.vec<-vector()
impute.start<-vector()
impute.end<-vector()
impute.tot<-vector()

#Parameter Means
HR.mean.vec<-vector()
Resp.mean.vec<-vector()

#Parameter sd's 
HR.SD.vec<-vector()
Resp.SD.vec<-vector() 

#Parameter R-hats
HR.Rhat.vec<-vector()
Resp.Rhat.vec<-vector()

#Parameter N_eff
HR.N_eff.vec<-vector()
Resp.N_eff.vec<-vector()

#Getting ML Convergence test
a1.ML.conv.vec<-vector()
a2.ML.conv.vec<-vector()
a3.ML.conv.vec<-vector()
r1.ML.conv.vec<-vector()
r2.ML.conv.vec<-vector()
r3.ML.conv.vec<-vector()
r4.ML.conv.vec<-vector()
r5.ML.conv.vec<-vector()

#Creating folder structure to allow for quick saving of Files
out.folder<-paste0(wd, '/Output')
IBI.orig.folder<-paste0(out.folder, '/IBI_Files/Original')
IBI.MR.folder<-paste0(out.folder, '/IBI_Files/Mean_Replacement')
IBI.HD.folder<-paste0(out.folder, '/IBI_Files/Hot_Deck')
IBI.impute.folder<-paste0(out.folder, '/IBI_Files/Imputed')
Peak.summary.folder<-paste0(out.folder, '/Peak_Detection_Summaries')
graphics.folder<-paste0(out.folder, '/Graphics')
Bayes.graphics<-paste0(graphics.folder, '/Bayes_Graphics')
model.folder<-paste0(out.folder, '/Model_Summaries')
diagnostics.folder<-paste0(out.folder, '/Diagnostics')

#Main simulation program:
for(i in 1:n.sims){
  #browser()
  #Getting the specifics for each run: 
  impute.window.temp<-runif(1, impute.min, impute.max)
  row.select.temp<-sample(x=1:nrow(DF.segment), size = 1)
  
  #Obtaining temporary dataset
  dat.temp<-DF.segment[row.select.temp,]
  header.temp<-ifelse(dat.temp$Population=='child', dat.temp$Child_Header, dat.temp$Parent_Header)
  col.select.temp<-ifelse(dat.temp$Population=='child', dat.temp$Child_Column, dat.temp$Parent_Column)
  
  #Selecting imputation boundaries
  impute.LB.temp<-runif(1, dat.temp$time.min, dat.temp$time.max)
  impute.UB.temp<-impute.LB.temp+impute.window.temp
  
  #getting PPG data
  PPG<-read.table(paste0(data.folder, '/', dat.temp$ID, '.txt'), 
                  skip = header.temp, 
                  header = F, 
                  sep = '\t'
  )
  
  PPG<-PPG[,col.select.temp]
  
  PPG<-data.frame(PPG, 
                  Time = seq(from = 0, by = .0005, length.out = length(PPG))
  )
  
  #Selecting target window for the case
  PPG.temp<-PPG[PPG$Time>=dat.temp$Start.adj & PPG$Time<=dat.temp$Stop.adj,]
  
  #Cleaning signal - de-spiking and smoothing heart rate
  PPG.temp$PPG<-as.numeric(smooth(PPG.temp$PPG))
  PPG.temp$PPG<-smooth.spline(PPG.temp$PPG, nknots = 10000)$y
  PPG.temp$PPG<-PPG.temp$PPG-predict(lm(PPG~Time, data = PPG.temp))
  
  #Obtaining original IBI values for section:
  IBI.orig.temp<-iter.IBI(PPG.temp$PPG, ds=2000)
  
  write.table(round(IBI.orig.temp$IBI.done, digits = 4), 
              paste0(IBI.orig.folder, 
                     '/', 
                     paste(dat.temp$ID, 
                           dat.temp$Population, 
                           dat.temp$Start, 
                           dat.temp$Stop, 
                           round(impute.LB.temp, digits = 2),
                           sep = '_'), 
                     '.txt'), 
              row.names = F
  )
  
  write.table(round(head(IBI.orig.temp$Z, n=20), digits = 3), 
              paste0(Peak.summary.folder, 
                     '/', 
                     paste(dat.temp$ID,
                           dat.temp$Population,
                           'original',
                           dat.temp$Start, 
                           dat.temp$Stop, 
                           round(impute.LB.temp, digits = 2),
                           sep = '_'),
                     '.txt'
              ), 
              sep = '\t', 
              row.names = F
  )
  
  #Mean replacement strategy - Will replace the number of IBIs in the affected range with mean values
  #Step 1 - Restoring accurate timing file
  
  IBI.orig.time<-sum.rev(IBI.orig.temp$IBI.done)+dat.temp$Start.adj
  
  #Running very simple mean imputation (not accounting for time-series nature of the data)
  IBI.mean.replace<-IBI.orig.temp$IBI.done
  IBI.mean.replace[IBI.orig.time>impute.LB.temp & IBI.orig.time<impute.UB.temp]<-rep(mean(IBI.orig.temp$IBI.done))
  
  #removing first and last IBI value (to make more equivalent comparison with imputed data set)
  write.table(round(IBI.mean.replace[c(-1, -length(IBI.mean.replace))], digits = 4), 
              paste0(IBI.MR.folder, 
                     '/', 
                     paste(dat.temp$ID, 
                           dat.temp$Population,
                           'Mean_replace',
                           dat.temp$Start, 
                           dat.temp$Stop, 
                           round(impute.LB.temp, digits = 2),
                           sep = '_'), 
                     '.txt'), 
              row.names = F
  )
  
  #Running hot deck imputation (does not account for time-series nature of the data)
  IBI.sample.vals<-IBI.orig.temp$IBI.done[IBI.orig.time<impute.LB.temp | IBI.orig.time>impute.UB.temp]
  IBI.hotdeck.replace<-IBI.orig.temp$IBI.done
  N.IBIs<-length(IBI.hotdeck.replace[IBI.orig.time>impute.LB.temp & IBI.orig.time<impute.UB.temp])
  IBI.hotdeck.replace[IBI.orig.time>impute.LB.temp & IBI.orig.time<impute.UB.temp]<-sample(x=IBI.sample.vals, 
                                                                                           size = N.IBIs, 
                                                                                           replace = T)
  write.table(round(IBI.hotdeck.replace[c(-1, -length(IBI.hotdeck.replace))], digits = 4), 
              paste0(IBI.HD.folder, 
                     '/', 
                     paste(dat.temp$ID, 
                           dat.temp$Population,
                           'HotDeck_replace',
                           dat.temp$Start, 
                           dat.temp$Stop,
                           round(impute.LB.temp, digits = 2),
                           sep = '_'), 
                     '.txt'), 
              row.names = F
  )
  
  for(j in 1:length(chains)){
    chain.temp<-chains[j]
    for(k in 1:length(warmup)){
      warmup.temp<-warmup[k]
      n.iter.temp<-n.iter[k]
      for(l in 1:length(Data.sample.Hz)){
        Hz.temp<-Data.sample.Hz[l]
        for(m in 1:length(impute.sample.Hz)){
          Hz.impute.temp<-impute.sample.Hz[m]
          for(n in 1:length(GP.sample.Hz)){
            Hz.GP.temp<-GP.sample.Hz[n]
            for(o in 1:length(impute.fac)){
              impute.fac.temp<-impute.fac[o]
              
              #Downsampling (if appropriate)
              select.int<-2000/Hz.temp
              PPG.impute<-PPG.temp[seq(1, nrow(PPG.temp), select.int),]
              IBI.impute<-iter.IBI(PPG.impute$PPG, ds=Hz.temp)
              
              #Getting some priors for the model
              mu_HP<-mean(1/IBI.impute$IBI.done[c(-1,-length(IBI.impute$IBI.done))], na.rm = T)
              sigma_HP<-sd(1/IBI.impute$IBI.done[c(-1,-length(IBI.impute$IBI.done))], na.rm = T)
              
              #Estimating Respiration parameter & respiration SD if population age is child
              if(dat.temp$Population == 'child'){
                spec<-mvspec(PPG.temp, 
                             spans = c(7,7), 
                             taper=.1, 
                             demean = T, 
                             log='no', 
                             plot = F)
                min.R<-20/60/Hz.temp
                max.R<-30/60/Hz.temp
                
                spec.trunc<-data.frame(freq=spec$freq[spec$freq>=min.R&spec$freq<=max.R],
                                       spec=spec$spec[spec$freq>=min.R&spec$freq<=max.R])
                spec.trunc$prob<-spec.trunc$spec/sum(spec.trunc$spec)
                tmp.dist<-sample(spec.trunc$freq, size = 10000, replace = T, prob = spec.trunc$prob)*Hz.temp
                mu_R<-mean(tmp.dist)
                sigma_R<-sd(tmp.dist)
              }
              
              #Estimating Respiration parameter & respiration SD if population age is adult (frequency filter changes)
              if(dat.temp$Population == 'adult'){
                spec<-mvspec(PPG.temp, 
                             spans = c(7,7), 
                             taper=.1, 
                             demean = T, 
                             log='no', 
                             plot = F)
                min.R<-12/60/Hz.temp
                max.R<-20/60/Hz.temp
                
                spec.trunc<-data.frame(freq=spec$freq[spec$freq>=min.R&spec$freq<=max.R],
                                       spec=spec$spec[spec$freq>=min.R&spec$freq<=max.R])
                spec.trunc$prob<-spec.trunc$spec/sum(spec.trunc$spec)
                tmp.dist<-sample(spec.trunc$freq, size = 10000, replace = T, prob = spec.trunc$prob)*Hz.temp
                mu_R<-mean(tmp.dist)
                sigma_R<-sd(tmp.dist)
              }
              
              init.list<-list()
              for(x in 1:chain.temp){
                init.list[[x]]<-list(mu_HR = mu_HP, 
                                     sigma_HR = sigma_HP, 
                                     mu_R = mu_R, 
                                     sigma_R = sigma_R)
              }
              
              if(Hz.temp>=Hz.impute.temp){
                #First selecting a series of time values to impute back in at the appropriate sampling rate
                select.int2<-Hz.temp/Hz.impute.temp
                PPG.impute2<-PPG.impute[seq(1, nrow(PPG.impute), select.int2),]
                PPG.impute.pred<-PPG.impute2[PPG.impute2$Time>impute.LB.temp & PPG.impute2$Time<impute.UB.temp,]
                
                #Getting what will be Prediction values:
                Xp<-PPG.impute.pred$Time
                N2<-length(PPG.impute.pred$Time)
                
                #Acquiring "good data"
                min.TIME2<-min(Xp)
                max.TIME2<-max(Xp)
                time.span<-max.TIME2-min.TIME2
                Y.vals<-rbind(PPG.impute2[PPG.impute2$Time>min.TIME2-impute.fac.temp*time.span & PPG.impute2$Time<min.TIME2,],
                              PPG.impute2[PPG.impute2$Time>max.TIME2 & PPG.impute2$Time<max.TIME2+impute.fac.temp*time.span,])
                Y.vals<-na.omit(Y.vals)
                tot.Y.vals<-length(Y.vals[,1])
                sel.Y.vals<-round(seq(1, tot.Y.vals, length.out = round(tot.Y.vals/Hz.impute.temp*Hz.GP.temp)))
                sel.Y.vals<-unique(sel.Y.vals)
                Y<-Y.vals$PPG[sel.Y.vals]
                X<-Y.vals$Time[sel.Y.vals]
                N1<-length(X)
                
                #Running Model - there is no need for respiration priors
                dat<-list(N1=N1,
                          N2=N2,
                          Xp=Xp,
                          X=X,
                          Y=Y,
                          mu_HR=mu_HP,
                          sigma_HR=sigma_HP, 
                          mu_R = mu_R, 
                          sigma_R = sigma_R
                )
                
                #Obtaining estimates for hyper-parameters using ML approach
                opt_model<-stan_model(file=paste0(wd, '/Stan_code/GP_2_opt.stan'))
                
                #need to find some way to assess convergence...
                #Identifying values for hyperparameters (will not attempt to estimate a posterior for these values)
                #Identifying values for hyperparameters (will not attempt to estimate a posterior for these values)
                print('Obtaining first ML estimate for hyperparameters')
                opt_fit1<-NULL
                attempt<-1
                while(is.null(opt_fit1) & attempt <=15){
                  print(paste('Convergence attempt', attempt, 'out of 15'))
                  attempt<- attempt + 1
                  try(
                    opt_fit1<-optimizing(opt_model, 
                                         data=dat, 
                                         init=list(mu_HR = mu_HP, mu_R = mu_R), 
                                         seed=sample(1:5000, size = 1), iter = 10000)
                  )
                }
                
                print('Obtaining second ML estimate for hyperparameters')
                opt_fit2<-NULL
                attempt<-1
                while(is.null(opt_fit2) & attempt <=15){
                  print(paste('Convergence attempt', attempt, 'out of 15'))
                  attempt<- attempt + 1
                  try(
                    opt_fit2<-optimizing(opt_model, 
                                         data=dat, 
                                         init=list(mu_HR = mu_HP, mu_R = mu_R), 
                                         seed=sample(1:5000, size = 1), iter = 10000)
                  )
                }
                
                if(!is.null(opt_fit1)){
                  alpha1.1 <- opt_fit1$par['a1']
                  alpha2.1 <- opt_fit1$par['a2']
                  alpha3.1 <- opt_fit1$par['a3']
                  rho1.1 <- opt_fit1$par['r1']
                  rho2.1 <- opt_fit1$par['r2']
                  rho3.1 <- opt_fit1$par['r3']
                  rho4.1 <- opt_fit1$par['r4']
                  rho5.1 <- opt_fit1$par['r5']
                }
                
                if(!is.null(opt_fit2)){
                  alpha1.2 <- opt_fit2$par['a1']
                  alpha2.2 <- opt_fit2$par['a2']
                  alpha3.2 <- opt_fit2$par['a3']
                  rho1.2 <- opt_fit2$par['r1']
                  rho2.2 <- opt_fit2$par['r2']
                  rho3.2 <- opt_fit2$par['r3']
                  rho4.2 <- opt_fit2$par['r4']
                  rho5.2 <- opt_fit2$par['r5']
                }                
                
                if(!is.null(opt_fit1) & !is.null(opt_fit2)){
                  #Simple Alpha convergence check
                  a1.ML.conv.vec<-c(a1.ML.conv.vec,(abs(alpha1.1-alpha1.2)<.0001))
                  a2.ML.conv.vec<-c(a2.ML.conv.vec,(abs(alpha2.1-alpha2.2)<.0001))
                  a3.ML.conv.vec<-c(a3.ML.conv.vec,(abs(alpha3.1-alpha3.2)<.0001))
                  
                  #Simple rho covergence check
                  r1.ML.conv.vec<-c(r1.ML.conv.vec,(abs(rho1.1-rho1.2)<.0001))
                  r2.ML.conv.vec<-c(r2.ML.conv.vec,(abs(rho2.1-rho2.2)<.0001))
                  r3.ML.conv.vec<-c(r3.ML.conv.vec,(abs(rho3.1-rho3.2)<.0001))
                  r4.ML.conv.vec<-c(r4.ML.conv.vec,(abs(rho4.1-rho4.2)<.0001))
                  r5.ML.conv.vec<-c(r5.ML.conv.vec,(abs(rho5.1-rho5.2)<.0001))
                }
                
                else{
                  #Simple Alpha convergence check
                  a1.ML.conv.vec<-c(a1.ML.conv.vec, NA)
                  a2.ML.conv.vec<-c(a2.ML.conv.vec, NA)
                  a3.ML.conv.vec<-c(a3.ML.conv.vec, NA)
                  
                  #Simple rho covergence check
                  r1.ML.conv.vec<-c(r1.ML.conv.vec, NA)
                  r2.ML.conv.vec<-c(r2.ML.conv.vec, NA)
                  r3.ML.conv.vec<-c(r3.ML.conv.vec, NA)
                  r4.ML.conv.vec<-c(r4.ML.conv.vec, NA)
                  r5.ML.conv.vec<-c(r5.ML.conv.vec, NA)
                }
                
                if(is.null(opt_fit1) & !is.null(opt_fit2)){
                  alpha1.1<-alpha1.2
                  alpha2.1<-alpha2.2
                  alpha3.1<-alpha3.2
                  rho1.1<-rho1.2
                  rho2.1<-rho2.2
                  rho3.1<-rho3.2
                  rho4.1<-rho4.2
                  rho5.1<-rho5.2
                }
              
                pars.to.monitor<-c('HR', 'R', 'Ypred')
                start.time.temp<-Sys.time()
                
                start.time.temp<-Sys.time()
                
                dat.opt<-list(N1=N1,
                          N2=N2,
                          Xp=Xp,
                          X=X,
                          Y=Y,
                          mu_HR=mu_HP,
                          sigma_HR=sigma_HP, 
                          mu_R = mu_R, 
                          sigma_R = sigma_R,
                          a1 = alpha1.1, 
                          a2 = alpha2.1,
                          a3 = alpha3.1,
                          r1 = rho1.1, 
                          r2 = rho2.1, 
                          r3 = rho3.1,
                          r4 = rho4.1,
                          r5 = rho5.1
                )
                
                fit.stan<-stan(file=paste0(wd, '/Stan_code/GP_2.stan'),
                               data = dat.opt, 
                               init = init.list,
                               warmup = warmup.temp,
                               iter = n.iter.temp,
                               refresh=100,
                               chains = chain.temp,
                               pars = pars.to.monitor,
                               control = list(adapt_delta = .95, 
                                              max_treedepth = 15)
                )
                
                run_time<-round(Sys.time()-start.time.temp, digits = 3)
                units(run_time)<-'mins'
                run_time<-as.numeric(run_time)
                
                y_pred<-rstan::extract(fit.stan, 'Ypred')
                PPG.new<-colMeans(y_pred$Ypred)
                PPG.DF.new<-data.frame(PPG = PPG.new,
                                       Time = Xp)
                
                GP.summary<-as.data.frame(summary(fit.stan, pars = pars.to.monitor[-length(pars.to.monitor)],
                                                  probs = c(.1, .9))$summary)
                
                write.table(round(GP.summary, digits = 3), 
                            paste0(model.folder, 
                                   '/', 
                                   paste(dat.temp$ID, 
                                         dat.temp$Population,
                                         dat.temp$Segment,
                                         'GP_impute',
                                         round(impute.LB.temp, digits = 2), 
                                         round(impute.UB.temp, digits = 2),
                                         Hz.temp, 
                                         Hz.impute.temp, 
                                         Hz.GP.temp,
                                         impute.fac.temp,
                                         chain.temp, 
                                         warmup.temp, 
                                         n.iter.temp,
                                         round(impute.LB.temp, digits = 2),
                                         sep = '_'), 
                                   '.txt'), 
                            row.names = T, 
                            sep='\t'
                )
                
                g1<-ggplot()+
                  geom_line(data=PPG.DF.new, aes(x=Time, y = PPG), color = 'red', lty='dashed')+
                  geom_line(data=PPG.impute.pred, aes(x=Time, y=PPG), color = 'black')+
                  ggtitle(paste(dat.temp$ID, dat.temp$Population, dat.temp$Segment,
                                'File Hz =', Hz.temp, 'GP Dataset Hz =', Hz.impute.temp, 
                                'Model Sampling Rate =', Hz.GP.temp, 'Chains =', chain.temp, 
                                'Warmup =', warmup.temp, 'Iter =', n.iter.temp))
                
                ggsave(filename = paste0(graphics.folder, paste('/',dat.temp$ID, dat.temp$Population, dat.temp$Segment,
                                                                Hz.temp, Hz.impute.temp, Hz.GP.temp, impute.fac.temp, 
                                                                chain.temp, warmup.temp, n.iter.temp, round(impute.LB.temp, digits = 2),
                                                                sep = '_'), '.png'), 
                       plot = g1, 
                       width = 11, 
                       height = 8, 
                       units = 'in', 
                       dpi = 300, 
                       device = 'png')
                
                g2<-traceplot(fit.stan, pars='HR', inc_warmup=T)+
                  ggtitle(paste('HR param', dat.temp$ID, dat.temp$Population, dat.temp$Segment,
                                'File Hz =', Hz.temp, 'GP Dataset Hz =', Hz.impute.temp, 
                                'Model Sampling Rate =', Hz.GP.temp, '\n', 'Chains =', chain.temp, 
                                'Warmup =', warmup.temp, 'Iter =', n.iter.temp))
                
                
                ggsave(filename = paste0(Bayes.graphics, paste('/',dat.temp$ID, dat.temp$Population, dat.temp$Segment,
                                                               Hz.temp, Hz.impute.temp, Hz.GP.temp, impute.fac.temp, 
                                                               chain.temp, warmup.temp, n.iter.temp, round(impute.LB.temp, digits = 2),
                                                               sep = '_'), '.png'),
                       plot = g2, 
                       width = 11, 
                       height = 8, 
                       units = 'in', 
                       dpi = 300, 
                       device = 'png')
                
                PPG.impute.fin<-PPG.impute2
                PPG.impute.fin$PPG[PPG.impute.fin$Time>impute.LB.temp & PPG.impute.fin$Time<impute.UB.temp]<-PPG.new
                
                IBI.GP.impute.temp<-iter.IBI(PPG.impute.fin$PPG, ds=Hz.impute.temp)$IBI.done
                IBI.GP.impute.time<-sum.rev(IBI.GP.impute.temp)+dat.temp$Start.adj
                write.table(round(IBI.GP.impute.temp[c(-1, -length(IBI.GP.impute.temp))], digits = 4), 
                            paste0(IBI.impute.folder, 
                                   '/', 
                                   paste(dat.temp$ID, 
                                         dat.temp$Population,
                                         dat.temp$Segment,
                                         'GP_impute',
                                         round(impute.LB.temp, digits = 2), 
                                         round(impute.UB.temp, digits = 2),
                                         Hz.temp, 
                                         Hz.impute.temp, 
                                         Hz.GP.temp,
                                         impute.fac.temp,
                                         chain.temp, 
                                         warmup.temp, 
                                         n.iter.temp,
                                         round(impute.LB.temp, digits = 2),
                                         sep = '_'), 
                                   '.txt'), 
                            row.names = F
                )
                
                HR.est<-rstan::extract(fit.stan, 'HR')
                Resp.est<-rstan::extract(fit.stan, 'R')
      
                #Basic simulation vectors
                ID.vec<-c(ID.vec, dat.temp$ID)
                impute.ID<-c(impute.ID, paste('Segment', i, sep = '_'))
                CPU.vec<-c(CPU.vec, benchmarkme::get_cpu()$model_name)
                Cores.vec<-c(Cores.vec, parallel::detectCores())
                RAM.vec<-c(RAM.vec, paste(round(benchmarkme::get_ram()/1073741824), 'GB'))
                pop.vec<-c(pop.vec, dat.temp$Population)
                N.IBIs.vec<-c(N.IBIs.vec, N.IBIs)
                start.vec<-c(start.vec, dat.temp$Start)
                stop.vec<-c(stop.vec, dat.temp$Stop)
                IBI.GP.rmssd.15<-c(IBI.GP.rmssd.15, rmssd(IBI.GP.impute.temp[IBI.GP.impute.time>impute.LB.temp-7.5 & 
                                                                               IBI.GP.impute.time<impute.UB.temp+7.5]))
                IBI.GP.rmssd.30<-c(IBI.GP.rmssd.30, rmssd(IBI.GP.impute.temp[IBI.GP.impute.time>impute.LB.temp-15 & 
                                                                               IBI.GP.impute.time<impute.UB.temp+15]))
                IBI.GP.rmssd.45<-c(IBI.GP.rmssd.45, rmssd(IBI.GP.impute.temp[IBI.GP.impute.time>impute.LB.temp-22.5 & 
                                                                               IBI.GP.impute.time<impute.UB.temp+22.5]))
                IBI.GP.rmssd.60<-c(IBI.GP.rmssd.60, rmssd(IBI.GP.impute.temp[IBI.GP.impute.time>impute.LB.temp-30 & 
                                                                               IBI.GP.impute.time<impute.UB.temp+30]))
                runtime<-c(runtime, run_time)
                warmup.vec<-c(warmup.vec, warmup.temp)
                iter.vec<-c(iter.vec, n.iter.temp)
                Hz.vec<-c(Hz.vec, Hz.temp)
                impute.Hz.vec<-c(impute.Hz.vec, Hz.impute.temp)
                GP.Hz.vec<-c(GP.Hz.vec, Hz.GP.temp)
                impute.fac.vec<-c(impute.fac.vec, impute.fac.temp)
                impute.start<-c(impute.start, min.TIME2)
                impute.end<-c(impute.end, max.TIME2)
                impute.tot<-c(impute.tot, time.span)
                
                #Parameter Means
                HR.mean.vec<-c(HR.mean.vec, mean(HR.est$HR))
                Resp.mean.vec<-c(Resp.mean.vec, mean(Resp.est$R))
                
                #Parameter sd's 
                HR.SD.vec<-c(HR.SD.vec, sd(HR.est$HR))
                Resp.SD.vec<-c(Resp.SD.vec, sd(Resp.est$R))
                
                #Parameter R-hats
                HR.Rhat.vec<-c(HR.Rhat.vec, GP.summary$Rhat[rownames(GP.summary)=='HR'])
                Resp.Rhat.vec<-c(Resp.Rhat.vec, GP.summary$Rhat[rownames(GP.summary)=='R'])
                
                #Parameter N_eff
                HR.N_eff.vec<-c(HR.N_eff.vec, GP.summary$n_eff[rownames(GP.summary)=='HR'])
                Resp.N_eff.vec<-c(Resp.N_eff.vec, GP.summary$n_eff[rownames(GP.summary)=='R'])
                
                #Getting a vector of original value rmssd's
                IBI.rmssd.15<-c(IBI.rmssd.15, rmssd(IBI.orig.temp$IBI.done[IBI.orig.time>impute.LB.temp-7.5 & 
                                                                             IBI.orig.time<impute.UB.temp+7.5]))
                IBI.rmssd.30<-c(IBI.rmssd.30, rmssd(IBI.orig.temp$IBI.done[IBI.orig.time>impute.LB.temp-15 & 
                                                                             IBI.orig.time<impute.UB.temp+15]))
                IBI.rmssd.45<-c(IBI.rmssd.45, rmssd(IBI.orig.temp$IBI.done[IBI.orig.time>impute.LB.temp-22.5 & 
                                                                             IBI.orig.time<impute.UB.temp+22.5]))
                IBI.rmssd.60<-c(IBI.rmssd.60, rmssd(IBI.orig.temp$IBI.done[IBI.orig.time>impute.LB.temp-30 & 
                                                                             IBI.orig.time<impute.UB.temp+30]))
                
                #Getting a vector of mean replacement rmssd's 
                IBI.MR.rmssd.15<-c(IBI.MR.rmssd.15, rmssd(IBI.mean.replace[IBI.orig.time>impute.LB.temp-7.5 & 
                                                                             IBI.orig.time<impute.UB.temp+7.5]))
                IBI.MR.rmssd.30<-c(IBI.MR.rmssd.30, rmssd(IBI.mean.replace[IBI.orig.time>impute.LB.temp-15 & 
                                                                             IBI.orig.time<impute.UB.temp+15]))
                IBI.MR.rmssd.45<-c(IBI.MR.rmssd.45, rmssd(IBI.mean.replace[IBI.orig.time>impute.LB.temp-22.5 & 
                                                                             IBI.orig.time<impute.UB.temp+22.5]))
                IBI.MR.rmssd.60<-c(IBI.MR.rmssd.60, rmssd(IBI.mean.replace[IBI.orig.time>impute.LB.temp-30 & 
                                                                             IBI.orig.time<impute.UB.temp+30]))
                
                #Getting a vector of hotdeck replacement rmssd's
                IBI.HD.rmssd.15<-c(IBI.HD.rmssd.15, rmssd(IBI.hotdeck.replace[IBI.orig.time>impute.LB.temp-7.5 & 
                                                                                IBI.orig.time<impute.UB.temp+7.5]))
                IBI.HD.rmssd.30<-c(IBI.HD.rmssd.30, rmssd(IBI.hotdeck.replace[IBI.orig.time>impute.LB.temp-15 & 
                                                                                IBI.orig.time<impute.UB.temp+15]))
                IBI.HD.rmssd.45<-c(IBI.HD.rmssd.45, rmssd(IBI.hotdeck.replace[IBI.orig.time>impute.LB.temp-22.5 & 
                                                                                IBI.orig.time<impute.UB.temp+22.5]))
                IBI.HD.rmssd.60<-c(IBI.HD.rmssd.60, rmssd(IBI.hotdeck.replace[IBI.orig.time>impute.LB.temp-30 & 
                                                                                IBI.orig.time<impute.UB.temp+30]))
              }
            }
          }
        }
      }
    }
  }
}

save(list = ls(), 
     file = paste0(out.folder, '/', paste('GP_2', Sys.time(), n.sims, sep = '_'), '.RData'))
