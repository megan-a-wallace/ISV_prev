# Tom van der Most virus presence/absence data. Data wrangling and prevalence calculations

# Loading packages, functions and themes ---------------------------------------------------------

##packages
library(tidyverse)
library(viridis)
library(stringr)
library(lubridate)

##my colour pallete combined from viridis
pal_meg<-c("#340597FF","#7A02A8FF","#B42D8DFF","#DD5E66FF","#F1814DFF","#FCAA34FF","#FBD724FF","#ffff6d","linen", "#CFE11CFF","#9FDA3AFF","#4AC16DFF","#21908CFF","#365C8DFF","#08306B")

# Colour blind friendly palette (7 values):
cbPalette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

#Longer pallette (10 values)
cbPalette_long <- c("#000000","#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7","white")


MegTheme <- theme_bw() +
  theme(axis.title.x = element_text(vjust = -0.35,
                                    size = 16,
                                    colour = "black"),
        axis.title.y = element_text(vjust = 1.2,
                                    size = 16,
                                    colour = "black"),
        axis.text = element_text(size = 14, 
                                 colour = "black"),
        legend.text = element_text(size = 13, 
                                   colour = "black"),
        legend.title = element_text(size = 15, 
                                    colour = "black"),
        plot.title = element_text(size=20),
        strip.background = element_rect(fill = "white", colour = "dark grey"))

theme_set(MegTheme)

# InferPrevalence function 

##DJ0 with help from JJ Welch
##24Sept09, revised 11sep2013

# Calculating the 2DlogLik bounds on viral prevelance 

require(compiler)
### Core Likelihood function 
calc.log.lik<-cmpfun(function(p,hitFlies,missFlies){
  p.misses<-log((1-p)^missFlies)
  p.hits<-log(1-((1-p)^hitFlies))
  return(sum(c(p.hits,p.misses)))
})

## Convenience function that is used to identify bounds 
calc.bound<-cmpfun(function(p,LL,hitFlies,missFlies){
  p.misses<-log((1-p)^missFlies)
  p.hits<-log(1-((1-p)^hitFlies))
  return(
    (LL-sum(c(p.hits,p.misses)))^2
  )
})

######### End-user function that caclulates the estimated prevelance and if requested 
#########(1) interval-log-likelihood bounds, 
#########(2) LRT for for consistency between pooled and non-pooled samples, 
#########(3) Plots the LL surface

InferPrevalence<-function(nFlies,hits,bounds=FALSE,interval=2,test.consistant=FALSE, plot=FALSE){
  #strip out any NAs
  hits<-as.numeric(hits[!is.na(hits)])
  nFlies<-as.numeric(nFlies[!is.na(nFlies)])
  if((length(hits)==0)|(length(nFlies)==0)){return(NA)}
  #seperate out the pots into hits and misses
  hitFlies<-nFlies[which(as.logical(hits))]
  missFlies<-nFlies[which(!hits)]
  
  # maximise the LL (which is negative)
  optimise(f = calc.log.lik, interval = c(0,1), hitFlies, missFlies, maximum = TRUE,tol = .Machine$double.eps)->estimate	
  estimate$maximum->p
  estimate$objective->ML
  
  #Find the bounds (one at a time, seaching above and below the ML estimate
  if(bounds){
    optimise(f = calc.bound, interval = c(0,p), ML-interval, hitFlies, missFlies, maximum = FALSE,tol = .Machine$double.eps)$minimum->lower
    optimise(f = calc.bound, interval = c(p,1), ML-interval, hitFlies, missFlies, maximum = FALSE,tol = .Machine$double.eps)$minimum->upper
  }
  
  #If required to test for consistency between pooled and un-pooled samples
  if(test.consistant){
    #separate the data
    if((sum(nFlies>1)>0)&(sum(nFlies==1)>0)){
      hitFlies[hitFlies>1]->B.hit
      missFlies[missFlies>1]->B.miss
      hitFlies[hitFlies==1]->S.hit
      missFlies[missFlies==1]->S.miss
      #make the two estimates
      optimise(f = calc.log.lik, interval = c(0,1), B.hit, B.miss, maximum = TRUE,tol = .Machine$double.eps)->B.estimate
      optimise(f = calc.log.lik, interval = c(0,1), S.hit, S.miss, maximum = TRUE,tol = .Machine$double.eps)->S.estimate		
      c(B.estimate$maximum,S.estimate$maximum)->two.p
      #calculate their joint LL and do an LRT
      sum(B.estimate$objective,S.estimate$objective)->ML2
      2*(ML2-ML)->TwoDeltaLL
      pchisq(TwoDeltaLL, df =  1,lower.tail=FALSE)->p.value
    }else{
      test.consistant<-FALSE
    }
  }
  
  #Plot, if requested
  if(plot){
    if(!test.consistant){
      surface<-unlist(lapply(seq(0,1,0.0001),calc.log.lik,hitFlies,missFlies))
      plot(seq(0,1,0.0001),surface,type="l",ylim=c((ML-10),ML),xlab="Prevelance",ylab="log Likelihood")
      abline(v=p,col="red",lwd=4)	
      if(bounds){
        abline(h=(ML-interval))
        abline(v=lower,col="red",lty=3)
        abline(v=upper,col="red",lty=3)
      }
    }
    if(test.consistant){
      surface<-unlist(lapply(seq(0,1,0.0001),calc.log.lik,hitFlies,missFlies))
      surface1<-unlist(lapply(seq(0,1,0.0001),calc.log.lik,B.hit,B.miss))
      surface2<-unlist(lapply(seq(0,1,0.0001),calc.log.lik,S.hit,S.miss))
      plot(seq(0,1,0.0001),surface1,type="l",ylim=c(min(c(max(surface),max(surface1),max(surface2)))-10,max(c(surface,surface1,surface2))),xlab="Prevelance",ylab="log Likelihood",col="green")
      points(seq(0,1,0.0001),surface2,type="l",col="blue")
      points(seq(0,1,0.0001),surface,type="l",col="red")
      abline(v=c(p,two.p[1],two.p[2]),col=c("red","green","blue"),lwd=4)
    }
    
  }
  
  #contruct the return list
  result<-list()
  result$prevelance<-p
  result$log.liklihood<-ML
  if(bounds){result$bounds<-c(lower,upper)}
  if(test.consistant){
    result$alt.estimates<-two.p
    result$two.delta.LL<-TwoDeltaLL
    result$p.value<-p.value
    if(p.value<0.05){result$consistent<-FALSE}else{result$consistent<-TRUE}
  }
  return(result)
}

# Data formatting  --------------------------------------------------------

## Ugly code sorry...

read.csv('AMV_dat.csv',header=TRUE)->AMV_dat
read.csv('NEGV_dat.csv',header=TRUE)->NEGV_dat
read.csv('UMAV_dat.csv',header=TRUE)->UMAV_dat

AMV_vec<-(rep('AMV',132))
NEGV_vec<-(rep('NEGV',132))
UMAV_vec<-(rep('UMAV',128))

AMV_dat$virus<-AMV_vec
NEGV_dat$virus<-NEGV_vec
UMAV_dat$virus<-UMAV_vec

colnames(AMV_dat)[7]<-'PA'
colnames(NEGV_dat)[7]<-'PA'
colnames(UMAV_dat)[7]<-'PA'

#combined dataset with all viruses
bind_rows(AMV_dat,NEGV_dat,UMAV_dat)->PA_dat

#get variables into correct format
PA_dat$Pool_size<-as.numeric(as.character(PA_dat$Pool_size))
PA_dat$Location<-as.factor(PA_dat$Location)
##checking levels and removing whitespace for location 
levels(PA_dat$Location)<-gsub(" ","_",levels(PA_dat$Location))
PA_dat$Year<-as.factor(PA_dat$Year)#actually its all 2020...
PA_dat$virus<-as.factor(PA_dat$virus)
PA_dat$Week<-factor(PA_dat$Week, ordered = TRUE) 
#converting the weeks to dates by ISO week date standard (ISO-8601), eg. monday to sunday weeks - taking mondays as the day, this is for some reason not working quite correctly with lubridates method so having to minus 1 week, I think because it is or isn't counting the part week at the beginning of the year
as.Date(paste(2020, as.numeric(as.character(PA_dat$Week))-1, 1), format = "%Y %W %w")->PA_dat$date
#making the PA a binary variable 
PA_dat$PA<-as.factor(PA_dat$PA)
ifelse(PA_dat$PA == "Yes",1,0)->PA_dat$PA #this should be numeric for the inferPrevalence function I think

##Looking at the number of mosquitoes in each category of virus, week, location eg. the sampling size - I think these should probably be added to the plots at a later stage as there's quite some variation
PA_dat %>% group_by(virus,Location) %>% summarise(totalmosq = sum(Pool_size)) ->samplesize_by_virusloc
PA_dat %>% group_by(virus,Week) %>% summarise(totalmosq = sum(Pool_size)) ->samplesize_by_virusweek

# Global prevalence of each virus ----

###Loop which calculates a list of virus prevalence across the whole dataset
comb_virus_data<-PA_dat

viruses<-list(levels(PA_dat$virus))
prev_list<-list()
##to sum the log likelihoods with no separation by any co-variate
LL_global_vec<-as.numeric(vector(length = length(viruses[[1]])))

#cycling through the viruses in comb_viruses and calculating a data frame of prevalence using the InferPrevalence function for each one
for (i in 1:length(viruses[[1]])) {
  
  #print chosen virus
  print(paste("calculating prevalence for",viruses[[1]][i],sep = " "))
  
  #filtering the dataframe for the focal virus
  virus_data<-comb_virus_data[comb_virus_data$virus==viruses[[1]][i],]
  
  prev_list[[i]]<-virus_data %>%
    summarise(prevalence = as.numeric(InferPrevalence(nFlies=Pool_size,hits=PA,bounds=TRUE,interval=2,test.consistant=TRUE,plot=FALSE)$prevelance),
              lower_bound = as.numeric(InferPrevalence(nFlies=Pool_size,hits=PA,bounds=TRUE,interval=2,test.consistant=TRUE,plot=FALSE)$bounds[1]),
              upper_bound = as.numeric(InferPrevalence(nFlies=Pool_size,hits=PA,bounds=TRUE,interval=2,test.consistant=TRUE,plot=FALSE)$bounds[2]),
              log_likelihood = as.numeric(InferPrevalence(nFlies=Pool_size,hits=PA,bounds=TRUE,interval=2,test.consistant=TRUE,plot=FALSE)$log.liklihood)) %>% 
    as.data.frame() 
  
  names(prev_list)[i] <- viruses[[1]][i]
  LL_global_vec[i] <- prev_list[[i]]$log_likelihood
  
  #updating progress of loop
  print(paste(viruses[[1]][i],"prevalence calculated",sep = " "))
}

##Making and exporting a table from this list with basic prevalence across the whole dataset for each virus
data.frame(virus=c("AMV","NEGV","UMAV"),
             prevalence=c(prev_list[[1]][[1]],prev_list[[2]][[1]],prev_list[[3]][[1]]),
           lower_bound=c(prev_list[[1]][[2]],prev_list[[2]][[2]],prev_list[[3]][[2]]),
           upper_bound=c(prev_list[[1]][[3]],prev_list[[2]][[3]],prev_list[[3]][[3]]),
           log_likelihood=c(prev_list[[1]][[4]],prev_list[[2]][[4]],prev_list[[3]][[4]]),
           stringsAsFactors = FALSE)->prev_by_virus_df

write.csv(prev_by_virus_df,file = "TvdM_globalprev_byvirus.csv")

###Virus prevalence by location and log 10 prevalence plots -----

##prepping the data
comb_virus_data<-PA_dat

viruses<-list(levels(PA_dat$virus))
prev_by_site_list<-list()

#cycling through the viruses in comb_viruses and calculating a by location data frame of prevalence using the InferPrevalence function for each one
for (i in 1:length(viruses[[1]])) {
  
  #print chosen virus
  print(paste("creating prevalence by site df for",viruses[[1]][i],sep = " "))
  
  #filtering the dataframe for the focal virus
  virus_data<-comb_virus_data[comb_virus_data$virus==viruses[[1]][i],]
  
  prev_by_site_list[[i]]<-virus_data %>%
    group_by(Location) %>%
    summarise(prevalence = as.numeric(InferPrevalence(nFlies=Pool_size,hits=PA,bounds=TRUE,interval=2,test.consistant=TRUE,plot=FALSE)$prevelance),
              lower_bound = as.numeric(InferPrevalence(nFlies=Pool_size,hits=PA,bounds=TRUE,interval=2,test.consistant=TRUE,plot=FALSE)$bounds[1]),
              upper_bound = as.numeric(InferPrevalence(nFlies=Pool_size,hits=PA,bounds=TRUE,interval=2,test.consistant=TRUE,plot=FALSE)$bounds[2]),
              log_likelihood = as.numeric(InferPrevalence(nFlies=Pool_size,hits=PA,bounds=TRUE,interval=2,test.consistant=TRUE,plot=FALSE)$log.liklihood)) %>% 
    as.data.frame() 
  
  #ordering the data by longitude
  #prev_by_site_list[[i]]<-prev_by_site_list[[i]][order(prev_by_site_list[[i]]$lon),]
  
  names(prev_by_site_list)[i] <- viruses[[1]][i]
  
  #updating progress of loop
  print(paste(viruses[[1]][i],"prevalence by site df created",sep = " "))
}

#once this loop has been run - check that you have a list of prevalence tables, one for each virus split by location/site
glimpse(prev_by_site_list)

#exporting a data table of prevalence and bounds by Site and virus 
prev_by_site_df<-data.frame(Virus=c(rep.int("AMV",16),rep.int("NEGV",16),rep.int("UMAV",16)),
                            Location=c(levels(prev_by_site_list$AMV$Location),levels(prev_by_site_list$NEGV$Location),levels(prev_by_site_list$UMAV$Location)),
                            Prevalence=c(prev_by_site_list$AMV$prevalence,prev_by_site_list$NEGV$prevalence,prev_by_site_list$UMAV$prevalence),
                            Lower_bound=c(prev_by_site_list$AMV$lower_bound,prev_by_site_list$NEGV$lower_bound,prev_by_site_list$UMAV$lower_bound),
                            Upper_bound=c(prev_by_site_list$AMV$upper_bound,prev_by_site_list$NEGV$upper_bound,prev_by_site_list$UMAV$upper_bound),
                            Log_likelihood=c(prev_by_site_list$AMV$log_likelihood,prev_by_site_list$NEGV$log_likelihood,prev_by_site_list$UMAV$log_likelihood),
                            stringsAsFactors = FALSE)

write.csv(prev_by_site_df,file = "TvdM_prev_by_location.csv")

#Loop for creating a plot series of prevalence by location for each virus
##You should get a barplot for each of the three viruses after this, uncomment the pdf creating bits to have it automatically explort a pdf plot for each one

#Make list of viruses in the form you want the title to take...eg. capitalised
virus_titles<-c("Alphamesonivirus","Negevirus","Umatilla virus") #I guessed at these

multi_virus_data<-prev_by_site_list

for (i in 1:length(multi_virus_data)) {
  
  #isolating data for only one virus
  data<-multi_virus_data[[i]]
  
  #preparing the data
  #proportions -> percentages
  data[,2:4] <- data[,2:4]*100
  #changing any value < 0.01% to 0.01%
  index <-data[,2:4] < 0.01
  data[,2:4][index] <- 0.01
  prev_index<-index[,1] #for changing colours etc. of sp w prev < 0.01
  lb_index<-index[,2] #for changing the lower bounds the error bars when prev is < 0.01 to 0.001
  data[,3][lb_index]<-0.001 #changing those lower bounds
  
  ##if none of the non-zero lower bounds are less than 0.1...changing the limit of the graph to 0.1, if none are less than 1, changing to 1
  non_index<-data[,3:4]>0.01
  
  if (min(data[,3:4][non_index])>1) {
    
    data[,2][prev_index]<-1#change prevalence to 1% if <0.01 so that that't the bottom of the graph
    
    line_list<-c(1.5,2,3,5,7.5,10,15,20,30,50,75,100)
    
    #making some fake data so that the whole range is encompassed - doesn't matter that these are not the right names, there just needs to be 16
    fake_sites<-c("CR","SI","CS","BR","BG","CE","TR","HL","BN","LW","IB","DK","IN","VG","OX","CK")
    fake_prevalence<-as.numeric(seq(1,100,length.out = 16))
    fake_data<-data.frame(fake_sites,fake_prevalence)
    
    y_axis<-c(1,1.5,2,3,5,7.5,10,15,20,30,50,75,100)
    x_text_y<-par("usr")[4]-1.24
    
  } else if (min(data[,3:4][non_index])>0.1) {
    
    data[,2][prev_index]<-0.1#change prevalence to 0.1% if <0.01 so that that't the bottom of the graph
    
    line_list<-c(0.2,0.3,0.5,1,2,3,5,10,20,30,50,100)
    
    #making some fake data so that the whole range is encompassed
    fake_sites<-c("CR","SI","CS","BR","BG","CE","TR","HL","BN","LW","IB","DK","IN","VG","OX","CK")
    fake_prevalence<-as.numeric(seq(0.1,100,length.out = 16))
    fake_data<-data.frame(fake_sites,fake_prevalence)
    
    y_axis<-c(0.1,0.2,0.3,0.5,1,2,3,5,10,20,30,50,100)
    x_text_y<-par("usr")[4]-1.93
    
  } else {
    
    data=data #keep the same, w bottom at 0.01%
    
    line_list<-c(0.03,0.05,0.1,0.3,0.5,1,3,5,10,30,50,100)
    
    #making some fake data so that the whole range is encompassed
    fake_sites<-c("CR","SI","CS","BR","BG","CE","TR","HL","BN","LW","IB","DK","IN","VG","OX","CK")
    fake_prevalence<-as.numeric(seq(0.01,100,length.out = 16))
    fake_data<-data.frame(fake_sites,fake_prevalence)
    
    y_axis<-c(0.01,0.03,0.05,0.1,0.3,0.5,1,3,5,10,30,50,100)
    x_text_y<-par("usr")[4]-1.994
  }
  
  #pdf(file = paste("log_prevalence_by_site_plots/log_prevalence_by_site_", sub(" ","_",virus_titles[i]), ".pdf", sep = ""),width = 12,height = 6)
  #setting margins for plotting 
  par(mar = c(10.1,6.1, 2.6, 2.1), # change the margins
      lwd = 1.7,# increase the line thickness
      cex.axis = 1.3, # increase default axis label size
      cex.lab = 1.4)
  
  # starting to prepare the plot
  #initally making the empty plot using some fake data, so that the whole plot space needed is created
  barplot(fake_data$fake_prevalence~fake_data$fake_sites,axes = FALSE,xaxt='n',yaxt='n', border= NA, col = NA,main="", xlab="", ylab="",log = 'y') # invisible bars - only plot
  #using a loop to add custom lines on the log scale
  for(j in 1:length(line_list)){
    list<-line_list
    par(xpd = TRUE)
    lines(x = c(par("usr")[1],par("usr")[2]), y = c(list[j],list[j]),lty = 1,lwd=2,col="grey90")
  }
  par(xpd=FALSE)#REMEMBER TO RUN THIS AFTER
  #adding y axis
  axis(2, at = y_axis, tick = TRUE, labels = y_axis, lwd = 2,lwd.ticks = 2,las = 1,mgp=c(3, 0.75, 0))
  #adding y axis label
  mtext(side = 2, line = 4, "Virus prevalence (%)", cex = 1.6,padj = 0.6)
  #adding x axis
  vec<-seq(par("usr")[1]+0.8,par("usr")[2]-0.448,length.out = 16)
  axis(1, at = vec,
       tick = FALSE,
       labels = FALSE)
  #x coordinates for labels 
  xtext_cols<-rep.int("black",20)
  xtext_cols[prev_index]<-"grey60"
  text(x = vec,
       y = x_text_y,
       labels = c("Almen","Almere De Vaarten","Asten","Bergumermeer","Leiden Hortus Botanicus","Lelystad","Maastricht","Nijmegen Berg en Dal","Overdinkel","Reusel","Rotterdam museumpark","Utrecht Eendenkooi","Utrecht Griend","Utrecht Vleuterweide","Wageningen Campus","Zwarte Meer"),
       ## Rotate the labels by 35 degrees.
       xpd = NA,
       srt = 30,
       adj = 1,
       cex = 1.5,
       col = xtext_cols)
  #adding x axis label
  mtext(side = 1, line = 8,"Sites", cex = 1.6)
  #Setting the amount of space to leave before each bar
  bar_spacing<-c(-0.22,rep.int(0.27,15)) 
  
  #making colour vector fr bars
  xbar_cols<-rep.int("#009E73",16)
  xbar_cols[prev_index]<-rgb(199/255,199/255,199/255,0.4)#alter alpha [4] here to show low value bars
  xborder_cols<-rep.int("black",16)
  xborder_cols[prev_index]<-rgb(199/255,199/255,199/255,0.4)#alter alpha [4] here to show low value bar borders
  #plotting actual data on to the plot
  barplot(data$prevalence~data$Location,add=TRUE,yaxt ="n",xaxt="n",log = "y",col=xbar_cols,xpd=TRUE,border=xborder_cols,space=bar_spacing)
  #
  #adding error bars
  barCenters<-barplot(data$prevalence~data$Location,plot=FALSE,yaxt ="n",xaxt="n",log = "y",space=bar_spacing)
  #space=bar_spacing
  arrows(barCenters,data$lower_bound,barCenters,data$upper_bound, lwd=1.7, angle=90, code=3,length = 0.1,col = xtext_cols)
  #adding text to plot with virus name 
  mtext(side=3,line=1,paste(virus_titles[i]),cex=1.8,padj = 0.25)
  
  #dev.off()
}

### Virus prevalence by week/collection date and log10 prevalence plots ----

##prepping the data
comb_virus_data<-PA_dat

viruses<-list(levels(PA_dat$virus))
prev_by_week_list<-list()

#cycling through the viruses in comb_viruses and calculating a by week data frame of prevalence using the InferPrevalence function for each one
for (i in 1:length(viruses[[1]])) {
  
  #print chosen virus
  print(paste("creating prevalence by week df for",viruses[[1]][i],sep = " "))
  
  #filtering the dataframe for the focal virus
  virus_data<-comb_virus_data[comb_virus_data$virus==viruses[[1]][i],]
  
  prev_by_week_list[[i]]<-virus_data %>%
    group_by(Week) %>%
    summarise(prevalence = as.numeric(InferPrevalence(nFlies=Pool_size,hits=PA,bounds=TRUE,interval=2,test.consistant=TRUE,plot=FALSE)$prevelance),
              lower_bound = as.numeric(InferPrevalence(nFlies=Pool_size,hits=PA,bounds=TRUE,interval=2,test.consistant=TRUE,plot=FALSE)$bounds[1]),
              upper_bound = as.numeric(InferPrevalence(nFlies=Pool_size,hits=PA,bounds=TRUE,interval=2,test.consistant=TRUE,plot=FALSE)$bounds[2]),
              log_likelihood = as.numeric(InferPrevalence(nFlies=Pool_size,hits=PA,bounds=TRUE,interval=2,test.consistant=TRUE,plot=FALSE)$log.liklihood)) %>% 
    as.data.frame() 
  
  names(prev_by_week_list)[i] <- viruses[[1]][i]
  
  #updating progress of loop
  print(paste(viruses[[1]][i],"prevalence by week df created",sep = " "))
}

#once this loop has been run - check that you have a list of prevalence tables, one for each virus split by week
glimpse(prev_by_week_list)

#exporting a data table of prevalence and bounds by week and virus 
prev_by_week_df<-data.frame(Virus=c(rep.int("AMV",12),rep.int("NEGV",12),rep.int("UMAV",12)),
                            Week=c(levels(prev_by_week_list$AMV$Week),levels(prev_by_week_list$NEGV$Week),levels(prev_by_week_list$UMAV$Week)),
                            Prevalence=c(prev_by_week_list$AMV$prevalence,prev_by_week_list$NEGV$prevalence,"NA",prev_by_week_list$UMAV$prevalence),
                            Lower_bound=c(prev_by_week_list$AMV$lower_bound,prev_by_week_list$NEGV$lower_bound,"NA",prev_by_week_list$UMAV$lower_bound),
                            Upper_bound=c(prev_by_week_list$AMV$upper_bound,prev_by_week_list$NEGV$upper_bound,"NA",prev_by_week_list$UMAV$upper_bound),
                            Log_likelihood=c(prev_by_week_list$AMV$log_likelihood,prev_by_week_list$NEGV$log_likelihood,"NA",prev_by_week_list$UMAV$log_likelihood),
                            stringsAsFactors = FALSE) #NAs because the 27th week wasn't tested for UMAV

write.csv(prev_by_week_df,file = "TvdM_prev_by_week.csv")

#Loop for creating a plot series of prevalence by week/date for each virus
##You should get a barplot for each of the three viruses after this, uncomment the pdf creating bits to have it automatically explort a pdf plot for each one

#IMPORTANT - THIS CURRENTLY ISN'T WORKING FOR UMAV BECAUSE OF THE MISSING DATA IN WK 27...FIGURE OUT A FIX. 

#Make list of viruses in the form you want the title to take...eg. capitalised
virus_titles<-c("Alphamesonivirus","Negevirus","Umatilla virus") #I guessed at these

multi_virus_data<-prev_by_week_list

for (i in 1:length(multi_virus_data)) {
  
  #isolating data for only one virus
  data<-multi_virus_data[[i]]
  
  #preparing the data
  #proportions -> percentages
  data[,2:4] <- data[,2:4]*100
  #changing any value < 0.01% to 0.01%
  index <-data[,2:4] < 0.01
  data[,2:4][index] <- 0.01
  prev_index<-index[,1] #for changing colours etc. of sp w prev < 0.01
  lb_index<-index[,2] #for changing the lower bounds the error bars when prev is < 0.01 to 0.001
  data[,3][lb_index]<-0.001 #changing those lower bounds
  
  ##if none of the non-zero lower bounds are less than 0.1...changing the limit of the graph to 0.1, if none are less than 1, changing to 1
  non_index<-data[,3:4]>0.01
  
  if (min(data[,3:4][non_index])>1) {
    
    data[,2][prev_index]<-1#change prevalence to 1% if <0.01 so that that't the bottom of the graph
    
    line_list<-c(1.5,2,3,5,7.5,10,15,20,30,50,75,100)
    
    #making some fake data so that the whole range is encompassed - doesn't matter that these are not the right names, there just needs to be 12
    fake_weeks<-c("CR","SI","CS","BR","BG","CE","TR","HL","BN","LW","IB","CK")
    fake_prevalence<-as.numeric(seq(1,100,length.out = 12))
    fake_data<-data.frame(fake_weeks,fake_prevalence)
    
    y_axis<-c(1,1.5,2,3,5,7.5,10,15,20,30,50,75,100)
    x_text_y<-par("usr")[4]-1.24
    
  } else if (min(data[,3:4][non_index])>0.1) {
    
    data[,2][prev_index]<-0.1#change prevalence to 0.1% if <0.01 so that that't the bottom of the graph
    
    line_list<-c(0.2,0.3,0.5,1,2,3,5,10,20,30,50,100)
    
    #making some fake data so that the whole range is encompassed
    fake_weeks<-c("CR","SI","CS","BR","BG","CE","TR","HL","BN","LW","IB","CK")
    fake_prevalence<-as.numeric(seq(0.1,100,length.out = 12))
    fake_data<-data.frame(fake_weeks,fake_prevalence)
    
    y_axis<-c(0.1,0.2,0.3,0.5,1,2,3,5,10,20,30,50,100)
    x_text_y<-par("usr")[4]-1.93
    
  } else {
    
    data=data #keep the same, w bottom at 0.01%
    
    line_list<-c(0.03,0.05,0.1,0.3,0.5,1,3,5,10,30,50,100)
    
    #making some fake data so that the whole range is encompassed
    fake_weeks<-c("CR","SI","CS","BR","BG","CE","TR","HL","BN","LW","IB","CK")
    fake_prevalence<-as.numeric(seq(0.01,100,length.out = 12))
    fake_data<-data.frame(fake_weeks,fake_prevalence)
    
    y_axis<-c(0.01,0.03,0.05,0.1,0.3,0.5,1,3,5,10,30,50,100)
    x_text_y<-par("usr")[4]-1.994
  }
  
  #pdf(file = paste("log_prevalence_by_week_plots/log_prevalence_by_week_", sub(" ","_",virus_titles[i]), ".pdf", sep = ""),width = 12,height = 6)
  #setting margins for plotting 
  par(mar = c(6.1,6.1, 2.6, 2.1), # change the margins
      lwd = 1.7,# increase the line thickness
      cex.axis = 1.3, # increase default axis label size
      cex.lab = 1.4)
  
  # starting to prepare the plot
  #initally making the empty plot using some fake data, so that the whole plot space needed is created
  barplot(fake_data$fake_prevalence~fake_data$fake_weeks,axes = FALSE,xaxt='n',yaxt='n', border= NA, col = NA,main="", xlab="", ylab="",log = 'y') # invisible bars - only plot
  #using a loop to add custom lines on the log scale
  for(j in 1:length(line_list)){
    list<-line_list
    par(xpd = TRUE)
    lines(x = c(par("usr")[1],par("usr")[2]), y = c(list[j],list[j]),lty = 1,lwd=2,col="grey90")
  }
  par(xpd=FALSE)#REMEMBER TO RUN THIS AFTER
  #adding y axis
  axis(2, at = y_axis, tick = TRUE, labels = y_axis, lwd = 2,lwd.ticks = 2,las = 1,mgp=c(3, 0.75, 0))
  #adding y axis label
  mtext(side = 2, line = 4, "Virus prevalence (%)", cex = 1.6,padj = 0.6)
  #adding x axis
  vec<-seq(par("usr")[1]+0.8,par("usr")[2]-0.448,length.out = 12)
  axis(1, at = vec,
       tick = FALSE,
       labels = FALSE)
  #x coordinates for labels 
  xtext_cols<-rep.int("black",12)
  xtext_cols[prev_index]<-"grey60"
  text(x = vec,
       y = x_text_y,
       labels = c("29th Jun","6th Jul","13th Jul","20th Jul","27th Jul","3rd Aug","10th Aug","17th Aug","24th Aug","31st Aug","7th Sep","14th Sep"),
       ## Rotate the labels by 35 degrees.
       xpd = NA,
       srt = 30,
       adj = 1,
       cex = 1.5,
       col = xtext_cols)
  #adding x axis label
  mtext(side = 1, line = 5,"Collection Date (week beginning)", cex = 1.6)
  #Setting the amount of space to leave before each bar
  bar_spacing<-c(-0.22,rep.int(0.27,11)) 
  
  #making colour vector fr bars
  xbar_cols<-rep.int("#009E73",12)
  xbar_cols[prev_index]<-rgb(199/255,199/255,199/255,0.4)#alter alpha [4] here to show low value bars
  xborder_cols<-rep.int("black",12)
  xborder_cols[prev_index]<-rgb(199/255,199/255,199/255,0.4)#alter alpha [4] here to show low value bar borders
  #plotting actual data on to the plot
  barplot(data$prevalence~data$Week,add=TRUE,yaxt ="n",xaxt="n",log = "y",col=xbar_cols,xpd=TRUE,border=xborder_cols,space=bar_spacing)
  #
  #adding error bars
  barCenters<-barplot(data$prevalence~data$Week,plot=FALSE,yaxt ="n",xaxt="n",log = "y",space=bar_spacing)
  #space=bar_spacing
  arrows(barCenters,data$lower_bound,barCenters,data$upper_bound, lwd=1.7, angle=90, code=3,length = 0.1,col = xtext_cols)
  #adding text to plot with virus name 
  mtext(side=3,line=1,paste(virus_titles[i]),cex=1.8,padj = 0.25)
  
  #dev.off()
}
