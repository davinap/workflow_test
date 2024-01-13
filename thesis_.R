setwd("C:/Users/davin/Documents/PhD/Writing/THESIS/scripttest/")
#Load packages
library(tidyverse)
library(gtools)
library(reshape2)
library(diptest)
library(diveRsity)
library(gridExtra)
library(openxlsx)
library(readxl)
library(fitdistrplus)
library(flexmix)
library(tictoc)

#NOTES:
# can run each 'MEDIA' chunk by itself but before doing that:
# 1. change the directory where the csv files are located
# 2. change the script where the control file is indicated i.e. where it says 'B8' to the correct well

tic()
#####ALL FUNCTIONS#####

#Read in files
#read csvs in a folder, remove negative values, log transform and create giant table
read.files.in.log <- function(f){
  all50 <- data.frame(cell=1:50000) #1 to 50000 column to cbind later
  csvin <-read.csv(f, header=TRUE, stringsAsFactors = FALSE)
  csvpos <- filter(csvin, csvin[1]>0) #remove rows containing negative values which would break log transform
  logit <- log10(csvpos)
  event<-1:nrow(logit) #make a column with rowlength of each individual file
  bind <- cbind(event, logit) #add cells column to all before extracting the bl1h channel
  join <- left_join(all50, bind, by=c("cell"="event")) #cells--all channel columns
  return(join) #get the exported csv with a column appended that contains the event (cell) number and bl1h col only
}


#Scales output of fitdist so you can plot it as density rather than cell count
sdnorm <- function(x, mean, sd, lambda){lambda*dnorm(x, mean, sd)}


#ggplot characteristics for 'get.fd.plots' function
#modified with new colours
#stat_function 1 = plot normal distribution of the sample - red
#stat_function 2 = plot normal distribution of the negative control - grey
get.norm.plot<-function(logged, sumdf){ggplot(melt(logged), aes(x=value)) + 
    geom_density(color="#A0A0A0") +
    stat_function(fun=dnorm,
                  args=list(mean=sumdf[[1]],
                            sd=sumdf[[2]]),
                  color="#FF6666", fill="#FF6666",alpha=0.3, geom="polygon") +
    stat_function(fun=dnorm,
                  args=list(mean=mean(ctrl$value, na.rm=TRUE), 
                            sd=sd(ctrl$value, na.rm=TRUE)),
                  geom = "area",
                  color="#E0E0E0",
                  fill="grey", alpha=0.3) +
    theme_minimal() +
    scale_x_continuous(limits=c(2,6)) +
    scale_y_continuous(limits = c(0,3.75)) +
    xlab("")+
    ylab("")
}


#Get the values from fitdist for each csv file and plot normal distributions
#lines 2-4 = read csv, remove negative values and log10 transform them
#set.seed sets the random number generator so it's reproducible - set.seed(1)
#fitdist fits data to normal curve - uses dnorm function as its foundation
#returns a ggplot object
get.fd.plots <- function (f){
  file<-read.csv(f) 
  pos <- filter(file, file[1]>0)
  logged <-log10(pos)
  set.seed(1) 
  fd<-fitdist(logged[[1]],"norm") #use fitdist function to generate a fit to normal distribution
  sumdf<-data.frame(mu=fd$estimate[1],sigma=fd$estimate[2], row.names = "stat") #make a dataframe containing the fit values
  preplot<-get.norm.plot(logged, sumdf) #use plotting function to plot the fitdist results
  return(preplot) #show the plot
}


#Get mean and standard deviation from fitdist
#run fitdist and return a dataframe with mean and sigma
get.fd.stats <- function (f){
  file<-read.csv(f) 
  pos <- filter(file, file[1]>0)
  logged <-log10(pos)
  set.seed(1)
  fd<-fitdist(logged[[1]],"norm")
  sumdf<-data.frame(mu=fd$estimate[1],sigma=fd$estimate[2], row.names = "stat")
  return(sumdf)
}


#Tidy the dataframe returned from get.fd.stats
tidy.fd.stats <- function (fd){
  df<-as.data.frame(t(fd)) 
  rown<-tibble::rownames_to_column(df, "sample") #make column 0 (sample names) the row names
  cleancol <- rown[2:nrow(rown),] #remove first column containing 'NA'
  odd<-seq(1,nrow(cleancol),2) #numbers for extracting odd rows (odd=mu) - start at 1, until the last row, go by every 2 numbers
  even <- seq(2,nrow(cleancol),2) #numbers for extracting even rows (even=sigma)
  stats_e<-cleancol[even,] #extract means (odd rows) and standard deviations (even rows)
  stats_o<-cleancol[odd,]
  bind<- cbind(stats_o,stats_e[2]) #cbind odd and even together
  names(bind)[2] <- "geometric mean" 
  names(bind)[3] <- "sd"
  return(bind)
}


#FUNCTIONS FOR FLEXFIT MODELS

#Gaussian models required for flexmix
mo1 <- FLXMRglm(family = "gaussian")
mo2 <- FLXMRglm(family = "gaussian")
mo3 <- FLXMRglm(family = "gaussian")

#For scaled density plots from flexfit data
plot_mix_comps <- function(x, mu, sigma, lam) {lam * dnorm(x, mu, sigma)}


#Plot 2 COMPONENT flexmix plots
#with modified colours, set.seed(1)
#lines 2-4 = read csv, remove negative values, log transform
#ff = get flexmix model for 2 components using the gaussian models mo1 and mo2
#c1/2 = get mu, sigma for pop1 and pop2
#c1df = make c1/2 output a dataframe so ggplot can access the values
#lam = get proportions of cells in pop1 and pop2
#stat_function 1/2 = plot populations scaled for density (plot_mix_comps) - red and turquoise
#stat_function 3 = plot negative control in grey
get.flex.plots2<-function(f){
  csvin <- read.csv(f, header = TRUE, stringsAsFactors = FALSE)
  csvpos <- filter(csvin, csvin[1]>0) #remove rows containing negative values which would break log transform
  logged <- log10(csvpos)
  set.seed(1)
  ff<-flexmix(logged$`Comp.BL1.H....BL1.H`~1, k=2, model=list(mo1, mo2))
  c1 <- parameters(ff, component=1)[[1]]
  c2 <- parameters(ff, component=2)[[1]]
  c1df<- as.data.frame(c1)
  c2df <- as.data.frame(c2)
  lam <- table(clusters(ff))
  ggplot(logged, aes(x=`Comp.BL1.H....BL1.H`)) + 
    geom_density(color="#A0A0A0") +
    stat_function(data=logged,
                  fun = plot_mix_comps,
                  args = list(c1df[1,1], c1df[2,1],lam[1]/sum(lam)),
                  color="#FF6666", fill="#FF6666",alpha=0.3, geom="polygon") +
    stat_function(data=logged, 
                  fun = plot_mix_comps,
                  args = list(c2df[1,1], c2df[2,1], lam[2]/sum(lam)),
                  color="#009999",fill="#00CCCC",alpha=0.3,geom="polygon") +
    stat_function(fun=dnorm,
                  args=list(mean=mean(ctrl$value, na.rm=TRUE), 
                            sd=sd(ctrl$value, na.rm=TRUE)),
                  geom = "area",
                  color="#E0E0E0",
                  fill="grey", alpha=0.3) +
    theme_minimal() +
    scale_x_continuous(limits=c(2,6)) +
    scale_y_continuous(limits = c(0,3.75)) +
    xlab("")+
    ylab("")
}


#Plot 3 COMPONENT flexmix plots - same general structure as get.flex.plots2
get.flex.plots3<-function(f){
  csvin <- read.csv(f, header = TRUE, stringsAsFactors = FALSE)
  csvpos <- filter(csvin, csvin[1]>0) #remove rows containing negative values which would break log transform
  logged <- log10(csvpos)
  set.seed(1)
  ff<-flexmix(logged$`Comp.BL1.H....BL1.H`~1, k=3, model=list(mo1, mo2, mo3))
  c1 <- parameters(ff, component=1)[[1]]
  c2 <- parameters(ff, component=2)[[1]]
  c3 <-  parameters(ff, component=3)[[1]]
  c1df<- as.data.frame(c1)
  c2df <- as.data.frame(c2)
  c3df <- as.data.frame(c3)
  lam <- table(clusters(ff))
  ggplot(logged, aes(x=`Comp.BL1.H....BL1.H`)) + 
    geom_density(color="#A0A0A0") +
    stat_function(data=logged,
                  fun = plot_mix_comps,
                  args = list(c1df[1,1], c1df[2,1],lam[1]/sum(lam)),
                  color="#FF6666", fill="#FF6666",alpha=0.4, geom="polygon") +
    stat_function(data=logged, 
                  fun = plot_mix_comps,
                  args = list(c2df[1,1], c2df[2,1], lam[2]/sum(lam)),
                  color="#009999",fill="#00CCCC",alpha=0.35,geom="polygon") +
    stat_function(data=logged, 
                  fun = plot_mix_comps,
                  args = list(c3df[1,1], c3df[2,1], lam[3]/sum(lam)),
                  color="#9933FF",fill="#AB00FF",alpha=0.2,geom="polygon") +
    stat_function(fun=dnorm,
                  args=list(mean=mean(ctrl$value, na.rm=TRUE), 
                            sd=sd(ctrl$value, na.rm=TRUE)),
                  geom = "area",
                  color="#E0E0E0",
                  fill="grey", alpha=0.25) +
    theme_minimal() +
    scale_x_continuous(limits=c(2,6)) +
    scale_y_continuous(limits = c(0,3.75)) +
    xlab("")+
    ylab("")
}


#Get stats for 2 COMPONENT flexmix models
#same general structure as get.flex.plots2
#statsff = make a list of each population's stats - mu, sigma and lambda
get.flex.stats2<-function(f){
  csvin <- read.csv(f, header = TRUE, stringsAsFactors = FALSE)
  csvpos <- filter(csvin, csvin[1]>0) #remove rows containing negative values which would break log transform
  logged <- log10(csvpos)
  set.seed(1)
  ff<-flexmix(logged$`Comp.BL1.H....BL1.H`~1, k=2, model=list(mo1, mo2))
  c1 <- parameters(ff, component=1)[[1]]
  c2 <- parameters(ff, component=2)[[1]]
  c1df<- as.data.frame(c1)
  c2df <- as.data.frame(c2)
  lam <- table(clusters(ff))
  statsff<- list(mu1=c1df[1,1], mu2=c2df[1,1],
                 sd1=c1df[2,1],sd2=c2df[2,1],
                 lam1=lam[[1]],lam2=lam[[2]], 
                 total=(lam[[1]]+lam[[2]]))
  return(statsff)
}

#Tidy up list returned from get.flex.stats2 (2 COMPONENTS)
#creates a transposed dataframe, renames columns and removes empty placeholder row
tidy.flex2 <- function(ff){
  transdf <- as.data.frame(t(ff))
  names(transdf)[1] <- "mu1" 
  names(transdf)[2] <- "mu2"
  names(transdf)[3] <- "sd1"
  names(transdf)[4] <- "sd2"
  names(transdf)[5] <- "lam1"
  names(transdf)[6] <- "lam2"
  names(transdf)[7] <- "n"
  out <- transdf[2:nrow(transdf),]
  return(out)
}


#get.flex.stats for 3 COMPONENTS
get.flex.stats3<-function(f){
  csvin <- read.csv(f, header = TRUE, stringsAsFactors = FALSE)
  csvpos <- filter(csvin, csvin[1]>0) #remove rows containing negative values which would break log transform
  logged <- log10(csvpos)
  set.seed(1)
  ff<-flexmix(logged$`Comp.BL1.H....BL1.H`~1, k=3, model=list(mo1, mo2,mo3))
  c1 <- parameters(ff, component=1)[[1]]
  c2 <- parameters(ff, component=2)[[1]]
  c3 <- parameters(ff, component=3)[[1]]
  c1df<- as.data.frame(c1)
  c2df <- as.data.frame(c2)
  c3df <- as.data.frame(c3)
  lam <- table(clusters(ff))
  statsff<- list(mu1=c1df[1,1], mu2=c2df[1,1], mu3=c3df[1,1],
                 sd1=c1df[2,1],sd2=c2df[2,1],sd3=c3df[2,1],
                 lam1=lam[[1]],lam2=lam[[2]],lam3=lam[[3]],
                 total=(lam[[1]]+lam[[2]])+lam[[3]])
  return(statsff)
}

#Tidy up list returned from get.flex.stats3
tidy.flex3 <- function(ff){
  transdf <- as.data.frame(t(ff))
  names(transdf)[1] <- "mu1" 
  names(transdf)[2] <- "mu2"
  names(transdf)[3] <- "mu3"
  names(transdf)[4] <- "sd1"
  names(transdf)[5] <- "sd2"
  names(transdf)[6] <- "sd3"
  names(transdf)[7] <- "lam1"
  names(transdf)[8] <- "lam2"
  names(transdf)[9] <- "lam3"
  names(transdf)[10] <- "n"
  out <- transdf[2:nrow(transdf),]
  return(out)
}


#Get Bayesian information criteria - values for whether 2 or 3 components gives the better fit
#the lower the BIC, the better. Most BICs are -ve because of the large sample size
#is very slow ~ 10m runtime from start to end
#runs flexmix modelling for 2 to 3 components (k=2:3) with 2 repetitions (nrep=2) using the gaussian model (mo1)
#BIC(ex) returns the BIC values only
get.BIC <- function(f){
  csvin <- read.csv(f, header = TRUE, stringsAsFactors = FALSE)
  csvpos <- filter(csvin, csvin[1]>0) #remove rows containing negative values which would break log transform
  logged <- log10(csvpos)
  set.seed(1)
  ex <- initFlexmix(logged$Comp.BL1.H....BL1.H~1, k = 2:3, model = mo1, nrep = 2) 
  df <- BIC(ex)
  return(df)
}



####===================================================MEDIA 1 - Y1=======================================================####

#read in files in A1, A2... order
y1_in <- mixedsort(list.files(path ="C:/Users/davin/Documents/PhD/Writing/THESIS/scripttest/", 
                              full.names=TRUE, pattern = ".csv"))
y1_in
read.files.in.log(y1_in[[1]]) #test function

#Loop to combine columns
loop_y1 <- data.frame(cell_number=1:50000) #empty dataframe for for loop
for(c in y1_in){loop_y1[[c]]<-read.files.in.log(c)}

write.csv(loop_y1, "y1.csv", row.names = FALSE) #need to write.csv so headers are read
y1_csv <- read.csv("y1.csv", header=TRUE) #read csv back in
all_y1 <- y1_csv %>% dplyr::select(1, contains("BL1.H")) #remove duplicate columns

####Diptest####
#Diptest with sapply - only need 2:ncol because 1st column is cell number
dip_y1<-sapply(all_y1[2:ncol(all_y1)], function(x) dip.test(x)$p.value)
write.csv(dip_y1, "test_pvals.csv", row.names = TRUE)



####Fitdist stats####

#test function, should get a mu and sigma value
get.fd.stats(y1_in[[1]]) 

#for loop - runtime <10s
stats_y1 <- data.frame(sample=NA) #new dataframe to write loop into
for(f in y1_in){stats_y1[[f]]<-cbind(get.fd.stats(f))}

#tidy fd stats and write as excel file
y1_fd <- tidy.fd.stats(stats_y1)
write.xlsx(y1_fd,'y1_fd_test.xlsx') #save



####Plot fitdist####

#melt large table to extract the control - easiest way to get it with the headers labelled
# (the csv alone does not have the sample identifying header which is why the large df will be melted)
melt_y1 <- melt(all_y1, id.vars = "cell_number", na.rm = TRUE)
head(melt_y1)

#extract control sample
#CHANGE 'B8' DEPENDING ON WELL OF CONTROL SAMPLE
#IMPORTANT: TO KEEP THE GET.NORM.PLOT FUNCTION THE SAME, 'ctrl' WILL BE OVERWRITTEN FOR
#           EACH MEDIA CONDITION
ctrl <- filter(melt_y1, grepl("B8",variable))
#head(ctrl) 

#for loop
#y1np='y1 media normal plots' -list object to store ggplots inside (requirement for grid.arrange())
y1np<-list() 
for(f in y1_in){y1np[[f]]<-get.fd.plots(f)}

#plot in a grid with 5 columns
y1_norm_plots <- grid.arrange(grobs=y1np,ncol=5, nrow=5,top="y1_norm")  #'top' adds a title to all grid arranged plots



####Flexmix plots#####

#for loop - 2 COMPONENT FM
#will get errors but plots still work
y1ff_2<-list()
for(f in y1_in){y1ff_2[[f]]<-get.flex.plots2(f)}

#plot and SAVE - y1 flexmix 2 components
y1ff_plots2 <- grid.arrange(grobs=y1ff_2,ncol=5,nrow=5,top="y1_2components") 


#for loop - 3 COMPONENT FM
y1ff_3<-list()
for(f in y1_in){y1ff_3[[f]]<-get.flex.plots3(f)}

#plot and SAVE - y1 flexmix 2 components
y1ff_plots3 <- grid.arrange(grobs=y1ff_3,ncol=5,nrow=5,top="y1_3components") 


####Flexmix stats####
get.flex.stats2(y1_in[[13]]) #test function

#'ffs2' = flexfit stats 2 COMPONENT
#for loop
y1_ffs2 <- data.frame(c(1,2,3,4,5,6,7))
for(f in y1_in){y1_ffs2[[f]]<-cbind(get.flex.stats2(f))}

#tidy output
y1_ffs2_tidy <- tidy.flex2(y1_ffs2)
write.xlsx(y1_ffs2_tidy, 'y1_2comp_stats.xlsx', rowNames=TRUE) #save


#3 COMPONENTS
get.flex.stats3(y1_in[[13]]) #test function

#'ffs3' = flexfit stats 3 COMPONENT
#for loop
y1_ffs3 <- data.frame(c(1,2,3,4,5,6,7,8,9,10))
for(f in y1_in){y1_ffs3[[f]]<-cbind(get.flex.stats3(f))}

#tidy output
y1_ffs3_tidy <- tidy.flex3(y1_ffs3)
write.xlsx(y1_ffs3_tidy, 'y1_3comp_stats.xlsx', rowNames=TRUE) #save


####Flexmix 2v3 fit####
y1_bic <- data.frame(c(2,3)) #started with 2,3 so that output table is clearer and points to 2 and 3 components
for(f in y1_in){y1_bic[[f]]<-get.BIC(f)}

y1_bict<- t(y1_bic)
y1_bict
write.csv(y1_bict, 'y1_bic.csv')



toc()




###pretty plots for thesis:###
#original function
get.flex.plots3<-function(f){
  csvin <- read.csv(f, header = TRUE, stringsAsFactors = FALSE)
  csvpos <- filter(csvin, csvin[1]>0) #remove rows containing negative values which would break log transform
  logged <- log10(csvpos)
  set.seed(1)
  ff<-flexmix(logged$`Comp.BL1.H....BL1.H`~1, k=3, model=list(mo1, mo2, mo3))
  c1 <- parameters(ff, component=1)[[1]]
  c2 <- parameters(ff, component=2)[[1]]
  c3 <-  parameters(ff, component=3)[[1]]
  c1df<- as.data.frame(c1)
  c2df <- as.data.frame(c2)
  c3df <- as.data.frame(c3)
  lam <- table(clusters(ff))
  ggplot(logged, aes(x=`Comp.BL1.H....BL1.H`)) + 
    geom_density(color="#A0A0A0") +
    stat_function(data=logged,
                  fun = plot_mix_comps,
                  args = list(c1df[1,1], c1df[2,1],lam[1]/sum(lam)),
                  color="#FF6666", fill="#FF6666",alpha=0.4, geom="polygon") +
    stat_function(data=logged, 
                  fun = plot_mix_comps,
                  args = list(c2df[1,1], c2df[2,1], lam[2]/sum(lam)),
                  color="#009999",fill="#00CCCC",alpha=0.35,geom="polygon") +
    stat_function(data=logged, 
                  fun = plot_mix_comps,
                  args = list(c3df[1,1], c3df[2,1], lam[3]/sum(lam)),
                  color="#9933FF",fill="#AB00FF",alpha=0.2,geom="polygon") +
    stat_function(fun=dnorm,
                  args=list(mean=mean(ctrl$value, na.rm=TRUE), 
                            sd=sd(ctrl$value, na.rm=TRUE)),
                  geom = "area",
                  color="#E0E0E0",
                  fill="grey", alpha=0.25) +
    theme_minimal() +
    scale_x_continuous(limits=c(2,6)) +
    scale_y_continuous(limits = c(0,3.75)) +
    xlab("")+
    ylab("")
}

#nice colours and removing background grids function
tgg2 <-theme(axis.text.x = element_text(family="Calibri Light", size=7.5),
            axis.title.y = element_blank(),
            axis.title.x = element_blank(),
            axis.text.y = element_text(family="Calibri Light", size=7.5),
            axis.ticks.y = element_line(colour="gray", size=0.5),
            axis.ticks.x = element_line(colour="gray", size=0.5),
            panel.grid.minor.y = element_blank(),
            panel.grid.major.y = element_blank(),
            panel.grid.minor.x=element_blank(),
            panel.grid.major.x = element_blank(),
            plot.margin = unit(c(0, 5.5, 0, 5.5), "pt"))

get.flex.plots3.pretty<-function(f){
  csvin <- read.csv(f, header = TRUE, stringsAsFactors = FALSE)
  csvpos <- filter(csvin, csvin[1]>0) #remove rows containing negative values which would break log transform
  logged <- log10(csvpos)
  set.seed(1)
  ff<-flexmix(logged$`Comp.BL1.H....BL1.H`~1, k=3, model=list(mo1, mo2, mo3))
  c1 <- parameters(ff, component=1)[[1]]
  c2 <- parameters(ff, component=2)[[1]]
  c3 <-  parameters(ff, component=3)[[1]]
  c1df<- as.data.frame(c1)
  c2df <- as.data.frame(c2)
  c3df <- as.data.frame(c3)
  lam <- table(clusters(ff))
  ggplot(logged, aes(x=`Comp.BL1.H....BL1.H`)) + 
    stat_function(data=logged,
                  fun = plot_mix_comps,
                  args = list(c1df[1,1], c1df[2,1],lam[1]/sum(lam)),
                  color="#C96643", fill="#C96643",alpha=0.6, geom="polygon") + #ORANGE
    stat_function(data=logged, 
                  fun = plot_mix_comps,
                  args = list(c2df[1,1], c2df[2,1], lam[2]/sum(lam)),
                  color="#84ADAD",fill="#84ADAD",alpha=0.6,geom="polygon") + #BLUE
    stat_function(data=logged, 
                  fun = plot_mix_comps,
                  args = list(c3df[1,1], c3df[2,1], lam[3]/sum(lam)),
                  color="#3B0535",fill="#3B0535",alpha=0.5,geom="polygon") + #PURPLE
    stat_function(fun=dnorm,
                  args=list(mean=mean(ctrl$value, na.rm=TRUE), 
                            sd=sd(ctrl$value, na.rm=TRUE)),
                  geom = "area",
                  color="#E0E0E0",
                  fill="grey", alpha=0.25) +
    geom_density(color="#A0A0A0") +
    theme_minimal() +
    scale_x_continuous(limits=c(2,6), expand = expansion(mult = c(0, 0))) +
    scale_y_continuous(limits = c(0,3.75)) +
    xlab("")+
    ylab("")+
    geom_hline(yintercept = 0, colour = "gray") +
    geom_vline(xintercept = 2, colour = "gray") +
    tgg2
}


#plot
y1ff_3_pretty<-list()
for(f in y1_in){y1ff_3_pretty[[f]]<-get.flex.plots3.pretty(f)}

#plot and SAVE - y1 flexmix 2 components
y1ff_plots3_pretty <- grid.arrange(grobs=y1ff_3_pretty,ncol=5,nrow=5,top="y1_3components") 
