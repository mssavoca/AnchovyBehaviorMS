## Libraries
library(SDMTools)
library(jpeg)
library(raster)
library(spdep)
library(ggplot2)
library(gtools)
library(plyr)
library(dplyr)
library(data.table)
library(xlsx)
library(nlme)
library(predictmeans)
library(pbkrtest)
library(rethinking)
library(gridExtra)

################ 
## Read in screenshots and process
################
## Set working diretory to appropriate location
setwd(dir = "~/Google Drive/Anchovy feeding videos/Analysis/data/")

## Read summary table of videos
vt <- read.csv("Anchovy_Project_Video_Synopsis.csv", 
               header = T, 
               stringsAsFactors = F)

## Set up empty data.frame
output <- NULL

## List folders (days)
fdays <- list.files()

## For each folder
for (days in 1:length(fdays)){
  
  ## List folders
  fd <- list.files(paste("./", fdays[days], "/", sep = ""))
  
  ## For each folder
  for (folder in 1:length(fd)){
    
    ## List files
    fis <- list.files(paste("./", fdays[days], "/", fd[folder], "/Black and white screenshots", sep = ""))
    
    ## If more than 50 files in folder, proceed
    if(length(fis) > 50) {
      
      ## Remove thumbs.db for each folder
      fis.f <- fis[-length(fis)]
      
      ## Read in jpgs for each file in each folder
      jpgs <- lapply(fis.f, function(x){readJPEG(paste("./", fdays[days], "/", fd[folder], "/Black and white screenshots/", x, sep = ""), native = T)})
      
      ## Rows
      cols <- nrow(jpgs[[1]]) #extracts the dimensions for raster/matrix analysis
      
      ## Cols
      rows <- ncol(jpgs[[1]]) #extracts the dimensions for raster/matrix analysis
      
      ## Change to matrices
      ms <- lapply(jpgs, function(x){matrix(x, nrow = rows, ncol = cols)})
      
      ## Then change to rasters
      rs <- lapply(ms, function(x){raster(x)})
      
      ## All values >= 0 and < 1 become 1. So: here the "fish" are colored as patch 1
      m <- c(-Inf, -5000000, 1, -5000000, Inf, NA)
      rclmat <- matrix(m, ncol=3, byrow=TRUE)
      rs.rcl <- lapply(rs, function(x){reclassify(x, rclmat, right = T)})
      
      ## Convert to ascs files
      ascs <- lapply(rs.rcl, function(x){asc.from.raster(x)})
      
      ## Calculate some spatial statistics for each image
      sp.sts <- lapply(ascs, function(x){ClassStat(x)})
      
      ## Rbind spatial statistics
      sp.df <- do.call(rbind, sp.sts)
      
      ## Add column "pic" based on filename
      sp.df$pic <- fis.f
      
      ## Add time index
      sp.df$time <- 1:nrow(sp.df)
      
      ## Add treatment type
      sp.df$treatment <- vt[vt$Video.name  == fd[folder],]$Treatment
      
      ## Add day
      sp.df$day <- fdays[days]
      
      ## Rbind to master output df
      output <- rbind(output, sp.df)
      
    }
    
  }
  
}

## Save file
write.csv(output, file = "Anchovy_project_spatial_stats_master_output.csv", row.names = F)

################ 
## Read in the processed data and add some extra info
################ 

## Read in master data
output_master <- read.csv("Anchovy_project_spatial_stats_master_output.csv", 
                          header = T, 
                          stringsAsFactors = F)
output_master$trial_uni_id <- group_indices(output_master, treatment, pic)

## Unique trial ID within each treatment
output_master <- setDT(output_master)[, trial_id := rleid(pic), treatment][]

## Create new label based on time interval, pre-treatment (first 2 mins), post-treatment (next 2 mins), >2 mins post-treatment (last 10 mins)
labels <- data.table(label = c("pre", "post", "last"), time = c(0,24,52), key = "time")

## Change to data.table
output.dt <- data.table(output_master, key = "time") #CHANGE df FOR DIFFERENT DATA FRAMES

## Merge to give label
merge <- labels[output.dt, roll = Inf]

## Change label to ordered
merge$label <- ordered(merge$label, c("pre", "post", "last"))

## Rename
out_sumr <- merge; rm(merge)

## Calculate Moran's I over first two minutes to give the 'baseline'
out_moran <-
  out_sumr %>% 
  dplyr::group_by(pic) %>%
  do(baseline_mean = with(filter(., time < 25), mean(moran)),
     baseline_sd = with(filter(., time < 25), sd(moran)))

## Change to vector
out_moran$baseline_mean <- unlist(out_moran$baseline_mean)
out_moran$baseline_sd <- unlist(out_moran$baseline_sd)

## Merge with output
out_sumr <- plyr::join(out_sumr, out_moran, by = "pic")

## Subtract baseline from each time step
out_sumr$moran_zscore <- (out_sumr$moran - out_sumr$baseline_mean)/out_sumr$baseline_sd

## Read in video synopsis file
vt <- read.xlsx("Anchovy_Project_Video_Synopsis.xlsx",
                sheetIndex = 1, 
                stringsAsFactors = F)

## Sort
vt <- arrange(vt, Treatment, Video.name)

## Create batch column in vt. This is a unique number for each video in each treatment. This is just a sequential and unique number for each video within each treatment.
vt <- vt %>%
  dplyr::group_by(Treatment) %>%
  dplyr::mutate(Batch = 1:length(unique(Video.name)))

## Add batch column from vt to 
vt_short <- select(vt, Video.name, Batch)

## Need to rename video column in Anchovy Video Synopsis file for merging with out_sumr, since pic IS video file name in out_sumr
vt_short <-  dplyr::rename(vt_short, pic = Video.name) 

## Merge
out_sumr <- merge(out_sumr, vt_short, by = "pic")

## Create and add a vector with the second sequence repeated the correct number of times to fit the full data.frame
out_sumr$time2 <- with(out_sumr, 
                       ifelse(time<=25, (time-25)*5,
                              ifelse(time ==25, 0, (time-25)*5)))

## Save
write.csv(out_sumr, 
          "Anchovy_project_spatial_stats_master_output_full.csv",
          row.names = F)

################ 
## Read in the new data and rearrange
################ 

## Read in data
output_master <- read.csv("Anchovy_project_spatial_stats_master_output_full.csv", header = T, stringsAsFactors = F)

## Select columns and sort
output_abridged <- output_master %>%
  select(day,
         treatment,
         trial_id,
         trial_uni_id,
         batch,
         pic,
         time,
         time2,
         moran_zscore,
         label) %>%
  arrange(day,
          treatment,
          time)

## Remove empty traetment
output_abridged <- output_abridged[is.na(output_abridged$treatment) == F,]

## Save
write.csv(output_abridged, 
          "Anchovy_project_spatial_stats_master_output_abridged.csv",
          row.names = F)

################ 
## Analyses of aggregation data
################ 

## Read in data
output_master <- read.csv("Anchovy_project_spatial_stats_master_output_abridged.csv", 
                          header = T, 
                          stringsAsFactors = F)
vt <- read.csv("anchovy_project_video_synopsis.csv", 
               header = T, 
               stringsAsFactors = F)

# Subsetting
before_treat <- filter(output_master, time <= 24)
after_treat <- filter(output_master, time >= 25)

## Bootstrapped, repeated measures ANOVA
RM_Moran = aov(moran_zscore ~ (treatment*time.bin) + Error(trial_uni_id/day), data=after_treat)
summary(RM_Moran)

## Set correlation structure
before_treat <- before_treat %>%
  arrange(pic, time) 
cor.b <- with(before_treat, (cor(x = head(moran_zscore, -1), tail(moran_zscore, -1))))
cor.structure <- corCAR1(value = cor.b, form = ~time|day)

## Pre model
out_sumr$treatment <- as.factor(as.character(out_sumr$treatment))
out_sumr$treatment <- relevel(out_sumr$treatment, ref = "SW Control") 
out_sumr$label <- as.factor(as.character(out_sumr$label))
out_sumr$label <- relevel(out_sumr$label, ref = "pre") 

##GLMM
MoranTest = lmer(moran_zscore ~ treatment*label -1 + (1|day) + (1|Batch),
                 data = out_sumr)

## Permutation test of the GLMM
permMoranTest <- permmodels(model=MoranTest, data=out_sumr, block = "Batch", group= "day", nsim=10000) 

##############
#Time to peak analysis
##############

## Read in combined data
out_comb <- read.csv("Anchovy_project_spatial_stats_combined.csv", 
                     header = T, 
                     stringsAsFactors = F)

## Calculate time to peak response for each picture
out_comb <-
  out_comb %>% 
  group_by(pic, measure) %>% 
  mutate(peak_response = max(value, na.rm = T),
         max_response_time = time[which.max(value == max(value))])

#this makes one line per video if to analyze max response intensity and time to peak intensity
out_comb_4a = ddply(out_comb, "pic", head, 1)

#need to join this with the out_sumr table to have random effects by day
out_sumr_4a = ddply(out_sumr, "pic", head, 1)
out_comb_4a = merge(x = out_comb_4a, y = out_sumr_4a, by = "pic", all = TRUE)


################ 
## Prepare treatment labels for out_comb_4a
################ 

## Renaming and reordering treatments
out_comb_4a$treatment <- ifelse(grepl("SW Control", out_comb_4a$treatment.x), "Seawater Control", out_comb_4a$treatment.y) 

out_comb_4a$treatment <- ordered(out_comb_4a$treatment, levels = c("Seawater Control", "Clean Plastic (low)", "Clean Plastic (med)", "Clean Plastic (high)", "Krill odor (low)", "Krill odor (med)",  "Bf Plastic (high)", "Krill odor (high)", "Bf Plastic (low)", "Bf Plastic (med)", "Food"))

## Concentration label
out_comb_4a$concentration <- ifelse(grepl("high", out_comb_4a$treatment.x), "High", 
                                    ifelse(grepl("med", out_comb_4a$treatment.x), "Medium",
                                           ifelse(grepl("low", out_comb_4a$treatment.x), "Low", NA)))

## Order
out_comb_4a$concentration <- ordered(out_comb_4a$concentration, levels = c("High", "Medium", "Low"))

## General treatment type
out_comb_4a$treament_general <- ifelse(grepl("Clean", out_comb_4a$treatment.x), "Clean Plastic",     #Yes treatment is misspelled here
                                       ifelse(grepl("odor", out_comb_4a$treatment.x), "Food Odor",
                                              ifelse(grepl("Bf", out_comb_4a$treatment.x), "Biofouled Plastic", 
                                                     NA)))
View(out_comb_4a)

# GLMM for Time To Peak response, remember that the estimates are the intercept COMBINED with the estimate
MoranTTP = lmer(max_response_time ~ treatment -1 + (1|day) + (1|Batch),
                data = out_comb_4a)  #from Anchovy_plots_CTfinal
summary(MoranTTP)

# to get raw data for Time To Peak
TTP_out = out_comb_4a %>%
  group_by (treatment.y) %>%  
  summarize(MeanTTP=median(max_response_time), StdEr = SE(max_response_time))


# GLMM for Intensity of Peak response, remember that the estimates are the intercept COMBINED with the estimate
MoranIP = lmer(peak_response ~ treatment + (1|day) + (1|Batch),
               data = out_comb_4a)  #from Anchovy_plots_CTfinal
summary(MoranIP)


################ 
## Analyses of rheotactic data
################ 

# writing a function to quantify logistic output from binomial models
logistic = function(x){1/(1+exp(1)^-x)}

#create function for standard error
SE = function(x){sd(x)/sqrt(sum(!is.na(x)))}

# Read in data
vt <- read.csv("anchovy_project_video_synopsis.csv", 
               header = T, 
               stringsAsFactors = F)
RT_master <- read.csv("rheotaxisAnalysis_Master.csv",
                      header = T, 
                      stringsAsFactors = F)

## Clean up data strutures
RT_master <- na.omit(RT_master) #remove NA rows
Proportion.fish.non.rheotactic = as.numeric(as.character(RT_master$Proportion.fish.non.rheotactic)) # prop non-rheotactic being read in as a factor, need to change this to a number
RT_master$Proportion.fish.non.rheotactic <- Proportion.fish.non.rheotactic
X..Non.rheotactic = as.numeric(as.character(RT_master$X..Non.rheotactic)) # number non-rheotactic being read in as a factor, need to change this to a number
RT_master$X..Non.rheotactic <- X..Non.rheotactic
Total...fish = as.numeric(as.character(RT_master$Total...fish)) # total number fish being read in as a factor, need to change this to a number
RT_master$Total...fish <- Total...fish
RT_master$Treatment <- as.factor(as.character(RT_master$Treatment))
RT_master$Treatment <- relevel(RT_master$Treatment, ref = "SW Control")
RT_master$Date.recorded <- as.factor(RT_master$Date.recorded)
RT_master$Screenshot.time <- as.factor(RT_master$Screenshot.time)
str(RT_master)

# counting how many replicates of each treatment there are
RT_master_count = RT_master %>% 
  group_by(Video.name, Treatment) %>%
  summarise(Unique_Elements = n_distinct(Video.name, Treatment))
count(RT_master_count, Treatment)

# create a sequence of values for second mark in video, relative to injection
time.rel.inj <- seq(-80, 100, by = 5)

# counts the number of videos in the data frame to know how many times to repeat the above sequence
num.uni.vids <- length(unique(as.character(RT_master$Video.name)))

# create a vector with the second sequence repeated the correct number of times to fit the full dataframe
sec.rel.inj <- rep(time.rel.inj, num.uni.vids) 

# append vector to dataframe
RT_master$sec.rel.inj <- sec.rel.inj

RT_master$Pre.or.post <- as.factor(as.character(RT_master$Pre.or.post))
RT_master$Pre.or.post <- relevel(RT_master$Pre.or.post, ref = "pre") # relevels labal (time:pre,post,last) factors relative to pre-treatment

## Trying to zscore transform Proportion non-Rheotactic to address rank deficiency
## Calculate prop fish non-rheotactic before solution injection to give the 'pre'
out_rheo <-
  RT_master %>% 
  dplyr::group_by(Video.name) %>%
  do(pre_mean = with(filter(., Screenshot.order < 25), mean(Proportion.fish.non.rheotactic)),
     pre_sd = with(filter(., Screenshot.order < 25), sd(Proportion.fish.non.rheotactic)))

## Change to vector
out_rheo$pre_mean <- unlist(out_rheo$pre_mean)
out_rheo$pre_sd <- unlist(out_rheo$pre_sd)

## Merge with output
RT_master <- plyr::join(RT_master, out_rheo, by = "Video.name")

## Subtract baseline from each time step
RT_master$PropNonRheo_zscore <- (RT_master$Proportion.fish.non.rheotactic - RT_master$pre_mean)/RT_master$pre_sd

#GLMM and Permutation test
m1 = glmer(PropNonRheo_zscore ~ Treatment*Pre.or.post + (1|Date.recorded) + (1|Observer),
           data=RT_master)
permRheoTest <- permmodels(model=m1, data=RT_master, block = "Observer", group= "Date.recorded", nsim=10000)
summary(permRheoTest)



################################################ 
## Figures
################################################ 

################ 
## Read in the data 
################ 

## Read in data
out_sumr <- read.csv("Anchovy_project_spatial_stats_master_output_abridged.csv", 
                     header = T, 
                     stringsAsFactors = F)

################ 
## Prepare treatment labels
################ 

## Renaming and reordering treatments
out_sumr$treatment <- ifelse(grepl("SW Control", out_sumr$treatment), "Seawater Control", out_sumr$treatment) 

out_sumr$treatment <- ordered(out_sumr$treatment, levels = c("Seawater Control", "Clean Plastic (low)", "Clean Plastic (med)", "Clean Plastic (high)", "Krill odor (low)", "Krill odor (med)",  "Bf Plastic (high)", "Krill odor (high)", "Bf Plastic (low)", "Bf Plastic (med)", "Food"))

## Concentration label
out_sumr$concentration <- ifelse(grepl("high", out_sumr$treatment), "High", 
                                 ifelse(grepl("med", out_sumr$treatment), "Medium",
                                        ifelse(grepl("low", out_sumr$treatment), "Low", NA)))

## Order
out_sumr$concentration <- ordered(out_sumr$concentration, levels = c("High", "Medium", "Low"))

## General treatment type
out_sumr$treament_general <- ifelse(grepl("Clean", out_sumr$treatment), "Clean Plastic", 
                                    ifelse(grepl("odor", out_sumr$treatment), "Food Odor",
                                           ifelse(grepl("Bf", out_sumr$treatment), "Biofouled Plastic", 
                                                  NA)))

## Order
out_sumr$treament_general <- ordered(out_sumr$treament_general, levels = c("Clean Plastic", "Biofouled Plastic", "Food Odor"))

################ 
## Time series of clustering
################ 

## Graph of the time series for just treaments with concentrations
(ts.conc <- ggplot(data = out_sumr[is.na(out_sumr$treament_general) == F & out_sumr$time2 < 600,],
                   aes(x= time2, y = moran_zscore, colour = as.factor(treament_general))) +
   geom_point(alpha = 1/5) +
   geom_smooth(data = filter(out_sumr[is.na(out_sumr$treament_general) == F & out_sumr$time2 < 600,], time2 < 0), 
               se = T, span = 0.9, size = 0.5, colour = "black") +
   geom_smooth(data = filter(out_sumr[is.na(out_sumr$treament_general) == F & out_sumr$time2 < 600,], time2 > 0), 
               se = T, span = 0.9, size = 0.5, colour = "black") +
   geom_vline(xintercept = 0, linetype="dashed") +
   ylim(-5, 30) +
   scale_x_continuous(breaks = seq(-120, 600, by = 120)) + 
   #xlab("Time (sec.) relative to treatment presentation") +
   ylab("Standardized Aggregation \n Metric (Moran's I)") +
   scale_color_manual(name = "Treatment", values = c("goldenrod2", "seagreen3", "darksalmon")) +
   theme_bw() +
   theme(plot.title = element_text(face="bold",vjust=30, size = 18),
         axis.title.x = element_text(face="bold", vjust=-0.25, size=12),
         axis.title.y = element_text(face="bold", vjust=14, size=12),
         axis.text.y  = element_text(vjust=0.5, size=10),
         axis.text.x  = element_text(vjust=0.5, size=10, angle = 45),
         strip.text.x = element_text(size = 12),
         legend.position = "none") + 
   facet_grid(concentration ~ treament_general, scales = "free_x") +
   annotate("text", x = -100, y = 25, label = c("C", "D", "E", "", "", "", "", "", ""), fontface = 2))

## Graph of the time series for SW control and food
(ts.sw.food <- ggplot(data = out_sumr[out_sumr$treatment == c("Seawater Control", "Food") & out_sumr$time2 < 600,],
                      aes(x= time2, y = moran_zscore, colour = as.factor(treatment))) +
    geom_point(alpha = 1/5) +
    geom_smooth(data = filter(out_sumr[out_sumr$treatment == c("Seawater Control", "Food") & out_sumr$time2 < 600,], time2 < 0), 
                se = T, span = 0.9, size = 0.5, colour = "black") +
    geom_smooth(data = filter(out_sumr[out_sumr$treatment == c("Seawater Control", "Food") & out_sumr$time2 < 600,], time2 > 0), 
                se = T, span = 0.9, size = 0.5, colour = "black") +
    geom_vline(xintercept = 0, linetype="dashed") +
    ylim(-5, 30) +
    scale_x_continuous(breaks = seq(-120, 600, by = 120)) + 
    labs(x = NULL,
         y = "Standardized Aggregation \n  Metric (Moran's I)") +
    scale_color_manual(name = "Treatment", values = c("darkcyan", "coral4")) +
    theme_bw() +
    theme(plot.title = element_text(face="bold",vjust=30, size = 18),
          axis.title.x = element_text(face="bold", vjust=-0.25, size=12),
          axis.title.y = element_text(face="bold", vjust=14, size=12),
          axis.text.y  = element_text(vjust=0.5, size=10),
          axis.text.x  = element_text(vjust=0.5, size=10, angle = 45),
          strip.text.x = element_text(size = 12),
          legend.position = "none") +
    facet_grid(. ~ treatment, scales = "free_x", drop = T) +
    annotate("text", x = -100, y = 25, label = c("A", "B"), fontface = 2))

## Multi panel of time series
mp.fig <- grid.arrange(ts.sw.food, ts.conc, 
                       nrow = 2, heights=c(0.6, 1))

## Save
ggsave(filename = "/Users/Tyson/Google Drive/Anchovy feeding videos/Analysis/figures/moran_timeseries.pdf", 
       plot = mp.fig,
       width= 160,
       height=174,
       dpi = 300,
       units = "mm")

################ 
## Time series of rheotactic behavior
################ 

## Read in rheotaxis data
rh <- read.csv("RheotaxisAnalysis_Master.csv", header = T, stringsAsFactors = F)
rh <- rh[complete.cases(rh),]
rh$treatment <- rh$Treatment
rh$Treatment <- NULL
rh$time <- as.numeric(as.difftime(rh$Screenshot.time, format = "%M:%S", units = "secs"))
rh$time_adj <- rh$time - 120

## Organize data

## Renaming and reordering treatments
rh$treatment <- ifelse(grepl("SW Control", rh$treatment), "Seawater Control", rh$treatment) 

rh$treatment <- ordered(rh$treatment, levels = c("Seawater Control", "Clean Plastic (low)", "Clean Plastic (med)", "Clean Plastic (high)", "Krill odor (low)", "Krill odor (med)",  "Bf Plastic (high)", "Krill odor (high)", "Bf Plastic (low)", "Bf Plastic (med)", "Food"))
# out_sumr$treatment2 <- factor(out_sumr$treatment, labels = c("Seawater Control","Clean Plastic Odor", "Biofouled Plastic Odor", "Food Odor", "Food"))

## Concentration label
rh$concentration <- ifelse(grepl("high", rh$treatment), "High", 
                           ifelse(grepl("med", rh$treatment), "Medium",
                                  ifelse(grepl("low", rh$treatment), "Low", NA)))

## Order
rh$concentration <- ordered(rh$concentration, levels = c("High", "Medium", "Low"))

## General treatment type
rh$treament_general <- ifelse(grepl("Clean", rh$treatment), "Clean Plastic", 
                              ifelse(grepl("odor", rh$treatment), "Krill Odor",
                                     ifelse(grepl("Bf", rh$treatment), "Biofouled Plastic", NA)))

## Order
rh$treament_general <- ordered(rh$treament_general, levels = c("Clean Plastic", "Biofouled Plastic", "Krill Odor"))


times <- seq(-60, 90, by = 30)

## Read in rheotaxis data
rh <- read.csv("RheotaxisAnalysis_Master.csv", header = T, stringsAsFactors = F)
rh <- rh[complete.cases(rh),]
rh$treatment <- rh$Treatment
rh$Treatment <- NULL
rh$time <- as.numeric(as.difftime(rh$Screenshot.time, format = "%M:%S", units = "secs"))
rh$time_adj <- rh$time - 120

## Organize data

## Renaming and reordering treatments
rh$treatment <- ifelse(grepl("SW Control", rh$treatment), "Seawater Control", rh$treatment) 

rh$treatment <- ordered(rh$treatment, levels = c("Seawater Control", "Clean Plastic (low)", "Clean Plastic (med)", "Clean Plastic (high)", "Krill odor (low)", "Krill odor (med)",  "Bf Plastic (high)", "Krill odor (high)", "Bf Plastic (low)", "Bf Plastic (med)", "Food"))

## Concentration label
rh$concentration <- ifelse(grepl("high", rh$treatment), "High", 
                           ifelse(grepl("med", rh$treatment), "Medium",
                                  ifelse(grepl("low", rh$treatment), "Low", NA)))

## Order
rh$concentration <- ordered(rh$concentration, levels = c("High", "Medium", "Low"))

## General treatment type
rh$treament_general <- ifelse(grepl("Clean", rh$treatment), "Clean Plastic", 
                              ifelse(grepl("odor", rh$treatment), "Food Odor",
                                     ifelse(grepl("Bf", rh$treatment), "Biofouled Plastic", NA)))

## Order
rh$treament_general <- ordered(rh$treament_general, levels = c("Clean Plastic", "Biofouled Plastic", "Food Odor"))


times <- seq(-60, 90, by = 30)

## Plot of the time series of proportion of non-rheotactic fish for just treaments with concentrations
(rh.conc <- ggplot(data = rh[is.na(rh$treament_general) == F,],
                   aes(x= time_adj, y = 1-Proportion.fish.non.rheotactic, colour = as.factor(treament_general))) +
    geom_point(alpha = 1/5) +
    geom_smooth(data = filter(rh[is.na(rh$treament_general) == F,], time_adj >= 0), 
                se = T, span = 0.9, size = 0.5, colour = "black") +
    geom_smooth(data = filter(rh[is.na(rh$treament_general) == F,], time_adj <= 0), 
                se = T, span = 0.9, size = 0.5, colour = "black") +
    geom_vline(xintercept = 0, linetype="dashed") +
    ylim(0, 1) +
    labs(x = "Time (sec.) relative to treatment presentation",
         y = "Proportion of Fish Rheotactic") +
    scale_x_continuous(breaks = times) +
    scale_color_manual(name = "Treatment", values = c("goldenrod2", "seagreen3", "darksalmon")) +
    theme_bw() +
    theme(plot.title = element_text(face="bold",vjust=30, size = 18),
          axis.title.x = element_text(face="bold", vjust=-0.25, size=12),
          axis.title.y = element_text(face="bold", vjust=0, size=12),
          axis.text.y  = element_text(vjust=0.5, size=10),
          axis.text.x  = element_text(vjust=0.5, size=10, angle = 45),
          strip.text.x = element_text(size = 12),
          legend.position = "none") + 
    facet_grid(concentration ~ treament_general, scales = "free_x") + 
    annotate("text", x = -70, y = 0.07, label = c("C", "D", "E", "", "", "", "", "", ""), fontface = 2))


## Graph of the time series for SW control and food
(rh.sw.food <- ggplot(data = rh[rh$treatment == c("Seawater Control", "Food") ,],
                      aes(x= time_adj, y = 1-Proportion.fish.non.rheotactic, colour = as.factor(treatment))) +
    geom_point(alpha = 1/5) +
    geom_smooth(data = filter(rh[rh$treatment == c("Seawater Control", "Food"),], time_adj >= 0), 
                se = T, span = 0.9, size = 0.5, colour = "black") +
    geom_smooth(data = filter(rh[rh$treatment == c("Seawater Control", "Food"),], time_adj <= 0), 
                se = T, span = 0.9, size = 0.5, colour = "black") +
    geom_vline(xintercept = 0, linetype="dashed") +
    ylim(0, 1) +
    scale_x_continuous(breaks = times) +
    labs(x = NULL,
         y = "Proportion of Fish Rheotactic") +
    scale_color_manual(name = "Treatment", values = c("darkcyan", "coral4")) +
    theme_bw() +
    theme(plot.title = element_text(face="bold",vjust=30, size = 18),
          axis.title.x = element_text(face="bold", vjust=-0.25, size=12),
          axis.title.y = element_text(face="bold", vjust=0, size=12),
          axis.text.y  = element_text(vjust=0.5, size=10),
          axis.text.x  = element_text(vjust=0.5, size=10, angle = 45),
          strip.text.x = element_text(size = 12),
          legend.position = "none") +
    facet_grid(. ~ treatment, scales = "free_x", drop = T) +
    annotate("text", x = -70, y =0.05, label = c("A", "B"), fontface = 2))

## Multi panel of time series
(mp.fig <- grid.arrange(rh.sw.food, rh.conc,
                        nrow = 2, heights=c(0.6, 1)))

## Save
ggsave(filename = "/Users/Tyson/Google Drive/Anchovy feeding videos/Analysis/figures/moran_timeseries.pdf", 
       plot = mp.fig,
       width= 200,
       height=215,
       dpi = 300,
       units = "mm")

################ 
## Boxplots of time to peak clustering
################ 

## Read in combined data
out_comb <- read.csv("Anchovy_project_spatial_stats_combined.csv", 
                     header = T, 
                     stringsAsFactors = F)

## Calculate time to peak response for each picture
out_comb <-
  out_comb %>% 
  group_by(pic, measure) %>% 
  mutate(peak_response = max(value, na.rm = T),
         max_response_time = time[which.max(value == max(value))])

################ 
## Prepare treatment labels
################ 

## Renaming and reordering treatments
out_comb$treatment <- ifelse(grepl("SW Control", out_comb$treatment), "Seawater Control", out_comb$treatment) 

out_comb$treatment <- ordered(out_comb$treatment, levels = c("Seawater Control", "Clean Plastic (low)", "Clean Plastic (med)", "Clean Plastic (high)", "Krill odor (low)", "Krill odor (med)",  "Bf Plastic (high)", "Krill odor (high)", "Bf Plastic (low)", "Bf Plastic (med)", "Food"))
# out_sumr$treatment2 <- factor(out_sumr$treatment, labels = c("Seawater Control","Clean Plastic Odor", "Biofouled Plastic Odor", "Food Odor", "Food"))

## Concentration label
out_comb$concentration <- ifelse(grepl("high", out_comb$treatment), "High", 
                                 ifelse(grepl("med", out_comb$treatment), "Medium",
                                        ifelse(grepl("low", out_comb$treatment), "Low", NA)))

## Order
out_comb$concentration <- ordered(out_comb$concentration, levels = c("High", "Medium", "Low"))

## General treatment type
out_comb$treament_general <- ifelse(grepl("Clean", out_comb$treatment), "Clean Plastic", 
                                    ifelse(grepl("odor", out_comb$treatment), "Krill Odor",
                                           ifelse(grepl("Bf", out_comb$treatment), "Biofouled Plastic", 
                                                  NA)))

## Order
out_comb$treament_general <- ordered(out_comb$treament_general, levels = c("Clean Plastic", "Biofouled Plastic", "Krill Odor"))

## Only retain clustering
out_comb <- out_comb[out_comb$measure == "moran"]

## Boxplots for time to peak response for each treatment with concentrations
(time.to.peak.plot.conc <- ggplot(data = out_comb[is.na(out_comb$treament_general) == F,]) + 
    geom_boxplot(aes(x = treament_general, y = max_response_time, fill = treament_general)) +
    labs(x = NULL, 
         y = "Time to Peak \n Clustering (seconds)") +
    theme_bw() +  
    theme(axis.title.x = element_text(face="bold", vjust=-0.25, size=12),
          axis.text.y  = element_text(vjust=0.5, size=12),
          axis.title.y = element_text(face="bold", vjust=14, size=12),
          axis.text.y  = element_text(vjust=0.5, size=12),
          # axis.text.x  = element_text(vjust=0.5, size=8),
          strip.text.x = element_text(size = 12),
          strip.text.y = element_text(size = 12),
          legend.position="none",
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.title.x = element_blank(),
          plot.title = element_text(face="bold",vjust=30, size = 18)) +
    facet_grid(concentration ~ treament_general, scales = "free_x", drop = T) +
    scale_fill_manual(name = "Treatment", values = c("khaki4", "seagreen3", "darksalmon")) +
    annotate("text", x = 0.6, y = 535, label = c("C", "D", "E", "", "", "", "", "", ""), fontface = 2))


## Graph of the time to peak response for SW control and food
(time.to.peak.plot.sw.food <- ggplot(data = out_comb[out_comb$treatment %in% c("Seawater Control", "Food"),]) + 
    geom_boxplot(aes(x = treatment, y = max_response_time, fill = treatment)) +
    labs(x = NULL, y = "Time to Peak \n Clustering (seconds)") +
    theme_bw() +
    theme(axis.title.x = element_text(face="bold", vjust=-0.25, size=12),
          axis.text.y  = element_text(vjust=0.5, size=12),
          axis.title.y = element_text(face="bold", vjust=14, size=12),
          axis.text.y  = element_text(vjust=0.5, size=10),
          strip.text.x = element_text(size = 12),
          strip.text.y = element_text(size = 12),
          legend.position="none",
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.title.x = element_blank(),
          plot.title = element_text(face="bold",vjust=30, size = 18)) +
    facet_grid(. ~ treatment, scales = "free_x", drop = T) +
    scale_fill_manual(name = "Treatment", values = c("darkcyan", "coral4")) +
    annotate("text", x = 0.6, y = 535, label = c("A", "B"), fontface = 2))

## Multi panel of time to peak response
(time.peak.response.boxplot <- grid.arrange(time.to.peak.plot.sw.food, time.to.peak.plot.conc, 
                                            nrow = 2, 
                                            heights=c(0.6, 1)))

## Save
ggsave(filename = "/Users/Tyson/Google Drive/Anchovy feeding videos/Analysis/figures/time_peak_response_boxplot.pdf", 
       plot = time.peak.response.boxplot,
       width=174,
       height=174,
       dpi = 300,
       units = "mm")

## Violin and Boxplots for time to peak response for each treatment with concentrations
(time.to.peak.plot.conc <- ggplot(data = out_comb_4a[is.na(out_comb_4a$treament_general) == F,]) + 
    geom_violin(aes(x = treament_general, y = max_response_time, fill = treament_general), width = 0.6) +
    geom_boxplot(aes(x = treament_general, y = max_response_time, fill = treament_general), width = 0.2) +
    labs(x = NULL, 
         y = "Time to Peak \n Aggregation (sec.)") +
    theme_bw() +  
    theme(axis.title.x = element_text(face="bold", vjust=-0.25, size=12),
          axis.text.y  = element_text(vjust=0.5, size=12),
          axis.title.y = element_text(face="bold", vjust=14, size=12),
          axis.text.y  = element_text(vjust=0.5, size=8),
          # axis.text.x  = element_text(vjust=0.5, size=8),
          strip.text.x = element_text(size = 12),
          strip.text.y = element_text(size = 12),
          legend.position="none",
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.title.x = element_blank(),
          plot.title = element_text(face="bold",vjust=30, size = 18)) +
    facet_grid(concentration ~ treament_general, scales = "free_x", drop = T) +
    scale_fill_manual(name = "Treatment", values = c("goldenrod2", "seagreen3", "darksalmon")) +
    annotate("text", x = 0.6, y = 535, label = c("C", "D", "E", "", "", "", "", "", ""), fontface = 2))


## Graph of the time to peak response for SW control and food
(time.to.peak.plot.sw.food <- ggplot(data = out_comb_4a[out_comb_4a$treatment.y %in% c("Seawater Control", "Food"),]) + 
    geom_violin(aes(x = treatment.y, y = max_response_time, fill = treatment.y), width = 0.6) +
    geom_boxplot(aes(x = treatment.y, y = max_response_time, fill = treatment.y), width = 0.2) +
    labs(x = NULL, y = "Time to Peak \n Aggregation (sec.)") +
    theme_bw() +
    theme(axis.title.x = element_text(face="bold", vjust=-0.25, size=12),
          axis.text.y  = element_text(vjust=0.5, size=12),
          axis.title.y = element_text(face="bold", vjust=14, size=12),
          axis.text.y  = element_text(vjust=0.5, size=8),
          strip.text.x = element_text(size = 12),
          strip.text.y = element_text(size = 12),
          legend.position="none",
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.title.x = element_blank(),
          plot.title = element_text(face="bold",vjust=30, size = 18)) +
    facet_grid(. ~ treatment.y, scales = "free_x", drop = T) +
    scale_fill_manual(name = "Treatment", values = c("darkcyan", "coral4")) +
    annotate("text", x = 0.6, y = 535, label = c("A", "B"), fontface = 2))

## Multi panel of time to peak response
(time.peak.response.ViolinBoxplot <- grid.arrange(time.to.peak.plot.sw.food, time.to.peak.plot.conc, 
                                                  nrow = 2, 
                                                  heights=c(0.6, 1)))

## Save
ggsave(filename = "/Users/matthewsavoca/Documents/Research Data/Anchovy behavior project/Figures/moran_timetopeakVB.pdf", 
       plot = time.peak.response.ViolinBoxplot,
       width= 160,
       height=174,
       dpi = 300,
       units = "mm")
