## Graphing Q outputs from cluster analyses
# started 2024-02-21

# Setting up the environment -----

library(stringr)
library(tidyverse)
library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(microViz)

library(devtools)
install_github('royfrancis/pophelper')
library(pophelper)

# Importing data -----

# bring in vcf in case (thought I needed it to pulll individual labels)

Humber_vcf <- read.table("c:/Users/Katherine/OneDrive - University of Guelph/Post-doc - Mandeville/Cyprinid project/Humber_vcf.vcf")

# import Amanda's list of individual names 

Humber_indivs <- read.delim("c:/Users/Katherine/OneDrive - University of Guelph/Post-doc - Mandeville/Cyprinid project/names_humber.txt", sep = "\t", header = F)
names(Humber_indivs) <- "Mandeville_ID"

# import Humber dataset metadata

Humber_metadata <- read.csv("c:/Users/Katherine/OneDrive - University of Guelph/Post-doc - Mandeville/Cyprinid project/Humber_Leuciscid_Metadata_May2023_KD.csv", encoding="UTF-8")

# __bring in fastSTRUCTURE Q files ----

setwd("c:/Users/Katherine/OneDrive - University of Guelph/Post-doc - Mandeville/Cyprinid project/fastSTRUCTURE analysis/Humber dataset run")

# bringing in k=4 and k=7 as per chooseK.py results
Humber_fSQ4 <- read.delim("HumberfS_k4.4.meanQ", sep="", header=F)
Humber_fSQ7 <- read.delim("HumberfS_k7.7.meanQ", sep="", header=F)

# merging Q values, individuals, and metadata
Q4df <- cbind(Humber_fSQ4, Humber_indivs)
k4df <- merge(Q4df, Humber_metadata, by="Mandeville_ID")

Q7df <- cbind(Humber_fSQ7, Humber_indivs)
k7df <- merge(Q7df, Humber_metadata, by="Mandeville_ID")


# for certain applications, need the data in long format

k4df_long <- gather(k4df, Q_col, Q_score, V1:V2:V3:V4)
k7df_long <- gather(k7df, Q_col, Q_score, V1:V2:V3:V4)

# save these as .csv so I can just import them directly in the future

write.csv(k4df_long, paste0("Humber_fS_k4_longformat.csv"), row.names=FALSE)
write.csv(k7df_long, paste0("Humber_fS_k7_longformat.csv"), row.names=FALSE)


# __bring in ADMIXTURE Q files ----

setwd("c:/Users/Katherine/OneDrive - University of Guelph/Post-doc - Mandeville/Cyprinid project/ADMIXTURE analysis/Humber dataset run 2")

Humber_ADQ4 <- read.delim("Humbervcf_tobed.4.Q", sep="", header=F)

# merging Q values, individuals, and metadata (chose 4 initially because of FS results)
AQ4df <- cbind(Humber_ADQ4, Humber_indivs)
Ak4df <- merge(AQ4df, Humber_metadata, by="Mandeville_ID")

# for certain applications, need the data in long format

Ak4df_long <- gather(Ak4df, Q_col, Q_score, V1:V2:V3:V4)


# save these as .csv so I can just import them directly in the future

write.csv(Ak4df_long, paste0("Humber_ADM_k4_longformat.csv"), row.names=FALSE)

# _____choosing K

Humber_ADM_cv <- read.csv("Humber_ADM_run2_CV.csv")

Humber_ADM_cv2 <- gather(Humber_ADM_cv, key=CV_fold, value=CV_error, 2:3)

ggplot(Humber_ADM_cv2, aes(x=factor(K), y=CV_error, fill=CV_fold, color=CV_fold, group=CV_fold)) +
  geom_point() +
  geom_line() +
  theme_grey() +
  scale_x_discrete()

# Basic plotting ----

ggplot(k4df_long, aes(x=Mandeville_ID, y=Q_score, fill=Q_col))+
  geom_bar(stat="identity")

ggplot(k4df_long, aes(x=Mandeville_ID, y=Q_score, fill=Common_name))+
  geom_bar(stat="identity") # just to see what the species are

ggplot(Ak4df_long, aes(x=Mandeville_ID, y=Q_score, fill=Q_col))+
  geom_bar(stat="identity")

ggplot(Ak4df_long, aes(x=Mandeville_ID, y=Q_score, fill=Common_name))+
  geom_bar(stat="identity") # just to see what the species are

# Trying stuff to get the x-axis ordered the way we want ----

# rounding Q_score values to make this a little more straightforward


k4df_long_fixed <- k4df_long_fixed %>% mutate(Rounded=round(Q_score, 3))

ggplot(k4df_long_fixed, aes(x=reorder(Mandeville_ID, Rounded, sum), y=Rounded, fill=Q_col))+
  geom_bar(stat="identity")

# __trying Pophelper -----

# it requires the raw, original file, so need to set wd and then import
setwd("c:/Users/Katherine/OneDrive - University of Guelph/Post-doc - Mandeville/Cyprinid project/fastSTRUCTURE analysis/Humber dataset run")

# this works, but it doesn't have the labels
qlist_fS4 <- readQ(files="HumberfS_k4.4.meanQ")

# going to save my merged list with ind labels as a file, then re-import as Pophelper wants
Q4df_2 <- Q4df[,c(5,1,2,3,4)] # reorder the columns
rownames(Q4df_2) <- Q4df_2[,1] # make Mandeville_ID the rownames
Q4df_2$Mandeville_ID <- NULL # delete Mandeville_ID column

write.table(Q4df_2, "HumberfS_k4_indlabs.meanQ", sep=" ", row.names=TRUE, col.names = FALSE)

qlist_fS4_indlabs <- readQ(files="HumberfS_k4_indlabs.meanQ") # OKAY for some reason it's eating one of the columns???

test1 <- read.table("HumberfS_k4_indlabs.meanQ", row.names=1) # this looks fine, so I don't think it's the input file that is the problem ...
# posted on Roy Francis's github Issues page, will see if he responds



# _____plotting -----

# first try, on file without ind labels
plotQ(qlist_fS4, 
      sortind="all",
      height = 2,
      barbordersize = 0.5,
      barbordercolour = "white",
      showlegend = TRUE,
      showdiv = TRUE,
      divsize = 2,
      outputfilename = "HumberfSk4_1.1",
      imgtype = "png",
      exportpath=getwd())

plotQ(qlist_fS4_indlabs,
      sortind="all",
      outputfilename = "HumberfSk4_2.0",
      exportpath=getwd()) 


# not working yet ----

k4df_long$Q_col <- factor(k4df_long$Q_col,
                          levels=names(sort(colSums(k4df_long[,18:19]), decreasing=TRUE)))

# didn't work

group_k4df_long <- k4df_long %>% group_by(Q_col)
ggplot(group_k4df_long, aes(x=reorder(Mandeville_ID,Q_col), y=Q_score, fill=Q_col))+
  geom_bar(stat="identity")
