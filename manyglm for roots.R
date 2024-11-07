library(dplyr)
library(tidyr)
library(ggplot2)
library(car)
library(DHARMa)
library(vegan)
library(magrittr)
library(mvabund)
library(ggplot2)
library(moments)
library(emmeans)
library(writexl)
library(mvabund)

Roots <- read.csv("ASV_deseq2_roots_nmds_nooutliers.csv")
Roots.matrix <- Roots[4:5665 ] # just the community data
Roots.factors <- Roots[1:3] # just the factors
Roots.factors$Root <- as.factor(Roots.factors$Root)
Roots.factors$Salinity <- as.factor(Roots.factors$Salinity)

Roots.mva=mvabund(Roots.matrix) # makes an mvabund object with the bacterial data (the transposed matrix of RAW COUNTS)

# Plot the data - look at the differences between ASVs between samples 
plot(Roots.mva ~ Roots.factors$Root, cex.axis=0.55, type="bx", legend=F) # Box plot

# Lets fit the model
nb.1 = manyglm(Roots.mva ~ Root*Salinity, data=Roots.factors)

# Check the model fitting plotting the residuals
plot(nb.1,n.vars=15) # looks like shit but oh well the model still runs 

# Let run the statistical test
#dat.aov = anova(nb.1, nBoot=1000, p.uni="adjusted") # used nBoot = 100 for trial,# 1000 for publication
#saveRDS(dat.aov, file= "C:/Users/z3532738/OneDrive - UNSW/UNSW ra/Student Data/Lizzy/Bioinformatics/Post_processing_outputs/R graphs/roots.dat.aov.rds")
dat.aov <- readRDS("roots.dat.aov.rds")

# Get the output of the test
dat.aov$table # multivariate table
#marginal interaction, effect of salinity, marginal effect of root. but we can stick wth the primer results

# We can look at test statistics and p values here. 
dat.aov$uni.test # univariate test statistics

dat.aov$uni.p # univariate adjusted p values, note we are interested in row 2 and 3

# We can use code to as: Is any of these significantly different? (i.e. p<0.05)
#we are interested in row 2 (roots) and row 3 (salinity)
any(dat.aov$uni.p[3,]<0.05) 
# Yes, some bacterial classes differed between species

# We can visualise this by plotting those p-values from smaller to larger
plot(sort(dat.aov$uni.p[2,]), main="Sorted univariate P values")
abline(h=0.05, lty=2, col="red") # line at the 0.05 cut-off

# Around 5-6 classes around p = 0.05. Who are these?

diffs=which(dat.aov$uni.p[2,]<0.05) # you could also ask for < 0.1 so we catch the group slightly above 0.05
diffs

stat_dat_aov <- t(dat.aov$uni.test) # univariate test statistics
p_dat_aov <- t(dat.aov$uni.p) # univariate test statistics; p-values

# Convert the matrix to a dataframe
stat_dat_aov_df <- as.data.frame(stat_dat_aov)
p_dat_aov_df <- as.data.frame(p_dat_aov)

# Save the dataframe as a CSV file
#write.csv(stat_dat_aov_df, "C:/Users/z3532738/OneDrive - UNSW/UNSW ra/Student Data/Lizzy/Bioinformatics/Post_processing_outputs/R graphs/p_roots.pls.csv")
write.csv(p_dat_aov_df, "C:/Users/z3532738/OneDrive - UNSW/UNSW ra/Student Data/Lizzy/Bioinformatics/Post_processing_outputs/R graphs/names.csv")
#weirdly CSV isnt working so using XL
write_xlsx(p_dat_aov_df, col_names=TRUE, "C:/Users/z3532738/OneDrive - UNSW/UNSW ra/Student Data/Lizzy/Bioinformatics/Post_processing_outputs/R graphs/p_roots.xlsx")
write_xlsx(stat_dat_aov_df, col_names=TRUE, "C:/Users/z3532738/OneDrive - UNSW/UNSW ra/Student Data/Lizzy/Bioinformatics/Post_processing_outputs/R graphs/stat_roots.xlsx")
#combine these two excel files into one, and then sort by "p value". 
#copy paste the ASVs for significant values (p<0.05) into a table and put that table in your thesis 

# We can then plot the abundances of only these groups.ugly graph below 
palette('default') # mvabund ruins your colour palette, so reset
plot(Roots.mva ~ Roots.factors$Root, var.subset=diffs, type="bx", cex.axis=0.55, legend=F)

#make nice ggplot with edited data 
#import data with ASVS as rows 
#make sure this data is RELATIVE ABUNDANCE asv table 
roots_manyglm <- read.csv("roots_manyglm.csv", header = TRUE, stringsAsFactors = FALSE)

#Change asv table from wide to long format with column names Sequence, ASV, SampleID and Abundance (in that order)
asv_table_ra_long <- roots_manyglm %>% 
  pivot_longer(cols = -1, # Pivot all columns except the first one (ID column = ASV)
               names_to = "Sample", # New column for the names of the original columns
               values_to = "Relabun")
write.csv(asv_table_ra_long,"C:/Users/z3532738/OneDrive - UNSW/UNSW ra/Student Data/Lizzy/Bioinformatics/Post_processing_outputs/R graphs/roots_long.csv")
#in excel add in the meta data (unhealthy vs healthy) and upload again below 

roots_manyglm <- read.csv("roots_long_edited.csv")
mean_relabun <- roots_manyglm %>%
  group_by(SampleID, Salinity) %>%
  summarise(Mean_Relabun = mean(Relabun, na.rm = TRUE),SEM = sd(Relabun, na.rm = TRUE) / sqrt(n()))
write.csv(mean_relabun,"C:/Users/z3532738/OneDrive - UNSW/UNSW ra/Student Data/Lizzy/Bioinformatics/Post_processing_outputs/R graphs/roots_abundance_salinity.csv")
mean_relabun <- roots_manyglm %>%
  group_by(SampleID, Root) %>%
  summarise(Mean_Relabun = mean(Relabun, na.rm = TRUE),SEM = sd(Relabun, na.rm = TRUE) / sqrt(n()))
write.csv(mean_relabun,"C:/Users/z3532738/OneDrive - UNSW/UNSW ra/Student Data/Lizzy/Bioinformatics/Post_processing_outputs/R graphs/roots_abundance_root.csv")

#in excel search for the significant ASVs and add them to their own excel file, and upload below as a plot file 

#take the ASVS you want to plot in excel and turn into csv file
all.abun.plot <- read.csv("all.abundance.plot.root.csv")
all.abun.plot$Root <- factor(all.abun.plot$Root, levels = c("Disrupted", "Intact"), labels = c("Disrupted", "Intact")) 

all.abun.plot.S <- read.csv("all.abundance.plot.salinity.csv")
all.abun.plot.S$Salinity <- factor(all.abun.plot.S$Salinity) 


# Setting colors
mycols2 <- c("Intact" = "#6B8E23",  # Soft mauve
             "Disrupted" = "#BDA0CB")  # Olive green
mycols3 <- c("Ambient" = "#000080", "Low" = "#a8c9e3", "Medium" = "#4292c6")

#roots
root2 <- ggplot(all.abun.plot, aes(x = SampleID, y = Mean_Relabun, fill = Root)) +
  geom_bar(stat = "identity", width = 0.75, position = position_dodge(width = -0.8), alpha = 0.8) +
  xlab("") + ylab("") +
  scale_fill_manual(values = mycols2) +
  geom_errorbar(aes(ymin = Mean_Relabun - SEM, ymax = Mean_Relabun + SEM), width = 0.4,
                position = position_dodge(width = -0.8)) +
  labs(title = "A)", fill = "Root treatment") +  # Change legend title here
  theme(
    axis.text.y = element_text(size = 12), 
        axis.title.y = element_text(size = 12, face = 'bold'),
        axis.text.x = element_text(size = 14), 
        axis.title.x = element_text(size = 12, face = 'bold'),
        legend.position = "none", 
        legend.direction = "horizontal",
        legend.title =  element_text(size = 12),  # Style the legend title, 
        legend.text = element_text(size = 11),
        panel.background = element_rect(fill = "white", colour = "black"),
        strip.text = element_text(size = 15, face = "bold")) +
  coord_flip()  # Flip coordinates for horizontal bars
root2

ggarrange(root2/root1)

#salinity
###input this graph into a nice format with ggplot2. 
sal1<- ggplot(all.abun.plot.S, aes(x = SampleID, y = Mean_Relabun, fill = Salinity)) +
  geom_bar(stat = "identity", width = 0.9, position = position_dodge(width = -0.8), alpha = 0.8) +
  xlab("") + ylab("") +
  scale_fill_manual(values = mycols3) +
  geom_errorbar(aes(ymin = Mean_Relabun - SEM, ymax = Mean_Relabun + SEM), width = 0.4,
                position = position_dodge(width = -0.8)) +
  labs(title = "A)",fill = "Salinity treatment") +  # Change legend title here
  theme(axis.text.y = element_text(size = 12), 
        axis.title.y = element_text(size = 12, face = 'bold'),
        legend.position = "none", 
        legend.direction = "horizontal",
        legend.title =  element_text(size = 12),  # Style the legend title, 
        legend.text = element_text(size = 11),
        panel.background = element_rect(fill = "white", colour = "black"),
        strip.text = element_text(size = 15, face = "bold")) +
  coord_flip()  # Flip coordinates for horizontal bars
sal1
