metata_table <- read.csv('../../HLI/Metadata/HLI-pull.csv', stringsAsFactors = F)

load("../../HLI/Metadata/HLI_available_pairs_dis_table.RData")


avail_table <- merge(metata_table, Available_paired_table, by.x = "pair_id", by.y = "BMT")

Age_table <- avail_table[, c(1:4, 12, 22, 27, 28, 31,32)]
library(eeptools)

Age_table$RAge <- sapply(1:dim(Age_table)[1], function(x) age_calc(as.Date(Age_table$Rbdte[x]), as.Date(Age_table$transplant_dte[x]), units = "years"))

Age_table$DAge <- sapply(1:dim(Age_table)[1], function(x) age_calc(as.Date(Age_table$Dbdte[x]), as.Date(Age_table$transplant_dte[x]), units = "years"))

Age_table$AgeDiff <- abs(Age_table$RAge - Age_table$DAge)

Age_table$RtoD_AgeDiff <- Age_table$RAge - Age_table$DAge

Age_table$RAge_range <- sapply(1:dim(Age_table)[1], function(x) if (Age_table$RAge[x] <= 9) "0~9" else if (Age_table$RAge[x] <= 19) "18~19" else if(Age_table$RAge[x] <= 29) "20~29" 
                               else if(Age_table$RAge[x] <= 39) "30~39" else if(Age_table$RAge[x] <= 49) "40~49" else if (Age_table$RAge[x] <= 59) "50~59" else "60+")

Age_table$DAge_range <- sapply(1:dim(Age_table)[1], function(x) if (Age_table$DAge[x] <= 19) "18~19" else if(Age_table$DAge[x] <= 29) "20~29" 
                               else if(Age_table$DAge[x] <= 39) "30~39" else if(Age_table$DAge[x] <= 49) "40~49" else if (Age_table$DAge[x] <= 59) "50~59" else "60+")


# write.csv(Age_table, file = "../../HLI/Metadata/Available_pair_ageTable.csv", row.names = F)
## 
Age_table <- read.csv("../../HLI/Metadata/Available_pair_ageTable.csv", stringsAsFactors = F)
## MiHA table
MiHA_count_table <- read.csv("../../HLI/Metadata/HLI_known_MiHA_table_032018.csv", stringsAsFactors = F)

#### combined table

MiHA_Age_table <- merge(x = MiHA_count_table, y = Age_table, by.x = "groupID", by.y = "GroupID")

MiHA_Age_table$RAge_range <- as.factor(MiHA_Age_table$RAge_range)
MiHA_Age_table$DAge_range <- as.factor(MiHA_Age_table$DAge_range)

MiHA_Age_table$AgeDiff_range <- sapply(1:dim(Age_table)[1], function(x) if (Age_table$AgeDiff[x] <= 9) "0~9" 
                                       else if(Age_table$DAge[x] <= 19) "10~19" 
                                       else if(Age_table$DAge[x] <= 29) "20~29" 
                                       else if(Age_table$DAge[x] <= 39) "30~39" 
                                       else if(Age_table$DAge[x] <= 49) "40~49"  else "50+")
MiHA_Age_table$AgeDiff_range <- as.factor(MiHA_Age_table$AgeDiff_range)
MiHA_Age_table$GVHD <- as.factor(MiHA_Age_table$GVHD)

## plot
library(ggplot2)
library(beeswarm)
library(plyr)

# color by RecipientAgeRange color - 7
# Recipient_colors = c("#ffffcc", "#c7e9b4", "#7fcdbb", "#41b6c4", "#1d91c0", "#225ea8", "#0c2c84") # Blue - 7
# Recipient_colors = c("#f0f9e8", "#ccebc5", "#ccebc5", "#7bccc4", "#4eb3d3", "#2b8cbe", "#08589e") # LightBlue - 7
# Recipient_colors = c("#edf8e9", "#c7e9c0", "#a1d99b", "#74c476", "#41ab5d", "#238b45", "#005a32") # Green - 7
# Recipient_colors = c("#ffffb2", "#fed976", "#feb24c", "#fd8d3c", "#fc4e2a", "#e31a1c", "#b10026") # Red - 7
# Recipient_colors = c("#f2f0f7", "#dadaeb", "#bcbddc", "#9e9ac8", "#807dba", "#6a51a3", "#4a1486") # Purple - 7
Recipient_colors = c("#d7191c", "#fdae61", "#5e3c99", "#b2182b", "#b35806", "#1b7837", "#2c7bb6")

# color by DonorAgeRange color - 5
Donor_colors = c("#d7191c", "#008837", "#0571b0", "#5e3c99", "#e66101") # Diverging- 5

# color by AgeDiffRagne - 6
AgeDiff_colors = c("#ca0020", "#008837", "#7b3294", "#d01c8b", "#0571b0", "#fdae61") # - 6
  
##### Restricted MiHA


# by Recipient Age range
beeswarm_ResMiHA_all_R <- beeswarm(HLA.Restricted ~ GVHD, data = MiHA_Age_table, 
                                 method = 'swarm', pwcol = RAge_range,
                                 spacing = 0.3)[, c(1,2, 4, 6)]
colnames(beeswarm_ResMiHA_all_R) <- c("x", "y", "RecipientAgeRange", "GVHD") 

ggplot(beeswarm_ResMiHA_all_R, aes(x, y)) +
  xlab("") +
  scale_y_continuous(expression("#Mismatched HLA Restricted Known MiHA SNPs")) + 
  #geom_point(aes(colour = RecipientAgeRange), size = 2, alpha = 0.5) +
  geom_jitter(width = 0.1, aes(color = RecipientAgeRange), size = 4, alpha = 0.6) +
  scale_colour_manual(values = Recipient_colors) + 
  scale_x_continuous(breaks = c(1:2), 
                     labels = c("acute GVHD", "non-GVHD"), expand = c(0, 0.05)) +
  geom_boxplot(aes(x, y, group = round_any(beeswarm_ResMiHA_all_R$x, 1, round)), outlier.shape = NA, alpha = 0) # +
  #theme(legend.position = "none")
Rage_table <- table(MiHA_Age_table$RAge_range, MiHA_Age_table$GVHD)
chisq.test(Rage_table) 
# Pearson's Chi-squared test
# 
# data:  Rage_table
# X-squared = 8.0333, df = 6, p-value = 0.2357  

Restricted_Rage_table <- table(MiHA_Age_table$RAge_range, MiHA_Age_table$HLA.Restricted)
chisq.test(Restricted_Rage_table) 

# Pearson's Chi-squared test
# 
# data:  Restricted_Rage_table
# X-squared = 56.543, df = 66, p-value = 0.7903

# by donor age range
beeswarm_ResMiHA_all_D <- beeswarm(HLA.Restricted ~ GVHD, data = MiHA_Age_table, 
                                 method = 'swarm', pwcol = DAge_range,
                                 spacing = 0.3)[, c(1,2, 4, 6)]
colnames(beeswarm_ResMiHA_all_D) <- c("x", "y", "DonorAgeRange", "GVHD") 

ggplot(beeswarm_ResMiHA_all_D, aes(x, y)) +
  xlab("") +
  scale_y_continuous(expression("#Mismatched HLA Restricted Known MiHA SNPs")) + 
  #geom_point(aes(colour = RecipientAgeRange), size = 2, alpha = 0.5) +
  geom_jitter(width = 0.1, aes(color = DonorAgeRange), size = 4, alpha = 0.7) +
  scale_colour_manual(values = AgeDiff_colors) + 
  scale_x_continuous(breaks = c(1:2), 
                     labels = c("acute GVHD", "non-GVHD"), expand = c(0, 0.05)) +
  geom_boxplot(aes(x, y, group = round_any(beeswarm_ResMiHA_all_D$x, 1, round)), outlier.shape = NA, alpha = 0) # +
#theme(legend.position = "none")
Dage_table <- table(MiHA_Age_table$DAge_range, MiHA_Age_table$GVHD)
chisq.test(Dage_table) 
# Pearson's Chi-squared test
# 
# data:  Dage_table
# X-squared = 6.8459, df = 4, p-value = 0.1443

Restricted_Dage_table <- table(MiHA_Age_table$DAge_range, MiHA_Age_table$HLA.Restricted)
chisq.test(Restricted_Dage_table) 

# Pearson's Chi-squared test
# 
# data:  Restricted_Dage_table
# X-squared = 45.735, df = 44, p-value = 0.3999

# by Age diff range
beeswarm_ResMiHA_all_AD <- beeswarm(HLA.Restricted ~ GVHD, data = MiHA_Age_table, 
                                   method = 'swarm', pwcol = AgeDiff_range,
                                   spacing = 0.3)[, c(1,2, 4, 6)]
colnames(beeswarm_ResMiHA_all_AD) <- c("x", "y", "AgeDifference", "GVHD") 

ggplot(beeswarm_ResMiHA_all_AD, aes(x, y)) +
  xlab("") +
  scale_y_continuous(expression("#Mismatched HLA Restricted Known MiHA SNPs")) + 
  #geom_point(aes(colour = RecipientAgeRange), size = 2, alpha = 0.5) +
  geom_jitter(width = 0.1, aes(color = AgeDifference), size = 4, alpha = 0.6) +
  scale_colour_manual(values = AgeDiff_colors) + 
  scale_x_continuous(breaks = c(1:2), 
                     labels = c("acute GVHD", "non-GVHD"), expand = c(0, 0.05)) +
  geom_boxplot(aes(x, y, group = round_any(beeswarm_ResMiHA_all_AD$x, 1, round)), outlier.shape = NA, alpha = 0) # +
#theme(legend.position = "none")

ADiff_table <- table(MiHA_Age_table$AgeDiff_range, MiHA_Age_table$GVHD)
chisq.test(ADiff_table) 

# Pearson's Chi-squared test
# 
# data:  ADiff_table
# X-squared = 5.6353, df = 5, p-value = 0.3433


Restricted_ADiff_table <- table(MiHA_Age_table$AgeDiff_range, MiHA_Age_table$HLA.Restricted)
chisq.test(Restricted_ADiff_table) 

# Pearson's Chi-squared test
# 
# data:  Restricted_ADiff_table
# X-squared = 65.55, df = 55, p-value = 0.1561

###### Unrestricted MiHA
# by recipient age range
beeswarm_UnResMiHA_all_R <- beeswarm(HLA.Unrestricted ~ GVHD, data = MiHA_Age_table,
                                   method = 'swarm', pwcol = RAge_range,
                                   spacing = 0.5)[, c(1,2,4, 6)]
colnames(beeswarm_UnResMiHA_all_R) <- c("x", "y", "RecipientAgeRange", "GVHD") 

ggplot(beeswarm_UnResMiHA_all_R, aes(x, y)) +
  xlab("") +
  scale_y_continuous(expression("#Mismatched HLA Unrestricted Known MiHA SNPs")) + 
  #geom_point(aes(colour = GVHD), size = 3, alpha = 0.7) +
  geom_jitter(width = 0.1, aes(color = RecipientAgeRange), size = 4, alpha = 0.7) +
  scale_colour_manual(values = Recipient_colors) + 
  scale_x_continuous(breaks = c(1:2), 
                     labels = c("acute GVHD", "non-GVHD"), expand = c(0, 0.05)) +
  geom_boxplot(aes(x, y, group = round_any(beeswarm_UnResMiHA_all_R$x, 1, round)), outlier.shape = NA, alpha = 0) 

Unrestricted_RAge_table <- table(MiHA_Age_table$RAge_range, MiHA_Age_table$HLA.Unrestricted)
chisq.test(Unrestricted_RAge_table)
# Pearson's Chi-squared test
# 
# data:  Unrestricted_RAge_table
# X-squared = 102.37, df = 120, p-value = 0.8761

## by donor age range
beeswarm_UnResMiHA_all_D <- beeswarm(HLA.Unrestricted ~ GVHD, data = MiHA_Age_table,
                                     method = 'swarm', pwcol = DAge_range,
                                     spacing = 0.5)[, c(1,2,4, 6)]
colnames(beeswarm_UnResMiHA_all_D) <- c("x", "y", "DonorAgeRange", "GVHD") 

ggplot(beeswarm_UnResMiHA_all_D, aes(x, y)) +
  xlab("") +
  scale_y_continuous(expression("#Mismatched HLA Unrestricted Known MiHA SNPs")) + 
  #geom_point(aes(colour = GVHD), size = 3, alpha = 0.7) +
  geom_jitter(width = 0.1, aes(color = DonorAgeRange), size = 4, alpha = 0.7) +
  scale_colour_manual(values = Donor_colors) + 
  scale_x_continuous(breaks = c(1:2), 
                     labels = c("acute GVHD", "non-GVHD"), expand = c(0, 0.05)) +
  geom_boxplot(aes(x, y, group = round_any(beeswarm_UnResMiHA_all_D$x, 1, round)), outlier.shape = NA, alpha = 0) 

Unrestricted_DAge_table <- table(MiHA_Age_table$DAge_range, MiHA_Age_table$HLA.Unrestricted)
chisq.test(Unrestricted_DAge_table)
# Pearson's Chi-squared test
# 
# data:  Unrestricted_DAge_table
# X-squared = 84.931, df = 80, p-value = 0.3319

## by age diff range
beeswarm_UnResMiHA_all_AD <- beeswarm(HLA.Unrestricted ~ GVHD, data = MiHA_Age_table,
                                     method = 'swarm', pwcol = AgeDiff_range,
                                     spacing = 0.5)[, c(1,2,4, 6)]
colnames(beeswarm_UnResMiHA_all_AD) <- c("x", "y", "AgeDifference", "GVHD") 

ggplot(beeswarm_UnResMiHA_all_AD, aes(x, y)) +
  xlab("") +
  scale_y_continuous(expression("#Mismatched HLA Unrestricted Known MiHA SNPs")) + 
  #geom_point(aes(colour = GVHD), size = 3, alpha = 0.7) +
  geom_jitter(width = 0.1, aes(color = AgeDifference), size = 4, alpha = 0.7) +
  scale_colour_manual(values = AgeDiff_colors) + 
  scale_x_continuous(breaks = c(1:2), 
                     labels = c("acute GVHD", "non-GVHD"), expand = c(0, 0.05)) +
  geom_boxplot(aes(x, y, group = round_any(beeswarm_UnResMiHA_all_AD$x, 1, round)), outlier.shape = NA, alpha = 0) 

Unrestricted_AgeDiff_table <- table(MiHA_Age_table$AgeDiff_range, MiHA_Age_table$HLA.Unrestricted)
chisq.test(Unrestricted_AgeDiff_table)
# Pearson's Chi-squared test
# 
# data:  Unrestricted_AgeDiff_table
# X-squared = 111.92, df = 100, p-value = 0.1954

##### Number of Missense Mismatch

# by recipient age range
beeswarm_MissenseMis_R <-  beeswarm(NumVar ~ GVHD, data = MiHA_Age_table,
                                  method = 'swarm',pwcol = RAge_range,
                                  spacing = 0.4)[, c(1,2, 4, 6)]
colnames(beeswarm_MissenseMis_R) <- c("x", "y", "RecipientAgeRange", "GVHD") 

ggplot(beeswarm_MissenseMis_R, aes(x, y)) +
  xlab("") +
  scale_y_continuous(expression("#Missense Mismathced SNPs")) + 
  #geom_point(aes(colour = GVHD), size = 3, alpha = 0.6) +
  geom_jitter(width = 0.1, aes(color = RecipientAgeRange), size = 3, alpha = 0.7) +
  scale_colour_manual(values = Recipient_colors) + 
  scale_x_continuous(breaks = c(1:2), 
                     labels = c("acute GVHD", "non-GVHD"), expand = c(0, 0.05)) +
  geom_boxplot(aes(x, y, group = round_any(beeswarm_MissenseMis_R$x, 1, round)), outlier.shape = NA, alpha = 0)

Miss_RAge_table <- table(MiHA_Age_table$RAge_range, MiHA_Age_table$NumVar)
chisq.test(Miss_RAge_table)
# Pearson's Chi-squared test
# 
# data:  Miss_RAge_table
# X-squared = 1085.6, df = 1032, p-value = 0.1203

# by donor age range
beeswarm_MissenseMis_D <-  beeswarm(NumVar ~ GVHD, data = MiHA_Age_table,
                                    method = 'swarm',pwcol = DAge_range,
                                    spacing = 0.4)[, c(1,2, 4, 6)]
colnames(beeswarm_MissenseMis_D) <- c("x", "y", "DonorAgeRange", "GVHD") 

ggplot(beeswarm_MissenseMis_D, aes(x, y)) +
  xlab("") +
  scale_y_continuous(expression("#Missense Mismathced SNPs")) + 
  #geom_point(aes(colour = GVHD), size = 3, alpha = 0.6) +
  geom_jitter(width = 0.1, aes(color = DonorAgeRange), size = 3, alpha = 0.7) +
  scale_colour_manual(values = Donor_colors) + 
  scale_x_continuous(breaks = c(1:2), 
                     labels = c("acute GVHD", "non-GVHD"), expand = c(0, 0.05)) +
  geom_boxplot(aes(x, y, group = round_any(beeswarm_MissenseMis_D$x, 1, round)), outlier.shape = NA, alpha = 0)

Miss_DAge_table <- table(MiHA_Age_table$DAge_range, MiHA_Age_table$NumVar)
chisq.test(Miss_DAge_table)
# Pearson's Chi-squared test
# 
# data:  Miss_DAge_table
# X-squared = 722.76, df = 688, p-value = 0.1737

# by Age Difference range
beeswarm_MissenseMis_AD <-  beeswarm(NumVar ~ GVHD, data = MiHA_Age_table,
                                    method = 'swarm',pwcol = AgeDiff_range,
                                    spacing = 0.4)[, c(1,2, 4, 6)]
colnames(beeswarm_MissenseMis_AD) <- c("x", "y", "AgeDifferenceRange", "GVHD") 

ggplot(beeswarm_MissenseMis_AD, aes(x, y)) +
  xlab("") +
  scale_y_continuous(expression("#Missense Mismathced SNPs")) + 
  #geom_point(aes(colour = GVHD), size = 3, alpha = 0.6) +
  geom_jitter(width = 0.1, aes(color = AgeDifferenceRange), size = 3, alpha = 0.7) +
  scale_colour_manual(values = AgeDiff_colors) + 
  scale_x_continuous(breaks = c(1:2), 
                     labels = c("acute GVHD", "non-GVHD"), expand = c(0, 0.05)) +
  geom_boxplot(aes(x, y, group = round_any(beeswarm_MissenseMis_AD$x, 1, round)), outlier.shape = NA, alpha = 0)

Miss_Age_table <- table(MiHA_Age_table$AgeDiff_range, MiHA_Age_table$NumVar)
chisq.test(Miss_Age_table)
# #	Pearson's Chi-squared test
# 
# data:  Miss_Age_table
# X-squared = 879.69, df = 860, p-value = 0.313


###########
# Restricted MiHAs By AGE
GVHD_color <- c("#ca0020", "#0571b0")#, "#008837")
## by recipient age range
beeswarm_ResMiHA_all_R <- beeswarm(HLA.Restricted ~ RAge_range, data = MiHA_Age_table, 
                                   method = 'swarm', pwcol = as.factor(GVHD),
                                   spacing = 0.3)[, c(1,2, 4, 6)]
colnames(beeswarm_ResMiHA_all_R) <- c("x", "y", "GVHD", "RAge_range") 

ggplot(beeswarm_ResMiHA_all_R, aes(x, y)) +
  xlab("Recipient age range") +
  scale_y_continuous(expression("#Mismatched HLA Restricted Known MiHA SNPs")) + 
  #geom_point(aes(colour = RecipientAgeRange), size = 2, alpha = 0.5) +
  geom_jitter(width = 0.05, aes(color = GVHD), size = 4, alpha = 0.6) +
  scale_colour_manual(values = GVHD_color) + 
  scale_x_continuous(breaks = c(1:length(levels(MiHA_Age_table$RAge_range))), 
                     labels = levels(MiHA_Age_table$RAge_range), expand = c(0, 0.05)) +
  geom_boxplot(aes(x, y, group = round_any(beeswarm_ResMiHA_all_R$x, 1, round)), outlier.shape = NA, alpha = 0) # +
#theme(legend.position = "none")

## by donor age range
beeswarm_ResMiHA_all_D <- beeswarm(HLA.Restricted ~ DAge_range, data = MiHA_Age_table, 
                                   method = 'swarm', pwcol = as.factor(GVHD),
                                   spacing = 0.3)[, c(1,2, 4, 6)]
colnames(beeswarm_ResMiHA_all_D) <- c("x", "y", "GVHD", "DAge_range") 

ggplot(beeswarm_ResMiHA_all_D, aes(x, y)) +
  xlab("Donor age range") +
  scale_y_continuous(expression("#Mismatched HLA Restricted Known MiHA SNPs")) + 
  #geom_point(aes(colour = RecipientAgeRange), size = 2, alpha = 0.5) +
  geom_jitter(width = 0.05, aes(color = GVHD), size = 4, alpha = 0.6) +
  scale_colour_manual(values = GVHD_color) + 
  scale_x_continuous(breaks = c(1:length(levels(MiHA_Age_table$DAge_range))), 
                     labels = levels(MiHA_Age_table$DAge_range), expand = c(0, 0.05)) +
  geom_boxplot(aes(x, y, group = round_any(beeswarm_ResMiHA_all_D$x, 1, round)), outlier.shape = NA, alpha = 0) # +
#theme(legend.position = "none")

## by Age difference
beeswarm_ResMiHA_all_AD <- beeswarm(HLA.Restricted ~ AgeDiff_range, data = MiHA_Age_table, 
                                   method = 'swarm', pwcol = as.factor(GVHD),
                                   spacing = 0.3)[, c(1,2, 4, 6)]
colnames(beeswarm_ResMiHA_all_AD) <- c("x", "y", "GVHD", "AgeDiff_range") 

ggplot(beeswarm_ResMiHA_all_AD, aes(x, y)) +
  xlab("Age difference") +
  scale_y_continuous(expression("#Mismatched HLA Restricted Known MiHA SNPs")) + 
  #geom_point(aes(colour = RecipientAgeRange), size = 2, alpha = 0.5) +
  geom_jitter(width = 0.05, aes(color = GVHD), size = 4, alpha = 0.6) +
  scale_colour_manual(values = GVHD_color) + 
  scale_x_continuous(breaks = c(1:length(levels(MiHA_Age_table$AgeDiff_range))), 
                     labels = levels(MiHA_Age_table$AgeDiff_range), expand = c(0, 0.05)) +
  geom_boxplot(aes(x, y, group = round_any(beeswarm_ResMiHA_all_AD$x, 1, round)), outlier.shape = NA, alpha = 0) # +
#theme(legend.position = "none")

####
# Unrestricted MiHAs By AGE
GVHD_color <- c("#ca0020", "#0571b0")#, "#008837")
## by recipient age range
beeswarm_UnResMiHA_all_R <- beeswarm(HLA.Unrestricted ~ RAge_range, data = MiHA_Age_table, 
                                   method = 'swarm', pwcol = as.factor(GVHD),
                                   spacing = 0.3)[, c(1,2, 4, 6)]
colnames(beeswarm_UnResMiHA_all_R) <- c("x", "y", "GVHD", "RAge_range") 

ggplot(beeswarm_UnResMiHA_all_R, aes(x, y)) +
  xlab("Recipient age range") +
  scale_y_continuous(expression("#Mismatched HLA Unrestricted Known MiHA SNPs")) + 
  #geom_point(aes(colour = RecipientAgeRange), size = 2, alpha = 0.5) +
  geom_jitter(width = 0.05, aes(color = GVHD), size = 4, alpha = 0.6) +
  scale_colour_manual(values = GVHD_color) + 
  scale_x_continuous(breaks = c(1:length(levels(MiHA_Age_table$RAge_range))), 
                     labels = levels(MiHA_Age_table$RAge_range), expand = c(0, 0.05)) +
  geom_boxplot(aes(x, y, group = round_any(beeswarm_UnResMiHA_all_R$x, 1, round)), outlier.shape = NA, alpha = 0) # +
#theme(legend.position = "none")

## by donor age range
beeswarm_UnResMiHA_all_D <- beeswarm(HLA.Unrestricted ~ DAge_range, data = MiHA_Age_table, 
                                   method = 'swarm', pwcol = as.factor(GVHD),
                                   spacing = 0.3)[, c(1,2, 4, 6)]
colnames(beeswarm_UnResMiHA_all_D) <- c("x", "y", "GVHD", "DAge_range") 

ggplot(beeswarm_UnResMiHA_all_D, aes(x, y)) +
  xlab("Donor age range") +
  scale_y_continuous(expression("#Mismatched HLA Unrestricted Known MiHA SNPs")) + 
  #geom_point(aes(colour = RecipientAgeRange), size = 2, alpha = 0.5) +
  geom_jitter(width = 0.05, aes(color = GVHD), size = 4, alpha = 0.6) +
  scale_colour_manual(values = GVHD_color) + 
  scale_x_continuous(breaks = c(1:length(levels(MiHA_Age_table$DAge_range))), 
                     labels = levels(MiHA_Age_table$DAge_range), expand = c(0, 0.05)) +
  geom_boxplot(aes(x, y, group = round_any(beeswarm_UnResMiHA_all_D$x, 1, round)), outlier.shape = NA, alpha = 0) # +
#theme(legend.position = "none")

## by Age difference
beeswarm_UnResMiHA_all_AD <- beeswarm(HLA.Unrestricted ~ AgeDiff_range, data = MiHA_Age_table, 
                                    method = 'swarm', pwcol = as.factor(GVHD),
                                    spacing = 0.3)[, c(1,2, 4, 6)]
colnames(beeswarm_UnResMiHA_all_AD) <- c("x", "y", "GVHD", "AgeDiff_range") 

ggplot(beeswarm_UnResMiHA_all_AD, aes(x, y)) +
  xlab("Age difference") +
  scale_y_continuous(expression("#Mismatched HLA Unrestricted Known MiHA SNPs")) + 
  #geom_point(aes(colour = RecipientAgeRange), size = 2, alpha = 0.5) +
  geom_jitter(width = 0.05, aes(color = GVHD), size = 4, alpha = 0.6) +
  scale_colour_manual(values = GVHD_color) + 
  scale_x_continuous(breaks = c(1:length(levels(MiHA_Age_table$AgeDiff_range))), 
                     labels = levels(MiHA_Age_table$AgeDiff_range), expand = c(0, 0.05)) +
  geom_boxplot(aes(x, y, group = round_any(beeswarm_UnResMiHA_all_AD$x, 1, round)), outlier.shape = NA, alpha = 0) # +
#theme(legend.position = "none")

########
####
# Missense mismatches By AGE
GVHD_color <- c("#ca0020", "#0571b0")#, "#008837")
## by recipient age range
beeswarm_MisMismatch_all_R <- beeswarm(NumVar ~ RAge_range, data = MiHA_Age_table, 
                                     method = 'swarm', pwcol = as.factor(GVHD),
                                     spacing = 0.3)[, c(1,2, 4, 6)]
colnames(beeswarm_MisMismatch_all_R) <- c("x", "y", "GVHD", "RAge_range") 

ggplot(beeswarm_MisMismatch_all_R, aes(x, y)) +
  xlab("Recipient age range") +
  scale_y_continuous(expression("#Missense Mismathced SNPs")) + 
  #geom_point(aes(colour = RecipientAgeRange), size = 2, alpha = 0.5) +
  geom_jitter(width = 0.05, aes(color = GVHD), size = 3, alpha = 0.6) +
  scale_colour_manual(values = GVHD_color) + 
  scale_x_continuous(breaks = c(1:length(levels(MiHA_Age_table$RAge_range))), 
                     labels = levels(MiHA_Age_table$RAge_range), expand = c(0, 0.05)) +
  geom_boxplot(aes(x, y, group = round_any(beeswarm_MisMismatch_all_R$x, 1, round)), outlier.shape = NA, alpha = 0) # +
#theme(legend.position = "none")

## by donor age range
beeswarm_MisMismatch_all_D <- beeswarm(NumVar ~ DAge_range, data = MiHA_Age_table, 
                                     method = 'swarm', pwcol = as.factor(GVHD),
                                     spacing = 0.3)[, c(1,2, 4, 6)]
colnames(beeswarm_MisMismatch_all_D) <- c("x", "y", "GVHD", "DAge_range") 

ggplot(beeswarm_MisMismatch_all_D, aes(x, y)) +
  xlab("Donor age range") +
  scale_y_continuous(expression("#Missense Mismathced SNPs")) + 
  #geom_point(aes(colour = RecipientAgeRange), size = 2, alpha = 0.5) +
  geom_jitter(width = 0.05, aes(color = GVHD), size = 4, alpha = 0.6) +
  scale_colour_manual(values = GVHD_color) + 
  scale_x_continuous(breaks = c(1:length(levels(MiHA_Age_table$DAge_range))), 
                     labels = levels(MiHA_Age_table$DAge_range), expand = c(0, 0.05)) +
  geom_boxplot(aes(x, y, group = round_any(beeswarm_MisMismatch_all_D$x, 1, round)), outlier.shape = NA, alpha = 0) # +
#theme(legend.position = "none")

## by Age difference
beeswarm_MisMismatch_all_AD <- beeswarm(NumVar ~ AgeDiff_range, data = MiHA_Age_table, 
                                      method = 'swarm', pwcol = as.factor(GVHD),
                                      spacing = 0.3)[, c(1,2, 4, 6)]
colnames(beeswarm_MisMismatch_all_AD) <- c("x", "y", "GVHD", "AgeDiff_range") 

ggplot(beeswarm_MisMismatch_all_AD, aes(x, y)) +
  xlab("Age difference") +
  scale_y_continuous(expression("#Missense Mismathced SNPs")) + 
  #geom_point(aes(colour = RecipientAgeRange), size = 2, alpha = 0.5) +
  geom_jitter(width = 0.05, aes(color = GVHD), size = 4, alpha = 0.6) +
  scale_colour_manual(values = GVHD_color) + 
  scale_x_continuous(breaks = c(1:length(levels(MiHA_Age_table$AgeDiff_range))), 
                     labels = levels(MiHA_Age_table$AgeDiff_range), expand = c(0, 0.05)) +
  geom_boxplot(aes(x, y, group = round_any(beeswarm_MisMismatch_all_AD$x, 1, round)), outlier.shape = NA, alpha = 0) # +
#theme(legend.position = "none")