setwd("C:Lab/MethylFlash")
# Clean up the environment and remove all the sample object files
rm(list = ls())
# make sure you are using the most recent version of R
install.packages("installr")
library(installr)
updateR()

install.packages("devtools")
library(devtools)

#Here is how you should set it up. 

Methylation_Percent = c(-0.273446569, 0.079291239, 0.119883874, 0.074253442,
                        -0.172644879, -0.076518095, -0.378956542, 
                        -0.016590583, -0.351684992, -0.375687625)


Methylation_Percent
Group = c("Untreated", "Untreated", "Untreated", "Untreated",
          "DMSO 0.02%", "DMSO 0.02%", "DMSO 0.02%",
          "DMSO 0.05%", "DMSO 0.05%", "DMSO 0.05%")


Group  
#Now link them all together
My_data = data.frame(Group, Methylation_Percent);

#Check your work
My_data$Group

library(dplyr)
group_by(My_data, Group) %>%
  summarise(
    count = n(),
    mean = mean(Methylation_Percent, na.rm = TRUE),
    sd = sd(Methylation_Percent, na.rm = TRUE)
  )

install.packages("ggpubr")

library("ggpubr")
ggboxplot(My_data, x = "Group", y = "Methylation_Percent", add ="jitter",
          main = "Global DNA Methylation Percentage
Following DMSO Treatment",
          fill = "Group",
          palette  = "BuPu", 
          order = c("Untreated", "DMSO 0.02%", "DMSO 0.05%"),
          ylab = expression(Delta~5~mC~"% Compared to Untreated Average"), xlab = "Treatment Conditions") +
  theme(plot.title = element_text(hjust = 0.5, size= 18, face="bold"), text = element_text(size = 15), legend.position = "none")

#Compute the one-way ANOVA test
res.aov <- aov(Methylation_Percent ~ Group, data = My_data)

#Output the summary of the ANOVA
summary(res.aov)

#Okay this tells us there is no statistical difference between the groups we can see the pairwise comparison between each group using the post hoc turkey test.
TukeyHSD(res.aov)

# Switching over to running a Kruskal-Wallis test because I don't think my data is normally distributed. Do I really know if my data is normally distributed by looking at it?....No. But there is a test for that too, the Shapiro Wilk test. 
# Lets start with tests for normalcy 

My_data
ggdensity(My_data$Methylation_Percent, 
          main = "Density Plot of Methylation Percent",
          xlab = expression(Delta~5~mC~"% Compared to Untreated Average"),
          ylab = "Density",
          ylim = c(1,2)) 
#Oh that does not look normal. Lets dig even more into this and get a Q-Q plot.

ggqqplot(My_data$Methylation_Percent,
         main = "Q-Q Plot of Methylation Percent")
shapiro.test(My_data$Methylation_Percent)

#Okay to the ggqplot and shapiro test indicate that the data is normal. However, these tests often will pass data sets with small sample sizes so I say that the because it failed the visual density test I belive the data to not be normal and I should move forward with the Kruskal-Wallis test.
My_data$Group <- ordered(My_data$Group,
                         levels = c("Untreated", "DMSO 0.02%", "DMSO 0.05%"))

group_by(My_data, Group) %>%
  summarise(
    count = n(),
    mean = mean(Methylation_Percent),
    sd = sd(Methylation_Percent),
    median = median(Methylation_Percent),
    IQR = IQR(Methylation_Percent)
  )
My_data

kruskal.test(Methylation_Percent ~ Group, data = My_data)
a = pairwise.wilcox.test(My_data$Methylation_Percent, My_data$Group,
                         p.adjust.method = "BH")
a
b = as.matrix(a$p.value)
install.packages("corrplot")
library(corrplot)
corrplot(b, method = "shade", type = 'lower', addCoef.col = 'black')

corrplot(b, p.mat = b , method = 'color', type = 'lower',
         sig.level = c(0.001, 0.01, 0.05), pch.cex = 0.9, 
         insig = 'label_sig', pch.col = 'grey20') 