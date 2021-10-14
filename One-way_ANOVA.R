setwd("C:Lab/MethylFlash")
# Clean up the environment and remove all the sample object files
rm(list = ls())
# make sure you are using the most recent version of R
install.packages("installr")
library(installr)
updateR()
Methylation_Percent = c(-0.273446569, 0.079291239, 0.119883874, 0.074253442,
                        -0.237819952, -0.15191891, -0.06592903, 
                        -0.310569367, -0.188723062, -0.493662425)


Methylation_Percent
Group = c("Untreated", "Untreated", "Untreated", "Untreated",
          "EtOH 0.02%", "EtOH 0.02%", "EtOH 0.02%", 
          "EtOH 0.05%", "EtOH 0.05%", "EtOH 0.05%")


Group  
#Now link them all together
My_data = data.frame(Group, Methylation_Percent);

#Check your work
My_data$Group

#Check the data to determine the shape of the distribution
ggdensity(My_data$Methylation_Percent, 
          main = "Density Plot of Methylation Percent",
          xlab = expression(Delta~5~mC~"% Compared to Untreated Average"),
          ylab = "Density",
          ylim = c(0.5,2)) 
#Cool that looks normal. Lets dig even more into this and get a Q-Q plot and do a shapiro test for data normality.

ggqqplot(My_data$Methylation_Percent,
         main = "Q-Q Plot of Methylation Percent")
shapiro.test(My_data$Methylation_Percent)
#Great, both the Q-Q plot looks normal with no points outside of the theoretical diagonal. The Shapiro test has a value higher than 0.05 so we can not reject the null hypothesis and indicates the data is normal. 
#Calculate the mean and sd using the dplyr package.
library(dplyr)
group_by(My_data, Group) %>%
  summarise(
    count = n(),
    mean = mean(Methylation_Percent, na.rm = TRUE),
    sd = sd(Methylation_Percent, na.rm = TRUE)
  )

#Plot the data as a box plot
install.packages("ggpubr")
library("ggpubr")
ggboxplot(My_data, x = "Group", y = "Methylation_Percent", add ="jitter",
          main = "Global DNA Methylation Percentage
Following EtOH Treatment",
          fill = "Group",
          palette  = "BuPu", 
          order = c("Untreated", "EtOH 0.02%", "EtOH 0.05%"),
          ylab = expression(Delta~5~mC~"% Compared to Untreated Average"), xlab = "Treatment Conditions") +
  theme(plot.title = element_text(hjust = 0.5, size= 18, face="bold"), text = element_text(size = 15), legend.position = "none")


#Compute the one-way ANOVA test
res.aov <- aov(Methylation_Percent ~ Group, data = My_data)

#Output the summary of the ANOVA
summary(res.aov)

#Okay this tells us there is significant difference but doesn't tell us between which groups so we need to a Tukey test between multiple groups
TukeyHSD(res.aov)
#Looks like there is no difference between the different groups of control and vehicle which is what we would expect.
