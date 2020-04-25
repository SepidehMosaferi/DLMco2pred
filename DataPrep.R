# ============================================================================================
# author: Sepideh Mosaferi
# Dynamic Linear Models
# code name: Data Preparation 
# version: 0.1 April 2019
# ============================================================================================

# Required Libraries -------------------------------------------------------------------------

library("car")
library("ggplot2")

setwd("/home/sepidehmosaferi/DLM")

Corn <- read.csv("Corn.csv" , header = TRUE)

Split1 <- strsplit(as.character(Corn$TIMESTAMP) , split = "/" , fixed = TRUE)

Corn$Month <- sapply(1:nrow(Corn) , function(i){as.numeric(Split1[[i]][1])})

Corn$Day <- sapply(1:nrow(Corn) , function(i){as.numeric(Split1[[i]][2])})

YearTime <- sapply(1:nrow(Corn) , function(i){Split1[[i]][3]})

Corn$Year <- as.numeric(substr(YearTime , 1 , 2))

Split2 <- sapply(1:nrow(Corn) , function(i){strsplit(YearTime[[i]] , split = " " , fixed = TRUE)})

Corn$Time <- sapply(1:nrow(Corn) , function(i){Split2[[i]][2]})

summary(Corn)

# Remove the missing values from the data ----------------------------------------------------

Corn <- Corn[complete.cases(Corn) , ]

Corn$Season <- ifelse(Corn$Month %in% c(1 , 2 , 12) , "Winter" , 
                      ifelse(Corn$Month %in% c(3 , 4 , 5) , "Spring" , 
                             ifelse(Corn$Month %in% c(6 , 7 , 8) , "Summer" , "Fall")))

Corn <- Corn[order(Corn$Year , Corn$Month , Corn$Day , Corn$Night.Day.Flag , Corn$TIMESTAMP) , ]

# Subsetting the Data set --------------------------------------------------------------------
# Only Summer 2015 be considered in the data -------------------------------------------------

CornSub <- subset(Corn , Corn$Season == "Summer" & Corn$Year == 15)

range(CornSub$Fc_wpl..mg.m.2.s1.)

DNSummer <- suppressWarnings(cbind(CornSub$Fc_wpl..mg.m.2.s1.[CornSub$Night.Day.Flag == 1] ,
                                   CornSub$Fc_wpl..mg.m.2.s1.[CornSub$Night.Day.Flag == 0]))
colnames(DNSummer) <- c("Day" , "Night")

range(CornSub$Fc_wpl..mg.m.2.s1.[CornSub$Night.Day.Flag == 1])
range(CornSub$Fc_wpl..mg.m.2.s1.[CornSub$Night.Day.Flag == 0])

p1 <- ggplot(CornSub , aes(x = CornSub$Fc_wpl..mg.m.2.s1.)) +
  geom_density(color = "black" , fill = "gray") + theme_bw() + scale_colour_hue(l = 25) +
  ggtitle("Kernel Density Plot of CO2 for Summer 2015") +
  theme(plot.title = element_text(hjust = 0.5)) +
    labs(x = "CO2 Flux", y = "Density" , lwd = 1.5) +
  annotate("text" , x = -2 , y = 0.65 , label = "range:[-3.06,1.28]" , size = 5)

Group <- ifelse(CornSub$Night.Day.Flag == 1 , "Day" , "Night")

df <- data.frame(group = Group , Value = CornSub$Fc_wpl..mg.m.2.s1.)

p2 <- ggplot(df , aes(x = group , y = Value)) + geom_boxplot() + theme_bw( ) +
  ggtitle("Boxplot of CO2 for Summer 2015")

cowplot::plot_grid(p1 , p2 , labels = c("A" , "B") , align = "v")

CornSub$Time <- paste(CornSub$Month , CornSub$Day , CornSub$Night.Day.Flag)
uniqueTime <- unique(CornSub$Time)
OBS <- sapply(1:length(unique(CornSub$Time)) ,
              function(i){mean(CornSub$Fc_wpl..mg.m.2.s1.[CornSub$Time == unique(uniqueTime)[i]])})
NewData <- as.data.frame(cbind(1:length(uniqueTime) , OBS))
colnames(NewData) <- c("Time" , "OBS")

par(mfrow = c(2 , 2))
ts.plot(NewData$OBS , col = "black" , lwd = 1.5 , ylab = "CO2 Flux" ,
        main = "Time-Series Plot of Aggregated \n Values for CO2 Flux")
acf(NewData$OBS , ylim = c(-1 , 1) , main = "Autocorrelation Plot")
pacf(NewData$OBS , ylim = c(-1 , 1) , main = "Partial Autocorrelation Plot")
qqnorm(NewData$OBS , ylim = c(-2 , 2))
qqline(NewData$OBS , col = "blue")
shapiro.test(NewData$OBS)
text(-0.75 , 1.5 , paste("Shapiro-Wilk Statistic: 0.933 \n P-value: 4.384e-05"))

# Box-Cox transformation BC of negative values -----------------------------------------------

NewData$OBSpos <- NewData$OBS - min(NewData$OBS[which(NewData$OBS < 0)]) + 1
LAMBDA <- as.vector(coef(powerTransform(NewData$OBSpos)))  # optimal lambda is "-0.023".
NewData$OBS <- ((NewData$OBSpos^LAMBDA) - 1) / LAMBDA
FinalData <- as.data.frame(cbind(NewData$Time , NewData$OBS))
colnames(FinalData) <- c("Time" , "OBS")
summary(lm(FinalData$OBS ~ FinalData$Time))

# END ========================================================================================

