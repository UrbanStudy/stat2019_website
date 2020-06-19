#Load libraries
library(ggplot2)
library(reshape)
library(randomForest)
library(dplyr)

#import data
data<-read.csv('Willamette_habitat_features.csv')

#melt dataset with habitat area and RKM to plot by feature
features_melt<-melt(data, id = c("Habitat_area", "RKM_2008"))

#plot by each feature
ggplot(features_melt, aes(Habitat_area, value))+
  geom_point()+
  facet_wrap(~variable, scales = "free")

#######################
##begin random forest##
#######################

#replace NA with 0
data[is.na(data)]<-0

#set new data source to manipulate
data_model<-data

#Set seed
set.seed(123)

#run with out of the box parameters
model_low<-randomForest(Habitat_area ~ ., data = data_model, proximity = FALSE)

#Summarize/view model results
model_low
varImpPlot(model_low)
plot(model_low)

##Remove less important datasets to check for improved model performance
#set new dataframe
data_model_sub<-data_model
#RKM appears relatively important, but not a physical variable
data_model_sub$RKM_2008<-NULL
data_model_sub$InverseAlcove_Area<-NULL
data_model_sub$Bedrock_Area<-NULL
data_model_sub$DisconnectedWater_Area<-NULL

#rerun model
model_low_sub<-randomForest(Habitat_area ~ ., data = data_model_sub, proximity = FALSE)

#summarize/view model results
model_low_sub
varImpPlot(model_low_sub)
plot(model_low_sub)


##Try with classifiction by breaking area into factors##
#create new df
data_model_fac<-data_model_sub

#Seperate data into n quartiles
n<-3
data_model_fac$quartile<-as.factor(ntile(data_model$Habitat_area, n))

#remove habitat area data to prevent model estimating from it
data_model_fac$Habitat_area<-NULL

#run model
model_fac<-randomForest(quartile ~., data = data_model_fac, importance = TRUE)
model_fac
varImpPlot(model_fac)
plot(model_fac)


write.csv(cbind(model_low$predicted,model_low$y),"prediction_random_Forest.csv", row.names = F)
write.csv(cbind(model_low_sub$predicted,model_low_sub$y),"prediction_random_Forest_sub.csv", row.names = F)
write.csv(cbind(model_fac$predicted,model_fac$y),"prediction_random_Forest_fac.csv", row.names = F)


