rm(list = ls())
library(tidyverse)
# Import data
table_habitat <- read.csv("Willamette_habitat_features.csv")

# fix some wrong values
table_habitat[42,]$Slope <- 0.000732
table_habitat[42,]$Floodplain_elevation<- 4.214
table_habitat[41,]$Slope <- 0.000586

# Add tiny value for NA
table_habitat[is.na(table_habitat)] <- 1e-8
# Reasign the index
table_habitat <-table_habitat[order(table_habitat$RKM_2008,decreasing = F),]
table_habitat$RKM_2008 <- 1:178

# Add new columns
table_habitat <- table_habitat%>%mutate(ConnectedWet_area=AllWetArea-DisconnectedWater_Area)%>% #Creat ConnectedWet_area
  mutate(perc_1_2m=perc_2m-perc_1m)%>% #Creat pure Area_2m
  mutate(Habitat_level=as.factor(ntile(table_habitat$Habitat_area, 3))) #Creat Habitat Area level

# Change to short names
original_name <- names(table_habitat) 
names(table_habitat) <- c("No","H_A","D1_A","D2_A","D1_P","D2_P","W_m_A","W_s_A","W_a_A","L_b_A","L_v_A","W_d_A","W_ia_A","W_r_A","W_L","W_A","W_m_L","S","FE","A","W_c_A","D12_P","H_A_L")
table <- cbind(original_name[1:23],names(table_habitat)[1:23])



# Transfer the Area to proportion
table_habitat_perc<- table_habitat[,-c(3,4)]
table_habitat_perc[,c(2,5:12,14,18,19)]<- table_habitat[,c(2,7:14,16,20,21)]/table_habitat[,21]
names(table_habitat_perc) <- c("No","H_P","D1_P","D2_P","W_m_P","W_s_P","W_a_P","L_b_P","L_v_P","W_d_P","W_ia_P","W_r_P","W_L","W_P","W_m_L","S","FE","P","W_c_P","D12_P","H_A_L")
glimpse(table_habitat_perc)


write.csv(table_habitat_perc,"table_habitat_perc.csv", row.names =T)
write.csv(table_habitat,"table_habitat.csv", row.names =T)

# Normalize the variables
# table_habitat[,3:22] <- scale(table_habitat[,3:22], center = T, scale = F)