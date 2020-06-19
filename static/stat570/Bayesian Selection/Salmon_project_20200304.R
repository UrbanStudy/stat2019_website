
rm(list = ls())
# options(digits=8)
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
  mutate(Habitat_level=as.integer(ntile(table_habitat$Habitat_area, 3))) #Creat Habitat Area level

# Change to short names
#original_name <- names(table_habitat) 
#names(table_habitat) <- c("No","H_A","D1_A","D2_A","D1_P","D2_P","W_m_A","W_s_A","W_a_A","L_b_A","L_v_A","W_d_A","W_ia_A","W_r_A","W_L","W_A","W_m_L","S","FE","A","W_c_A","D12_P","H_A_L")
# table <- rbind(original_name[1:8],names(table_habitat)[1:8],
#               original_name[9:16],names(table_habitat)[9:16],
#               original_name[17:24],names(table_habitat)[17:24])

# Transfer the Area to proportion
table_habitat_perc<- table_habitat[,-c(3,4)]
table_habitat_perc[,c(2,5:12,14,18,19)]<- table_habitat[,c(2,7:14,16,20,21)]/table_habitat[,21]
names(table_habitat_perc) <- c("No","H_P","D1_P","D2_P","W_m_P","W_s_P","W_a_P","L_b_P","L_v_P","W_d_P","W_ia_P","W_r_P","W_L","W_P","W_m_L","S","FE","P","W_c_P","D12_P","H_A_L")
glimpse(table_habitat_perc)

# Normalize the variables
# table_habitat[,3:22] <- scale(table_habitat[,3:22], center = T, scale = T)
# table_habitat_perc[,3:20] <- scale(table_habitat_perc[,3:20], center = T, scale = F)




####################################################
ind.insample <- sample(1:178,120)
X <- data.frame(table_habitat[-c(2,23)])
# y <- table_habitat[23] # Use descrete response
y <- table_habitat[2]  # Use continual response
####################################################
# proportion data
y_perc <- table_habitat_perc[2]
X_perc <- data.frame(table_habitat_perc[-c(2,21)])


###################################################

source("VarSelectHC.R")
source("summaryout.R")

########## 7 candidates after AIC selection #######

# vtest <- c("W_c_A","W_s_A","No","S","W_m_A","W_L","A","W_a_A","FE")
vtest <- c('RKM_2008','area_1m','perc_1_2m','MainChannel_Area','SideChannel_Area','Alcove_Area','BareBar_Area','VegetatedBar_Area','InverseAlcove_Area','Bedrock_Area','AllWetLength','MainChannelLength','Slope','Floodplain_elevation')
####################################################

#----------------------------------------------------------------
#with habitat area as a response
datain <- data.frame(y=y[ind.insample,],X[ind.insample,vtest])  # c(vbase,vtest)
data.holdout <- data.frame(y=y[-ind.insample,],X[-ind.insample,vtest])  # c(vbase,vtest)
modpriorvec=c("HOP","HIP","HUP")

# baseformula <- as.formula(paste(".~ ",paste0(vbase,collapse="+")))
theformula <- as.formula(paste("y ~",paste0(vtest,collapse="+"))) # c(vbase,vtest)

res=VarSelectHC(full.formula=theformula,
                data=datain,
                base.formula=as.formula(. ~ 1),#baseformula,#
                maxdeg=2,
                nodes.to.remove=NULL,
                SH = T,
                model.prior.type=modpriorvec,
                model.prior.pars = "children",
                beta.prior.type = "IP",
                beta.prior.pars = list(alpha=1,nu=1),
                niter=5000)

summary.res <- summaryout(mcmc.out=res,insampledata=datain,modelprior.nams=modpriorvec,
                          shr.adj=T,outsampledata=data.holdout,respnam="y",top.ave=10,betaprtype="IP",
                          parsprbeta=list(alpha=1,nu=1))

#----------------------------------------------------------------

vtest <- c("W_c_P","W_s_P","No","S","W_m_P","W_L","P","W_a_P","FE")
#with proportion of habitat area and other variables
datain.prop <- data.frame(y=y_perc[ind.insample,],X_perc[ind.insample,vtest]) # c(vbase,vtest)
data.holdout.prop <- data.frame(y=y_perc[-ind.insample,],X_perc[-ind.insample,vtest]) # c(vbase,vtest)

theformula <- as.formula(paste("y ~",paste0(vtest,collapse="+"))) # c(vbase,vtest)

res.prop=VarSelectHC(full.formula=theformula,
                 data=datain.prop,
                 base.formula=as.formula(. ~ 1),#baseformula,#
                 maxdeg=2,
                 nodes.to.remove=NULL,
                 model.prior.type=modpriorvec,
                 model.prior.pars = "children",
                 beta.prior.type = "IP",
                 beta.prior.pars = list(alpha=1,nu=1),
                 niter=5000)

summary.res.prop <- summaryout(mcmc.out=res.prop,insampledata=datain.prop,modelprior.nams=modpriorvec,
                               shr.adj=T,outsampledata=data.holdout.prop,respnam="y",top.ave=10,betaprtype="IP",
                               parsprbeta=list(alpha=1,nu=1))

# save(file="7plan.RData",
#     list=c("res","summary.res","res.prop","summary.res.prop"))

# model_bayes_perc<- lm(H_P ~ W_s_P+W_m_P+P+W_s_P^2+W_s_P*W_m_P+W_s_P*P+W_m_P^2, data=table_habitat_perc) 
# summary(model_bayes_perc)

# model_bayes<- lm(H_A ~ W_c_A+W_s_A+W_L+W_s_A*W_L, data=table_habitat) 
# ols_regress(model_bayes)

# model_4<- lm(H_A ~ W_m_A+W_s_A+W_a_A+FE+W_m_A*W_a_A+W_s_A*FE, data = table_habitat)
# ols_regress(model_4)
