############## repeat the lecture example

library(Rchoice)
data(Health)
summary(Health$hhinc)
data.Health <- Health %>%subset(hhinc>0) %>%mutate(doc=ifelse(docvis==0,0,1))# %>%subset(hhinc<=10000&year==1986)

data.Health$hhinc <- (data.Health$hhinc)/10000

probit.visit <- Rchoice(doc ~ female + age + hhinc + married + hhkids+ educ+
                          female:married+ female:age + hhinc:age,
                        data = data.Health,
                        family = binomial("probit"))
summary(probit.visit) 

Rchoice(doc ~  hhinc*age ,
        data = data.Health,
        family = binomial("probit"))%>%summary() 

Rchoice(doc ~  hhinc*age + I(age^2) + I(hhinc^2),
        data = data.Health,
        family = binomial("probit"))%>%summary() 

################ The relationship of probabilities, age, income, and their interaction.
kd <- data.frame(probit.visit$mf$hhinc ,
                  probit.visit$mf$age,
                  probit.visit$mf$hhinc * probit.visit$mf$age,
                  probit.visit$probabilities,
                  probit.visit$mf$doc)
names(kd) <- c("inc","age","inc.age","p","doc")
GGally::ggpairs(kd,aes(color =as.factor(doc),alpha=0.1),
                columns=c("inc","age","inc.age","p"))+theme_light()

################# 3D plot of probabilities, age, and income
library(plotly)
kd1 <- kd[kd$p>0.5,]
plot_ly (x=kd1$inc ,y=kd1$age ,z =kd1$p, size=.01, opacity=0.8)#,color = as.factor(kd2$newhsat) #%>%

################# 3D plot of Kernel Density of income v.s. age
# define a 3D grid 
kd2 <- MASS::kde2d(kd1$inc,kd1$age, n = 25)
# Plot the 
plot_ly (type='surface',x=kd2$x,y=kd2$y,z =scale(kd2$z,center=F,scale=T), opacity=0.95)%>% 
 layout(scene = list(yaxis=list(
 nticks = 8,
  range = c(25,40)
 )))

