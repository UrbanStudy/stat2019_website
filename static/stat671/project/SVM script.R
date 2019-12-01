library(e1071)
iris <- read.csv("C:/Users/shiwei/Google Drive/LAPTOP Backup/STAT 610 Machine Learning I-III/III/iris.csv")[,2:6]
iris <- read.csv("~/Google Drive File Stream/My Drive/LAPTOP Backup/STAT 610 Machine Learning I-III/III/iris.csv")[,2:6]

iris <- subset(iris, iris$Species!="I. versicolor")
iris$Species <- factor(iris$Species)
iris.train <- iris[c(1:40,51:90,101:140),]
iris.test <- iris[c(41:50,91:100,141:150),]

plot(iris.train$Petal.length,iris.train$Petal.width, col=iris.train$Species)
col <- c("Sepal.length","Petal.width","Species")

classifier <- svm(formula = Species~ Petal.width + Petal.length, data = iris.train, 
                  type = 'C-classification', kernel = 'linear')

test_pred <- predict(classifier, type = 'response', newdata = iris.test[-5])

# Making Confusion Matrix
table(iris.test[,5], test_pred)

plot(classifier, iris.train, Petal.width ~ Petal.length,
     slice = list(Sepal.Width = 3, Sepal.Length = 4))

attr(test_pred, "decision.values")
classifier$decision.values
