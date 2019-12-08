library(e1071)
data(iris)
# iris <- read.csv("C:/Users/shiwei/Google Drive/LAPTOP Backup/STAT 610 Machine Learning I-III/III/iris.csv")[,2:6]
# iris <- read.csv("~/Google Drive File Stream/My Drive/LAPTOP Backup/STAT 610 Machine Learning I-III/III/iris.csv")[,2:6]

# iris <- subset(iris, iris$Species!="I. versicolor")
# iris$Species <- factor(iris$Species)
iris.train <- iris[c(1:40,51:90,101:140),]
iris.test <- iris[c(41:50,91:100,141:150),]

plot(iris.train$Petal.Length,iris.train$Petal.Width, col=iris.train$Species)
col <- c("Sepal.Length","Petal.Width","Species")

classifier <- svm(formula = Species~ Petal.Width + Petal.Length, data = iris.train, 
                  type = 'C-classification', kernel = 'linear')

test_pred <- predict(classifier, type = 'response', newdata = iris.test[-5])

# Making Confusion Matrix
table(iris.test[,5], test_pred)

plot(classifier, iris.train, Petal.Width ~ Petal.Length,
     slice = list(Sepal.Width = 3, Sepal.Length = 4))

attr(test_pred, "decision.values")
classifier$decision.values
