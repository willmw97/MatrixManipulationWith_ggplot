library(ggplot2)
library(tidyverse)
library(lattice)
library(ash)
library(plotly)

# Compares Matrix A and shows its symmetric when A=A

A <- matrix(c(9,-2,-2,6), nrow = 2, ncol = 2, byrow = T)
A

At <- t(A)
At

#An example of using eigenanalysis of A
eigA <-eigen(A)
eigA

P = eigA$vectors  # Form P whose columns are the eigenvectors of A
P
t(P)  #symmetric

solve(P)  #The inverse of P is the transpose of P

options(digits=2)
#Matrix transformation P'*P
t(P)%*%P # P and its transpose/inverse are orthogonal matrices

#P*P
P%*%t(P) 

Lam = diag(eigA$values) 
Lam

P%*%Lam%*%t(P)  # P(Lam)P' 

# Check that A = P(Lam)
A   

#verfies that PP' = P'P = I using P
P%*%t(P)    #P*P'
t(P)%*%P    #P'*P

#Shows the inverse A^-1
Ainv = solve(A)
Ainv


#Finds egienvalues and egienvectors of A^-1
eigAinv <- eigen(Ainv)
eigAinv

#Demonstrates the spectral decomposition

Pinv <- eigAinv$vectors  # Form P whose columns are the eigenvectors of Ainv
Pinv
t(Pinv)   #symmetric

solve(Pinv)  #inverse of P is the transpose of P

options(digits = 2)
#compute lambda
Laminv = diag(1/eigAinv$values)
Laminv

#Pinv*Lambda-1*Pinv'
P%*%Laminv%*%t(P)

Ainv #check that A^-1 =P*Lam^(-1)*P' 


#Find the matrices A^(1/2) and A^(-1/2).  Verify that A^(1/2) A^(1/2)=A and A^(-1/2) A^(-1/2)=A^(-1)
#Find A^1/2

#Lam^1/2
sqrtLam = diag(sqrt(eigA$values))
sqrtLam

#Lam^-1/2
sqrtLaminv = diag(1/sqrt(eigA$values))
sqrtLaminv

#A^(1/2) = PLam^(1/2)P'
sqrtA = P%*%sqrtLam%*%t(P)
sqrtA
#check the sqaure root matrix condition A = A^(1/2)*A^(1/2)
sqrtA%*%sqrtA
A  #condition checks

#Performs Operation
#A^(-1/2) = PLam^(-1/2)P'
sqrtAinv = P%*%sqrtLaminv%*%t(P)
sqrtAinv
#check the sqaure root matrix condition A^(-1) = A^(-1/2)*A^(-1/2)
sqrtAinv%*%sqrtAinv
Ainv  #condition checks


#Doiung another 3 by 2 matrix
# A = [[4,8,8],[3,6,-9]]

#Matrix declared
#B = UDV'

B = matrix(c(4,8,8,3,6,-9), nrow = 2, ncol = 3, byrow = T)
B #B
t(B) #B'

#B*B'
BBt = B%*%t(B) 
BBt

#B'*B
BtB = t(B)%*%B 
BtB

BBteigen = eigen(BBt)
BBteigen

BtBeigen = eigen(BtB)
BtBeigen


# Obtaining singular value de
#Code obtains SVD of B
Bsvd = svd(B)
Bsvd

U = Bsvd$u
U

D = diag(Bsvd$d)
D

V = Bsvd$v
V

#UDV' = B

#Check that UDV' = B
U%*%D%*%t(V)
B #checks transformation

#Thus when conducting a SVD in R, 
#use of svd function in R rather than 
#doing eigenanalyses of XX′ and X′X directly

#Calculating Eigen Vectors with a 2 by 3 matrix

C = matrix(c(1,1,2,-2,2,2), nrow = 3, ncol = 2, byrow = T)
C

#calculates egienvalues and egienvectors

#Finds matrix manipulation C'
t(C)
#Finds transformation of C*C'
CCt = C%*%t(C)
CCt
#Gets eigen values and eigen vectors
CCteigen = eigen(CCt)
CCteigen

#confirms nonzero eigenvalues are the same

#Batrix manipulation of C'*C
CtC = t(C)%*%C
CtC
#Gets eigen values and eigen vectors
CtCeigen = eigen(CtC)
CtCeigen

#Obtains singular value decomposition

#Gets SVD of C
Csvd = svd(C)
Csvd

U = Csvd$u
U
#finds diagonal bc of constructing diagonal of the matrix
D = diag(Csvd$d)
D

V = Csvd$v
V

#UDV' = B

#Check that UDV' = C
U%*%D%*%t(V)
C #it checks

#Researchers were interested in characterizing fatty acids in olive oils based
#on the region they were grown in Italy.

#These plots help vizualize the distributuin based on the 9 gorwing regions 
# Contstructs distribution plots of all types of acids 
a = ggplot(Olives,aes(Area.Name,Palmitic))
a + geom_point(fill="lightblue") + xlab("Area Name") +
  ylab("Palmitic Acid") + ggtitle("Palmitic Acid vs. Growing Area")

a = ggplot(Olives,aes(Area.Name,Palmitoleic))
a + geom_point(fill="lightblue") + xlab("Area Name") +
  ylab("Palmitoleic Acid") + ggtitle("Palmitoleic Acid vs. Growing Area")

a = ggplot(Olives,aes(Area.Name,Strearic))
a + geom_point(fill="lightblue") + xlab("Area Name") +
  ylab("Strearic Acid") + ggtitle("Strearic Acid vs. Growing Area")

a = ggplot(Olives,aes(Area.Name,Oleic))
a + geom_point(fill="lightblue") + xlab("Area Name") +
  ylab("Oleic Acid") + ggtitle("Oleic Acid vs. Growing Area")

a = ggplot(Olives,aes(Area.Name,Linoleic))
a + geom_point(fill="lightblue") + xlab("Area Name") +
  ylab("Linoleic Acid") + ggtitle("Linoleic Acid vs. Growing Area")

a = ggplot(Olives,aes(Area.Name,Linolenic))
a + geom_point(fill="lightblue") + xlab("Area Name") +
  ylab("Linolenic Acid") + ggtitle("Linolenic Acid vs. Growing Area")

a = ggplot(Olives,aes(Area.Name,Eicosanoic))
a + geom_point(fill="lightblue") + xlab("Area Name") +
  ylab("Eicosanoic Acid") + ggtitle("Eicosanoic Acid vs. Growing Area")

a = ggplot(Olives,aes(Area.Name,Eicosenoic))
a + geom_point(fill="lightblue") + xlab("Area Name") +
  ylab("Eicosenoic Acid") + ggtitle("Eicosenoic Acid vs. Growing Area")


#Constructs a 2D plot of Olive Acids
pairs(Olives[,2:9], pch = 19)

olive.mat = Olives[,c(3,4,5,6,7,8,9)]
pairs(olive.mat)

# Constructs a semi-3D plot Olive Acids with a color code based on the region their from
splom(~Olives[,2:9],groups=Area.Name,auto.key=list(columns=7),data=Olives)

pairs.persp <- function(x) {
  par(bg="sky blue")
  pairs(x,panel=function(x,y) {
    foo <- bin2(cbind(x,y),nbin=c(8500,8500))
    foo <- ash2(foo,m=c(8,8))
    par(new=T)
    persp(foo,xlab="",ylab="",theta=-45,phi=35,col="royal blue",
          shade=.75,box=F,scale=F,border=NA,expand=.9)
  })
  par(new=F,bg="white")
}
pairs.persp(Olives[,2:9])


#Observations show Palmitic acid seems to doverge the most.
#Olic acid in Umbia is also noteworthy becasue of the little variation it has
#This show by in the graph becasue of the "cigar" like shape shown

#When looking at which region is the most homogenous, Olive oils in Sicily would be based on its large variation

pairs(Olives[2:9,], pch = 19)

olive.mat = Olives[c(3,4,5,6,7,8,9),]
pairs(olive.mat)
a = ggplot(Olives,aes(Area.Name,Olives[,2:9]))
a + geom_boxplot(fill="lightblue") + xlab("Area Name") +
  ylab(" Acid") + ggtitle("Oleic Acid vs. Growing Area")

pairs.trendsd <- function(data,...) { 
  pairs(data,lower.panel=panel.cor,upper.panel=panel.trendsd,
        diag.panel=panel.hist,...)}

pairs.trendsd(Olives[,c(6,5,7)])


#The next follwoing plots use ggplot to help better visualize the distribution of fatty acids

#This was just testing the histogram on city data
City = read.csv("City77.csv")
a = ggplot(City,aes(medinc))
a + geom_histogram(binwidth=1000,aes(y=..density..),fill="blue",color="black") + 
  geom_density(linetype=2,linewidth=2,col="red") + 
  xlab("Median Income") + 
  ylab("Density") +
  ggtitle("Histogram of Median Income") +
  theme_classic()

#Creates a box plot based on the olive data
Olives = read.csv("Olives.txt")
head(Olives)
a = ggplot(Olives,aes(Area.Name,Oleic))
a + stat_boxplot(geom ='errorbar') +
  geom_boxplot(fill="white") + xlab("Area Name") + 
  ylab("Oleic Acid") + ggtitle("Oleic Acid vs. Growing Area") +
  theme_classic()


#Creates a violin plot based on the Olive Data
Olives = read.csv("Olives.txt")
names(Olives)

ggplot(Olives, aes(x = Area.Name, y = Oleic)) +
  geom_violin(trim = FALSE, width = 2) + 
  stat_summary(Olives.Oleic=mean_sdl,  Olives.Area.Name = list(mult=1), 
               geom="pointrange", color = "blue") +
  theme_classic()



#Does a pie chart with bar graph from admisssions data

#Pie charts and bar graphs

data("UCBAdmissions")
UCBAdmissions
temp = data.frame(UCBAdmissions)
temp

GenderAdmit = margin.table(UCBAdmissions,c(1,2))
gender = data.frame(GenderAdmit) 
gender

admitted_rejected <- data.frame(Admit = c("Admitted", "Rejected"), Male = c(1198, 1493),
                                Female = c(557, 1278))
head(admitted_rejected)

#Pie chart for males
ggplot(admitted_rejected, aes(x = "", y = Male, fill = Admit))+ 
  geom_bar(width = 3, stat = "identity", color = "white") +
  coord_polar("y", start = 1.5) + 
  theme_void() + 
  labs(title = "Males") +
  scale_fill_brewer("Males") 

#Pie chart for females
ggplot(admitted_rejected, aes(x = "", y = Female, fill = Admit)) + 
  geom_bar(width = 3, stat = "identity", color = "white") + 
  coord_polar("y", start = 1.5) + 
  theme_void() + 
  ggtitle("Female") + 
  scale_fill_brewer("Female")


#bar graph 1: stacked

ggplot(gender, aes(x = gender$Gender,y = gender$Freq, fill=gender$Admit)) + 
  geom_bar(stat="identity",color='black',width = 0.4,position=position_stack()) + 
  guides(fill = "none") + 
  labs(x = "Gender", 
       y = "Frequency") +
  scale_fill_manual(values = c("lightgrey","black")) +
  theme_classic()

#bar graph 2: side by side
ggplot(gender, aes(x = gender$Gender, y = gender$Freq, fill=gender$Admit)) + 
  geom_bar(stat="identity",color='black',width = 0.4, position="dodge") + 
  guides(fill = "none")+xlab('Gender')+ylab('Frequency') + 
  scale_fill_manual(values = c("black","lightgrey")) +
  theme_classic()

#Looks the built in Iris data and creates a facet based on Petal and Sepal length
library(lattice)
data(iris)
names(iris)

t <-ggplot(iris, aes(x=Petal.Length, y=Sepal.Length)) +
  geom_point(shape=1, color="blue") + facet_grid(.~iris$Species) 

t + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_rect("lightpink"),
          axis.line = element_line(colour = "black"), strip.background = element_rect(fill = "yellow")
) +
  labs(x="Petal Length", y="Sepal Length")


#Creates a 2-D density plot based Swiss data showing the density of the bill diagnoal length and 
#its height
library(ash)
Swiss = read.csv("Swiss.csv")

m <- ggplot(Swiss, aes(x = diagon, y = right)) +
  geom_point(pch=as.character(Swiss$genu),cex=.4)+
  xlim(137, 143) +
  ylim(128.5, 131.5) 
m + stat_density_2d(aes(fill = stat(density)), geom = 'raster', contour = FALSE) +       
  scale_fill_gradient2(low = "red", mid = "orange", high = "yellow", midpoint = .2)+
  coord_cartesian(expand = FALSE) +
  geom_density2d(color="black")+
  geom_point(shape = '.', col = 'black')

#Creates a bubble plot based on the amount of goals for anf against with
#the circle size shwoing the amount of wins
NHL = read.csv("NHL.csv")
names(NHL)
NHL$W
ggplot(data=NHL, aes(x=GF, y=GA,label=TEAM),guide=FALSE)+
  geom_point(aes(size = NHL$W),shape = 21,color="black",fill="lightblue") + scale_size_continuous(range = c(2, 6)) +
  scale_size(range = c(1,12)) +
  scale_x_continuous(name="Goals For", limits=c(150,300))+
  scale_y_continuous(name="Goals Against", limits=c(150,300))+
  geom_text(label =NHL$TEAM,size=2,nudge_x = -2,nudge_y = -2)+
  theme_classic()+
  theme(legend.title=element_blank())




