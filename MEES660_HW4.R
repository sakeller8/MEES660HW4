library(vegan)
library(MASS)

fish <- read.csv("C:/Users/socce/OneDrive/Desktop/Fall 2022/MEES660/Quantitative Exercises/MEES660.HW4/fish_comm.csv",
                 header=TRUE,sep=",",row.names=1) 
#import fish assemblage data
#Note sample IDs read in as labels rather than as actual data

class(fish) #imported correctly as a data.frame
dim(fish) #check dimensions of the data.frame

head(rownames(fish)) #lists first 6 row names

head(colnames(fish)) #lists first 6 column names

#import metadata on site information
metadata <- read.csv("C:/Users/socce/OneDrive/Desktop/Fall 2022/MEES660/Quantitative Exercises/MEES660.HW4/site_data.csv",header=TRUE,sep=",",row.names=1)

head(rownames(metadata)) #looks good

head(colnames(metadata)) #looks good

#verify our fish data and our metadata match in terms of row names
all.equal(rownames(fish),rownames(metadata))

boxplot(specnumber(fish) ~ metadata$Habitat, ylab = "Species richness")
boxplot(diversity(fish) ~ metadata$Habitat, ylab = "Shannon index")

#compare species richness and Shannon diversity between estuaries
op = par(mfrow=c(1,2)) #requests 2 panel boxplot space
boxplot(specnumber(fish) ~ metadata$Habitat, ylab = "Species richness")
boxplot(specnumber(fish) ~ metadata$Season, ylab = "Species richness")

boxplot(diversity(fish) ~ metadata$Habitat, ylab = "Shannon diversity H'")
boxplot(diversity(fish) ~ metadata$Season, ylab = "Shannon diversity H'")

#### I just recently learned how to reorder classes in figures such as the 
# Season classes above. You can use the following line of code to specify
# the order, then replotting:
metadata$Season = factor(metadata$Season, levels=c("Spring", "Early summer", "Summer"))

par(op) # Restores plot parameters 

#use anova to test for a differences in species richness
R_test <- aov(specnumber(fish) ~ metadata$Habitat + metadata$Season +
                  metadata$Habitat*metadata$Season)
drop1(R_test,~.,test="F") #look at output: Main effects and interaction not significant

R_test2 <- aov(specnumber(fish) ~ metadata$Habitat + metadata$Season)
drop1(R_test2,~.,test="F") #look at output:Season and Habitat are significant

TukeyHSD(R_test2) #test for pairwise differences between estuaries and among seasons

#let's check our ANOVA assumptions by visually inspecting our residuals
op = par(mfrow=c(2,2)) 
boxplot(residuals(R_test2) ~ metadata$Habitat, ylab = "Richness residuals")
boxplot(residuals(R_test2) ~ metadata$Season, ylab = "Richness residuals")
hist(residuals(R_test2))
qqnorm(residuals(R_test2))
par(op) # Restore plot parameters 

#use anova to test for a differences in Shannon diversity
H_test <- aov(diversity(fish) ~ metadata$Habitat + metadata$Season +
                  metadata$Habitat*metadata$Season)
drop1(H_test,~.,test="F") #look at output: Main effects and interaction not significant

set.seed(321)
lnfish <- log1p(fish) #log(x+1) transformation of fish assemblage data
fish.mds <- metaMDS(lnfish, trace = FALSE) 
#wrapper function that calculates dissimilarities, ordinations, stress
fish.mds #call for results

op = par(mfrow=c(1,1)) 
plot(fish.mds, type = "t") #plotting our 2-D nMDS ordination

sh_name <- make.cepnames(names(lnfish)) #shorten names
sh_name[1:5] #check on shorter names by looking at first 5 species names

abund <- colSums(lnfish) #sums the columns of the lnfish dataset
plot (fish.mds, display = "sites", type = "p") #plots nmds of sites only
sel <- orditorp(fish.mds, dis="sp", lab=sh_name, priority= abund, pcol = "gray", pch="+") #this overlays our species data points with labeled short names

ef <- envfit(fish.mds, metadata, permutations = 999) #runs envfit analysis
ef #call for the results

plot(fish.mds, display = "sites", type = "p")
with(metadata, ordiellipse(fish.mds, Season, kind = "se", conf = 0.95))
with(metadata, ordispider(fish.mds, Season, col = "blue", label= TRUE))
with(metadata, ordihull(fish.mds, Season, col="blue", lty=2))

plot(fish.mds, display = "sites", type = "p")
with(metadata, ordihull(fish.mds, Season, col="blue", lty=2))
sel <- orditorp(fish.mds, dis="sp", lab=sh_name, priority= abund, pcol = "gray", pch="+") #this overlays our data points with labeled short names
