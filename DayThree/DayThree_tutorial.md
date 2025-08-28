### Day Three
Today, we'll again be working with the Qfly dataset. We'll use three key methods to determine the significance of predictor climatic variables 
in explaining allele frequency patterns in the native and introduced ranges of the Queensland fruit fly:

1. RDA
3. LFMM
4. Gradient Forest
   
Let's get into it!

### Open R and load in required packages
First, set you working directory to something logical (e.g., a 'DayThree' folder on your desktop), and run all code from within that directory.
IMPORTANT - ONLY LOAD THESE PACKAGES FOR NOW as there will be problems later otherwise.
```
library(geodata)
library(raster)
library(ozmaps)
```
### Extract climate variables
First, we need to extract relevant climate variables for all sites. We'll extract these from the WorldClim database, as follows:
```
raster_path_folder <- tempdir() 
raster_stack <- worldclim_country(country = "AUS",path = raster_path_folder ,version="2.1",res=0.5, var="bio")
```
Now, let's extract the 19 bioclimate variables to individual raster files:
```
bio_1 <- raster(raster_stack$wc2.1_30s_bio_1)
bio_2 <- raster(raster_stack$wc2.1_30s_bio_2)
bio_3 <- raster(raster_stack$wc2.1_30s_bio_3)
bio_4 <- raster(raster_stack$wc2.1_30s_bio_4)
bio_5 <- raster(raster_stack$wc2.1_30s_bio_5)
bio_6 <- raster(raster_stack$wc2.1_30s_bio_6)
bio_7 <- raster(raster_stack$wc2.1_30s_bio_7)
bio_8 <- raster(raster_stack$wc2.1_30s_bio_8)
bio_9 <- raster(raster_stack$wc2.1_30s_bio_9)
bio_10 <- raster(raster_stack$wc2.1_30s_bio_10)
bio_11 <- raster(raster_stack$wc2.1_30s_bio_11)
bio_12 <- raster(raster_stack$wc2.1_30s_bio_12)
bio_13 <- raster(raster_stack$wc2.1_30s_bio_13)
bio_14 <- raster(raster_stack$wc2.1_30s_bio_14)
bio_15 <- raster(raster_stack$wc2.1_30s_bio_15)
bio_16 <- raster(raster_stack$wc2.1_30s_bio_16)
bio_17 <- raster(raster_stack$wc2.1_30s_bio_17)
bio_18 <- raster(raster_stack$wc2.1_30s_bio_18)
bio_19 <- raster(raster_stack$wc2.1_30s_bio_19)
```
You can see the definition for each of the bio variables here: https://www.worldclim.org/data/bioclim.html.
Try plotting a few of these to see what they look like:
```
plot(bio_1)
plot(bio_18)
```
Now, let's pull out just the locations we're interested in:
```
places <- read.delim("Qfly_coordinates.txt", header = TRUE)
coords<-data.frame(lon=places[,2], lat=places[,3])
coords$lon <- as.numeric(coords$lon)
coords$lat <- as.numeric(coords$lat)
coordinates(coords) <- c("lon","lat")
```
Double check the distribution of points on the map (except for pop 28, which is out of scale):
```
ozmap(x = "country")
points(coords, pch=16)
```
Next, use the "extract" function to get climate variables just for our 28 locations:
```
val_bio_1 <- extract(x=bio_1, y=coords)
val_bio_2 <- extract(x=bio_2, y=coords)
val_bio_3 <- extract(x=bio_3, y=coords)
val_bio_4 <- extract(x=bio_4, y=coords)
val_bio_5 <- extract(x=bio_5, y=coords)
val_bio_6 <- extract(x=bio_6, y=coords)
val_bio_7 <- extract(x=bio_7, y=coords)
val_bio_8 <- extract(x=bio_8, y=coords)
val_bio_9 <- extract(x=bio_9, y=coords)
val_bio_10 <- extract(x=bio_10, y=coords)
val_bio_11 <- extract(x=bio_11, y=coords)
val_bio_12 <- extract(x=bio_12, y=coords)
val_bio_13 <- extract(x=bio_13, y=coords)
val_bio_14 <- extract(x=bio_14, y=coords)
val_bio_15 <- extract(x=bio_15, y=coords)
val_bio_16 <- extract(x=bio_16, y=coords)
val_bio_17 <- extract(x=bio_17, y=coords)
val_bio_18 <- extract(x=bio_18, y=coords)
val_bio_19 <- extract(x=bio_19, y=coords)
df <- data.frame(val_bio_1,val_bio_2,val_bio_3,val_bio_4,val_bio_5,val_bio_6,val_bio_7,val_bio_8,val_bio_9,val_bio_10,val_bio_11,val_bio_12,val_bio_13,val_bio_14,val_bio_15,val_bio_16,val_bio_17,val_bio_18,val_bio_19)
write.table(df, "Qfly_AllBioclim.txt", sep = "\t")
```
For later analysis, we'll need to remove highly correlated (|r| > 0.7) bioclimatic variables, so let's identify these now:
```
library(psych)
library(caret)
my_cor <- cor.plot(df)
hc  <- findCorrelation(my_cor, cutoff=0.70)
hc  <- sort(hc)
reduced_Data = df[,-c(hc)]
print(reduced_Data)
```
You should find that five variables are uncorrelated enough to use to represent the full dataset: val_bio_4, val_bio_5, val_bio_9, val_bio_17, val_bio_18. Remember these for later!
You'll also maybe spot that there is missing data for the three populations that are out of scale. Don't worry, we'll provide the full dataset to read in later.

### Generate required allele frequency data
```
library(robust)
library(qvalue)
library(dplyr)
library(vegan)
library(vcfR)
library(adegenet)
```
Read in the VCF and convert it to a matrix (individuals in rows and SNPs in columns).
SNPs will also need to be converted to genotype format with 0, 1, and 2: 0 when the individual is homozygous for the major allele, 1 when the individual is heterozygous, and 2 when the individual is homozygous for the second (or alternative) allele.
```
data <- read.vcfR("Qfly.vcf")
data
geno <- extract.gt(data)
dim(geno)

G <- geno
G[geno %in% c("0/0")] <- 0  
G[geno  %in% c("0/1")] <- 1  
G[geno %in% c("1/1")] <- 2  
```
Transpose the genotype matrix:
```
gen <- t(G)   
dim(gen)
```
Look for any missing data:
```
sum(is.na(gen))  
```
We can see there are missing data here (62935 out of a total of 2018807 = ~3%). We will impute them below to fill in these gaps.
Next, we need to assign population information to the genotype matrix (we use the nomatch part of the code to make sure all the individuals in the popmap file are present in our genotype file).
```
Genotypes <- as.data.frame(gen)
popmap <- read.table("Qfly_popmap.txt")
Genotypes <- Genotypes[match(popmap$V1, row.names(Genotypes), nomatch = 0),]
write.table(Genotypes,"Genotypes.txt", sep = "\t")
```
We also need to do a bit of work on the Genotype file, as later steps don't like the presence of the ':', which is in the genotype names (e.g., NC_052499.1:84).
I've replaced the ':' with an underscore using sed in the command line as the file is too big to open with something like notepad.
Sed works just like a search and replace: sed 's/:/_/g' Genotypes.txt > Genotypes2.txt.
Now, read in the correctly formatted file:
```
Genotypes2 <- read.delim("Genotypes2.txt")
```
We will use population allele frequencies (ranging from 0 to 1) instead of individual genotypes in the subsequent analyses.
This is because we have several individual genotypes at each sampling site, meaning they all experienced the exact same climatic conditions. 
Using population allele frequencies will also help with any uneven sample sizes across populations.
Estimating population allele frequencies:
```
AllFreq <- aggregate(Genotypes2, by = list(popmap$V2), function(x) mean(x, na.rm = T)/2)
row.names(AllFreq) <- as.character(AllFreq$Group.1)
```
Now, impute the missing genotypes using the median of the locus allele frequencies across all populations:
```
for (i in 1:ncol(AllFreq))
{
  AllFreq[which(is.na(AllFreq[,i])),i] <- median(AllFreq[-which(is.na(AllFreq[,i])),i], na.rm=TRUE)
}
```
And, filter on minor allele frequency (MAF) again (although this was done during the original SNP calling to produce the Qfly.vcf file, we may have now imputed some sites with small MAF, which we want to remove:
```
freq_mean <- colMeans(AllFreq[,-1])
```
We now have a final genetic matrix that includes allele frequncy data for 28 populations and 6,527 SNPs (as we indeed lost some sites with the MAF cutoff):
```
AllFreq <- AllFreq[,-which(freq_mean>=0.95 | freq_mean<=0.05)] 
write.table(AllFreq, "Qfly_alleleFreq.txt", sep = "\t") 
```
### RDA analysis
Partial RDA (pRDA) analysis can be used to identify the contribution of climate, population structure, and geography in explaining genetic variation. 
We will therefore use the 5 bioclimate variables (environmental data); population scores along the first three axes of a genetic PCA (genetic data); 
and population longitude and latitude coordinates (geographic variation) in our analysis.
First, read in the environmental and genetic data:
```
env_data <- read.delim("Qfly_AllBioclim_Final.txt", header = TRUE)
AllFreq <- read.table("Qfly_alleleFreq.txt")
```
Now, standardise the env_data to avoid discrepancy in mean and standard deviation among variables and ensure 
that the variable units are comparable (e.g., precipitation in mm, and temperature in °C).
```
env_data_std <- scale(env_data[,-1:-3], center=TRUE, scale=TRUE)
```
Recover scaling coefficients:
```
scale_env <- attr(env_data_std, 'scaled:scale')
center_env <- attr(env_data_std, 'scaled:center')
```
Generate a climatic table:
```
env_data_std <- as.data.frame(env_data_std)
row.names(env_data_std) <- c(env_data$pop)
```
Infer pop structure from allele frequnency data:
```
pca <- rda(AllFreq[,-1], scale=T)
screeplot(pca, type = "barplot", npcs=10, main="PCA Eigenvalues")  
```
We'll use the first three PCs:
```
PCs <- scores(pca, choices=c(1:3), display="sites", scaling=0)
PopStruct <- data.frame(Population = AllFreq[,1], PCs)
colnames(PopStruct) <- c("Population", "PC1", "PC2","PC3")
```
Now merge the PC data with the coordinates data (places, from above), and environmental data:
```
Variables <- data.frame(places, PopStruct[,-1], env_data) 
```
Subset this data to retain our uncorrelated bio variables from above:
```
env_subset <- subset(Variables, select = c(longitude,latitude,PC1,PC2,PC3,bio_4,bio_5,bio_9,bio_17,bio_18))
```
Look at the distribution of these remaining variables:
```
pairs.panels(env_subset[,6:10], scale=T)
```
Looks pretty good! We're ready to run the full model:
```
pRDAfull <- rda(AllFreq[,-1] ~ PC1 + PC2 + PC3 + longitude + latitude + bio_4+bio_5+bio_9+bio_17+bio_18, Variables)
RsquareAdj(pRDAfull) 
anova(pRDAfull)
```
Based on the RsquareAdj value, we can see that the full model explains about 15% of variation.
Let's now run models based on genetic and geography matrices (where env + lat/long, and env + genetic are factored out, respectively):
```
pRDAstruct <- rda(AllFreq[,-1] ~ PC1 + PC2 + PC3 + Condition(longitude + latitude +  bio_4,bio_5,bio_9,bio_17,bio_18),  Variables)
RsquareAdj(pRDAstruct)  
anova(pRDAstruct)

pRDAgeog <- rda(AllFreq[,-1] ~ longitude + latitude + Condition(bio_4+bio_5+bio_9+bio_17+bio_18+ PC1 + PC2 + PC3),  Variables)
RsquareAdj(pRDAgeog) #
anova(pRDAgeog)
```
You should find that the genetic model is significant and explains about 16% of variation, while the geographic one is non-significant and explains <1% of variation.
The variation explained is not very high for the full and genetic models, but this low explanatory power is not surprising given that we expect most of the SNPs in our dataset will not show a relationship with the environmental predictors (e.g., most SNPs will be neutral).
Let's look at the full model in more detail now:
```
summary(eigenvals(pRDAfull, model = "constrained"))
screeplot(pRDAfull)
```
Here, we can see that the first three constrained axes explain ~52% of the total variance. 
We can also check our RDA model for significance using formal tests, assessing both the full model and each constrained axis using F-statistics. 
The null hypothesis is that no linear relationship exists between the SNP data and the environmental predictors.
```
signif.full <- anova.cca(pRDAfull, parallel=getOption("mc.cores")) 
signif.full
```
The full model is significant, as we already saw above. For the next tests, each constrained axis is tested using all previous constrained axes as conditions. 
The purpose is to determine which constrained axes we should investigate for candidate loci:
```
signif.axis <- anova.cca(pRDAfull, by="axis", parallel=getOption("mc.cores"))
signif.axis
```
The first two axes are significant and the third is marginally non-significant, consistent with the scree plot from earlier.

Finally, vegan has a simple function for checking Variance Inflation Factors for the predictor variables used in the model:
```
vif.cca(pRDAfull)
```
All values are below 10, which means that multicolinearity among these predictors shouldn’t be a problem for the model.

Let's now plot the PCA, first for axes 1 and 2; then for 1 and 3:
```
plot(pRDAfull, scaling=3)     
plot(pRDAfull, choices = c(1, 3), scaling=3)
```
In these plots, the SNPs are in red (in the centre of each plot), and the individuals are the black circles. 
The blue vectors are the environmental predictors. The relative arrangement of these items in the ordination space reflects their relationship with the ordination axes, which are linear combinations of the predictor variables.
Note that these plot are not very useful for our dataset; you can look at the RDA tutorial to a nicer visual example.

Next, we want to identify the candidate outliers.
To do this, we use the loadings of the SNPs in the ordination space to determine which SNPs are candidates for local adaptation.
The SNP loadings are stored as 'species' in the RDA object. We’ll extract the SNP loadings from the first significant constrained axis:
```
load.rda <- scores(pRDAfull, choices=c(1), display="species") 
hist(load.rda[,1], main="Loadings on RDA1")
```
The histogram of the loadings on each RDA axis should show relatively normal distributions. 
SNPs loading at the centre of the distribution have no relationship with the environmental predictors; those loading in the tails are more likely to be under selection as a function of those predictors (or some other predictor correlated with them).
Let's run the function (called outliers, where x is the vector of loadings and z is the number of standard deviations to use) provided by the tutorial now:
```
outliers <- function(x,z){
  lims <- mean(x) + c(-1, 1) * z * sd(x)     # find loadings +/-z sd from mean loading     
  x[x < lims[1] | x > lims[2]]               # locus names in these tails
}
```
Now apply the function to axis 1. We’ll use a 3 standard deviation cutoff (two-tailed p-value = 0.0027), but you can also change this to other values (e.g., 2.5 corresponds to a p-value of 0.012).
```
cand1 <- outliers(load.rda[,1],3)   
length(cand1) 
cand1 <- cbind.data.frame(rep(1,times=length(cand1)), names(cand1), unname(cand1))
colnames(cand1) <- c("axis","snp","loading")
write.table(cand1, "RDAoutliers.txt")
```
If you change the sd value, you should find you get 24 candidates on axis 1 at sd=3; 98 at sd=2.5; and 843 at sd=1.5.
Export the values at sd=3 and keep this list for later comparison to other methods.

### LFMM
LFMM tests for association between loci and environmental variables, taking population structure into account by introducing hidden latent factors into the model and identifying non-random associations between SNPs and environmental variables. 
As a first step, the optimal number of genetic clusters in the dataset should be determined to set the number of latent factors in the model. 
We used the sNMF function on Day Two to show that K=3 was optimal, so we'll use that today.
Load required packages:
``` 
library(qvalue)   
library(lfmm)
```
Read in the population allele frequency data from above:
```
gen <- read.table("Qfly_alleleFreq.txt")[,-1]
```
Read in the environmental data, and subset to the five bio variables of interest:
```
All_env <- read.delim("Qfly_AllBioclim_Final.txt")
env_subset <- All_env[,c("bio_4","bio_5","bio_9","bio_17","bio_18")]
```
Run the lfmm:
```
Qfly.lfmm <- lfmm_ridge(Y = gen, X = env_subset, K = 3)
```
Perform association testing using the fitted model:
```
pv <- lfmm_test(Y = gen, X = env_subset, lfmm = Qfly.lfmm, calibrate="gif")
```
Plot QQ-plot:
```
pvalues <- pv$calibrated.pvalue 
qqplot(rexp(length(pvalues), rate = log(10)),
       -log10(pvalues), xlab = "Expected quantile",
       pch = 19, cex = .4)
abline(0,1)
```
QQ-plot looks good! Also check genomic inflation factor (GIF) to get a sense of how well the model has accounted for confounding factors in the data:
```
pv$gif
```
An appropriately calibrated set of tests will have a GIF of around 1, so these values look good.
Next, we can check to see how applying the GIF to the pvalues changes the pvalue distribution for each env variable. 
In each case, we should see GIF-adjusted histogram looking generally flatter, with a peak for smaller pvalues (at the left hand side of the plot, near zero):
```
hist(pv$pvalue[,1], main="Unadjusted p-values bio_4")        
hist(pv$calibrated.pvalue[,1], main="GIF-adjusted p-values bio_4") 
hist(pv$pvalue[,2], main="Unadjusted p-values bio_5")        
hist(pv$calibrated.pvalue[,2], main="GIF-adjusted p-values bio_5") 
hist(pv$pvalue[,3], main="Unadjusted p-values bio_9")        
hist(pv$calibrated.pvalue[,3], main="GIF-adjusted p-values bio_9") 
hist(pv$pvalue[,4], main="Unadjusted p-values bio_17")        
hist(pv$calibrated.pvalue[,4], main="GIF-adjusted p-values bio_17") 
hist(pv$pvalue[,5], main="Unadjusted p-values bio_18")        
hist(pv$calibrated.pvalue[,5], main="GIF-adjusted p-values bio_18") 
```
Let's now estimate adjusted p-values:
```
pvalues <- pv$pvalue 
zs <- pv$score  
summary(zs) 
```
The summary function shows the z scores for each of the env variables that were used in the analysis.

Finally, we need to convert the adjusted p-values to q-values, which provide a measure of each SNP’s significance, while taking into account the fact that thousands are simultaneously being tested. We'll then use an FDR threshold to control the number of false positive detections (given that our p-value distribution is “well-behaved”).
This code will then print out the candidate loci for each bio variable:
```
scaffolds <- row.names(pvalues)
calibratedpvalues <- as.data.frame(pv$calibrated.pvalue)
qv_bio_4 <- as.data.frame(qvalue(calibratedpvalues$bio_4)$qvalues)
candidates_bio_4 <- scaffolds[which(qv_bio_4 < 0.05)]
candidates_bio_4
qv_bio_5 <- as.data.frame(qvalue(calibratedpvalues$bio_5)$qvalues)
candidates_bio_5 <- scaffolds[which(qv_bio_5 < 0.05)]
candidates_bio_5
qv_bio_9 <- as.data.frame(qvalue(calibratedpvalues$bio_9)$qvalues)
candidates_bio_9 <- scaffolds[which(qv_bio_9 < 0.05)]
candidates_bio_9
qv_bio_17 <- as.data.frame(qvalue(calibratedpvalues$bio_17)$qvalues)
candidates_bio_17 <- scaffolds[which(qv_bio_17 < 0.05)]
candidates_bio_17
qv_bio_18 <- as.data.frame(qvalue(calibratedpvalues$bio_18)$qvalues)
candidates_bio_18 <- scaffolds[which(qv_bio_18 < 0.05)]
candidates_bio_18
```
You should find that the number of outlier SNPs (i.e., those that show a signficant association with environmental variables) is 2, 3, 8, 4, and 3 for bio_4, 5, 9, 17, and 18, respectively.
Let's combine all outliers into a list and print to file:
```
allsnps = list(bio_4 = candidates_bio_4, bio_5 = candidates_bio_5, bio_9 = candidates_bio_9, bio_17 = candidates_bio_17, bio_18 = candidates_bio_18)
capture.output(allsnps, file = "LFMMoutliers.txt")
```
Overall using K=3, the default GIF correction, and an FDR threshold of 0.05, we have detected 20 outliers under putative selection in response to 5 bioclimatic variables. 
17 of these SNPs are unique and can be considered as the final LFMM candidate SNP list.

To wrap up this section, let's plot those outliers:
```
library(ggplot2)
library(dplyr)
```
Read in the chromosome and SNP position for the 6,526 SNPs that we retained after our MAF adjustment above:
```
scaffolds2 <- read.table("LFMM_allSNPs_CHRinfo6526.txt", header = T)
```
Let's now plot the candidates for bio_4. 
Note that you will need to read in the provided LFMMoutliers_bio4_final.txt file, where I have converted the original outlier list to one that shows chromosome information for plotting.
```
pvals_df <- as.data.frame(pvalues)
pval_bio_4 <- as.data.frame(-log10(pvals_df$bio_4))
pval_bio_4 <- cbind(pval_bio_4,scaffolds2)
pval_bio_4 <- rename(pval_bio_4, log10pval=`-log10(pvals_df$bio_4)`)
outliers_bio_4 <- read.table("LFMMoutliers_bio4_final.txt")
outliers_bio_4_vec <- as.vector(outliers_bio_4$V2)

selected_snps_bio4 <- pval_bio_4[pval_bio_4$SNP %in% outliers_bio_4_vec, ] #selecting candidates SNPs under bio_4

ggplot(pval_bio_4, aes(x=MRK, y = abs(log10pval), color = CHR)) + ##abs <- absolute value of plog
  geom_point(show.legend = FALSE, alpha = 0.8, size = 3) +
  scale_color_gradient(low="#D19C1D", high="#472C1B") +
  geom_point(selected_snps_bio4, mapping = aes(x=MRK, y = abs(log10pval)), color = "#EB4511", size = 4) +
  theme_minimal() +
  ggtitle("bio_4")
```
Now, adapting the above code, plot the outlier SNPs for the other four bioclim variables. 
For example, you'll change the first line to read: 'pval_bio_5 <- as.data.frame(-log10(pvals_df$bio_5))' and make similar adjustments all the way down to 'ggtitle("bio_5")'.
The required LFMMoutliers_bio5_final.txt file (and those for the other bio variables) are provided in the Day3/RequiredFiles/ folder.

### Gradient Forest
Finally, we'll explore the gradient forest method.
First, load in the required packages.
```
install.packages("gradientForest", repos="http://R-Forge.R-project.org")
install.packages("extendedForest", repos="http://R-Forge.R-project.org")
library(gradientForest)
library(extendedForest)
library(tidyr)
library(sp)
require(raster)
library(tidyr)
library(sp)
library(vegan)
```
Load in the environmental data from above, retaining only the data for each population and then extracting out PCNM spatial variables.
Note that PCNMs are principal coordinates of neighbor matrices - a bit more sophisticated than just using lat/long data.
```
clim.layer <- All_env
clim.points <- env_data
places <- read.delim("Qfly_coordinates.txt", header = TRUE)
pcnm <- pcnm(dist(places)) 
keep <- round(length(which(pcnm$value > 0))/2)
pcnm.keep <- scores(pcnm)[,1:keep] 
pcnm.keep
```
Now, run the gf analysis to model the associations of spatial and climate variables with allele frequencies (genotypes) of individuals.
First, we need to create an env.fg object that includes the climate and PCNM spatial variables. 
```
env.gf <- cbind(clim.points[ , c("bio_4", "bio_5", "bio_9", "bio_17", "bio_18") ], pcnm.keep)
```
A maximum number of 'splits' can be used to evaluate the gf model, as per the developers suggestion:
```
maxLevel <- log2(0.368*nrow(env.gf)/2)
```
And now run the model. The input is the combined climate, PCNM, and SNP data.
The other parts of the command define the predictor and response variables, and other parameters that are set to developer recommendations.
You can ignore any errors that might show up in red.
```
gf <- gradientForest(cbind(env.gf, AllFreq), predictor.vars=colnames(env.gf), response.vars=colnames(AllFreq), ntree=500, maxLevel=maxLevel, trace=T, corr.threshold=0.50)
```
Let's plot bar graphs depicting the importance of each spatial and climate variable:
```
plot(gf, plot.type = "O")
```
Identify the important variables:
```
gf
```
You should find: PCNM1  bio_18 PCNM2  bio_4  bio_17.
How do these results compare to the Qfly paper (Figure 4)?  Note that we used a few different bio variables to that paper, but you can compare the R2 weighted importance results for
bio5, bio9 and PCNM1 and PCNM2. Our results highlight the importance of the spatial (PCNM1 and 2), precipitation  (bio_17, bio_18), and temperature (bio_4) variables.

Let's now determine which variables are most important in the gf model by plotting the 'turnover functions' that show how allelic composition changes along the spatial or environmental gradients.
These plots are nonlinear and large jumps show steep genetic changes along certain portions of the environmental gradient. The height that the function acheives on the right side of the plot is the total importance and should match the R2 weighted importance barplot. First, organise the variables by importance and then plot:
```
by.importance <- names(importance(gf))
plot(gf, plot.type = "C", imp.vars = by.importance, show.species = F, common.scale = T, cex.axis = 1, cex.lab = 1.2, line.ylab = 1, par.args = list(mgp = c(1.5, 0.5, 0), mar = c(2.5, 2, 2, 2), omi = c(0.2, 0.3, 0.2, 0.4)))
```
What patterns do you see? Note areas where genetic variation changes abruptly and/or acheives high cumulative importance.
We can also make plots of turnover functions for individual loci. In this case, each line within each panel represents allelic change at a single SNP. Notice that in each panel some SNPs show very steep changes 
along the environmental gradient. These SNPs might be especially good candidates for local adaptation along the gradient and are highlighted in the legend. 
```
plot(gf, plot.type = "C", imp.vars = by.importance, show.overall = F, legend = T, leg.posn = "topleft", leg.nspecies = 5, cex.lab = 0.7, cex.legend = 0.8, cex.axis = 0.6, ylim = c(0, 0.5), line.ylab = 0.9, par.args = list(mgp = c(1.5, 0.5, 0), mar = c(2.5, 1, 0.1, 0.5), omi = c(0, 0.3, 0, 0)))
```
### Final steps
Have a look at your outlier lists for each of the three methods, i.e., compare files:
1. RDAoutliers.txt
2. LFMMoutliers.txt
3. Look at the legends in the last plots generated for the GF analysis to see some outliers there.
Are any of the outliers common across the different methods? 

### Bonus step
You can also rerun the model on only the adaptive SNPs

### Resources from today
https://popgen.nescent.org/2018-03-27_RDA_GEA.html

https://github.com/Elahep/LFMM

https://github.com/pgugger/LandscapeGenomics/blob/master/2019/Exercise4.md
