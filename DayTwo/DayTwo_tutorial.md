### Day Two
Today, we'll be working again with the Qfly dataset. We'll use two key methods to identify potential outliers that may be under selection in the different populations:

1. PCAdapt
2. BayPass

Note that we will have to run some of the steps outside of R, as BayPass is a standalone program.
   
Let's get into it!

### Load in the required libraries
First, set your working directory to something sensible, like 'DayTwo', then load in the required libraries.
Note, you may have to install them first with either install.packages(x) or BiocManager::install("x")!
```
library("pcadapt")
library("ggplot2")
library(qvalue)
library(dplyr)
library(topGO)
```
Read in the VCF file:
```
Qfly_vcf <- "Qfly.vcf"
```
Run PCadapt:
```
Qfly_pcadapt <- read.pcadapt(Qfly_vcf, type = "vcf")
```
Ignore the error about converting from vcf to pcadapt - it does still work currently; in future you might need to first convert your VCF file into BED file (e.g., using PLINK software).
Now, run the pcadapt function and plot the resulting screen plot:
```
Qfly_pcadapt_kplot <- pcadapt(input = Qfly_pcadapt, K = 20)
plot(Qfly_pcadapt_kplot, option = "screeplot")
```
You're looking for the point after which the curve starts to flatten out, so a K-value of 3 is most appropriate here. Now, run a PCA using K=3:
```
Qfly_pcadapt_pca <- pcadapt(Qfly_pcadapt, K = 3)
summary(Qfly_pcadapt_pca)
```
Now check to see whether there is actually additional population structure on other PCs. First, read in the population data, then plot PCAs coloured by population for axes 1 vs 2 and axes 4 vs 5:
```
pop.data2 <- read.table("Qfly_popmap.txt", header = F)
print(pop.data2)
plot(Qfly_pcadapt_kplot, option = "scores", i = 1, j = 2, pop = pop.data2$V2)
plot(Qfly_pcadapt_kplot, option = "scores", i = 4, j = 5, pop = pop.data2$V2)
```
It looks like we capture the structure well with K=3.
Now, let's plot a manhattan plot and a qq-plot:
```
plot(Qfly_pcadapt_pca, option = "manhattan")
plot(Qfly_pcadapt_pca, option = "qqplot")
```
Let's examine the p-value frequency data and perform a Bonferroni correction for alpha=0.1:
```
Qfly_pcadapt_pvalues <- as.data.frame(Qfly_pcadapt_pca$pvalues)
hist(Qfly_pcadapt_pca$pvalues, xlab = "p-values", main = NULL, breaks = 50, col = "grey")

Qfly_pcadapt_padj <- p.adjust(Qfly_pcadapt_pca$pvalues,method="bonferroni")
alpha <- 0.1
outliers <- which(Qfly_pcadapt_padj < alpha)
length(outliers)
```
This should identify 18 outliers. View them with:
```
outliers
```
You can also print them to file:
```
write.table(outliers, file="Qfly_pcadapt_outliers.txt")
```
That's it for PCAdapt! We'll come back to the results later and see how they compare to BayPass.

### BayPass
Now, you need to move outside of R. The first step is to install BayPass. To do this, follow the instructions here for your system: https://forge.inrae.fr/mathieu.gautier/baypass_public
Here's an example that worked on my linux machine, typed in a terminal window:
```
git clone https://forge.inrae.fr/mathieu.gautier/baypass_public.git
cd baypass_public/sources
make clean all FC=gfortran
```
Now, copy the resulting g_baypass file from the baypass_public/sources directory to your working directory. Note, your version of BayPass  may be named differently (e.g., ifx_baypass or i_baypass).
Now you're ready to run BayPass. In a terminal window, type:
Note that I've provided the Qfly.bp (which I created using: python ./reshaper_baypass.py Qfly.vcf Qfly_popmap.txt Qfly.bp; the python script is available here: https://gitlab.com/YDorant/Toolbox/-/tree/master?ref_type=heads).
To run BayPass you also need to have specific constrast files set up. BayPass performs pairwise contrasts, so here is an example for how to run it using all native populations vs the Alice population. Basically, you identify this in BayPass by using a '1' to indicate the native populations, a '0' to indicate populations that should be ignored in the current analysis; and a '-1' to indicate the population for contrast (in this case, Alice). The ecotypes files you will need are all provided in the DayTwo folder.
```
./g_baypass -gfile Qfly.bp -contrastfile Qfly_NativevsAlice.ecotype -efile Qfly_NativevsAlice.ecotype -outprefix Qfly_NativevsAlice -nthreads 6
```
Note that you might need to change the number of threads to suit your computer (i.e., adjust -nthreads 6). This will take anywhere from ~5 mins to 30 mins or more, depending on the power of your computer.
This example will produce a series of 9 output files with the prefix 'Qfly_NativevsAlice'. Of these, we are most interested in the 'Qfly_NativevsAlice_summary_contrast.out' file, which we will now read into R.

#### Return to R
Read in the contrast file and rename column 4 to be a bit more useful:
```
Alice.C2 <- read.table("Qfly_NativevsAlice_summary_contrast.out",h=T) 
colnames(Alice.C2)[4]  <- "Log10pvals"
```
Now, import the provided (DayTwo folder) scaffold list, which tells us how the loci are located onto scaffolds. The was achieved by running some commands in bash on the Qfly.vcf file: cat Qfly.vcf | grep -v "#" | cut -f3 > Qfly_scaffold_list.txt
Import this file into R, and append the column to the Alice.C2 dataframe:
```
scaffolds <- read.table("Qfly_scaffold_list.txt")
Alice.C2 <- as.data.frame(cbind(Alice.C2, scaffolds))
colnames(Alice.C2)[5] <- "Scaffold"
```
Now, let's check the behaviour of the p-values associated with the C2 statistic:
```
hist(10**(-1*Alice.C2$Log10pvals),freq=F,breaks=40)
abline(h=1)
```
We need to convert these p-values from log10 to normal and convert to q-values:
```
pvalues <- as.data.frame(10^-Alice.C2[,4]) 
Alice.C2 <- cbind(Alice.C2,pvalues) 
colnames(Alice.C2)[6]  <- "pvals"
pvalues <- as.vector(Alice.C2$pvals)
qval <- qvalue(p = pvalues)
plot(qval$qvalues)
abline(h=0.05)

qvalues <- as.data.frame(qval$qvalues)
Alice.C2 <- cbind(Alice.C2,qvalues)
colnames(Alice.C2)[7]  <- "Qvals"
```
Note that q-values < 0.05 are highly significant. Let's extract the SNPs that have q-values < 0.05:
```
selected_SNPs <- Alice.C2[Alice.C2$Qvals < 0.05, ]
write.table(selected_SNPs,"BPoutliers_NativevsAlice_FDR5.txt", sep = "\t")
```
There should be 62 outliers. Let's plot them out:
```
ggplot(Alice.C2, aes(x=MRK, y = Log10pvals)) + 
  geom_point(show.legend = FALSE, alpha = 1, size = 2) +
  scale_color_gradient(low="#D19C1D", high="#472C1B") +
  geom_point(selected_SNPs, mapping = aes(x=MRK, y = Log10pvals), color = "#EB4511", size = 2) +
  theme_classic()
```
We can also replot these by chromosome, using the provided file, where scaffold names are updated to chromosome numbers (with the number 6 used for unplaced scaffolds):
```
scaffolds2 <- read.table("Qfly_chromosomes-chrNumbers.txt", header = TRUE) 
Alice.C2_2 <- cbind(Alice.C2, scaffolds2)
ggplot(Alice.C2_2, aes(x=MRK, y = Log10pvals, color = CHR)) + 
  geom_point(show.legend = FALSE, alpha = 1, size = 2) +
  scale_color_gradient(low="#D19C1D", high="#472C1B") +
  geom_point(selected_SNPs, mapping = aes(x=MRK, y = Log10pvals), color = "#EB4511", size = 2) +
  theme_classic()
```

### Next steps
Re-run the above steps for the other two BayPass comparisons. In BayPass (outside of R), run:
```
./g_bapyass -gfile Qfly.bp -contrastfile Qfly_NativevsExpanded.ecotype -efile Qly_NativevsExpanded.ecotype -outprefix Qfly_NativevsExpanded -nthreads 6
./g_bapyass -gfile Qfly.bp -contrastfile Qfly_NativevsIslands.ecotype -efile Qly_NativevsIslands.ecotype -outprefix Qfly_NativevsIslands -nthreads 6
```
Then follow the rest of the code above to generate outlier lists for these two additional comparisons, creating these two files:
1. Qfly_NativevsExpanded_summary_contrast.out
2. Qfly_NativevsIslands_summary_contrast.out
   
Then, back in R, repeat the rest of the code to generate files:
1. BPoutliers_NativevsExpanded_FDR5.txt
2. BPoutliers_NativevsIslands_FDR5.txt
Now, make a venn diagram to check for common SNPs between the three BayPass comparisons:
```
Native <- read.table("BPoutliers_NativevsAlice_FDR5.txt")
Expanded <- read.table("BPoutliers_NativevsExpanded_FDR5.txt")
Islands <- read.table("BPoutliers_NativevsIslands_FDR5.txt")

venn1 <- list(Native=Native$V2, Expanded=Expanded$V2, Islands=Islands$V2)
ggvenn(venn1, fill_color = c("#F8D210", "#593F1F", "#DEE3CA"), fill_alpha = 0.7, stroke_size = 0.2, set_name_size = 4, stroke_color = "black", show_percentage = FALSE)
```
Let's also check to see how the BayPass result compares with PCAdapt. We already generated the table of PCAdapt outliers above, so:
```
PCadapt <- read.table("Qfly_pcadapt_outliers_a10.txt")
venn2 <- list(Native=Native$V2, Expanded=Expanded$V2, Islands=Islands$V2, PCAdapt=PCadapt$x)
ggvenn(venn2, fill_color = c("#F8D210", "#593F1F", "#DEE3CA", "grey"), fill_alpha = 0.7, stroke_size = 0.2, set_name_size = 4, stroke_color = "black", show_percentage = FALSE)
```
Your results should look like those provided in the DayTwo folder (BayPassVenn.jpeg and FullVenn.jpeg).

### Outlier annotation
Functional annotation of candidate outlier genes is a big task and we don't have time to complete the full pipeline today, so we will skip straight to the GO analysis (gene ontology) analysis step in R, where we'll use the topGO package. However, to give you an overview, here's an outline of the steps that we're not doing today:
1. Extract the outlier SNPs from the VCF file, using VCFTools - this would generate a new file (NativevsAlice_outliers.vcf) so we can annotate the outliers that define Native vs Invasive Alice population: vcftools --vcf Qfly.vcf --snps NativevsAlice_candidates.txt --recode --recode-INFO-all; then mv out.recode.vcf NativevsAlice_outliers.vcf
2. Convert vcf file to bed file using BEDOPS: vcf2bed < NativevsAlice_outliers.vcf > NativevsAlice_outliers.bed
3. Generate a fai.fa file with samtools using the fasta genome file for Qfly: samtools faidx Qfly.fa
4. Create a bedtools genome format file: awk -v OFS='\t' {'print $1,$2'} Qfly.fa.fai > Qfly.txt
5. Extract flanking regions (10kb downstream and upstream of the outlier SNP) using BEDTOOLS: bedtools flank -i NativevsAlice_outliers.bed -g Qfly.txt -b 10000 | cut -f1-3 > NativevsAlice_outliers_10kb_flank.bed
6. Upload NativevsAlice_outliers_10kb_flank.bed file to UCSC Table Browswer (https://genome.ucsc.edu/) to get transcript IDs associated with these sequences
7. Peform Go terms analysis to assess the gene functions of the putatively adaptive loci: (a) Get the transcript ID for all loci using snpEff:
   a. snpEff -c snpEff.config mygenome Qfly.vcf > Qfly_FullData.anno.vcf
   b. cat snpEff_genes.txt | grep -v "#" | cut -f3 | uniq > Qfly_FullData_transIDfromsnpeff.txt
8. Run in Interpro scan:
   a. interproscan.sh -i Qfly_FullData_batchentrez.fasta -t n --goterms
9. Remove any transcripts that lack GOterm IDs:
   a. grep -w "GO" Qfly_FullData_batchentrez.fasta.tsv | cut -f1,14 > Qfly_FullData_GOlist.txt

Now, we're ready to move into R to run the GO terms analysis. For more information on the above steps, see the full tutorial here: https://github.com/Elahep/B.tryoni_PopGenomics/tree/main/4-GOterms

#### Back in R
library(topGO)

First, import the tab delimited file of all genes and their GO terms:
```
geneID2GO <- readMappings("Qfly_FullData_GOlist.txt")  
geneID2GO$XM_040099783.1   #check the GO terms for some of the transcript IDs
str(head(geneID2GO))
geneNames <- names(geneID2GO)
head(geneNames)
```
Next, import transcript IDs for the outlier SNPs:
```
interesting_genes = read.table("NativevsAlice_outliers_transcriptIDs.txt")
myInterestingGenes <- as.vector(interesting_genes$V1)
geneList <- factor(as.integer(geneNames %in% myInterestingGenes))
names(geneList) <- geneNames
str(geneList)
```
Now, create a topGo object for the "Biological Processes" ontology terms (BP), and get the list of significant genes:
```
GOdata_BP <- new("topGOdata", ontology = "BP", allGenes = geneList, annot = annFUN.gene2GO, gene2GO = geneID2GO) 
sig_genes = sigGenes(GOdata_BP) 
```
Generate the GO graph:
```
resultFisher <- runTest(GOdata_BP, algorithm="weight01", statistic="fisher")
```
See how many GO terms are significant and write the restuls to file: 
```
resultFisher  
allRes <- GenTable(GOdata_BP, raw.p.value = resultFisher, classicFisher = resultFisher, ranksOf = "classicFisher", Fis = resultFisher, topNodes = length(resultFisher@score)) 
allRes
write.table(allRes, "NativevsAlice_GOresults_BP.txt", sep = "\t") 
```
Now, repeat using the "Molecular Function" (MF) ontology:
```
GOdata_MF <- new("topGOdata", ontology = "MF", allGenes = geneList, annot = annFUN.gene2GO, gene2GO = geneID2GO)
GOdata_MF
resultFisher_MF <- runTest(GOdata_MF, algorithm="weight01", statistic="fisher")
resultFisher_MF  #this shows how many GO terms are significant
allRes_MF <- GenTable(GOdata_MF, raw.p.value = resultFisher_MF, classicFisher = resultFisher_MF, ranksOf = "classicFisher", Fis = resultFisher_MF, topNodes = length(resultFisher_MF@score)) 
allRes_MF
write.table(allRes_MF, "NativevsAlice_GOresults_MF.txt", sep = "\t")
```
And the "Cellular Component" terms (CC):
```
GOdata_CC <- new("topGOdata", ontology = "CC", allGenes = geneList, annot = annFUN.gene2GO, gene2GO = geneID2GO)
GOdata_CC
resultFisher_CC <- runTest(GOdata_CC, algorithm="weight01", statistic="fisher")
resultFisher_CC  #this shows how many GO terms are significant
allRes_CC <- GenTable(GOdata_CC, raw.p.value = resultFisher_CC, classicFisher = resultFisher_CC, ranksOf = "classicFisher", Fis = resultFisher_CC, topNodes = length(resultFisher_CC@score)) 
allRes_CC
write.table(allRes_CC, "NativevsAlice_GOresults_CC.txt", sep = "\t")
```
Explore these tables to highlight xxx














































### Resources for today:
https://bcm-uga.github.io/pcadapt/articles/pcadapt.html
