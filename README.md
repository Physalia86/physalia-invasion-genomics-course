# Physalia Invasion Genomics Course
## Invasion genomics in (mostly) R

https://www.physalia-courses.org/courses-workshops/invasion-genomics/

### Dr Ang McGaughran
### 8th, 10th and 11th September 2025

### COURSE OVERVIEW
This course on invasion genomics offers a comprehensive three-day program designed to introduce participants to key concepts and techniques in the field. Through a combination of recorded lectures and self-paced tutorials, attendees will delve into three main topics, including movement dynamics, adaptation dynamics, and environment dynamics relevant to invasive species. On Day One, participants will learn a variety of methods for understanding population structure and introduction pathways. Day Two will include an introduction to selection outlier scanning methods. Day Three will focus on the role of environmental variables in determining adaptive outcomes.

Participants will gain fundamental knowledge (video presentations) and practical experience (self-paced code tutorials) that advances their research in this rapidly evolving field. Participants will also be able to interact with the instructors regarding any issues encountered or discussion points raised.

### TARGET AUDIENCE
The course is aimed at graduate students, researchers, and professionals interested in invasion genomics. Participants should have a basic understanding of genetics/genomics and evolutionary biology, as well as familiarity with R. To ensure you can follow the course effectively, we recommend that participants without prior experience in R complete this R tutorial (https://r-coder.com/r-introduction/r-basics/) before attending.

### BACKGROUND READING
•    McGaughran A, Liggins L, Marske KA, Dawson MD, Schiebelhut LM, Lavery S, Knowles L, Moritz C, Riginos C (2022). Comparative phylogeography in the genomic age: opportunities and challenges. Journal of Biogegoraphy, 49, 2130-2144. https://doi.org/10.1111/jbi.14481.

•    North HL, McGaughran A, Jiggins C (2021). Insights into invasive species from whole-genome resequencing. Molecular Ecology, 30, 6289-6308. https://doi.org/10.1111/mec.15999.

### SESSION CONTENT
#### Monday - 8th September (DayOne)
Movement dynamics: inferring invasion sources, dispersal pathways, and population structure.
Recorded lecture: Introduction to invasion genomics and the key kinds of research questions and methods we use to address them. We’ll also cover analytical challenges associated with invasive species projects.
Tutorial 1: Analysis of population structure: FST, PCA, sNMF, population assignment.

#### Wednesday - 10th September (DayTwo)
Adaptation dynamics: inferring selection processes in the invaded range.
Recorded lecture: Selection analyses and how we can use them to infer adaptive process in invasive species/populations.
Tutorial 2: Analysis using outlier detection methods: PCAdapt, BayPass

#### Thursday - 11th September (DayThree)
Environment dynamics: inferring adaptation to local conditions.
Recorded lecture: Introduction to environmental association analysis and its potential use for understanding invasive species adaptation.
Tutorial 3: Analysis of three genotype-environment association (GEA) packages: RDA, LFMM, GF.

### COURSE PREPARATION
#### R
The code has been tested on version 4.4.2, so I would suggest that you have at least version 4.2 (and preferably version 4.3 or higher) installed prior to the course start date.
Note that R and RStudio are two different things: updating RStudio will be insufficient, you'll also need to update R.
To check what version of R you have installed, you can run:
```
version
```
The process of updating R will vary depending on your operating system, but you should be able to find instructions for this online.

We will make use of several R packages throughout the course that you'll need to have installed:
- adegenet
- caret
- corrplot
- dplyr
- extendedForest
- geodata
- ggplot2
- ggvenn
- gradientForest
- LEA
- lfmm
- ozmaps
- pcadapt
- qvalue
- raster
- remotes
- reshape2
- robust
- pophelper
- psych
- qvalue
- sp
- StAMPP
- tidyr
- UpSetR
- vegan
- vcfR

You will need to install these packages using a few different ways. Some are Bioconductor packages, which are installed with BiocManager:
```
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("LEA")
BiocManager::install("qvalue")
```
For the package remotes, you'll need to run:
```
library(remotes)
remotes::install_github('royfrancis/pophelper')
```
The rest should work with install.packages in R, for example:
```
install.packages("vcfR")
install.packages("adegenet")
install.packages("StAMPP")
```
Don't underestimate how long it may take to get all of these packages installed! My recommendation is that you install each package one by one, and then confirm the installation has worked by running:
```
library(packagename)
```
For example:
```
library(vcfR)
```
Once you have achieved this, the tutorials will be able to be run. During each tutorial, you will run the library(packagename) command to load the package as required.

#### BayPass
There is one section of the course that is run outside of R,  in the program BayPass. So, you will also need to install BayPass on your system as part of the course preparation.
You should follow instructions here for your system:
https://forge.inrae.fr/mathieu.gautier/baypass_public.
In my case, I am working on a linux machine and used the following code in a terminal window:
```
git clone https://forge.inrae.fr/mathieu.gautier/baypass_public.git
cd baypass_public/sources
make clean all FC=gfortran
```
You will then need to copy the g_baypass file (or your version of BayPass may be called ifx_baypass or i_baypass) to your working directory (from baypass_public/sources).
You can copy and paste the file across, or do it via the terminal/shell (replacing 'unname' with whatever your's is called):
```
cp /home/uname/baypass_public/sources/g_baypass ~/Desktop/Physalia/DayTwo/.
