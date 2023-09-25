if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager") BiocManager::install("affy")
library(affy) #Loads the affy package.
setwd("C:/Users/HP/OneDrive/Documents/Su_CELs/")
mydata = ReadAffy()
#Reads all *.CEL (*.cel) files in your current working directory and stores them into the AffyBatch object 'mydata'.
#BiocManager::install("tkWidgets")
library(tkWidgets)
mydata = ReadAffy(widget=TRUE) #Opens file browser to select specific CEL files.
eset = rma(mydata)
eset = mas5(mydata)
#BiocManager::install("gcrma")
library(gcrma)  
eset_gcrma = gcrma(mydata)
#Use this command to acquire gcrma data. The 'library(gcrma)' needs to be loaded first.
#BiocManager::install("plier")
library(plier)
eset_plier = justPlier(mydata)
#Use this command to generate plier data. The 'library(plier)' needs to be loaded first.
eset_PMA = mas5calls(mydata)
#Generates MAS 5.0 P/M/A calls. The command creates ExpressionSet with P/M/A calls in the 'exprs' slot and the wilcoxon p-values in the 'se.exprs' slot. To access them see below.
eset = expresso(mydata, normalize.method="invariantset", bg.correct=FALSE, pmcorrect.method="pmonly", summary.method="liwong")
# Generates expression calls similar to dChip (MBEI) method from Li and Wong.
#BiocManager::install("affycoretools")
library(affycoretools)
affystart(plot=T, express="mas5")
write.exprs(eset, file="mydata.txt")
#Writes expression values to text file in working directory.
x = data.frame(exprs(eset), exprs(eset_PMA), assayDataElement(eset_PMA, "se.exprs"))
x = x[,sort(names(x))];
write.table(x, file="mydata_PMA.xls", quote=F, col.names = NA, sep="\t") #Writes expression values, PMA values and wilcoxon p-values to text file in working directory. Remember: the old command for accessing the wilcoxon pvalues in BioC versions.
mypm = pm(mydata)
#Retrieves PM intensity values for single probes
mymm = mm(mydata)
#Retrieves MM intensity values for single probes
myaffyids = probeNames(mydata) #Retrieves Affy IDs
result = data.frame(myaffyids, mypm, mymm) #Combines the above information in data frame
eset
pData(eset)
#Provides summary information of ExpressionSet object 'eset' and lists the analyzed file names.
exprs(eset)[1:2,1:4]
exprs(eset)[c("244901_at","244902_at"),1:4]
#Retrieves specific rows and fields of ExpressionSet object. To learn more about this format class, consult ExpressionSet manual with command '?ExpressionSet'
test = as.data.frame(exprs(eset));
eset2 = new("ExpressionSet", exprs = as.matrix(test), annotation="ath1121501");
eset2
#Example for creating an ExpressionSet object from a data frame. To create the object from an external file, use the read.delim() function first and then convert it accordingly.
data.frame(eset) #Prints content of 'eset' as data frame to STDOUT.
exprs(eset_PMA)[1:2,1:2];
assayDataElement(eset_PMA, "se.exprs")[1:2,1:2]
#Returns from ExpressionSet of 'mas5calls(mydata)' the PMA values from its 'exprs' slot and the p-values from its
# 'se.exprs' slot.
#BiocManager::install("ath1121501.db")
library(ath1121501.db)
#Opens library with annotation data.
library(help=ath1121501.db)
#Shows availability and syntax for annotation data.
ath1121501()
#Provides a summary about the available annotation data sets of an annotation library.
#BiocManager::install("ath1121501cdf")
library(ath1121501cdf);
ls(ath1121501cdf)
#Retrieves all Affy IDs for a chip in vector format.
x = c("245265_at", "260744_at", "259561_at", "254759_at", "267181_at")
#Generates sample data set of Affy ID numbers.
contents(ath1121501ACCNUM)[1:40] #Retrieves locus ID numbers for Affy IDs.
myAnnot = data.frame(ACCNUM=sapply(contents(ath1121501ACCNUM), paste, collapse=", "), SYMBOL=sapply(contents(ath1121501SYMBOL), paste, collapse=", "), DESC=sapply(contents(ath1121501GENENAME), paste, collapse=", ")) 
myAnnot[x, ]
#Organizes full set of annotation features in a data frame, here: ACCNUM, SYMBOL and GENENAME. The annotations from probe sets with several mappings are collapsed into single strings.
mget(x, ath1121501GO, ifnotfound=NA) # Retrieves GO information for Affy IDs.
#BiocManager::install("ath1121501probe")
library(ath1121501probe) #Opens library with probe sequence data.
print.data.frame(ath1121501probe[1:22,])
#Prints probe sequences and their positions for first two Affy IDs.
library(Biostrings)
pm = DNAStringSet(ath1121501probe$sequence)
names(pm) = ath1121501probe$Probe.Set.Name
# Stores probe sequences as DNAStringSet object. See HT-Seq manual for more details.
#Generates a comprehensive QC report for the AffyBatch object 'mydata' in PDF format. See affyQCReport for details.
deg = AffyRNAdeg(mydata); summaryAffyRNAdeg(deg); plotAffyRNAdeg(deg)
#Performs RNA degradation analysis. It averages on each chip the probes relative to the 5'/3' position on the target genes. A summary list and a plot are returned.
image(mydata[ ,1])
#Reconstructs image with log intensities of first chip.
hist(mydata[ ,1:2])
#Plots histogram of PM intensities for 1st and 2nd array.
hist(log2(pm(mydata[,1])), breaks=100, col="blue")
#Plots bar histogram of the PM ('pm') or MM ('mm') log intensities of 1st array.
boxplot(mydata,col="red")
#Generates a box plot of un-normalized log intensity values.
boxplot(data.frame(exprs(eset)), col="blue", main="Normalized Data") #Generates a box plot of normalized log intensity values.
mva.pairs(pm(mydata)[,c(1,4)])
#Creates MA-plot for un-normalized data. A MA-plot is a plot of log-intensity ratios (M- values) versus log-intensity averages (A-values) between selected chips (here '[1,4]').
mva.pairs(exprs(eset)[,c(1,4)]) #Creates MA-plot for normalized data.
	  		
		  		
		  		
		  		
		  		
		  		
		  		