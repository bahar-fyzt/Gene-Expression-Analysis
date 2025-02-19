#### Project reference: 

Zambelli et al 2022, Picardi et al 2015 and D'Erchia et al 2014.

#### Project description:

We are going to work on a human gene expression data, profiling RNA-seq dataset samples from 6 human organs/tissues. These are obtained from 3 donors: S7, S12, S13.

We have two files:

-Counts2025.csv: a table containing gene expression values (read counts) for the 28.391 human genes in the 18 replicates.

-Metadata2025.csv: a metadata (design) table containing the experimental design. It has the following columns: Tissue (name of the tissue), Donor (name of the sample/donor).

We previously installed the package edgeR available through Bioconductor in Rstudio.

```{r}
library('edgeR')
```


*edgeR* is an R package used for analyzing gene expression data, particularly RNA-Seq count data. It helps identifying genes that are differentially expressed between experimental conditions.

In brief edgeR can: compute normalization and filter data, test for differential expression of genes, create meaningful plots etc..


We are going to perform the following steps:


### 1. Functions

#### 1.1 Check that names match and are aligned between the two tables

#### 1.2 Filter genes

#### 1.3 Compute normalization: ComputeCPM()

#### 1.4 Compute Average: ComputeAvg()
 
### 2. Data import

### 3. Common part

### 4. Second part: 

# 1. Functions 
We declare all the functions that we are going to use in this gene expression analysis.

## 1.1 Check that names match and are aligned
This function takes as input a count table and a design table. Then checks that the column names of the counts table are consistent with the row names of the design table.

```{r}
checkNames<-function(cTable,dTable)
{
  check<-all(colnames(cTable)==rownames(dTable))   # Use all() to check if all the names are the same and in the same order.
  return(check)
}
```

## 1.2 Filter genes
This function checks, for every gene, if that gene is expressed above a threshold and for at least a number of times. The input is a counts dataframe. 
ThrCounts is used to set the minimum counts threshold (we want ThrCounts as 10).
ThrTimes is used to set the minimum number of times threshold (we want ThrTimes as 1).

```{r}
FiltGenes<-function(cTable,ThrCounts=10,ThrTimes=1)
{
  ExpressedGenes<-cTable>=ThrCounts     # Check if counts are >= than a threshold, the result will be a matrix of logical values.
  HowManyTimes<-apply(ExpressedGenes,1,sum)   # Number of times expressed is equal to the number of TRUE values. 
  FiltCounts<-HowManyTimes>=ThrTimes          # Check that number of times is >= thrTimes.
  return(FiltCounts)                    # The result will be a vector of logical values with TRUE or FALSE for every gene.
}
```
## 1.3 Compute normalization: ComputeCPM()
This function takes a counts table in input and for every column:

1. It computes the total number of million reads

2. It divides that column for that number
```{r}
ComputeCPM<-function(cTable)
{
  for (i in 1:ncol(cTable))       # For every column.
  { 
    scale<-sum(cTable[,i])/1e6    # Compute the scale factor.
    cTable[,i]<-cTable[,i]/scale  # Divide that column by the scale factor.
    }
  return(cTable)
}
```
## 1.4 Compute Average: ComputeAvg()
This function computes the average of gene expression across all the replicates of each tissue. 
```{r}
computeAvg<-function(cTable,conds)
{
  uniqueConditions<-unique(conds)   # Identify unique conditions.
  avgTable<-matrix(ncol=length(uniqueConditions), 
                   nrow=nrow(cTable))    # Create an empty matrix to store results.
  rownames(avgTable)<-rownames(cTable)   # Set row names according to the content (genes).
  colnames(avgTable)<-uniqueConditions   # Set col names according to the content (structures/tissues).
  for (condition in uniqueConditions) 
  {
    selection<-conds==condition   # Select the condition.
    avgTable[,condition]<-apply(cTable[,selection],1,mean)   # Compute average on the rows (1st margin) using the apply().
  }
  return(avgTable)
}
```

# 2. Data import
We load the data. We use the *read.table()* function. 
We load the counts table: a table containing gene expression values (read counts) for the 28.391 human genes in the 18 replicates.
We load the data and we check the dimension of our data using *dim()*. 
```{r}
counts<-read.table('Counts2025.csv',
                   sep='\t',      # In our file values are delineated by tabulations.
                   header=TRUE,   # The file has a header (first line).
                   row.names=1)   # The file has a column (1st) with row names.
dim(counts)   # Check that we have data.
head(counts)  # Show the first 6 rows of the counts table.
```

We load the metadata: a table containing the experimental design. It has the following columns: Tissue (name of the tissue), Donor (name of the sample/donor).
```{r}
design<-read.table('Metadata2025.csv',
                   sep='\t',      # Also this file is 'tab' delineated.
                   header=TRUE,   # The file has a header (first line).
                   row.names=1)   # The file has a column (1st) with row names.
dim(design) # Check that we have data.
head(design) # Show the first 6 rows of the design table.
```
# 3. Common part
We have to take into account that in the counts matrix sample names are indicated as: *<Tissue_Donor>*.
So there is a concatenation between Tissue and Donor_id (separator='_').

For consistency we will add a row.names to the metadata table according to the same syntax.
```{r}
rownames(design)<-paste(design$Tissue, design$Donor, sep='_')   # Use paste to concatenate Tissue and Donor and save the resulting vector as the rownames of design.
head(design)   # We check with head() if the desired rownames are in place.
```

Now, as we proceed, first and foremost we must check if samples have the same name and same order. We must have consistent sample names between the two tables (counts and design) in order to perform gene analysis.
```{r}
checkNames(counts,design)   # Check if samples have the same name and same order.
```
The output is FALSE, which indicates that either the order of the samples is different or that the samples themselves differ between the counts and design table. We check if at least we have the same samples independently from the order. 
```{r}
all(colnames(counts) %in% rownames(design))   # %in% is an operator used to check if elements of one vector are present in another vector. The output is a logical vector.
```
The output is TRUE. This means that we have the same samples independently from the order.
The counts and design table must contain the samples in the same order. We order the counts table according to the design table.
```{r}
counts<-counts[,rownames(design)] # Order counts based on design.
checkNames(counts,design)         # Now it is TRUE.
head(counts)                      # Check the head of counts after ordering.
```
Up to here we did data ingestion.

Create an object of type DGEList. 
The main components of a DGEList object are: a 'counts' table, a data.frame containing the integer counts, a data.frame called 'samples' containing information about the samples. The data.frame 'samples' must also contain a column 'group' so as to define the experimental groups or conditions for differential expression analysis.

```{r}
data<-DGEList(counts=counts,
              group=design$Tissue,
              sample=design)
head(data)

```
We have to keep only the genes that are likely to be expressed. We have to take genes that have more than 10 reads in at least 1 replicate. We use the previously declared function.

```{r}
filter<-FiltGenes(data$counts,ThrCounts=10,ThrTimes=1)   # Filter the genes that are likely to be expressed. The output is a logical vector.
FiltData<-data[filter,]   # These are the genes corresponding to each TRUE value.
head(FiltData)
```

Normalization. We use the function *calcNormFactors()*.
Highly expressed and highly variable genes are excluded from the computation.
The default method for computing these scale factors uses a trimmed mean of m-values (TMM) between each pair of samples.
In the TMM normalization, one (any) of our sample is set as a baseline and very highly expressed genes are excluded. Then, ratio of expression of baseline w.r.t. other samples is computed. The genes that change the most are excluded. Then, the scale factor is determined, which is the median value of the ratio of expression of genes in baseline compared to another condition. Finally, in the TMM normalization process, the raw counts in each sample are multiplied by their respective scale factors. This adjustment ensures that differences in sequencing depth and RNA composition between samples are accounted for, enabling more accurate comparisons of gene expression across conditions.

```{r}
NormFiltData<-calcNormFactors(FiltData) 
head(NormFiltData)   #  We see that now the column called 'norm.factors' contains the corresponding values. Before, it was composed by all 1.
```

Now, we perform a MDS (Multidimensional Scaling), similar to PCA (Principal Component Analysis), plot of the data.
The *plotMDS()* function can be used to perform PCA-like analysis of our data. This type of plot helps visualize how samples group, based on overall expression differences, aiding in understanding biological variability and patterns across tissues or conditions.
```{r}
colMDS<-rainbow(length(levels(as.factor(design$Tissue))))   # We set the colors of the different tissues.

plotMDS(NormFiltData, col=colMDS[as.factor(design$Tissue)], # Colors are selected based on the different tissues.
        pch=16, # Chosen symbol.
       cex=1.5, # Dimension of the chosen symbol.
       main='MDS plot') 
# We use the legend() function to add the legend.
legend(4.5,    # Set the position of the legend in the plot (x-axis). 
       -1,     # Set the position of the legend in the plot (y-axis). 
       fill=colMDS, # We assign colors to the legend items based on colMDS.
       legend=levels(as.factor(design$Tissue)))   # The names of the legend corresponding colors.
      
```

#### Interpretation of the MDS Plot:
The x-axis represents the first principal dimension of variation (Leading logFC dim 1), accounting for 32% of the total variability in the dataset. The y-axis represents the second principal dimension of variation (Leading logFC dim 2), explaining 20% of the variability.
Each labelled point corresponds to a sample, and the clustering of points indicates the similarity in gene expression profiles among the samples. For instance, samples belonging to the same tissue type (e.g., "brain", "lung", etc.) tend to cluster together, which suggests that gene expression is more similar within tissues than across tissues.
The separation along the axes reflects differences in gene expression, with larger distances between clusters indicating greater differences. 
The proximity of clusters (e.g., "lung" and "kidney") may suggest some overlap or shared patterns in their gene expression profiles.

Up to here we did the quality control.

#### DE analysis between two tissues

We select two different (and meaningful) tissues and we perform differential gene expression analysis between those two conditions.
The differential gene expression analysis is a three step process. 

Step1: We use the NB (negative binomial) model to estimate gene expression levels (both variability and average) with *estimateDisp()*. 

```{r}
NormFiltDataDisp<-estimateDisp(NormFiltData)   # Compute average expression and variability for every gene.
```
Step2: Compare two conditions. 
Perform a differential expression analysis between those 2 conditions.

We select 2 different and meaningful tissues. We select brain and muscle. Indeed each sample of the same tissue cluster together and they are farther apart from the samples of the other tissue, in term of logFC, and this suggests differential expression.
The groups are in NormFiltDataDisp$samples$group. 

```{r}
table(NormFiltDataDisp$samples$group)   # Need to provide to the test conditions to be compared.
```
To test for differential expression we use the *exactTest()* function. This function can be used only to perform pairwise tests.
The object that is returned is basically a list containing: 
a dataframe called 'table' which holds the results of the DE test, as well as logFC and expression values, a vector called 'comparison', which contains the names of the conditions that were compared, and an optional data.frame called 'genes' which contains the annotations of the genes.

```{r}
BrainVsMuscle<-exactTest(NormFiltDataDisp,
                  pair=c('brain','muscle'))   # The 'pair' parameter can be used to specify the 2 conditions that are to be compared (in the form of a vector).
BrainVsMuscle   # Print results to screen.
```
Step3: Apply FDR (False Discovery Rate) to correct for multiple tests. Indeed we have to take into account the fact that we tested thousands of genes and we do not want FALSE results (FDR). Differentially expressed genes are genes that have low FDR (corrected p-value).
We use the *topTags()* function. An object of class topTags will be returned. 
This one is another list which contains:
our table data.frame to which the additional FDR column has been added, a character string specifying the type of correction that was applied, called 'correction', another string specifying the name of the test that was applied and finally a vector with the names of the conditions that were compared.

The function adds the FDR as a column.
By default only the 10 most significant genes are reported. We need to report all the genes or at least all the genes that are DE.

```{r}
# Create top tags 
BrainVsMuscle_FDR<-topTags(BrainVsMuscle, n=nrow(NormFiltDataDisp$counts))   # We use 'n=nrow()' to force the topTags() function to include all the genes.
head(BrainVsMuscle_FDR)

BrainVsMuscle_FDR_table<-BrainVsMuscle_FDR$table   # The main results of TopTags are in $table, we make a copy of the table with the main results.
head(BrainVsMuscle_FDR_table)

DEGs<-rownames(BrainVsMuscle_FDR_table)[BrainVsMuscle_FDR_table$FDR<=0.01]   # We set a threshold of FDR<=0.01 because DE genes are genes that have low FDR. We obtain a list of DEGs.
head(DEGs)
```
Assign all the genes into one of the 4 possible classes as follow:
```{r}
DE_UPgenes<-rownames(BrainVsMuscle_FDR)[BrainVsMuscle_FDR_table$FDR<=0.01 & BrainVsMuscle_FDR_table$logFC>0] 
# We extract the genes that are DE_UP. Select genes with significant differential expression (FDR<=0.01) and log fold change (logFC>0), meaning they are up regulated in brain compared to muscle.
head(DE_UPgenes)

DE_DOWNgenes<-rownames(BrainVsMuscle_FDR)[BrainVsMuscle_FDR_table$FDR<=0.01 & BrainVsMuscle_FDR_table$logFC<0] 
# We extract the genes that are DE_DOWN. Select genes with significant differential expression (FDR<=0.01) and log fold change (logFC<0), meaning they are down regulated in brain compared to muscle.
head(DE_DOWNgenes)

notDE_UPgenes<-rownames(BrainVsMuscle_FDR)[BrainVsMuscle_FDR_table$FDR>0.01 & BrainVsMuscle_FDR_table$logFC>0] 
# We extract the genes that are notDE_UP. Select genes with not significant differential expression (FDR>0.01) and log fold change (logFC>0), meaning they are not significantly up regulated.
head(notDE_UPgenes) 

notDE_DOWNgenes<-rownames(BrainVsMuscle_FDR)[BrainVsMuscle_FDR_table$FDR>0.01 & BrainVsMuscle_FDR_table$logFC<0] 
# We extract the genes that are notDE_DOWN. Select genes with not significant differential expression (FDR>0.01) and log fold change (logFC<0), meaning they are not significantly down regulated.
head(notDE_DOWNgenes)
```

We want to create boxplot based on logFC:
```{r}
# Extract logFC in all the classes.
DE_UP<-BrainVsMuscle_FDR_table[DE_UPgenes, 'logFC'] 

DE_DOWN<-BrainVsMuscle_FDR_table[DE_DOWNgenes, 'logFC']

notDE_UP<-BrainVsMuscle_FDR_table[notDE_UPgenes, 'logFC']

notDE_DOWN<-BrainVsMuscle_FDR_table[notDE_DOWNgenes, 'logFC']

```

We make a boxplot of the logFC of the genes belonging to each class. 
We can use the *boxplot()* function. To make the plot more easy to visualize we use the 'outline' parameter and we set it to FALSE.
Indeed, when we have too many outliers it might be difficult to compare distributions. So, we set the outline parameter accordingly.
```{r}
boxplot(DE_DOWN, notDE_DOWN, DE_UP, notDE_UP,
        outline=FALSE, names=c("DE_DOWN","notDE_DOWN","DE_UP","notDE_UP"), 
        lwd=2,  # Thickness of the lines used to draw the boxplot elements.
        col=rainbow(4), 
        main='BrainVsMuscle DEGs and not DEGs', ylab='logFC')
```

#### Interpretation of the boxplot : 
The x-axis has the log FC distribution for 4 categories of genes (DE_DOWN, notDE_DOWN, DE_UP, notDE_UP) during a DGE analysis between brain and muscle.
The y-axis represents the log fold change of gene expression. 
LogFC measures the magnitude and direction of the expression difference between the two tissues (e.g. brain and muscle). Positive values indicate up regulation in one tissue and negative values indicate down regulation.

As we said before:

DE_DOWN: Differentially expressed genes with significant down regulation (lower expression in brain compared to muscle).

notDE_DOWN: Genes with logFC<0 but not statistically significant.

DE_UP: Differentially expressed genes with significant upregulation (higher expression in brain compared to muscle).

notDE_UP: Genes with logFC>0 but that are not statistical significant.

From the graph we can observe that:

DE_DOWN: The median logFC is negative and the range is spread into negative values. This confirms that these genes are significantly down regulated in brain compared to muscle.

notDE_DOWN: The median logFC is close to 0, with less variability compared to DE_DOWN. These genes show a slight tendency toward down regulation but are not significant.

DE_UP: The median logFC is positive, with a wide spread toward higher positive values. These genes are significantly up regulated in brain compared to muscle.

notDE_UP: The median logFC is around 0, similar to notDE_DOWN, with minimal spread. These genes show a slight trend toward up regulation but are not significant.

In conclusion we can say that the DE_DOWN and DE_UP categories represent the most biologically interesting genes. Indeed they are significantly differentially expressed between brain and muscle.
DE_UP may include brain-specific genes related to metabolic functions of brain.
DE_DOWN may represent muscle-specific genes involved in contraction for example.
The notDE_DOWN an notDE_UP genes do not show significant differences in expression and represent genes with common expression across tissues.

To sum up:

DE_UP have bigger fold change, they change more than notDE.

DE_DOWN have lower fold change, they change more than notDE.

So, the main difference between DE and notDE is that DE genes (DEGs) have a bigger change in expression (FC) compared to genes that are not DE.


# 4. Second Part: 
We have to compute the average expression across tissues.
RNAseq libraries have a different sequencing depth for every sample/replicate.
The first normalization that we perform aims to make the number of reads comparable (normalization per library size).
We scale the counts by the total number (per million) number of counts.
To do that we use the *computeCPM()* function.

```{r}
CPM<-ComputeCPM(NormFiltDataDisp$counts)   # Computing Counts Per Million normalization.
```

Then, we compute the average by condition using the previously written function.
```{r}
AvgExpression<-computeAvg(CPM, design$Tissue)
head(AvgExpression)
```
We compare the average expression across the 6 tissues and we assign the genes to either of the following classes: tissue specific genes, tissue elevated genes and not-elevated genes.

```{r}
# Create a list to store elevated and specific genes for each tissue.
tissue_specific <- list()
tissue_elevated <- list()
not_elevated <- c()

for (i in 1:nrow(AvgExpression)) {        # For every gene in the table.
  
  gene<-rownames (AvgExpression) [i]      # Take the gene.
  values<-as.numeric(AvgExpression [i,])  # Take the values from the entire row corresponding to the current gene.
  
max_value <- max(values)                  # Get the maximum value among 'values'.
second_max <- sort (values, decreasing = TRUE) [2] # Take the second highest value among 'values'.

avg_other_tissues<-(sum(values)-max_value)/(length (values)-1)  # Find the mean of the values excluding 'max_value'.

# Tissue-specific: max value is at least 4x second-highest value.
if (max_value>=4*second_max) {
tissue<-colnames(AvgExpression) [which.max(values)] # Take the tissue to which the 'max_value' for that gene corresponds.
tissue_specific[[tissue]]<-c(tissue_specific[[tissue]], gene) # Add 'gene' to the list for the specified 'tissue'.

# Tissue-elevated: max value is at least 2x the average of all other tissues.
 } else if (max_value>=2*avg_other_tissues) { 
tissue<-colnames (AvgExpression) [which.max(values)]  # Take the tissue to which the 'max_value' for that gene corresponds.
tissue_elevated[[tissue]]<-c(tissue_elevated[[tissue]], gene) # Add 'gene' to the list for the specified 'tissue'.

# Not-elevated: genes that do not fit other categories.
} else {
not_elevated<-c(not_elevated, gene)
  }
}
str (tissue_specific) # Show the results.
str (tissue_elevated)
head (not_elevated)
```

Now, we count the number of genes in each class. 
We will save the results in a matrix with 2 rows and 6 columns for better visualization. Rows will hold the 2 categories of genes: tissue specific, tissue elevated, and the 6 columns will correspond to every tissue.

```{r}
finalResults<-matrix(ncol=6,nrow=2)   # Create a matrix to store the final results.
colnames(finalResults)<-colnames(AvgExpression)
rownames(finalResults)<-c("specific","elevated")

# For every tissue
for (tissue in colnames(finalResults))
{
  specific_genes<-tissue_specific[[tissue]]   # Get specific genes for the tissues.
  elevated_genes<-tissue_elevated[[tissue]]   # Get elevated genes for the tissues.
  
  # We compute the number of specific and elevated genes using the length() function and we put it in the 'finalResults' matrix.
  finalResults['specific',tissue]<-length(specific_genes) 
  finalResults['elevated',tissue]<-length(elevated_genes)
}
finalResults # We display the total number of tissue specific genes and the total number of tissue elevated genes.
```


We display the results in a barplot. 

```{r}
barplot((finalResults), main="Classification of gene expression", # We use 'main' to set the title of the barplot.
        
        beside = TRUE, # Setting beside = TRUE, groups bars side by side for each column instead of stacking them.
        
        col=c("purple","pink"), # Sets the colors of the bars
        
        density = 45, # Density of the pattern inside the bars. Lower values make the lines more spaced out.
        angle=45,    # Angle at which the lines are drawn.
        legend=TRUE, # We want to add a legend to the plot.
        horiz=TRUE,  # Bars are drawn horizontally instead of vertically (the default).
        las=1, # Orientation of axis labels: 1 means always horizontal.
        xlim=c(0,4000), 
        xlab="ClassifiedGenes")  # Label for the x-axis.
```


#### Results interpretation:
These results show the classification of genes based on their expression levels across 6 tissues: brain, heart, kidney, liver, lung, muscle. The categories are: Specific and Elevated.

Specific: These are the genes whose expression is highly specific in one tissue compared to the others.
We can see that the brain has the highest number of tissue-specific genes (2006). The lung, kidney, liver, heart and muscle have, in a decreasing way, relatively fewer specific genes, suggesting comparatively less tissue-specific expression patterns in these tissues.

Elevated: These genes have an expression level elevated (likely 2x or more) in a specific tissue compared to others.
The brain, again, has the highest number of elevated genes (3266), followed by lung, kidney, heart, liver and muscle. This suggests that the brain not only has many highly specific genes but also a wide range of genes with elevated expression.
The muscle shows the lowest number of elevated genes (397), which could imply that it shares many genes with similar expression levels in other tissues.

Biological Relevance:
These results align with expectations, as the brain’s specialized functions require many unique genes. Conversely, muscle tissue, with its more generalized role across the body, shows a more uniform gene expression profile.



#### DE analysis for genes in brain and muscle.
Now we have to consider the DEGs that we obtained in the 'common part' and the two tissues that were compared, we have to find:
the number of DEGs that are 'elevated' in one or the other tissue and the number of DEGs that are 'specific' in one or the other tissue.

We use the *intersect()* function to find the all the genes that are elevated/specific in Brain/Muscle that are also DEGs between Brain and Muscle. We obtain their number using *length()*.

```{r}
head(DEGs) # These are all the previously calculated DE genes (DE_UP and DE_DOWN between Brain and Muscle).

tissues<-c("brain", "muscle") # Selected tissues.
gene_types<-list(ElevatedGenes=tissue_elevated, SpecificGenes=tissue_specific) 
results<-matrix(nrow=2, ncol=2)   # Create a matrix that will store the results.
rownames(results)<-names(gene_types) 
colnames(results)<-tissues
for (gene_type in names(gene_types)) {       # Iteration on ElevatedGenes and SpecificGenes.
  for (tissue in tissues) {                  # Iteration on brain and muscle.
    
intersect_genes <- intersect(gene_types[[gene_type]][[tissue]], DEGs) # Identify genes that are both in the specified gene type for the given tissue and in the DEGs list.

results[gene_type, tissue]<-length(intersect_genes) # We store the result in the matrix.
  }
}
results
```

We create a barplot showing the results.

```{r}
barplot(results, beside= TRUE,   # Bars plotted side by side.
        col = c("purple", "pink"),
        legend = rownames(results), 
        main = "Significant DEGs by Tissue and Gene Type", 
        ylim=c(0,2000),
        xlab = "Tissue", ylab = "Count of Significant DEGs")
```

#### Results interpretation:
The brain shows a higher overall count of significant DEGs compared to muscle, with both the categories of genes being more abundant in the brain.

#### Creating a function

Write a function that can take as the input the name of the gene and report:

1. the average expression of that gene in all the structures

2. the list of tissues (if any) where the gene was 'elevated'

3. the list of tissues (if any) where the gene was 'specific'

4. the output should also include a barplot of the average expression levels of all the replicates in all the structures.


```{r}
GeneAnalysis <- function(GeneName) {
  if (!(GeneName %in% rownames(CPM))) {     # Check if the gene exists in the counts table.
    stop('Gene not found among the expressed genes!')
  }
  geneExpression<-CPM[GeneName, ]         # Extract expression values for the gene
  
  levs<-unique(design$Tissue)             # Identify unique tissues

  
  # Create an empty matrix to store results.
  avgTable<-matrix(ncol=length(levs), nrow = 1)  
  rownames(avgTable)<-GeneName
  colnames(avgTable)<-levs
  
  # Compute average expression per tissue.
  for (tissue in levs) {
    selection<-design$Tissue==tissue  # Select samples from the current tissue.
    avgTable[, tissue]<-mean(geneExpression[selection], na.rm = TRUE)
  }
  
  # Initialize vectors to store elevated and specific tissues.
  elevated_tissues<-c()
  specific_tissues<-c()
  
  # Get the gene's expression values across tissues.
  values<-as.numeric(avgTable)

  max_value<-max(values, na.rm = TRUE) 
  second_max<-sort(values, decreasing = TRUE)[2]
  avg_other_tissues<-(sum(values, na.rm = TRUE)-max_value)/(length(values)-1)
  
  # Identify elevated and specific tissues.
  max_tissue<-colnames(avgTable)[which.max(values)]
  
  if (max_value >= (4 * second_max)) {
    specific_tissues<-c(specific_tissues, max_tissue)
  } else if (max_value >= (2 * avg_other_tissues)) {
    elevated_tissues<-c(elevated_tissues, max_tissue)
  }
  
  # Barplot for visualization.
  barplot(values, 
          names.arg = colnames(avgTable),
          col = "pink",
          main = paste('Average Expression of gene:', GeneName), 
          ylab = 'Expression Level', 
          xlab = 'Tissues', 
          las = 1)
  
  # Prepare result as a list.
  result<-list(
    AverageExpression=avgTable,
    ElevatedTissues=ifelse(length(elevated_tissues) > 0, elevated_tissues, "None"),
    SpecificTissues=ifelse(length(specific_tissues) > 0, specific_tissues, "None")
  )
  
  return(result)
}

# Example usage
GeneAnalysis ('A2M') # We test the function with a gene that we know is expressed in the normalized counts table. 
```
A2M is predominantly expressed in the lung, as indicated by the much higher expression value compared to other tissues.
It is not broadly elevated across multiple tissues, but rather specifically expressed in the lung.
This suggests A2M may have a lung-specific function, potentially related to lung physiology or pathology.



