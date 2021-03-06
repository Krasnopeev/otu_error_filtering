# otu_error_filtering is an algorithm for estimate of OTU's representation relative error. 
This pipeline using for analysis of high throughput sequencing data, where is final result is an abundance table. For example common used OTU table. With this tool you can filter your dataset based on relative error of OTU representation in abundance table.
# installing dependent packages
You need to install `plyr` and `tidyverse` from R repo
```R
install.packages(c("plyr", "tidyverse"))

```
At first we need to load our function to the R

```R
source("C:\\path_to_file\\main_fin.r")
```

Input *.shared file from mothur

```R
data <- read.table("shared.shared", header=TRUE, sep="\t")
```

little modify of shared file from mothur
for further work you just need an OTU table where rows are OTUs and columns are Samples

```R
rownames(data) <- data$Group

data <- as.data.frame(t(data[,-c(1:3)]))
```

Generating a subsampled dataset with `cutoff` - trimm single an doubletons (you can icrease trimming `cutoff` threshold for your purposes), `value` - subsample of dominant OTUs by total count (1 equal 100%, 0.95 is a 95% pool etc.), `subsample` is depending on `subsize`. Set `subsample=TRUE` to get euqal size subsamples for each sample in your data with a given `subsize` value. You can set `subsize` depending of smallest sample size in your dataset (for ex. smallest sample in pool is 15000 sequences).

```R
ndt <- getSub(data, value=1, cutoff=2, subsample=FALSE, subsize=15000)
```

generating OTUs list for each sample

```R
otu_list <- list()

for (i in names(ndt))
  {
    otu_list[[i]] <- getOtuSample(ndt[i], trimmedby=0); print(i); flush.console()
  }
```

generating OTUs replicas for each sample via bootstrapping of OTU vectors from `otu_list`

```R
rep_list <- vector("list", length(names(otu_list)))

names(rep_list) <- names(otu_list)

for (j in names(otu_list))
  {  
    print(j); flush.console()
    rep_list[[j]] <- getReplica(otu_list[[j]], n=1000, scaling="int")
  }
```

Estimate relative error for OTUs. Feel free to try different `rate` values. The deafault error `rate` value is `0.2`.
For `repError` function In case if `filterErorr=TRUE` all error values higher than `rate` will be dropped and missed values replaced with `1`. You can skip filtering step by using `filterErorr=FALSE` and take a look at final result without error filtering.
Moreover you can look at `VmaxTable` and `VminTable` which are contain maximum and minimum values of OTU abundance after bootstrapping procedure. Values in these tables are converted to percent.

```R

ErrorsList <- vector("list", length(names(otu_list)))
VmaxList <- vector("list", length(names(otu_list)))
VminList <- vector("list", length(names(otu_list)))

names(ErrorsList) <- names(otu_list)
names(VmaxList) <- names(otu_list)
names(VminList) <- names(otu_list)

for(i in names(rep_list))
  {
    ErrorsList[[i]] <- as.data.frame((repError(rep_list[i], rate=0.2, filterErorr=FALSE)$error))
    VmaxList[[i]] <- as.data.frame((repError(rep_list[i], rate=0.2, filterErorr=FALSE)$Vmax))
    VminList[[i]] <- as.data.frame((repError(rep_list[i], rate=0.2, filterErorr=FALSE)$Vmin))
  }

OtuErrorTable <- reduce(ErrorsList, full_join, by = "otuid") %>% replace(., is.na(.), 1);
VmaxTable <- reduce(VmaxList, full_join, by = "otuid") %>% replace(., is.na(.), 1);
VminTable <- reduce(VminList, full_join, by = "otuid") %>% replace(., is.na(.), 1);

```

Generating new OTU table with reliable values. OTU values which estimated error rates were higher than values in `OtuErrorTable` will be replaced with `0`.

```R

NewOtuList <- vector("list", length(names(otu_list)))

names(NewOtuList) <- names(otu_list)

for(i in names(ErrorsList))
  {
    OtuNames <-rownames(ErrorsList[[i]])
    NewOtuList[[i]] <- data.frame(otuid=OtuNames, ndt[OtuNames, i, drop=FALSE])
  }

newOtuTable <- reduce(NewOtuList, full_join, by = "otuid") %>% replace(., is.na(.), 0);

```

Write to disk new OTU table

```R
write.table(newOtuTable, "newOtuTable.txt", sep='\t', col.names=TRUE, quote=FALSE)
```

