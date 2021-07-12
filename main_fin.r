if (!require("plyr")) install.packages("plyr")

library(plyr)

source("http://peterhaschke.com/Code/multiplot.R")

#getSub - функция позволяет получить заданную процентную выборку OTUs.
#На вход подаётся таблица otu. Строки - otu, колонки - samples
#getSub function generates subsample by percentile or by sample size.

getSub <- function(data, value=0.95, cutoff=1, subsample=FALSE, subsize=NULL) {
	#getSub - функция позволяет получить заданную процентную выборку OTUs.
	#На вход подаётся shared файл из мотура.
    data <- data[rowSums(data) > cutoff,]
    if (subsample)
    {
        print(sprintf("using subsample size: %i", subsize))

        #заводим список с первым элементом, который содержит имена всех известных OTUs
        out <-list()
        out[["origin"]] <- data.frame(Var1=rownames(data))

        for (i in names(data))
        {
            X <- data[i]
            
            print(i); flush.console()
            
            sumOtus <- rowSums(X)

            S<-c()
            for (n in 1:length(sumOtus))
            {
                S <- c(S, rep(names(sumOtus)[n], sumOtus[[n]]))
            }

            DF <- data.frame(table(sample(S, size=subsize, replace = FALSE)), stringsAsFactors=FALSE)
            names(DF) <- c("Var1", i)
            out[[i]] <- DF
        }
        res <- join_all(out, by = "Var1", type = 'left', match="first")
        res[is.na(res)]<-0
        rownames(res)<-res$Var1
        res <- res[-1]
        print("Number of Otus for each sample:")
        print(colSums(res !=0))
        print("Number of seqs for each sample:")
        print(colSums(res))
        
    }
    
    else
    {
        print(paste("Selected sample size: ", value*100))
        print(paste("Selected cutoff: ", cutoff))
        flush.console()
        #удаляем синглтоны
        #datan <- data[rowSums(data) > cutoff,]
        #сортируем по увеичению размеры выборки
        ndf <- data[order(colSums(data), decreasing=F)]
        #Теперь первый столбец это самоя малочисленныя проба
        #суммируем otus
        sumOtus <- rowSums(ndf)
        #считаем общее число последовательностей
        total <- sum(sumOtus)
        #создаём реплики имён OTU чтоб вычислить процентную часть сообщества
        S<-c()
        for (n in 1:length(sumOtus))
         {
           S <- c(S, rep(names(sumOtus)[n], sumOtus[[n]])) 
         }
        #например, выделение 95% выборки
        #находим значение колиства последовательностей для выборки в 95% от общего числа последовательностей
        samSize <- as.integer(value * total)
        #получеам процентную выборку нашей самой маленькой пробы
        percSample <- S[1:samSize] 
        #преобразуем OTUs в таблицу 
        # item_table <- table(percSample) 
        # a <- names(item_table[length(item_table)]) 
        # создаём выборку исходного набора данных, зная имя крайней OTU нашей % выборки
        res <- subset(data, rownames(data)<=percSample[length(percSample)])        
        print("Number of Otus for each sample:")
        print(colSums(res !=0))
        print("Number of seqs for each sample:")
        print(colSums(res))
    }
    
	return(res);
}


#создаём реплики Otu заданной выборки. На вход подаётся один образец. Например getSample(nd[1])
#generates a vector of OTU sample. 
#Example: getOtuSample(nd[1])

getOtuSample <- function(data, trimmedby = 2)
{
  trimmedData <- data[data > trimmedby,, drop=FALSE]
  sumOtus <- rowSums(trimmedData)
  total <- sum(sumOtus)
  S<-c()
	for (n in 1:length(sumOtus))
	 {
	   S <- c(S, rep(names(sumOtus)[n], sumOtus[[n]])) 
	 }
  return(S);
}

#generates OTU replicas. 
#Example: getReplica(otu_list[[j]], n=1000, scaling="int")

getReplica <- function(OtuVector, n = 1000, scaling="int")
{
    otusNames <- unique(OtuVector) #OtuVector is character vector of otus
    itemSize<-length(OtuVector)
    res <- data.frame(Var1 = otusNames)
    convToDf <- function(Var1, itemSize)
    {
        repOtusNames <- unique(Var1)
        DF <- data.frame(table(sample(Var1, size=itemSize, replace = TRUE)))
        return(DF)
    }
    XX <- replicate(n+1, convToDf(OtuVector), simplify = FALSE)
    XX[[1]] <- data.frame(Var1 = otusNames)
    XXres <- join_all(XX, by = 'Var1', type = 'left')
    XXres[is.na(XXres)]<-0
    rownames(XXres) <- otusNames
    XXres <- XXres[-1]
    colnames(XXres) <- seq(1, length(XXres))
    if (scaling == "percent")
        { XXres <- prop.table(as.table(as.matrix(XXres)), 2) }
    else if (scaling == "log")
        { XXres <- log(XXres+1) }
    XXres <- t(XXres)
    return(XXres)
}


#вычисляем относительную ошибку для наших выборок
#на вход элемент списка rep_list[i]
#estimate relative error. Enter rep_list[i]
repError <- function(data, rate = 0.2)
{
	setName <- names(data)
	x <- data[[1]]
	V <- data.frame()
	quantileMatrix <- apply(x, 2, quantile, probs = c(0.025,0.975))
	for (i in 1:ncol(quantileMatrix))
		{
			D <- quantileMatrix[, i, drop=FALSE]
			V <- rbind(V, data.frame(set=((max(D)-min(D))/max(D)), row.names=colnames(D)))
		}
	out <- V[V[,1] <= rate,,drop=FALSE]
	names(out)<-setName
    out$otuid <- rownames(out)
    out <- out[c("otuid", setName)]
	return(out)
}


