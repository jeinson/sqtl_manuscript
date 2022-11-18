meanDiff = function(data, func){
        data = func(data)
        m1 <- median(data[data$V2=='sQTL', "totalData"])
        m2 <- median(data[data$V2=='non_sQTL', "totalData"])
        return(m1-m2)
}

inv_norm <- function(x) qnorm((rank(x, na.last='keep', ties.method = "random") - 0.5)/sum(!is.na(x)))

my_bootstrap = function(x, R, func, func1){
        cores=detectCores()
        cl <- makeCluster(cores[1]-1) #not to overload your computer
        registerDoParallel(cl)
        dist<- foreach(
                i = 1:R, 
                .combine = c,
                .packages = 'foreach'
        ) %dopar% (func(x, func1))
        parallel::stopCluster(cl)
        return(dist)
}

select_samples = function(x){
        d = x[x[,1] %in% sample(x[,1], replace = TRUE),]
        colnames(d) = c("totalData", "V2")
        print(nrow(d))
        return(d)
}


ks_stat = function(x, y){
        x = x[!is.na(x)]
        y = y[!is.na(y)]
        pval =  ks.test(x, y)$p.value
        totalData = c(x, y)
        totalData = cbind(totalData, c(replicate(length(x), "sQTL"), replicate(length(y), "non_sQTL")))
        totalData[,1] = inv_norm(as.numeric(totalData[,1]))
        totalData = data.frame(totalData)
        totalData = totalData[!is.na(totalData[,1]), ]
        totalData[,1] = as.numeric(totalData[,1])
        print(mean(totalData[totalData$V2 == "sQTL", 'totalData']) - mean(totalData[totalData$V2 == "non_sQTL", 'totalData']))
        
        totalBoot = my_bootstrap(totalData, R=1000, meanDiff, select_samples)
        totalBootCI = quantile(totalBoot, probs=c(0.05, 0.5, 0.95))
        conf_low = totalBootCI[1]
        conf_top = totalBootCI[3]
        # stat = mean(totalData[totalData$V2 == "sQTL", 'totalData']) - mean(totalData[totalData$V2 == "non_sQTL", 'totalData'])
        stat = totalBootCI[2]
        print(c(pval, stat, conf_low, conf_top))
        return(c(pval, stat, conf_low, conf_top))
}