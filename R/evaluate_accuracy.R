LR_evaluation<-function(partial_EBV,  #验证群个体的EBV， 根据partial数据集计算
						whole_EBV,     #完整数据得到的验证群EBV， 根据完整数据计算
						candidate_id=NULL,
						type=c("bias","disperation","accuracy","ratio_of_accuracy","ratio_of_reliability"),
						partial_sigma_a2=NULL, # 利用 partial 数据估计得到的方差组分
						K=NULL,                   # relationship matrix 
						F=NULL,
						f2=NULL
						){
cat("The format of user-provided EBV should contain two columns: id  and EBV! \n")
if(is.matrix(partial_EBV)){partial_EBV=data.frame(partial_EBV,stringsAsFactors=F)}
colnames(partial_EBV)=c("Id","partial_Effect")

result=partial_EBV
result$whole_Effect=whole_EBV[match(result[,1],whole_EBV[,1]),2]

#candidate
result=result[result[,1]%in%candidate_id,]
if(nrow(result)==0){stop("Couldn't find candidate individual in provided EBV data!")}
bias=disperation=accuracy=ratio_of_accuracy=ratio_of_reliability=NULL
# statistics
if("bias"%in%type){bias=mean(result$partial_Effect)-mean(result$whole_Effect)}

if("disperation"%in%type){disperation=as.numeric(lm(result$whole_Effect~result$partial_Effect)[[1]][2])} # cov(w,p)/var(p)

if("accuracy"%in%type){
if(is.null(K)&is.null(F)&is.null(f2)){stop("Please provide Kinship matrix or inbreeding coefficients !")}
if(!is.null(K)){
cat("Please make sure the colnames of Kinship matrix are individual name! \n")
if(is.null(partial_sigma_a2)){stop("Please provide partial_sigma_a2!")}
IND_name=rownames(K)
K=K[match(as.character(candidate_id),IND_name),match(as.character(candidate_id),IND_name)];gc();
F=mean(diag(K)-1)
f2=mean(K[upper.tri(K)])
rm(K);gc();
}
accuracy=sqrt(cov(result$partial_Effect,result$whole_Effect)/((1+F-f2)*partial_sigma_a2))
}

if("ratio_of_accuracy"%in%type){ratio_of_accuracy=cor(result$partial_Effect,result$whole_Effect)}

if("ratio_of_reliability"%in%type){ratio_of_reliability=cov(result$whole_Effect,result$partial_Effect)/var(result$whole_Effect)}

return(list(bias=bias,disperation=disperation,accuracy=accuracy,
			ratio_of_accuracy=ratio_of_accuracy,ratio_of_reliability=ratio_of_reliability))
}



#相关系数检验
Hotelling_test <- function(gen_value,phe,two_tailed=TRUE){
	
   cat("Please make sure you've been installed package:psych! \n")	
   cat("The format of gen_value is: id,ebv1,ebv2...... \n")   
   cat("The format of phe is: id,phe \n")

   gen_value=data.frame(gen_value,stringsAsFactors=F)
   phe=data.frame(phe,stringsAsFactors=F)
   if(ncol(gen_value)<=2){stop("The column of gen_value should >=3 !")}
   if(!ncol(phe)==2){stop("The column of phe should ==2 !")}
   
   
   intersect_id=intersect(gen_value[,1],phe[,1])
   if(length(intersect_id)==0){stop("There are no identical id between gen_value and phe!")}
   gen_value=gen_value[match(intersect_id,gen_value[,1]),]
   phe=phe[match(intersect_id,phe[,1]),]
   gen_value=gen_value[,-1]
   phe=phe[,-1]
   
cat("Off-digonal is p-value \n")   
   gen_name <- colnames(gen_value)
   results <- matrix(ncol=ncol(gen_value),nrow=ncol(gen_value))
   rownames(results)=colnames(results)=gen_name
   for (i in 1:ncol(gen_value)){
     for (j in 1:ncol(gen_value)){
       y <- gen_value[,i]
       z <- gen_value[,j]
       if(i==j){next}
       xy <- cor(phe,y)
       xz <- cor(phe,z)
       yz <- cor(y,z)
       n <- length(y)
       results [i,j] <- psych::paired.r(xy, xz, yz, n,twotailed=two_tailed)$p
     }
   }
   return(results)
 }	
