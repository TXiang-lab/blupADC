CVdatasets = R6Class("CVdatasets",

	#inherit = ADCmodel,
     public = list(
	
	iters=1,
	train_sets=NULL,
	test_sets=NULL,
	ebv_sets=NULL,

	initialize=function(train_sets=NULL,test_sets=NULL,ebv_sets=NULL){


		self$train_sets=train_sets 
		self$test_sets=test_sets			
		self$ebv_sets=ebv_sets
		
		
	},
	print=function(...){
	
		if(!is.null(self$train_sets)){
		
			iters=length(self$train_sets)
		
			cat(paste0("<CVdatasets::",iters,"> \n"))
			
			
			
			for(i in 1:iters){
				n_train_col=ncol(self$train_sets[[i]])
				n_train_row=nrow(self$train_sets[[i]])		
				n_test_col=ncol(self$test_sets[[i]])
				n_test_row=nrow(self$test_sets[[i]])	
				
				if(i>=3){cat("...... \n");cat("...... \n");break}
				#for train sets
				cat(paste0("***<Train set",i,"(",n_train_row,"x",n_train_col,")> \n"))	
				if(n_train_col>=10){
					print(head(self$train_sets[[i]][,1:10]))
				}else{
					print(head(self$train_sets[[i]]))
				}	
				
				#for test sets
				cat(paste0("***<Test set",i,"(",n_test_row,"x",n_test_col,")> \n"))	
				if(n_test_col>=10){
					print(head(self$test_sets[[i]][,1:10]))
				}else{
					print(head(self$test_sets[[i]]))
				}	

				#for predicted sets
				if(!is.null(self$predict_sets)){
					n_predicted_row=nrow(self$predict_sets[[i]])
					n_predicted_col=ncol(self$predict_sets[[i]])
					cat(paste0("***<Predicted set",i,"(",n_predicted_row,"x",n_predicted_col,")> \n"))	
					if(n_predicted_col>=10){
						print(head(self$predict_sets[[i]][,1:10]))
					}else{
						print(head(self$test_sets[[i]]))
					}	
				}
			}
		
		}
	
	}
	
	)
)

CVpredi = R6Class("CVpredi",

	#inherit = ADCmodel,
     public = list(
	sets=NULL,
	seed=NULL,	
	task=NULL,train_rows=NULL,test_rows=NULL,
	partion_type=NULL,
	kfold=NULL,
	p_out=NULL,
	ratio=NULL,
	dmu_prior_file=NULL,
	dmu_module=NULL,
	order_date=NULL,split_date=NULL,birth_pos=NULL,
	estimators=NULL,
	iters=1,
	bias=NULL,
	accuracy=NULL,
	datasets=NULL,
	output_path=NULL,
	initialize=function(task=NULL,train_rows=NULL,test_rows=NULL,datasets=CVdatasets$new(),
						partion_type="holdout",# holdout, kfold, lpout(leave p out), birthout(need to provide birth date information), boostrap(next step)  using future package for multiple process
						kfold=5,
						p_out=1,
						ratio=0.8, # train/total
						order_date=FALSE,split_date=NULL,birth_pos=NULL,
						output_path=NULL,
						estimators="classical", #LR method or normal
						dmu_prior_file=NULL,
						dmu_module=NULL,
						seed=9498){ #DMU=NULL,BLUPF90=NULL,Bayes=NULL){				


		self$seed=seed 
		self$partion_type=partion_type
		self$kfold=kfold
		self$ratio=self$ratio
		self$order_date=order_date
		self$split_date=split_date
		self$birth_pos=birth_pos
		self$estimators=estimators
		self$task=task
		self$train_rows=train_rows
		self$test_rows=test_rows
		self$datasets=datasets
		self$dmu_module=dmu_module
		self$dmu_prior_file=dmu_prior_file
		
		if(!is.null(output_path)){
			self$output_path=output_path
		}
		
		#split dataset, this requires all task must have phenotype			

		if("BLUP"%in%class(task)){
			
			if(is.null(task$pars$phe))task$pars$read_phe()
			
			phe=task$pars$phe	
		
		}else if("Bayes"%in%class(task)){
		
	
		}

		train_sets=NULL 
		test_sets=NULL
		n_row=nrow(phe)
		test_rows=NULL
		
		
		#split dataset
		if("kfold"%in%partion_type){
	
			self$iters=kfold
			n_select=ceiling(n_row*1/kfold)
			set.seed(seed)
			test_rows=list(sample(1:n_row,n_select,replace=FALSE))
			train_rows=list(setdiff(1:n_row,test_rows[[1]]))
			
			for(i in 2:kfold){
				
				set.seed(seed)
				
				if(i!=kfold){
				
					i_test_rows=sample(setdiff(1:n_row,do.call(base::c,test_rows)),n_select,replace=FALSE)
				}else{
				
					i_test_rows=setdiff(1:n_row,do.call(base::c,test_rows))
				}		

				test_rows=c(test_rows,list(i_test_rows))	
				train_rows=c(train_rows,list(setdiff(1:n_row,i_test_rows)))
					
					#i_train_rows=setdiff(1:n_row,i_test_rows)
					#i_test_phe=phe[i_test_rows,]
					#i_train_phe=phe[i_train_rows,]
					#train_sets=c(train_sets,list(i_train_phe))
					#test_sets=c(test_sets,list(i_test_phe))
		
			}
					
			train_sets=lapply(train_rows,function(x)phe[x,])
			test_sets=lapply(test_rows,function(x)phe[x,])

		}else if("holdout"%in%partion_type){
			self$iters=1
			if(!is.null(train_rows)){
			
				i_train_phe=phe[train_rows,]
				i_test_phe=phe[setdiff(1:n_row,train_rows),]
			
			}else{
			
				set.seed(seed)
				i_test_rows=sample(1:n_row,ceiling(n_row*(1-ratio)),replace=FALSE)
				i_train_rows=setdiff(1:n_row,i_test_rows)
					
				i_test_phe=phe[i_test_rows,]
				i_train_phe=phe[i_train_rows,]			
			}
			
			train_sets=list(i_train_phe)
			test_sets=list(i_test_phe)
			
			
		}else if("birthout"%in%partion_type){		
				self$iters=1
				cat("Please make sure phenotype contains birth date column! \n")
				
				if(is.null(birth_pos))stop("Please specify the position of birth data in phenotype!")
				
				 if(!is.null(split_date)){
				
					 cat(paste0("Individuals with birth date earlies than ",split_date," would be divided into train set, the rest would be divided into test set \n"))
					
					birth_order=which(as.numeric(phe[,birth_pos]) >=as.numeric(split_date))
					i_test_phe=phe[birth_order,]
					i_train_phe=phe[setdiff(1:n_row,birth_order),]
					
					cat(paste0("The ratio of train set on whole dataset is: ",round(nrow(i_train_phe)/n_row,3)," \n"))
					
				 }else{
				
					 tmp_split_data=sort(as.numeric(phe[,birth_pos]))[ceiling(n_row*ratio)]
					
					 birth_order=which(as.numeric(phe[,birth_pos])>=as.numeric(tmp_split_data))
					 i_test_phe=phe[birth_order,]
					 i_train_phe=phe[setdiff(1:n_row,birth_order),]
					 cat(paste0("Based on provided ratio:",ratio,",individuals with birth date earlies than ",
								tmp_split_data," would be divided into train set, the rest would be divided into test set \n"))
									
			
				 }
				
				
				 train_sets=list(i_train_phe)
				 test_sets=list(i_test_phe)

		
			} 
		self$datasets$train_sets=train_sets
		self$datasets$test_sets=test_sets
				
		#define the predict ebv set 
		self$datasets$ebv_sets=vector(mode = "list", length = self$iters)
		
		self$bias=rep(NA, self$iters)
		self$accuracy=rep(NA, self$iters)
		},
		
		
	get_set=function(i=NULL){
	
			if(!is.null(i)){
				return(list(train=self$datasets$train_sets[[i]],test=self$datasets$test_sets[[i]]))
			}else{
				return(list(train=self$datasets$train_sets,test=self$datasets$test_sets))
			}

	},

	predict_sets=function(i=NULL){

		if(!is.null(self$dmu_module)){self$task$pars$dmu_module=self$dmu_module}
		if(!is.null(self$dmu_prior_file)){self$task$pars$dmu_prior_file=self$dmu_prior_file}
		
		if(!is.null(self$output_path)){
			self$task$pars$output_path=self$output_path
		}
		unchanged_path=self$task$pars$output_path
		
		sets=NULL
		

		
		
		#specify the train set and test set
		if(!is.null(i)){
		
			if(i>self$iters){
			
				stop(paste0("Index outbound of iters, please make sure the iters <= ",self$iters," !"))
			
			}
		
			iters_set=i
		}else{
		
			iters_set=1:self$iters
		}
		
		for(i_iters in iters_set){	
		
			cat(paste0("Start genetic evalution on dataset ",i_iters," ! \n"))
			
			self$task$pars$output_path=paste0(unchanged_path,"/CV",i_iters)
			self$task$pars$write_phe(phe=self$get_set(i_iters)$train,file_name=paste0("pheno_CV",i_iters,".txt"))
			self$task$pars$phe_colnames=colnames(self$get_set(i_iters)$train)
			
			self$task$run()			
			predict_set=self$task$ebv
			self$datasets$ebv_sets[[i_iters]]=predict_set #predict_set
			test_set=self$datasets$test_sets[[i_iters]]
			
			intersect_ids=intersect(test_set[,1],predict_set[,1])
			i_pos1=match(intersect_ids,test_set[,1])  #assume id in the first column
			i_pos2=match(intersect_ids,predict_set[,1]) #assume id in the first column
					
			test_set_ebv=as.numeric(test_set[i_pos1,6])
			predict_set_ebv=as.numeric(predict_set[i_pos2,2])
			
			accuracy=cor(test_set_ebv,predict_set_ebv)
			bias=as.numeric(lm(test_set_ebv~predict_set_ebv)[[1]][2])
			
			cat(paste0("Accuracy = ",round(accuracy,3),"; Bias = ",round(bias,3),"\n"))
		
			self$bias[i_iters]=bias
			self$accuracy[i_iters]=accuracy
		
		}
		
		self$task$pars$output_path=unchanged_path
		
	}	
	#perform analysis	
				
	
	
	),
	
	active=list(
		
		dynamic_sets=function(partion_type=self$partion_type){
			#split dataset
		if("kfold"%in%partion_type){
		
			self$iters=kfold
			n_select=ceiling(n_row*1/kfold)
			set.seed(seed)
			test_rows=list(sample(1:n_row,n_select,replace=FALSE))
			train_rows=list(setdiff(1:n_row,test_rows[[1]]))
			
			for(i in 2:kfold){
				
				set.seed(seed)
				
				if(i!=kfold){
				
					i_test_rows=sample(setdiff(1:n_row,do.call(base::c,test_rows)),n_select,replace=FALSE)
				}else{
				
					i_test_rows=setdiff(1:n_row,do.call(base::c,test_rows))
				}		

				test_rows=c(test_rows,list(i_test_rows))	
				train_rows=c(train_rows,list(setdiff(1:n_row,i_test_rows)))
					
					#i_train_rows=setdiff(1:n_row,i_test_rows)
					#i_test_phe=phe[i_test_rows,]
					#i_train_phe=phe[i_train_rows,]
					#train_sets=c(train_sets,list(i_train_phe))
					#test_sets=c(test_sets,list(i_test_phe))
		
			}
					
			train_sets=lapply(train_rows,function(x)phe[x,])
			test_sets=lapply(test_rows,function(x)phe[x,])

		}
		self$datasets$train_sets=train_sets
		self$datasets$test_sets=test_sets
		}
	
	
	)
	
)














