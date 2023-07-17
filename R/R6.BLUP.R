library(R6)



#calculate the se of variance components
my_cal_se<-function(expr,vars_value,inv_ai=NULL){  #expr( R expresson ), vars_value(variance components)
			ai_name=colnames(inv_ai)  #random effect names of ai matrix
			for (i in 1:length(ai_name))
				assign(ai_name[i], vars_value[i])
				grad <- t(as.numeric(attr(eval(deriv(expr, ai_name)), "gradient")))
				h2<-eval(parse(text=as.character(expr)[-1]))
				return(list(value=h2,value_se=ifelse(diag(grad %*% inv_ai %*% t(grad))<=0,0,sqrt(diag(grad%*% inv_ai %*% t(grad))))))
}

vars_se <- R6Class("vars_se",
  public = list(
		invAI_mat=NULL,
		vars_mat=NULL,
		r2=NULL,
		h2=NULL,

initialize=function(invAI_mat=NULL,
					vars_mat=NULL,
					r2=NULL,
					h2=NULL
					){
		
	self$invAI_mat=invAI_mat
	self$vars_mat=vars_mat 
	self$h2=h2
	self$r2=r2
	},
	
	cal_se=function(expr,mat=NULL,invAI_mat=NULL){
		
		if(is.null(mat)&is.null(invAI_mat)){
		return(my_cal_se(expr,self$vars_mat[,1],self$invAI_mat))
		}else{
		return(my_cal_se(expr,mat,invAI_mat))		
		}
	}
	
)
)


	 
BLUPpar = R6Class("BLUPpar",

     public = list(
	 dmu_ped_inbred_number=2,phe_file=NULL,phe_colnames=NULL,analysis_model=NULL,genetic_effect="Id",add_pe=FALSE, add_domi=FALSE, missing_value= -9999,
	 iteration_criteria=1.0e-7,kinship_file=NULL,residual_cov_trait=NULL,
	 output_path=base::getwd(),output_name=NULL,show_message=TRUE,
	 # corrected phenotype	
	 residual_average=TRUE,cal_debv=FALSE, debv_pedigree_file=NULL,cal_reliability=FALSE,debv_id=NULL,
	 dmu_software_path=NULL,blupf90_software_path=NULL,IND_geno_file_name="IND_geno.txt",IND_geno=NULL,SSBLUP_omega=0.05,plot_variance=FALSE,
	 
	 #for dmu 
	 dmu_module="dmuai",dmu_algorithm_code=NULL,dmu_prior_file=NULL,  dmu_integer_n=NULL, blupf90_algorithm="aireml",
	 
	 #genotype data
	 geno=NULL,
	 
	 
	 #phenotype data
	 phe=NULL,
	 
	 #result 
	 result=NULL,
	 
    initialize = function(				   
					     dmu_ped_inbred_number=2, #wheter consider inbreeding in constructing A
						phe=NULL,ped=NULL,geno=GenoData$new(), #user only need to provide the data, then DMU can convert these files directly, and to the subsquent analysis
						phe_file=NULL,phe_colnames=NULL,analysis_model="PBLUP_A",genetic_effect="Id",add_pe=FALSE, add_domi=FALSE, missing_value= -9999,
					     iteration_criteria=1.0e-7,kinship_file=NULL,dmu_module="dmuai",dmu_algorithm_code=NULL,dmu_prior_file=NULL, residual_cov_trait=NULL, 
					     dmu_integer_n=5,  #number of integer variable
						output_path=base::getwd(),output_name=NULL,show_message=TRUE,
						blupf90_algorithm="aireml",
						#task_type="DMU",
					     # corrected phenotype	
					     residual_average=TRUE,cal_debv=FALSE, debv_pedigree_file=NULL, #pedigree for calculating DEBV		
					     cal_reliability=FALSE,debv_id=NULL,
						dmu_software_path=blupADC_software_path(),
						blupf90_software_path=blupADC_software_path(),						
					     IND_geno_file_name="IND_geno.txt",IND_geno=NULL,SSBLUP_omega=0.05,
					     plot_variance=FALSE){
						 
                     self$phe_colnames=phe_colnames
				 self$phe=phe
                     self$phe_file=phe_file
                     self$analysis_model=analysis_model 
                     self$dmu_ped_inbred_number=dmu_ped_inbred_number
                     self$genetic_effect=genetic_effect
                     self$add_pe=add_pe 					 
                     self$add_domi=add_domi
                     self$missing_value=missing_value
                     self$iteration_criteria=iteration_criteria 					 	 
                     self$kinship_file=kinship_file
                     self$dmu_module=dmu_module
                     self$dmu_algorithm_code=dmu_algorithm_code 						 
                     self$dmu_prior_file=dmu_prior_file
                     self$residual_cov_trait=residual_cov_trait
                     self$dmu_integer_n=dmu_integer_n 						 
                     self$output_path=output_path
                     self$output_name=output_name 							 				 
				  # corrected phenotype	 
                     self$residual_average=residual_average
                     self$cal_debv=cal_debv
                     self$debv_pedigree_file=debv_pedigree_file 
                     self$cal_reliability=cal_reliability
                     self$debv_id=debv_id
                     self$dmu_software_path=dmu_software_path 	
				 self$blupf90_software_path=blupf90_software_path=NULL	
                     self$IND_geno_file_name=IND_geno_file_name
                     self$IND_geno=IND_geno
                     self$SSBLUP_omega=SSBLUP_omega 					 	 
                     self$plot_variance=plot_variance	
				 #self$task_type=task_type
				 self$show_message=show_message
    },
	
	write_phe=function(phe=NULL,file_name=NULL){
		if(is.null(file_name))file_name="temp_pheno.txt"
		if(is.null(phe))phe=self$phe
		#if(!is.null(phe)&is.null(self$phe_file)){	 #phenotype
		if(!is.null(phe)){	 #phenotype
			if(is.null(colnames(phe)))stop("Phenotype should have colnames!")
			self$phe_colnames=colnames(self$phe)
			
			if(!file.exists(self$output_path))dir.create(self$output_path,recursive=TRUE)
			
			fwrite(phe,paste0(self$output_path,"/",file_name),quote=F,row.names=F,col.names=F,sep=" ")
			self$phe_file=paste0(self$output_path,"/",file_name)
		}
	},
	
	read_phe=function(){
	
		if(is.null(self$phe)&!is.null(self$phe_file)){	 #phenotype

			if(is.null(self$phe_colnames))stop("Parameter: phe_colnames shouldn't be NULL!")
			
			self$phe=fread(self$phe_file,data.table=F)
			colnames(self$phe)=self$phe_colnames
		}
	 invisible(self)	
	}
	
  )
)

BLUP = R6Class("BLUP",

	inherit = ADCmodel,

     public = list(
		
		pars=NULL,
		task_type=NULL,
		ebv=NULL,	
		r2=NULL,h2=NULL,
		vars_se=vars_se$new(),
		
	initialize=function(fixed=NULL,covariate=NULL,random=NULL,polyno=NULL,id_name="id",dam_name="dam",pe_name=NULL,formulas=ADCformula$new(),
						task_type="DMU",pars=NULL){#pars=DMUpar$new(),blupf90_pars=NULL){
	
				
	
				super$initialize(fixed=fixed,covariate=covariate,random=random,polyno=polyno,id_name=id_name,dam_name=dam_name,pe_name=pe_name,formulas=formulas)
				
				# pars is the same for DMU and BLUPF90
				# if(task_type=="DMU"&"BLUPF90par"%in%class(pars)){
				
					# stop("The class of task_type and pars is incompatible!")

				# }else if(task_type=="BLUPF90"&"DMUpar"%in%class(pars)){
				
					# stop("The class of task_type and pars is incompatible!")
				
				# }
				
				# if(is.null(pars)){
					
					# if("DMU"%in%task_type){
					
						# pars=DMUpar$new()
					
					# }else if("BLUPF90"%in%task_type){
						
						# pars=BLUPF90par$new()
					
					# }
					
				# }	
				
				self$pars=pars
				self$task_type=task_type
	},
	
	
	
	run=function(dmu_ped_inbred_number=NULL,phe=NULL,ped=NULL,geno=NULL,phe_file=NULL,phe_colnames=NULL,analysis_model=NULL,
	               genetic_effect="Id",add_pe=FALSE, add_domi=FALSE, missing_value=NULL,
				iteration_criteria=NULL,kinship_file=NULL,dmu_module=NULL,dmu_algorithm_code=NULL,dmu_prior_file=NULL,residual_cov_trait=NULL, 
				dmu_integer_n=NULL,output_path=NULL,output_name=NULL,
				residual_average=TRUE,cal_debv=FALSE, debv_pedigree_file=NULL, 
				cal_reliability=FALSE,debv_id=NULL,dmu_software_path=NULL,blupf90_software_path=NULL,show_message=NULL,
				IND_geno_file_name=NULL,IND_geno=NULL,SSBLUP_omega=NULL,plot_variance=FALSE,task_type=NULL){
	
				if(!is.null(dmu_ped_inbred_number)){self$pars$dmu_ped_inbred_number=dmu_ped_inbred_number}
				if(!is.null(phe)){self$pars$phe=phe}
				if(!is.null(ped)){self$pars$ped=ped}
				if(!is.null(geno)){self$pars$geno=geno}
				if(!is.null(phe_file)){self$pars$phe_file=phe_file}
				if(!is.null(phe_colnames)){self$pars$phe_colnames=phe_colnames}
				if(!is.null(analysis_model)){self$pars$analysis_model=analysis_model}
				if(!is.null(genetic_effect)){self$pars$genetic_effect=genetic_effect}
				if(!is.null(add_pe)){self$pars$add_pe=add_pe}
				if(!is.null(add_domi)){self$pars$add_domi=add_domi}
				if(!is.null(missing_value)){self$pars$missing_value=missing_value}
				if(!is.null(iteration_criteria)){self$pars$iteration_criteria=iteration_criteria}
				if(!is.null(kinship_file)){self$pars$kinship_file=kinship_file}
				if(!is.null(dmu_module)){self$pars$dmu_module=dmu_module}			
				if(!is.null(dmu_algorithm_code)){self$pars$dmu_algorithm_code=dmu_algorithm_code}
				if(!is.null(dmu_prior_file)){self$pars$dmu_prior_file=dmu_prior_file}
				if(!is.null(residual_cov_trait)){self$pars$residual_cov_trait=residual_cov_trait}
				if(!is.null(dmu_integer_n)){self$pars$dmu_integer_n=dmu_integer_n}
				if(!is.null(output_path)){self$pars$output_path=output_path}
				if(!is.null(output_name)){self$pars$output_name=output_name}
				if(!is.null(residual_average)){self$pars$residual_average=residual_average}
				if(!is.null(cal_debv)){self$pars$cal_debv=cal_debv}
				if(!is.null(debv_pedigree_file)){self$pars$debv_pedigree_file=debv_pedigree_file}
				if(!is.null(cal_reliability)){self$pars$cal_reliability=cal_reliability}
				if(!is.null(debv_id)){self$pars$debv_id=debv_id}
				if(!is.null(dmu_software_path)){self$pars$dmu_software_path=dmu_software_path}
				if(!is.null(blupf90_software_path)){self$pars$blupf90_software_path=blupf90_software_path}
				if(!is.null(IND_geno_file_name)){self$pars$IND_geno_file_name=IND_geno_file_name}
				if(!is.null(IND_geno)){self$pars$IND_geno=IND_geno}
				if(!is.null(SSBLUP_omega)){self$pars$SSBLUP_omega=SSBLUP_omega}
				if(!is.null(plot_variance)){self$pars$plot_variance=plot_variance}
				if(!is.null(task_type)){self$task_type=task_type}
				if(!is.null(show_message)){self$pars$show_message=show_message}
				
		final_effect_name=self$get_effect_list()
		
		if(self$task_type=="DMU"){
		
	       result=run_DMU(
					genetic_effect_name=self$pars$genetic_effect,   
					target_trait_name=final_effect_name$trait_list,
					fixed_effect_name=final_effect_name$fixed_list,  #list 
					random_effect_name=final_effect_name$random_list, #list		   
					covariate_effect_name=final_effect_name$covariate_list, #list
					random_regression_effect_name=NULL, #list, e.g. list(c("",""),c("",""))
					maternal_effect_name=NULL,    #list		   
					include_social_effect=FALSE, #whether including social effect in genetic evaluation
					group_effect_name=NULL,       #for social genetic effect evaluation 											
					integer_group_names=NULL,    #integer name for genetic group
					real_group_names=NULL,       #real name for genetic group
					
					phe_file=self$pars$phe_file,
					kinship_file=self$pars$kinship_file,
					prior_file=self$pars$dmu_prior_file,
					phe_col_names=self$pars$phe_colnames,					
					ped_inbred_number=self$pars$dmu_ped_inbred_number,         #whether consider inbreeding 
					analysis_model=self$pars$analysis_model,				
					included_permanent_effect=self$pars$add_pe, 
					included_dominance_effect=self$pars$add_domi,
					missing_value= self$pars$missing_value,
					iteration_criteria=self$pars$iteration_criteria,
					dmu_module=self$pars$dmu_module,
					dmu_algorithm_code=self$pars$dmu_algorithm_code,
					integer_n=self$pars$dmu_integer_n,  #number of integer varable in phenotype
					residual_cov_trait=self$pars$residual_cov_trait,  #resitricted-residual trait 
					output_result_path=self$pars$output_path,
					#cor_phe
					#genetic_effect_number=NULL,
					residual_average=self$pars$residual_average,			
					cal_debv=self$pars$cal_debv,     
					cal_reliability=self$pars$cal_reliability,
					debv_id=NULL,   
					#output_ebv_path=self$pars$output_path,
					#output_ebv_name=self$pars$output_name,			   
					DMU_software_path=self$pars$dmu_software_path,
					IND_geno_file_name="IND_geno.txt", #
					IND_geno=NULL,
					SSBLUP_omega=self$pars$SSBLUP_omega,
					plot_variance=self$pars$plot_variance,
					return_result=TRUE,
					show_message=self$pars$show_message
				  )
		
		self$ebv=result[[1]]
		dmu_vars=result$vars
		
		self$vars_se$invAI_mat=dmu_vars$inv_ai_mat
		self$vars_se$vars_mat=dmu_vars$vars_mat
		self$vars_se$h2=dmu_vars$h2
		self$vars_se$r2=list(Cor=dmu_vars$gen_cor,SE=dmu_vars$gen_cor_se)

		
		
		}else{ #for BLUPF90
			run_BLUPF90(
					genetic_effect_name=self$pars$genetic_effect,   
					target_trait_name=final_effect_name$trait_list,
					fixed_effect_name=final_effect_name$fixed_list,  #list 
					random_effect_name=final_effect_name$random_list, #list		   
					covariate_effect_name=final_effect_name$covariate_list, #list
					random_regression_effect_name=NULL, #list, e.g. list(c("",""),c("",""))
	
					phe_file=self$pars$phe_file,
					kinship_file=self$pars$kinship_file,					
					phe_col_names=self$pars$phe_colnames,										
					analysis_model=self$pars$analysis_model,				
					included_permanent_effect=self$pars$add_pe, 
					
					missing_value= self$pars$missing_value,
					
					#cor_phe
					#genetic_effect_number=NULL,			
					#output_ebv_path=self$pars$output_path,
					#output_ebv_name=self$pars$output_name,
					output_result_path=self$pars$output_path,	
					BLUPF90_software_path=self$pars$blupf90_software_path,
					plot_variance=self$pars$plot_variance
				  )
			
		
		}	
				
		self$r2=self$vars_se$r2
		self$h2=self$vars_se$h2
	
	},
	cal_se=function(expr,mat=NULL,invAI_mat=NULL){
		
		if(is.null(mat)&is.null(invAI_mat)){
		return(my_cal_se(expr,self$vars_se$vars_mat[,1],self$vars_se$invAI_mat))
		}else{
		return(my_cal_se(expr,mat,invAI_mat))		
		}
	},	
	
	print=function(...){
		
		 if(identical(self$task_type,"DMU")){
			cat(paste0("<BLUP::DMU>"," \n"))
		 }else{
			cat(paste0("<BLUP::BLUPF90>"," \n"))
		 }
		 super$print();
	
	}
	
	)
)	





