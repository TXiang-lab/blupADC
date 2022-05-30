create_DMU <- R6Class("DMU",
  public = list(
	parameters_list=NULL,
	trait_model=paste0(1," trait model"),
	analysis_model="PBLUP_A",
	EBV=NULL,
	DIR=NULL,
	h2=NULL,
	get_EBV=function(){
	 temp=private$result[[1]]	
	 temp[temp=="-9999"]=NA
	 self$EBV=temp 			
	},
	get_h2=function(){
	 self$h2=data.frame(private$result[[2]],stringsAsFactors=F) 	
	},
	get_DIR=function(){
	self$DIR=private$result[[3]]
	},
    initialize = function(phe_col_names=NULL,
					   target_trait_name=NULL,
					   fixed_effect_name=NULL, #列表
					   random_effect_name=NULL, #列表，不包括永久环境效应			   
					   covariate_effect_name=NULL, #列表
					   random_regression_effect_name=NULL, #列表，随机回归效应   list(c("",""),c("",""))
					   maternal_effect_name=NULL,    #列表		   
					   include_social_effect=FALSE,    #是否评估social effect
					   group_effect_name=NULL,  #评估social effect时, group的名称,自动生成group_phe,												
					   integer_group_names=NULL,        #整型group的名称
					   real_group_names=NULL,            #实型group的名称					   
					   ped_inbred_number=2, #是否考虑近交构建A
					   provided_effect_file_path=NULL, #各个性状的效应文件
					   provided_effect_file_name=NULL,
					   provided_DIR_file_path=NULL,
					   provided_DIR_file_name=NULL,
					   phe_path=NULL,
					   phe_name=NULL,
					   analysis_model=NULL,
					   genetic_effect_name="Id",
					   included_permanent_effect=FALSE, 
					   included_dominance_effect=FALSE,
					   missing_value= -9999,
					   iteration_criteria=1.0e-7,
					   relationship_name=NULL,
					   relationship_path=NULL,
					   dmu_module="dmuai",
					   dmu_algorithm_code=NULL,
					   provided_prior_file_path=NULL,
					   provided_prior_file_name=NULL,
					   integer_n=NULL,  #整型数目
					   residual_cov_trait=NULL,  #限定残差协方差的性状
					   output_result_path=base::getwd(),
					   #cor_phe
					   genetic_effect_number=NULL,
					   residual_average=TRUE,
					   debv_pedigree_path=NULL, #计算debv需要提供系谱
					   debv_pedigree_name=NULL, #计算debv需要提供系谱						
					   cal_debv=FALSE,      #是否计算debv
					   cal_reliability=FALSE,
					   debv_id=NULL,    #计算这些个体的校正表型
					   output_ebv_path=NULL,
					   output_ebv_name=NULL,			   
					   DMU_software_path=ifelse(as.character(Sys.info()["sysname"])=="Linux",system.file("extdata/bin", package = "blupADC"),system.file("extdata/bin_windows", package = "blupADC")),
					   IND_geno_file_name="IND_geno.txt", #
					   IND_geno=NULL,
					   SSBLUP_omega=0.05,
					   plot_variance=FALSE,
					   return_result=TRUE) {
                        self$parameters_list$phe_col_names=phe_col_names
                        self$parameters_list$target_trait_name=phe_col_names
                        self$parameters_list$fixed_effect_name=fixed_effect_name #列表
                        self$parameters_list$random_effect_name=random_effect_name #列表，不包括永久环境效应			   
                        self$parameters_list$covariate_effect_name=covariate_effect_name #列表
                        self$parameters_list$random_regression_effect_name=random_regression_effect_name #列表，随机回归效应   list(c("",""),c("",""))
                        self$parameters_list$maternal_effect_name=maternal_effect_name    #列表		   
                        self$parameters_list$include_social_effect=include_social_effect    #是否评估social effect
                        self$parameters_list$group_effect_name=group_effect_name  #评估social effect时, group的名称,自动生成group_phe,												
                        self$parameters_list$integer_group_names=integer_group_names        #整型group的名称
                        self$parameters_list$real_group_names=real_group_names            #实型group的名称					   
                        self$parameters_list$ped_inbred_number=ped_inbred_number #是否考虑近交构建A
                        self$parameters_list$provided_effect_file_path=provided_effect_file_path #各个性状的效应文件
                        self$parameters_list$provided_effect_file_name=provided_effect_file_name
                        self$parameters_list$provided_DIR_file_path=provided_DIR_file_path
                        self$parameters_list$provided_DIR_file_name=provided_DIR_file_name
                        self$parameters_list$phe_path=phe_path
                        self$parameters_list$phe_name=phe_name
                        self$parameters_list$analysis_model=analysis_model
                        self$parameters_list$genetic_effect_name=genetic_effect_name
                        self$parameters_list$included_permanent_effect=included_permanent_effect 
                        self$parameters_list$included_dominance_effect=included_dominance_effect
                        self$parameters_list$missing_value= missing_value
                        self$parameters_list$iteration_criteria=iteration_criteria
                        self$parameters_list$relationship_name=relationship_name
                        self$parameters_list$relationship_path=relationship_path
                        self$parameters_list$dmu_module=dmu_module
                        self$parameters_list$dmu_algorithm_code=dmu_algorithm_code
                        self$parameters_list$provided_prior_file_path=provided_prior_file_path
                        self$parameters_list$provided_prior_file_name=provided_prior_file_name
                        self$parameters_list$integer_n=integer_n  #整型数目
                        self$parameters_list$residual_cov_trait=residual_cov_trait  #限定残差协方差的性状
                        self$parameters_list$output_result_path=output_result_path
                        self$parameters_list$genetic_effect_number=genetic_effect_number
                        self$parameters_list$residual_average=residual_average
                        self$parameters_list$debv_pedigree_path=debv_pedigree_path #计算debv需要提供系谱
                        self$parameters_list$debv_pedigree_name=debv_pedigree_name #计算debv需要提供系谱						
                        self$parameters_list$cal_debv=cal_debv      #是否计算debv
                        self$parameters_list$cal_reliability=cal_reliability
                        self$parameters_list$debv_id=debv_id    #计算这些个体的校正表型
                        self$parameters_list$output_ebv_path=output_ebv_path
                        self$parameters_list$output_ebv_name=output_ebv_name			   
                        self$parameters_list$DMU_software_path=DMU_software_path
                        self$parameters_list$IND_geno_file_name=IND_geno_file_name #
                        self$parameters_list$IND_geno=IND_geno
                        self$parameters_list$SSBLUP_omega=SSBLUP_omega
                        self$parameters_list$plot_variance=plot_variance
                        self$parameters_list$return_result=return_result
    },	
	
	
	do_analysis =function(){
		self$trait_model=paste0(length(self$parameters_list$target_trait_name)," trait model")
		self$analysis_model=self$parameters_list$analysis_model
		private$result=run_DMU(
                        phe_col_names=self$parameters_list$phe_col_names,
                        target_trait_name=self$parameters_list$target_trait_name,
                        fixed_effect_name=self$parameters_list$fixed_effect_name, #列表
                        random_effect_name=self$parameters_list$random_effect_name, #列表，不包括永久环境效应			   
                        covariate_effect_name=self$parameters_list$covariate_effect_name, #列表
                        random_regression_effect_name=self$parameters_list$random_regression_effect_name, #列表，随机回归效应   list(c("",""),c("",""))
                        maternal_effect_name=self$parameters_list$maternal_effect_name,    #列表		   
                        include_social_effect=self$parameters_list$include_social_effect,    #是否评估social effect
                        group_effect_name=self$parameters_list$group_effect_name,  #评估social effect时, group的名称,自动生成group_phe,												
                        integer_group_names=self$parameters_list$integer_group_names,        #整型group的名称
                        real_group_names=self$parameters_list$real_group_names,            #实型group的名称					   
                        ped_inbred_number=self$parameters_list$ped_inbred_number, #是否考虑近交构建A
                        provided_effect_file_path=self$parameters_list$provided_effect_file_path, #各个性状的效应文件
                        provided_effect_file_name=self$parameters_list$provided_effect_file_name,
                        provided_DIR_file_path=self$parameters_list$provided_DIR_file_path,
                        provided_DIR_file_name=self$parameters_list$provided_DIR_file_name,
                        phe_path=self$parameters_list$phe_path,
                        phe_name=self$parameters_list$phe_name,
                        analysis_model=self$parameters_list$analysis_model,
                        genetic_effect_name=self$parameters_list$genetic_effect_name,
                        included_permanent_effect=self$parameters_list$included_permanent_effect, 
                        included_dominance_effect=self$parameters_list$included_dominance_effect,
                        missing_value= self$parameters_list$missing_value,
                        iteration_criteria=self$parameters_list$iteration_criteria,
                        relationship_name=self$parameters_list$relationship_name,
                        relationship_path=self$parameters_list$relationship_path,
                        dmu_module=self$parameters_list$dmu_module,
                        dmu_algorithm_code=self$parameters_list$dmu_algorithm_code,
                        provided_prior_file_path=self$parameters_list$provided_prior_file_path,
                        provided_prior_file_name=self$parameters_list$provided_prior_file_name,
                        integer_n=self$parameters_list$integer_n,  #整型数目
                        residual_cov_trait=self$parameters_list$residual_cov_trait,  #限定残差协方差的性状
                        output_result_path=self$parameters_list$output_result_path,
                        genetic_effect_number=self$parameters_list$genetic_effect_number,
                        residual_average=self$parameters_list$residual_average,
                        debv_pedigree_path=self$parameters_list$debv_pedigree_path, #计算debv需要提供系谱
                        debv_pedigree_name=self$parameters_list$debv_pedigree_name, #计算debv需要提供系谱						
                        cal_debv=self$parameters_list$cal_debv,      #是否计算debv
                        cal_reliability=self$parameters_list$cal_reliability,
                        debv_id=self$parameters_list$debv_id,    #计算这些个体的校正表型
                        output_ebv_path=self$parameters_lis$toutput_ebv_path,
                        output_ebv_name=self$parameters_list$output_ebv_name,			   
                        DMU_software_path=self$parameters_list$DMU_software_path,
                        IND_geno_file_name=self$parameters_list$IND_geno_file_name, #
                        IND_geno=self$parameters_list$IND_geno,
                        SSBLUP_omega=self$parameters_list$SSBLUP_omega,
                        plot_variance=self$parameters_list$plot_variance,
                       return_result=self$parameters_list$return_result
        )
	
    }
  ),
    private = list(
    result=NULL
  )
)
