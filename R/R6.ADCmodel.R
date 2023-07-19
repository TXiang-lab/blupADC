library(R6)
library(data.table)
ADCformula = R6Class("ADCformula",
    public = list(
		fixed=NULL,
		covariate=NULL,
		random=NULL,
		polyno=NULL,
	initialize=function(fixed=NULL,
						covariate=NULL,
						random=NULL,
						polyno=NULL){
		
				self$fixed=fixed
				self$covariate=covariate 
				self$random=random
				self$polyno=polyno	

		  },
		  
	print=function(...){}
		  
		  ),
	active=list(

			trait=function(fixed=self$fixed)all.vars(fixed)[1]  #dynamic changing with the value of self$fixed
														         #user can't modify trait directly
		    )		   
)

ADCmodel= R6Class("ADCmodel",

     public = list(
		formulas=NULL,
		id_name=NULL,
		dam_name=NULL,
		pe_name=NULL,
	initialize=function(fixed=NULL,
						covariate=NULL,
						random=NULL,
						polyno=NULL,
						id_name="id",
						dam_name="dam",
						pe_name=NULL,
						formulas=ADCformula$new()
						){
						
	self$id_name=id_name
	self$dam_name=dam_name	
	

	if("ADCformula"%in%class(fixed)){formulas=fixed;fixed=NULL};

	if("ADCformula"%in%class(formulas)&!is.null(formulas$fixed)){ # for null ADCformula object  
		if(is.null(formulas$fixed)&is.null(fixed)){
			self$formulas=NULL
		}else{
			self$formulas=c(NULL,formulas)	#make sure formulas is a list 
			names(self$formulas)=sapply(self$formulas, function(x)x$trait)

		}
	
	}else{ 
	
		#formulas=
		if(!is.null(fixed)){
			formulas$fixed=fixed
			formulas$covariate=covariate 
			formulas$random=random
			formulas$polyno=polyno
		}
		print(formulas)
		self$formulas=c(NULL,formulas)	#make sure formulas is a list 

		names(self$formulas)=sapply(self$formulas, function(x)x$trait)
	}

	},
		  
	 add_formula=function(fixed=NULL,
						  covariate=NULL,
					       random=NULL,
					       polyno=NULL){
	 
	 if("ADCformula"%in%class(fixed)){
		i_formula=fixed         # R6 didn't support function polymorphism, thus we used the first parameter as the default parameter
	 }else{
		i_formula=ADCformula$new(fixed,covariate,random,polyno)
	 }

	if("ADCformula"%in%class(self$formulas)&is.null(self$formulas$fixed)){ #consider that vars is null lmt_vars object
	 
		self$formulas=c(NULL,i_formula)
	 
	 }else{

		traits_name=sapply(self$formulas,function(x)x$trait) 	
	 
		 if(i_formula$trait%in%traits_name){ 
			cat(paste0("Trait:",i_formula$trait," was already exists in formulas, the old one will be coverd by the new one!","\n"))
			self$formulas[[match(i_formula$trait,traits_name)]]=i_formula
		 }else{
			self$formulas=c(self$formulas,i_formula)
		 }	
	}	
		names(self$formulas)=sapply(self$formulas, function(x)x$trait)
	
	    invisible(self)
	} ,

	 add_model=function(fixed=NULL,
						  covariate=NULL,
					       random=NULL,
					       polyno=NULL){
	 
	 if("ADCformula"%in%class(fixed)){
		i_formula=fixed         # R6 didn't support function polymorphism, thus we used the first parameter as the default parameter
	 }else{
		i_formula=ADCformula$new(fixed,covariate,random,polyno)
	 }

	if("ADCformula"%in%class(self$formulas)&is.null(self$formulas$fixed)){ #consider that vars is null lmt_vars object
	 
		self$formulas=c(NULL,i_formula)
	 
	 }else{

		traits_name=sapply(self$formulas,function(x)x$trait) 	
	 
		 if(i_formula$trait%in%traits_name){ 
			cat(paste0("Trait:",i_formula$trait," was already exists in formulas, the old one will be coverd by the new one!","\n"))
			self$formulas[[match(i_formula$trait,traits_name)]]=i_formula
		 }else{
			self$formulas=c(self$formulas,i_formula)
		 }	
	}	
		names(self$formulas)=sapply(self$formulas, function(x)x$trait)
	
	    invisible(self)
	} ,
	 rm_model=function(trait=NULL){
	 
		traits_name=sapply(self$formulas,function(x)x$trait) 	
	 
		 if(trait%in%traits_name){ 
			self$formulas=self$formulas[-c(match(trait,traits_name))]
		 }else{
		 	cat(paste0("Trait:",trait," is not exists in formulas,please check your input!","\n"))
		 }	
		names(self$formulas)=sapply(self$formulas, function(x)x$trait)
		invisible(self)
	 },
	 rm_formula=function(trait=NULL){
	 
		traits_name=sapply(self$formulas,function(x)x$trait) 	
	 
		 if(trait%in%traits_name){ 
			self$formulas=self$formulas[-c(match(trait,traits_name))]
		 }else{
		 	cat(paste0("Trait:",trait," is not exists in formulas,please check your input!","\n"))
		 }	
		names(self$formulas)=sapply(self$formulas, function(x)x$trait)
		invisible(self)
	 },


	get_effect_list=function(...){
			nTraits=length(self$formulas)
			trait_list=NULL
			fixed_list=NULL
			random_list=NULL
			covariate_list=NULL
			
			for(i in 1:nTraits){
				
				 i_self_formula=self$formulas[[i]]
				 trait_list=c(trait_list,list(i_self_formula$trait))

				 i_fixed=list(get_nested(i_self_formula$fixed)$effect_name)
				 i_random=list(get_nested(i_self_formula$random)$effect_name) 
				 i_covariate=list(get_nested(i_self_formula$covariate)$effect_name) 
				 
				 
				 
				 i_print=NULL
				 fixed_list=c(fixed_list,(i_fixed))
				 random_list=c(random_list,(i_random))
				 covariate_list=c(covariate_list,(i_covariate))
				 
				 # i_fixed=get_nested(i_self_formula$fixed)$effect_name 
				 # i_random=get_nested(i_self_formula$random)$effect_name 
				 # i_covariate=get_nested(i_self_formula$covariate)$effect_name 
				 
				 
				 
				 # i_print=NULL
				 # if(!is.null(i_fixed)){
				 
					# fixed_list=c(fixed_list,list(i_fixed))
					
				 # }

				 # if(!is.null(i_random)){
				 
					# random_list=c(random_list,list(i_random))
					
				 # }

				 # if(!is.null(i_covariate)){
				 
					# covariate_list=c(covariate_list,list(i_covariate))
					
				 # }		
				 
			 }

		return(list(trait_list=trait_list,fixed_list=fixed_list,random_list=random_list,covariate_list=covariate_list))

	},


	
	get_hiblup_format=function(phe_colnames){ #need users to provided colnames of phenotype
			nTraits=length(self$formulas)
			
			pos_trait_list=NULL
			pos_fixed_list=NULL
			pos_random_list=NULL
			pos_covariate_list=NULL
			
			for(i in 1:nTraits){
				
				 i_self_formula=self$formulas[[i]]
				 i_fixed=get_nested(i_self_formula$fixed)$effect_name 
				 i_random=get_nested(i_self_formula$random)$effect_name 
				 i_random=setdiff(i_random,"id") #hiblup will use id as the default animal effect
				 i_covariate=get_nested(i_self_formula$covariate)$effect_name 
				 
				 i_pos_trait=match(i_self_formula$trait,phe_colnames)
				 
				 if(is.na(i_pos_trait)){  #for trait name
					stop(paste0("Can't find trait:",i_self_formula$trait," in the colnames of phenotype data, please check your data carefully!"))
				 }
							 
				 if(!is.null(i_fixed)){  # for fixed effect 
				      i_pos_fixed=match(i_self_formula$fixed,phe_colnames)
					 if(sum(is.na(i_pos_fixed))>=1){
						stop(paste0("Can't find fixed effect:",paste(i_self_formula$fixed[is.na(i_pos_fixed)],collapse="&")," in the colnames of phenotype data, please check your data carefully!"))
					 }
				 }


				 if(!is.null(i_random)){  # for random effect 
				      i_pos_random=match(i_self_formula$random,phe_colnames)
					 if(sum(is.na(i_pos_random))>=1){
						stop(paste0("Can't find random effect:",paste(i_self_formula$random[is.na(i_pos_random)],collapse="&")," in the colnames of phenotype data, please check your data carefully!"))
					 }
				 }

				 if(!is.null(i_covariate)){  # for covariate effect 
				      i_pos_covariate=match(i_self_formula$covariate,phe_colnames)
					 if(sum(is.na(i_pos_covariate))>=1){
						stop(paste0("Can't find covariate effect:",paste(i_self_formula$covariate[is.na(i_pos_covariate)],collapse="&")," in the colnames of phenotype data, please check your data carefully!"))
					 }
				 }

				 
				 pos_trait_list=c(pos_trait_list,i_pos_trait)				 
				 pos_fixed_list=c(pos_fixed_list,paste(i_pos_fixed,collapse=","))
				 pos_random_list=c(pos_random_list,paste(i_pos_random,collapse=","))
				 pos_covariate_list=c(pos_covariate_list,paste(i_pos_covariate,collapse=","))	
				 
			 }

			 pos_trait_list=paste(pos_trait_list,collapse=" ")
			 pos_fixed_list=paste(pos_fixed_list,collapse=" ")
			 pos_random_list=paste(pos_random_list,collapse=" ")
			 pos_covariate_list=paste(pos_covariate_list,collapse=" ")			 

		return(list(pos_trait_list=pos_trait_list,pos_fixed_list=pos_fixed_list,pos_random_list=pos_random_list,pos_covariate_list=pos_covariate_list))

	},
		
	print=function(...){
	
	
		 nTraits=length(self$formulas)
		
		 if(nTraits>=1){

			 cat(paste0("<ADCmodel::",nTraits," trait",ifelse(nTraits>1,"s",""),"> \n"))
			
			#format(formula)
			
			 for(i in 1:nTraits){
				 i_self_formula=self$formulas[[i]]
				 
				 i_fixed=get_nested(i_self_formula$fixed)$effect_name 
				 i_random=get_nested(i_self_formula$random)$effect_name 
				 i_covariate=get_nested(i_self_formula$covariate)$effect_name 
				 
				 
				 
				 i_print=NULL
				 if(!is.null(i_fixed)){
				 
					i_print=c(i_print,paste0(i_fixed,"(fixed)"))
					
				 }

				 if(!is.null(i_random)){
				 
					i_print=c(i_print,paste0(i_random,"(random)"))
					
				 }

				 if(!is.null(i_covariate)){
				 
					i_print=c(i_print,paste0(i_covariate,"(covariate)"))
					
				 }		

				 #using + join differnt part(exclude trait)
				 i_print=paste0(i_print,collapse="+")
				 
				 i_final_print=paste0("----Trait:",i_self_formula$trait,"~",i_print) #add trait~
				 
				 cat(paste0(i_final_print,"\n"))
				# formula_dataframe[i,1]=i_self_formula$fixed  #formual can't convert characters directly
				# formula_dataframe[i,2]=i_self_formula$random 
				# formula_dataframe[i,3]=i_self_formula$covariate 
			 }
		
		
		
		 }
	 }
	 
)
)


 # m1= ADCformula$new(fixed=tr1~f11+f12+f13,covariate=~c1c1,random=~id)
 # m2= ADCformula$new(fixed=tr2~f11+f12+f13,covariate=~c1c1,random=~id)

 # t1=ADCmodel$new()
 # t1$add_formula(m1)$add_formula(m2)
# t1$get_effect_list()

 # t2=ADCmodel$new(fixed=tr1~f11+f12+f13,covariate=~c1c1,random=~id)
# t2$get_effect_list()
