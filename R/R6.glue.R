#make sure the software path is correct no matter which version of blupADC is used
blupADC_software_path<-function(){

#for older version 
if(packageVersion("blupADC")<"1.1.0"){

software_path=ifelse(as.character(Sys.info()["sysname"])=="Linux",system.file("extdata/bin_linux", package = "blupADC"),
			  ifelse(as.character(Sys.info()["sysname"])=="Windows",system.file("extdata/bin_windows", package = "blupADC"),
													                   system.file("extdata/bin_mac", package = "blupADC")))
}else{
#for new version


software_path=ifelse(as.character(Sys.info()["sysname"])=="Linux",system.file("extdata/bin_linux", package = "blupSUP"),
			  ifelse(as.character(Sys.info()["sysname"])=="Windows",system.file("extdata/bin_windows", package = "blupSUP"),
													                   system.file("extdata/bin_mac", package = "blupSUP")))
}


return(software_path)

}


check_blupADC<-function(){

#for older version 
if(packageVersion("blupADC")>="1.1.0"){

	if(!requireNamespace("blupSUP", quietly = TRUE)){
	
		cat(paste0('For version of blupADC >= 1.1.0, User has to install package:blupSUP, install_github("TXiang-lab/blupSUP") !','\n'))
		stop("Please make sure to install the blupSUP package, which only needs to be installed once")
	
	}


}

}

#function for nested situation
get_nested<-function(effect){
	
	
	if(!is.null(effect)){
	
		effect_expression=attr(terms(effect), which = "term.labels") # get the name of effect (include expression format)	

		is_response=(attr(terms(effect),which = "response")==1) #if response variable is exist

		if(is_response){ # there is reponse variable in the left of ~
			effect_name=all.vars(effect)[-1] # get the name of effect(not include outside of  bracket)
		}else{
			effect_name=all.vars(effect) # get the name of effect(not include outside of  bracket)
		}
			
		in_nested_effect=NULL    #inside the bracket
		out_nested_effect=NULL   #outside the bracket

		if(!identical(effect_name,effect_expression)){ #  nested effect exists

			nested_pos=which(effect_name!=effect_expression)

			if(is_response){
				tmp=attr(terms(effect), which = "variables")[-c(1:2)]  # get all variable outside and inside the bracket
			}else{
				tmp=attr(terms(effect), which = "variables")[-1] # get all variable outside and inside the bracket
			}

			in_nested_effect=effect_name[nested_pos]
			
			for(pos in nested_pos){
			
				out_nested_effect=c(out_nested_effect,as.character(tmp[[pos]])[1]) #correspond to the nested_random_effect
			
			}

		}
	}else{
	
		effect_name=in_nested_effect=out_nested_effect=NULL
	
	}
	
	
	return(list(effect_name=effect_name,
				in_nested_effect=in_nested_effect,
				out_nested_effect=out_nested_effect))
}