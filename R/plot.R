#绘制方差组分的图
plot_dmu_blupf90_prior<-function(target_trait_name=NULL,random_effect_name=NULL,output_path=NULL,genetic_effect_name="Id"){
library(ggplot2)
library(RColorBrewer)
color_set=brewer.pal(9, "Set1")
color_set=color_set[c(1,3,2,4:9)]
unique_random=c(unique(do.call(c,random_effect_name)),"Residual")
group_color=color_set[1:length(unique_random)]
names(group_color)=unique_random

for(i in 1:length(target_trait_name)){

a=data.table::fread(paste0(target_trait_name[i],"_heritability_result.txt"),data.table=F)
a[,-1]=round(a[,-1],2)

a$prior_se_blanket=paste0(a$prior,"(",a$prior_se,")")
a$h2_se_blanket=paste0(a$h2,"(",a$h2_se,")")
a$h2=a$prior

colnames(a)[2]="Proportion of variance"
colnames(a)[4]="Variance"

b=reshape2::melt(a)
b=b[b[,"variable"]%in%c("Proportion of variance","Variance"),]
#b$variable=factor(b$variable,levels = c("Variance","Proportion of variance"))
b$prior_se_blanket=ifelse(b$variable=="Variance",b$prior_se_blanket,b$h2_se_blanket)

random_effect_order=c(genetic_effect_name,sort(setdiff(unique(b$Random_effect_name),genetic_effect_name)))
b$Random_effect_name=factor(b$Random_effect_name,levels =rev(random_effect_order))

p=ggplot2::ggplot(data = b, aes(x=variable, y=value, fill = Random_effect_name)) +
  geom_col(aes(fill = Random_effect_name), position = "dodge") +
  geom_text(aes(label =prior_se_blanket), position = position_dodge(0.9),size=2*length(target_trait_name),vjust = 0.5,hjust=-0.25,family = "serif")+
  labs(fill='Random effect')+  #修改图例
  scale_fill_manual(values =group_color)+  #设置颜色
  theme_bw() +ylab("Estimates(SE)")+
  coord_flip()+
  theme(axis.line = element_line(colour = "black"),
        axis.title.y=element_blank(), #去除x轴名称
        axis.text.y = element_text(hjust =0.5,face="bold",size=6.5*length(target_trait_name), family = "serif"),#调整x轴标签的内容显示
        axis.title.x=element_text(face="bold",size=6.5*length(target_trait_name), family = "serif"),
		axis.ticks.x=element_blank(), #去除y轴刻度线
	    axis.text.x=element_blank(),   #去除y轴刻度内容
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
	   legend.text=element_text(size=5*length(target_trait_name)), #调整图例文本字体大小
	   legend.title=element_text(size=5*length(target_trait_name)), #调整图例标题文本字体大小	   
        panel.background = element_blank())+ 
        scale_y_continuous(limits = c(0,max(b$value)*1.75))+
		theme(legend.position='bottom')+
		guides(fill = guide_legend(reverse = TRUE))
		

if(i==1){total=p}
if(i>=2){library(patchwork);total=total+p}
}

if(i>=2){
final=total+ patchwork::plot_annotation(tag_levels = 'A')&theme(plot.tag = element_text(size = 10*length(target_trait_name),face="bold"))
}else{
final=total
}


if(!is.null(output_path)){
if(length(target_trait_name)==1){
ggsave(final,filename = paste0(output_path,"/Variance_plot.png"),width=5*length(target_trait_name),height=2.7*length(target_trait_name))
}

if(length(target_trait_name)>=2){
ggsave(final,filename = paste0(output_path,"/Variance_plot.png"),width=6*length(target_trait_name),height=2*length(target_trait_name))
}
}else{
return(final)
}

}



detect_pca_plot<-function(input_data_numeric=NULL,Breed=NULL){
library(ggplot2)
pca<-data.frame(rownames(input_data_numeric),prcomp(input_data_numeric)$x[,1:2],stringsAsFactors=F);colnames(pca)=c("Id",paste0("PC",1:2))
PCA_Breed=data.frame(Breed,stringsAsFactors=F);colnames(PCA_Breed)=c("Id","Breed")
pca=merge(pca,PCA_Breed,by="Id")

if(nrow(pca)==0){stop("Can't match the records  between genotype data and provided breed information!!!")}
message("Start detect the breed of population by PCA......")

#detect outlier breed by kmeans
target_Breed_number=length(unique(PCA_Breed[,2]))
km=kmeans(pca[,2:3],target_Breed_number)
aaa <- data.frame(pca, km$cluster)

Data=NULL
for(i in 1:target_Breed_number){

tmp_Data <- pca[which(aaa$km.cluster == i),]
tmp_Data$Expeced_Breed=names(table(tmp_Data[,"Breed"]))[table(tmp_Data[,"Breed"]) == max(table(tmp_Data[,"Breed"]))]
Data=rbind(Data,tmp_Data)
}

outlier_Data=Data[Data[,"Breed"]!=Data[,"Expeced_Breed"],]

library(RColorBrewer)
color_set=brewer.pal(8, "Set1")
color_set=color_set[1:target_Breed_number]

origin_pca=ggplot2::ggplot(pca)+geom_point(aes(x=PC1,y=PC2,shape=Breed,color=Breed),size=5)+
				theme(legend.title = element_text(size=20),
				      legend.text = element_text(size = 15),
					  axis.title = element_text(size=20, color="black", face= "bold", vjust=0.5,hjust=0.5),
					  axis.text = element_text(size=15, color="black", face= "bold", vjust=0.5,hjust=0.5),
				     panel.grid.major = element_blank(),
                         panel.grid.minor = element_blank(),
                         panel.border = element_blank(),
	                    panel.background = element_blank(),
                        axis.line = element_line(colour = "black")						
					 )+
				   scale_shape_manual(values=c(18, 17,15))+ 
				    scale_colour_manual(values=color_set)  #设置颜色

kmeans_pca=ggplot2::ggplot(Data)+geom_point(aes(x=PC1,y=PC2,shape=Expeced_Breed,color=Expeced_Breed),size=5)+
				theme(legend.title = element_text(size=20),
				      legend.text = element_text(size = 15),
					  axis.title = element_text(size=20, color="black", face= "bold", vjust=0.5,hjust=0.5),
					  axis.text = element_text(size=15, color="black", face= "bold", vjust=0.5,hjust=0.5),
				     panel.grid.major = element_blank(),
                         panel.grid.minor = element_blank(),
                         panel.border = element_blank(),
	                    panel.background = element_blank(),
                        axis.line = element_line(colour = "black")						
					 )+ 
					 labs(shape='Expect Breed',color='Expect Breed')+			
					 scale_shape_manual(values=c(18, 17,15))+ 
				    scale_colour_manual(values=color_set)  #设置颜色


library(patchwork)
p=origin_pca+kmeans_pca+patchwork::plot_annotation(tag_levels = 'A')&theme(plot.tag = element_text(size = 20,face="bold"))
#ggsave(origin_pca+kmeans_pca+patchwork::plot_annotation(tag_levels = 'A')&theme(plot.tag = element_text(size = 20,face="bold")),
#  filename = paste0(output_plot_path,"/PCA_detect.png"),width=20,height=10)
return(list(outlier=outlier_Data[,-c(2:3)],plot=p))
}




test123<-function(){

#从这个例子就可以看出来 trace_pedigree 追溯的系谱代数计算方法了。DD18922905仍然被视作初代，尽管其和第一代是夫妻，但是由于其无父母记录，因此将其视为初代。
ped=data.frame(Offspring=c("DD20967904","DD18923611"),
                  Sire=c("DD18922905","DD17773903"),
			    Dam=c("DD18923611","DD17516512"))
				
				
ped2=data.frame(
 Offspring=c("DD20967904","DD18923611","DD17516512","DD16546520","DD12546527","DD10546538","DD10246533"),
 Sire=c("DD18922905","DD17773903","DD16526514","DD12546527","DD10546538","DD10246533","DD09246524"),
 Dam=c("DD18923611","DD17516512","DD16546520","DD13546529","DD11546519","DD10546579","DD08246517"))				


}



#计算谱系数据
get_family<-function(ped=NULL){



single_pedigree_result=plot_ped_single_pedigree(ped)
Ind=na.omit(unique(c(ped[,1],ped[,2],ped[,3])))
data_set=data.frame(Id=Ind,Family=NA)

Fn=single_pedigree_result[[3]]

get_related<-function(id,ped){
parents=as.vector(match_parents(id,ped)[[1]])
sibs=do.call(c,match_sibs(id,ped))
offsprings=match_offsprings(id,ped)
related_ids=na.omit(unique(c(id,parents,sibs,offsprings)))
ped=as.matrix(ped)
related_ids=na.omit(unique(as.vector(ped[ped[,1]%in%related_ids|ped[,2]%in%related_ids|ped[,3]%in%related_ids,])))
return(related_ids)
}

k=0
for(j in 1:length(Fn)){
id=Fn[j]
if(all(!is.na(data_set[,2]))){break} #当所有个体都有family记录时，退出for循环
if(!is.na(data_set[data_set[,1]%in%id,2])){next} #当前个体有family记录时，退出当前循环
k=k+1
related_ids=get_related(id,ped)
n=1
while(length(related_ids)>n){
n=length(related_ids)
for(i in 1:length(related_ids)){
related_ids=unique(na.omit(c(related_ids,get_related(related_ids[i],ped))))
}
}
data_set[data_set[,1]%in%related_ids,2]=k
}
return(data_set)
}


plot_ped_single_pedigree<-function(ped){   #针对单次的排序，挑选出初代,中代,终代， 

n=unique(na.omit(ped[,1]))

P=unique(na.omit(c(ped[,2],ped[,3])))   #假设缺失值是NA, 这些个体是父母，它们有后代

f1=P[!P %in% n] # 这些个体没有再往上的系谱了， 所以假定它们是初代

f2=P[P %in% n]  #这些个体可以往上追溯的父母，我们主要是将这些个体进行排序， 看谁是谁的曾祖代,祖代。。。。。。。

fn=n[!n %in% P] #这些个体不是父母，在父母那一栏里找不到这些个体， 所以这些个体是系谱最底层的那些个体

return(list(f1,f2,fn))

}

#获得系谱中所有个体的代数
get_generation<-function(ped){
ped[ped==0]=NA
single_pedigree_result=plot_ped_single_pedigree(ped)
Whole_ped=data.frame(Id=unique(do.call(c,single_pedigree_result)))
Whole_ped$Generation=NA
Whole_ped$X_score=NA  #父亲-1，母亲+1
Whole_ped$Y_score=NA  #亲本+1
Fn=single_pedigree_result[[3]] #无父母,  存在 个体无后代，但其全同胞有后代的情况，这些个体的代数有问题，后面将解决代数问题

#去除Fn中-全同胞有后代的个体 
fullsibs=NULL
for(i in Fn){
#将全同胞为父母的个体，从基础群中去除
if(length(match_sibs(i,ped)$full_sibs)>0){
if(sum(match_sibs(i,ped)$full_sibs%in%ped[,2])>0|sum(match_sibs(i,ped)$full_sibs%in%ped[,3])>0){
Fn=setdiff(Fn,i)
}
}
}
j=0
data=data.frame(Id=Fn,Generation=j,stringsAsFactors=F)

while(length(Fn)>=1&match_parents(Fn,ped)[[2]]=="Match"){

  j=j+1
  Fn=setdiff(na.omit(unique(as.vector(match_parents(Fn,ped)[[1]]))),Fn)
 
  for( tmp in Fn){
  Fn=unique(c(Fn,do.call(c,match_sibs(tmp,ped))))
  }
  if(length(Fn)==0){break}
  data=rbind(data,data.frame(Id=Fn,Generation=j))
}
return(data)
}



#获得个体及父母的 X轴，Y轴位置

get_X_Y_position_old<-function(X_score=0,Y_score=0,Fn=NULL,ped=NULL){
final_data=NULL
#根据nodes将 Fn排序
Fn=Fn[order(get_nodes(Fn,ped),decreasing = T)]

for(i in 1:length(Fn)){
Full_sibs_M=NULL
Full_sibs_F=NULL

id=Fn[i]
individual_data=data.frame(Id=id,
		                       X_score=X_score,
				             Y_score=Y_score,
						   X_score_end=X_score,
						   Y_score_end=Y_score+1,
						   Type="Individual")

if(i>1){
if(length(final_data[final_data[,1]%in%id,"Type"])>0){
if(final_data[final_data[,1]%in%id,"Type"]=="Sibs"){
next 
#已经在final_data里的个体，代表其位置信息已经确定，不需要再进行计算
}	
}
}					   
						   

if(match_parents(id,ped)[[2]]=="Complete_Unmatch"){
final_data=rbind(final_data,individual_data)
}else{

F=match_parents(id,ped)[[1]][1,1]  #父亲个体号
M=match_parents(id,ped)[[1]][1,2]  #母亲个体号

F_X_score=X_score-1
F_Y_score=Y_score+1

M_X_score=X_score+1
M_Y_score=Y_score+1

#判断个体有无全同胞，有全同胞的话，母亲X轴坐标向右移动(移动的长度等于全同胞的个数)
sibs_data=NULL
if(length(match_sibs(id,ped)$full_sibs)>0){
sibs_Id=match_sibs(id,ped)$full_sibs#全同胞的个体号
M_X_score=M_X_score+1*length(sibs_Id)

sibs_X_score=X_score+1:length(sibs_Id)
sibs_Y_score=Y_score

sibs_data=cbind(sibs_Id,sibs_X_score,sibs_Y_score)

sibs_data=data.frame(Id=sibs_Id,
		                X_score=sibs_X_score,
				       Y_score=sibs_Y_score,
					  X_score_end=sibs_X_score,
					  Y_score_end=sibs_Y_score+1,
					  Type="Sibs")

}

#从母亲出发，计算母亲上面分支的节点数
M_X_score=M_X_score+get_nodes(M,ped)*1 #母亲X轴向右平移
F_X_score=F_X_score-get_nodes(M,ped)*1 #父亲X轴向左平移


#从父亲出发，计算父亲上面分支的节点数
M_X_score=M_X_score+get_nodes(F,ped)*1 #母亲X轴向右平移
F_X_score=F_X_score-get_nodes(F,ped)*1 #父亲X轴向左平移

#统计父母的全同胞数目
Full_sibs_F=match_sibs(F,ped)$full_sibs
Full_sibs_M=match_sibs(M,ped)$full_sibs

parents_data=data.frame(Id=c(F,M),
		                   X_score=c(F_X_score,M_X_score),
				         Y_score=c(F_Y_score,M_Y_score),
					    X_score_end=c(M_X_score,M_X_score),
					    Y_score_end=c(M_Y_score,M_Y_score),
						Type=c("Father","Mother"))


final_data=rbind(final_data,sibs_data,parents_data,individual_data)
final_data=final_data[!is.na(final_data[,1]),]
#return(final_data)
}

X_score=max(final_data$X_score,final_data$X_score_end)+2*get_nodes(Fn[i+1],ped)
Y_score=0

#考虑母亲有多个全同胞，需要将下一个个体的X轴往右偏移
X_score=X_score+length(Full_sibs_M)

}

return(final_data)
}


get_X_Y_position<-function(X_score=0,Y_score=0,Fn=NULL,ped=NULL){
final_data=NULL
#根据nodes将 Fn排序

Fn=Fn[order(get_nodes(Fn,ped),decreasing = T)]

for(i in 1:length(Fn)){
Full_sibs_M=NULL
Full_sibs_F=NULL

id=Fn[i]
individual_data=data.frame(Id=id,
		                       X_score=X_score,
				             Y_score=Y_score,
						   X_score_end=X_score,
						   Y_score_end=Y_score+1,
						   Type="Individual")

if(i>1){
if(length(final_data[final_data[,1]%in%id,"Type"])>0){
if(final_data[final_data[,1]%in%id,"Type"]=="Sibs"){
next 
#已经在final_data里的个体，代表其位置信息已经确定，不需要再进行计算
}	
}
}					   
						   

if(match_parents(id,ped)[[2]]=="Complete_Unmatch"){
final_data=rbind(final_data,individual_data)
}else{

F=match_parents(id,ped)[[1]][1,1]  #父亲个体号
M=match_parents(id,ped)[[1]][1,2]  #母亲个体号

F_X_score=X_score-1
F_Y_score=Y_score+1

M_X_score=X_score+1
M_Y_score=Y_score+1

#判断个体有无全同胞，有全同胞的话，母亲X轴坐标向右移动(移动的长度等于全同胞的个数)
sibs_data=NULL
if(length(match_sibs(id,ped)$full_sibs)>0){
sibs_Id=match_sibs(id,ped)$full_sibs#全同胞的个体号


if(!is.na(id)&id%in%ped[,2]){
sibs_X_score=X_score-1:length(sibs_Id)  #当个体性别为公时，全同胞的X轴向左移动
F_X_score=F_X_score-1*length(sibs_Id)
}else {
sibs_X_score=X_score+1:length(sibs_Id)  #当个体性别为母时，全同胞的X轴向右移动
M_X_score=M_X_score+1*length(sibs_Id)
}

sibs_Y_score=Y_score

sibs_data=cbind(sibs_Id,sibs_X_score,sibs_Y_score)

sibs_data=data.frame(Id=sibs_Id,
		                X_score=sibs_X_score,
				       Y_score=sibs_Y_score,
					  X_score_end=sibs_X_score,
					  Y_score_end=sibs_Y_score+1,
					  Type="Sibs")

}

#从母亲出发，计算母亲上面分支的节点数
M_X_score=M_X_score+get_nodes(M,ped)*1 #母亲X轴向右平移
F_X_score=F_X_score-get_nodes(M,ped)*1 #父亲X轴向左平移


#从父亲出发，计算父亲上面分支的节点数
M_X_score=M_X_score+get_nodes(F,ped)*1 #母亲X轴向右平移
F_X_score=F_X_score-get_nodes(F,ped)*1 #父亲X轴向左平移

#统计父母的全同胞数目
Full_sibs_F=match_sibs(F,ped)$full_sibs
Full_sibs_M=match_sibs(M,ped)$full_sibs

parents_data=data.frame(Id=c(rep(F,length(F_X_score)),rep(M,length(M_X_score))),
		                   X_score=c(F_X_score,M_X_score),
				         Y_score=c(F_Y_score,M_Y_score),
					    X_score_end=c(M_X_score,M_X_score),
					    Y_score_end=c(M_Y_score,M_Y_score),
					    Type=c("Father","Mother"))


final_data=rbind(final_data,sibs_data,parents_data,individual_data)
final_data=final_data[!is.na(final_data[,1]),]
#return(final_data)
}

#X_score=max(final_data$X_score,final_data$X_score_end)+2*get_nodes(Fn[i+1],ped)
#id_max_X_score=final_data[which.max(final_data$X_score),1][1]

X_score=max(final_data$X_score,final_data$X_score_end)+4*get_nodes(Fn[i+1],ped)
#X_score=max(final_data$X_score,final_data$X_score_end)+2*get_nodes(id_max_X_score,ped)+2*get_nodes(Fn[i+1],ped)
if(get_nodes(Fn[i+1],ped)==0){X_score=X_score+1}
Y_score=0

#考虑母亲有多个全同胞，需要将下一个个体的X轴往右偏移
X_score=X_score+length(Full_sibs_M)

}

return(final_data)
}

#获得个体上面的代数(等价于系谱树的分支数)
get_nodes<-function(id,ped){
total_nodes=NULL
for(i in id){
nodes=0
while(match_parents(i,ped)[[2]]=="Match"){
nodes=nodes+1
i=na.omit(as.vector(match_parents(i,ped)[[1]]))
}
total_nodes=c(total_nodes,nodes)
}
return(total_nodes)
}

match_parents<-function(id,ped){ # ped 为三列向量（个体 父亲 母亲）， id 为一列向量 

parents=matrix(NA,ncol=2,nrow=length(id))

parents[,1]=ped[match(id,ped[,1]),2]
parents[,2]=ped[match(id,ped[,1]),3]
match_state=ifelse(sum(is.na(c(parents[,1],parents[,2])))==c(length(id)*2),"Complete_Unmatch","Match")
return(c(list(parents),list(match_state)))
}

#寻找个体的子代
match_offsprings<-function(id=NULL,ped=NULL){
offspring=NULL
if(!NA%in%id){
offspring=na.omit(unique(c(ped[ped[,2]%in%id,1],ped[ped[,3]%in%id,1])))
}
return(offspring)
}


#找寻个体的全同胞合半同胞
match_sibs<-function(id=NULL,ped=NULL){
id_F=ped[ped[,1]%in%id,2]
id_M=ped[ped[,1]%in%id,3]
if(length(id_F)==0){id_F=NA}
if(length(id_M)==0){id_M=NA}

sibs_F=NULL
sibs_M=NULL
if(!is.na(id_F)){sibs_F=ped[ped[,2]%in%id_F,1]}
if(!is.na(id_M)){sibs_M=ped[ped[,3]%in%id_M,1]}

full_sibs=NULL
half_sibs=NULL
sibs=unique(c(id,sibs_F,sibs_M))
sibs=setdiff(sibs,id)

if(!is.na(id_F)&!is.na(id_M)){ #当父母其中有一方为缺失值的时候，sibs必定为半同胞
for(i in sibs){
if(ped[ped[,1]%in%i,2]%in%id_F&ped[ped[,1]%in%i,3]%in%id_M){
full_sibs=c(full_sibs,i)
}else{
half_sibs=c(half_sibs,i)
}
}
}else{
half_sibs=c(half_sibs,sibs)
}
if(length(half_sibs)==0){half_sibs=NULL}
if(length(full_sibs)==0){full_sibs=NULL}

return(list(full_sibs=full_sibs,half_sibs=half_sibs))
}

#性别的形状

get_image<-function(shape_type=NULL,male_url=NULL,female_url=NULL){

if(!is.null(male_url)&!is.null(female_url)){
cat("Using user-provided figure to plot......")
male=male_url
female=female_url
}else{

if(shape_type==2){
female="https://qsmei-markdown.oss-cn-shanghai.aliyuncs.com/markdown-img/20210708210603.png"
male="https://qsmei-markdown.oss-cn-shanghai.aliyuncs.com/markdown-img/20210708210612.png"
}

if(shape_type==6){
male="https://qsmei-markdown.oss-cn-shanghai.aliyuncs.com/markdown-img/20210709155300.png"
female="https://qsmei-markdown.oss-cn-shanghai.aliyuncs.com/markdown-img/20210709155303.png"
}

if(shape_type==5){
female="https://qsmei-markdown.oss-cn-shanghai.aliyuncs.com/markdown-img/20210708213611.png"
male="https://qsmei-markdown.oss-cn-shanghai.aliyuncs.com/markdown-img/20210708213621.png"
}

if(shape_type==3){
male="https://qsmei-markdown.oss-cn-shanghai.aliyuncs.com/markdown-img/20210711221548.png"
female="https://qsmei-markdown.oss-cn-shanghai.aliyuncs.com/markdown-img/20210711221551.png"
}

if(shape_type==7){
male="https://qsmei-markdown.oss-cn-shanghai.aliyuncs.com/markdown-img/20210711220705.png"
female="https://qsmei-markdown.oss-cn-shanghai.aliyuncs.com/markdown-img/20210711220709.png"
}

if(shape_type==8){
male="https://qsmei-markdown.oss-cn-shanghai.aliyuncs.com/markdown-img/20210711223116.png"
female="https://qsmei-markdown.oss-cn-shanghai.aliyuncs.com/markdown-img/20210711223103.png"
}

if(shape_type==9){
male="https://qsmei-markdown.oss-cn-shanghai.aliyuncs.com/markdown-img/20210711223704.png"
female="https://qsmei-markdown.oss-cn-shanghai.aliyuncs.com/markdown-img/20210711223713.png"
}

}
return(list(male=male,female=female))
}

#ped=tmp_ped2
ggped<-function(input_pedigree=NULL,ind_text_size=4,ind_text_vjust=3,value_text_size=2,shape_size=8,plot_family=FALSE,trace_id=NULL,trace_sibs=FALSE,ind_sex=NULL,plot_inbred=FALSE,
				shape_type=4,male_url=NULL,female_url=NULL,show_curve=FALSE){

ped=input_pedigree

if(plot_inbred==TRUE){
K=cal_kinship(input_pedigree = ped,kinship_type = "P_A",return_result = TRUE)$P_A$A
Inbred=data.frame(Id=rownames(K),Inbred=diag(K)-1,stringsAsFactors=F)
rm(K);gc();
}

if(!is.null(trace_id)){
ped=trace_pedigree(input_pedigree=ped,trace_id=trace_id,trace_sibs=trace_sibs)[[1]][,1:3]
}

ped[ped==0]=NA
ped=ped[!(is.na(ped[,2])&is.na(ped[,3])),]
data_set=get_generation(ped)

Fn=data_set[data_set$Generation==0,1]

plot_pos=get_X_Y_position(X_score=0,Y_score=0,Fn=Fn,ped=ped) #起始位置

for(i in 1:nrow(data_set)){

Start_id=data_set[i,1]

if(Start_id %in% plot_pos[,1]){
if(match_parents(Start_id,ped)[[2]]=="Match"&all(plot_pos[plot_pos[,1]%in%Start_id,"Type"]!="Sibs")){
i_X_score=unique(plot_pos[plot_pos[,1]%in%Start_id,"X_score"])
i_Y_score=unique(plot_pos[plot_pos[,1]%in%Start_id,"Y_score"])
i_plot_pos=get_X_Y_position(X_score=i_X_score,Y_score=i_Y_score,Fn=Start_id,ped=ped)

plot_pos=rbind(plot_pos,i_plot_pos)
}

}else{
next  #已经在plot_pos里的个体，代表其位置信息已经确定，不需要再进行计算
}
plot_pos=plot_pos[!duplicated(paste0(plot_pos[,1],plot_pos[,4],plot_pos[,5])),]
}


if(!is.null(ind_sex)){

plot_pos$Sex=ind_sex[match(plot_pos[,1],ind_sex[,1]),2]
if(sum(is.na(plot_pos$Sex))>0){
stop("Please provide sex information for all individuals!")
}

}else{
plot_pos$Sex="Male"
plot_pos[plot_pos[,1]%in%ped[,3],"Sex"]="Female"
}

#统计同胞个体
Start_X_Y=paste0(plot_pos$X_score,"_",plot_pos$Y_score)
half_sibs_id=plot_pos$Id[duplicated(plot_pos$Id)&!duplicated(Start_X_Y)] #个体号相同，位置信息却不相同的id
half_sibs_id=unique(half_sibs_id)

#全同胞个体用虚线连接起来
total_half_sibs=NULL
if(length(half_sibs_id)>0){
for(i in half_sibs_id){
half_sibs_score=plot_pos[plot_pos$Id%in%i,]
half_sibs_score[1:(nrow(half_sibs_score)-1),4:5]=half_sibs_score[-1,2:3]
half_sibs_score=half_sibs_score[-nrow(half_sibs_score),]
total_half_sibs=rbind(total_half_sibs,half_sibs_score)
}

#check是否有重复的连线
status1=rep(NA,nrow(total_half_sibs))
status2=rep(NA,nrow(total_half_sibs))

for(i in 1:nrow(total_half_sibs)){
status1[i]=paste(total_half_sibs[i,1:5],collapse = "_")
status2[i]=paste(total_half_sibs[i,c(1,4,5,2,3)],collapse = "_")
}
total_half_sibs=total_half_sibs[!duplicated(status1%in%status2),]
}




library(ggplot2)
library(RColorBrewer)
color_set=brewer.pal(9, "Set1")

#将近交系数添加到数据集里
if(plot_inbred==TRUE){
plot_pos$Inbred_value=Inbred[match(plot_pos$Id,Inbred$Id),2]
}



if((shape_type!=1&shape_type!=4)|(!is.null(male_url)&!is.null(female_url))){
sex_url=get_image(shape_type,male_url=male_url,female_url=female_url)
male=sex_url$male
female=sex_url$female
plot_pos$Image=male
plot_pos[plot_pos$Sex=="Female","Image"]=female
}



if(plot_family==FALSE){
p=ggplot(plot_pos,aes(x=X_score,y=Y_score))+
        geom_segment(aes(x=X_score,
		              xend=X_score_end,
				    y=Y_score,
				    yend=Y_score_end))

if(show_curve==TRUE){
if(length(half_sibs_id)>=1){

p=p+geom_curve(aes(x =X_score,
                 y =Y_score,
                 xend=X_score_end, 
                 yend =Y_score_end),
				 data=total_half_sibs)

}}

if(shape_type==4){
p=p+geom_point(aes(x=X_score,y=Y_score,shape=Sex,color=Sex),size=shape_size)
}else if(shape_type==1){
p=p+geom_point(aes(x=X_score,y=Y_score,shape=Sex),size=shape_size)
}else{
p=p+ggimage::geom_image(size =shape_size/150,aes(image=Image))
}		
	
	p=p+geom_text(aes(label =Id),size=ind_text_size,vjust = ind_text_vjust,family = "serif")+
   theme(
        axis.title.y=element_blank(), #去除y轴名称
	   axis.title.x=element_blank(), #去除x轴名称
        axis.ticks.x=element_blank(), #去除x轴刻度线
	   axis.ticks.y=element_blank(), #去除y轴刻度线	
	   axis.text.x=element_blank(),   #去除y轴刻度内容
	   axis.text.y=element_blank(),   #去除y轴刻度内容
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
	   legend.position = "none", #去除所有图例内容   
        panel.background = element_blank())+		
		scale_y_continuous(limits = c(min(plot_pos$Y_score)-1/max(plot_pos$Y_score),max(plot_pos$Y_score)))+
	    scale_x_continuous(limits = c(min(plot_pos$X_score),max(plot_pos$X_score)))+
	   scale_shape_manual(values = c(16,15))+
	   scale_colour_manual(values =color_set)  #设置颜色
	   
if(plot_inbred==TRUE){
p=p+geom_text(aes(label =Inbred_value),size=value_text_size,vjust = 0,family = "serif")
}	   
	   
return(p)
}else{
family_set=get_family(ped)
plot_pos$Family=family_set[match(plot_pos$Id,family_set$Id),2]
total_half_sibs$Family=family_set[match(total_half_sibs$Id,family_set$Id),2]
Family_value=unique(plot_pos$Family)


tmp_plot_pos=plot_pos[plot_pos$Family%in%1,]
p=ggplot(tmp_plot_pos,aes(x=X_score,y=Y_score))+
        geom_segment(aes(x=X_score,
		              xend=X_score_end,
				    y=Y_score,
				    yend=Y_score_end))

if(show_curve==TRUE){
if(length(half_sibs_id)>=1){
tmp_total_half_sibs=total_half_sibs[total_half_sibs$Family%in%1,]
p=p+geom_curve(aes(x =X_score,
                 y =Y_score,
                 xend=X_score_end, 
                 yend =Y_score_end),
			data=tmp_total_half_sibs)
}}		

if(shape_type==4){
p=p+geom_point(aes(x=X_score,y=Y_score,shape=Sex,color=Sex),size=shape_size)
}else if(shape_type==1){
p=p+geom_point(aes(x=X_score,y=Y_score,shape=Sex),size=shape_size)
}else{
p=p+ggimage::geom_image(size =shape_size/150,aes(image=Image))
}
	p=p+geom_text(aes(label =Id),size=ind_text_size,vjust = ind_text_vjust,family = "serif")+
	  xlab("Family 1")+
   theme(
	   axis.title.x=element_text(face="bold",size=15, family = "serif"),
        axis.title.y=element_blank(), #去除y轴名称
	   #axis.title.x=element_blank(), #去除x轴名称
        axis.ticks.x=element_blank(), #去除x轴刻度线
	   axis.ticks.y=element_blank(), #去除y轴刻度线	
	   axis.text.x=element_blank(),   #去除x轴刻度内容
	   axis.text.y=element_blank(),   #去除y轴刻度内容
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
	   legend.position = "none", #去除所有图例内容   
        panel.background = element_blank())+		
		scale_y_continuous(limits = c(min(tmp_plot_pos$Y_score)-1/max(plot_pos$Y_score),max(tmp_plot_pos$Y_score)))+
	    scale_x_continuous(limits = c(min(tmp_plot_pos$X_score),max(tmp_plot_pos$X_score)))+
	   scale_shape_manual(values = c(16,15))+
	   scale_colour_manual(values =color_set)  #设置颜色

if(plot_inbred==TRUE){
p=p+geom_text(aes(label =Inbred_value),size=value_text_size,vjust = 0,family = "serif")
}
	   
sum_p=p	   
if(length(Family_value)>=2){

for(i in 2:length(Family_value)){

tmp_plot_pos=plot_pos[plot_pos$Family%in%i,]
p=ggplot(tmp_plot_pos,aes(x=X_score,y=Y_score))+
        geom_segment(aes(x=X_score,
		              xend=X_score_end,
				    y=Y_score,
				    yend=Y_score_end))

if(length(half_sibs_id)>=1){
tmp_total_half_sibs=total_half_sibs[total_half_sibs$Family%in%i,]
p=p+geom_curve(aes(x =X_score,
                 y =Y_score,
                 xend=X_score_end, 
                 yend =Y_score_end),
				 data=tmp_total_half_sibs)
}

if(shape_type==4){
p=p+geom_point(aes(x=X_score,y=Y_score,shape=Sex,color=Sex),size=shape_size)
}else if(shape_type==1){
p=p+geom_point(aes(x=X_score,y=Y_score,shape=Sex),size=shape_size)
}else{
p=p+ggimage::geom_image(size =shape_size/150,aes(image=Image))
}	
	p=p+geom_text(aes(label =Id),size=ind_text_size,vjust = ind_text_vjust,family = "serif")+
		xlab(paste0("Family ",i))+
   theme(
	   axis.title.x=element_text(face="bold",size=15, family = "serif"),
        axis.title.y=element_blank(), #去除y轴名称
	   #axis.title.x=element_blank(), #去除x轴名称
        axis.ticks.x=element_blank(), #去除x轴刻度线
	   axis.ticks.y=element_blank(), #去除y轴刻度线	
	   axis.text.x=element_blank(),   #去除x轴刻度内容
	   axis.text.y=element_blank(),   #去除y轴刻度内容
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
	   legend.position = "none", #去除所有图例内容   
        panel.background = element_blank())+		
		scale_y_continuous(limits = c(min(tmp_plot_pos$Y_score)-1/max(plot_pos$Y_score),max(tmp_plot_pos$Y_score)))+
	    scale_x_continuous(limits = c(min(tmp_plot_pos$X_score),max(tmp_plot_pos$X_score)))+
	   scale_shape_manual(values = c(16,15))+
	   scale_colour_manual(values =color_set)  #设置颜色
library(patchwork)	 
sum_p=sum_p+p
}}
sum_p=sum_p+ patchwork::plot_annotation(tag_levels = 'A')&theme(plot.tag = element_text(size = 20,face="bold",, family = "serif"))

return(sum_p)

}	   
}


show_ped_type<-function(){
library(ggplot2)
library(patchwork)
Theme=theme(
        axis.title.y=element_blank(), #去除y轴名称
	   axis.title.x=element_blank(), #去除x轴名称
        axis.ticks.x=element_blank(), #去除x轴刻度线
	   axis.ticks.y=element_blank(), #去除y轴刻度线	
	   axis.text.x=element_blank(),   #去除y轴刻度内容
	   axis.text.y=element_blank(),   #去除y轴刻度内容
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
	   legend.position = "none", #去除所有图例内容   
        panel.background = element_blank())	

plot_pos=data.frame(Id=c("Male","Female"),Sex=c("Male","Female"),X_score=c(1,1),Y_score=c(1,2.3))

four=ggplot(plot_pos,aes(x=X_score,y=Y_score))+
  geom_point(aes(x=X_score,y=Y_score,shape=Sex,color=Sex),size=10)+
  geom_text(aes(label =Id),size=4,vjust =3,family = "serif")+
  scale_shape_manual(values = c(16,15))+
	   scale_colour_manual(values =color_set)+Theme+
	   scale_y_continuous(limits = c(-1,3))
  
one=ggplot(plot_pos,aes(x=X_score,y=Y_score))+
  geom_point(aes(x=X_score,y=Y_score,shape=Sex),size=10)+
  geom_text(aes(label =Id),size=4,vjust =3,family = "serif")+
  scale_shape_manual(values = c(16,15))+
	   scale_colour_manual(values =color_set)+Theme+
	   scale_y_continuous(limits = c(-1,3))

for(shape_type in c(2:3)){
sex_url=get_image(shape_type=shape_type)
male=sex_url$male
female=sex_url$female
plot_pos$Image=male
plot_pos[plot_pos$Sex=="Female","Image"]=female  
p=ggplot(plot_pos,aes(x=X_score,y=Y_score))+
  geom_text(aes(label =Id),size=4,vjust =3,family = "serif")+
  ggimage::geom_image(size =.08,aes(image=Image))+Theme+
  scale_y_continuous(limits = c(-1,3))
one=one+p
}

one=one+four

for(shape_type in c(5:6)){
sex_url=get_image(shape_type=shape_type)
male=sex_url$male
female=sex_url$female
plot_pos$Image=male
plot_pos[plot_pos$Sex=="Female","Image"]=female  
p=ggplot(plot_pos,aes(x=X_score,y=Y_score))+
  geom_text(aes(label =Id),size=4,vjust =3,family = "serif")+
  ggimage::geom_image(size =.07,aes(image=Image))+Theme+
  scale_y_continuous(limits = c(-1,3))
one=one+p
} 

for(shape_type in c(7:9)){
sex_url=get_image(shape_type=shape_type)
male=sex_url$male
female=sex_url$female
plot_pos$Image=male
plot_pos[plot_pos$Sex=="Female","Image"]=female  
p=ggplot(plot_pos,aes(x=X_score,y=Y_score))+
  geom_text(aes(label =Id),size=4,vjust =3,family = "serif")+
  ggimage::geom_image(size =0.1,aes(image=Image))+Theme+
  scale_y_continuous(limits = c(-1,3))
one=one+p
} 
one=one+plot_annotation(tag_levels = list(paste0("Shape type=",1:9)))
return(one)
}

#ggsave(show_ped_type(),filename="c:/Users/26564/Desktop/test1.png",width = 10,height =7)		


#ped_plot=ggped(input_pedigree = blupADC::plot_pedigree,trace_id = c("121","125"),trace_sibs = T,shape_type =6,shape_size = 4)
#ggsave(ped_plot,filename="c:/Users/26564/Desktop/test-t.png",width = 10,height =7)



#ped_plot=ggped(input_pedigree = blupADC::plot_pedigree,trace_sibs = T,shape_type =4,
#               shape_size =8,trace_id = c("121","122","210","214","205","206","207"),plot_family=T)
#ggsave(ped_plot,filename="c:/Users/26564/Desktop/test-4.png",width =13,height =8)