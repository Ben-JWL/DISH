##############################################################################################
####
#### DISH (Direct Imputing Summary association statistics of HLA variants):
####
#### This R script is designed to impute summary association statistics of HLA variants 
#### from SNP summary association statistics based on linkage disequilibria 
#### in European and Asian populations. 
####
##############################################################################################
####
#### V1 was written by Jiwoo Lim and Kwangwoo Kim (kkim@khu.ac.kr) on July-16, 2018
#### V2 was modified by Jiwoo Lim (jiwooim0922@gmail.com) on July-30, 2018
####
##############################################################################################
####
#### usage: Rscript DiSH.r input_file input_type hg_version ethnicity MAF_threshold stat_type output (lambda)
####
##############################################################################################
####       Argument check start       ########################################################
start = Sys.time()
#Read arguments
print("reading arguments...")
args = commandArgs(TRUE)
argsLen <- length(args)
if (argsLen > 7) stop('error: You added too many arguments - usage: Rscript DISH.v2.r input_file MAF ethnicity(European/Asian) output')
if (argsLen < 6) stop('error: You need more auguments. check the usage - usage: Rscript DISH.v2.r input_file MAF ethnicity(European/Asian) output')

study_map_file = args[1]
if(!file.exists(study_map_file)){stop('error: The incorrect path of your input file')}

input_type = args[2]
if(input_type != "T" & input_type != "P") stop('error: The incorrect type of your input file')

hg_ver = args[3]
if(tolower(hg_ver)!="hg19" & tolower(hg_ver)!="hg18"){stop('error: The incorrect genomic version was used')}

ethnicity = args[4]
if(ethnicity != "EUR" & ethnicity != "ASN") stop('error: Ethnicity should be European or Asian')

maf = as.numeric(args[5])
if(is.na(maf) |!is.numeric(maf)) stop('error: MAF must be a number')   
if(maf <= 0 | maf >= 0.5) stop('error: MAF should be >0 and <0.5')

stat = args[6]
if(stat == "Z"){
	print("Z scores are used in this analysis") 
}else if(stat == "T"){ 
	print("T statistics are used in this analysis")
}else stop("error: Please define the statistics type")

output = args[7]
if(file.create(output) != TRUE) stop('error: The incorrect path of your output file or no permission to write an output in the specified path')

lambda = outfile <- if (argsLen == 7 ) 0.15 else as.numeric(args[8])
if(is.na(lambda) |!is.numeric(lambda)) stop("error: Lambda must be a number")

####       Argument check finished       ####################################################

####       Read inputs                   ####################################################
## Check & Read the reference file 
if(ethnicity=="EUR") {
  print("reading the European reference...")
  Ref = readRDS("European_hg18_hg19.Rds")
} else if(ethnicity == "ASN") {
  print("reading the Asian reference...")
  Ref = readRDS("Asian_hg18_hg19.Rds")
}
Ref_map = Ref[[2]]
Ref = Ref[[1]]

if(hg_ver == "hg18"){
	names(Ref_map)[names(Ref_map)=="POS18"] <- c("SNP_pos")
}else if(hg_ver == "hg19"){
	names(Ref_map)[names(Ref_map)=="POS19"] <- c("SNP_pos")
}

## Check & Read the user defined input file
print("reading your input file...")
if(input_type=="T"){
	study_map = complete.cases(read.table(study_map_file, header=T, sep="\t", stringsAsFactors=F))
}else if(input_type=="P"){
	print("Merge your .frq file with your statistics information")
	Ffile = read.table(pasete0(study_map_file,".frq",collapse=NULL), header=T,stringsAsFactors=F)
	Tfile = read.table(pasete0(study_map_file,".txt",collapse=NULL), header=T,stringsAsFactors=F)
	Mfile = complete.cases(merge(x=Ffile,y=Tfile,by='SNP'))
	Mfile = Mfile[order(Mfile$SNP_pos),]
	study_map = Mfile[,c("SNP","SNP_pos","A1","A2","STAT")]
	if(stat=="Z"){
		colnames(study_map) <- c("SNP_id","SNP_pos","Effect_allele","Non_effect_allele","Z")
	}else if(stat=="T"){
		colnames(study_map) <- c("SNP_id","SNP_pos","Effect_allele","Non_effect_allele","T")
	}
}


## Verify the statistics type and input file 
if(stat=="Z"){
	if(ncol(study_map) != 5) stop("error: incorrect column number of your input file; check if the input was tab-delimited")
	if(length(grep("TRUE", (colnames(study_map) == c("SNP_id", "SNP_pos", "Effect_allele", "Non_effect_allele", "Z")))) !=5) stop("error: incorrect field names in your input file")
}else if(stat=="T"){			## In T statistics method, we convert T statistics to Z score
	if(ncol(study_map) != 5) stop("error: incorrect column number of your input file; check if the input was tab-delimited")
	if(length(grep("TRUE", (colnames(study_map) == c("SNP_id", "SNP_pos", "Effect_allele", "Non_effect_allele", "T")))) !=5) stop("error: incorrect field names in your input file")
	for(i in 1:nrow(study_map)){
		if(study_map$T[i]>=0){
			study_map$Z_from_T[i] <- abs(qnorm((2*pt(-abs(study_map$T[i]), df=2, log.p=T))-log(2),lower.tail=F,log.p=T))
		}else if(study_map$T[i]<0){
			study_map$Z_from_T[i] <- -abs(qnorm((2*pt(-abs(study_map$T[i]), df=2, log.p=T))-log(2),lower.tail=F,log.p=T))
		}
	}	
}else{ 
	stop("Please define correct statistics and check your input file")
}

####      filter using a maf threshold    ###################################################
print(paste("filtering out SNPs with allele frequency < ",maf," or > ",1-maf," ...",sep=""))
filtered_idx = which(Ref_map$maf >= maf & Ref_map$maf <= 1 - maf)
sigma = Ref[filtered_idx,filtered_idx]
sigma = as.matrix(sigma)
Ref_map = Ref_map[filtered_idx,]


####      Partition typed and imputed markers through hg version  ###########################

temp <- merge(x=Ref_map, y=study_map, by='SNP_pos')
temp_pos <- temp[(temp$Major==temp$Effect_allele & temp$Minor==temp$Non_effect_allele) | (temp$Major==temp$Non_effect_allele & temp$Minor==temp$Effect_allele),"SNP_pos"]
temp_pos <- sort(temp_pos)
			
print("partitioning typed and imputed markers...")
typed_idx <- which(Ref_map$SNP_pos %in% temp_pos)
if(length(typed_idx)>=1){
	sigma_tt = sigma[typed_idx,typed_idx]
	sigma_tt = sigma_tt + lambda*diag(ncol(sigma_tt))
	sigma_it = sigma[-typed_idx,typed_idx]
	Ref_map_t = Ref_map[typed_idx,]
	Ref_map_i = Ref_map[-typed_idx,]
}else{
	stop("error: No matched SNP. Please check your input file's SNP position")
}

####      check the direction of Z scores based on alleles      ####
print("changing the direction of Z scores if effect alleles is not minor alleles ... ")
typed_in_both_idx <- which(study_map$SNP_pos %in% temp_pos)
if(length(typed_in_both_idx)>=1){
	study_map = study_map[typed_in_both_idx,]
	if(stat=="Z"){
		Zt = rep(NA,length(study_map$Z))
		for (i in 1:length(Zt)) {
			if(as.character(study_map$Effect_allele[i]) != as.character(Ref_map_t$Minor[which(Ref_map_t$SNP_pos == study_map$SNP_pos[i])]) ) {
				Zt[which(Ref_map_t$SNP_pos == study_map$SNP_pos[i])] = -study_map$Z[i]
			} else {
				Zt[which(Ref_map_t$SNP_pos == study_map$SNP_pos[i])] = study_map$Z[i]
			}
		}
	}else if(stat=="T"){
		Zt = rep(NA,length(study_map$Z_from_T))
		for (i in 1:length(Zt)) {
			if(as.character(study_map$Effect_allele[i]) != as.character(Ref_map_t$Minor[which(Ref_map_t$SNP_pos == study_map$SNP_pos[i])]) ) {
				Zt[which(Ref_map_t$SNP_pos == study_map$SNP_pos[i])] = -study_map$Z_from_T[i]
			} else {
				Zt[which(Ref_map_t$SNP_pos == study_map$SNP_pos[i])] = study_map$Z_from_T[i]
			}
		}
	}	
}
	
Zt = as.matrix(Zt)

        
print("calculating Z scores at untyped markers ... ")
print("it will take several minutes ... ")

# make an inverse matrix of sigma_tt
print(paste("det=",det(sigma_tt)))
sigma_tt_inv = solve(sigma_tt)
# impute Z scores at untyped markers
#sigma_it=as.matrix(sigma_it)
weight = sigma_it %*% sigma_tt_inv
Zi = weight %*% Zt

         
# calculate r2pred values ... 
print("calculating r2pred values ... ")
print("it will take several minutes ... ")
var = matrix(NA,nrow(weight),1)
for(i in 1:nrow(weight)) {
  var[i] = sum(weight[i,] %o% weight[i,] * sigma_tt )
}

         
#### print output ####
print("printing your output. Bye!")
results = cbind(as.character(Ref_map_i$SNP_id), as.character(Ref_map_i$SNP_pos), as.character(Ref_map_i$Minor), as.character(Ref_map_i$Majo), Zi, var, 2*pnorm(-abs(Zi)))
colnames(results) = c("Marker_id", "Marker_pos", "Effect_allele", "Non_effect_allele", "Imputed_Z",  "r2pred", "imputed_P")
write.table(results,output, row.names=F, col.names=T, quote=F, sep="\t")
end = Sys.time()
print(end)
print(end-start)
rm(list=ls())
