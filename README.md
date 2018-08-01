# DISH

# DISH
Direct Imputing Summary association statistics of HLA variants

This R script is designed to impute summary association statistics of HLA variants 
from SNP summary association statistics based on linkage disequilibria 
in European and Asian populations. 

!You can download the reference files from this link

https://drive.google.com/drive/folders/1DvyT4Ja_zRpS2cA0LkKIeXrXWQR2LlkC?usp=sharing

!Reference files and R script must be located in same directory

usage: Rscript DISH.r input_file input_type(T/P) hg_version(hg18/hg19) ethnicity(EUR/ASN) MAF_threshold stat_type(Z/T) output (lambda)

example) 
Rscript DISH.r EUR_sample.txt T hg18 EUR 0.005 T EUR_imputed.txt (lambda)
          
          or
         
Rscript DISH.r EUR_sample P hg18 EUR 0.005 T EUR_imputed.txt (lambda)



----------------------------------REQUIRED-------------------------------------

[input_file] must be a tab delimited file
-If [input_type] is T, the file extension must be included in [input_file]. If [input_type] is P, the file extension must be excluded.

[input_type] must be a "T" or "P"
-A T type input is a file that contains both SNP information and statistical values, 
 and a P type input means that both plink's .frq file and the user-defined file which include SNP with specific position and statistics values are used .

[hg_version] must be "hg18" or "hg19"

[ethnicity]  must be "EUR" or "ASN"
-EUR means European and ASN means Asian ethnicity

[MAF_threshold] must be a numeric value
-MAF_threshold must be >0 and <0.5

[stat_type] must be "Z" or "T"

[output] is an output file prefix

----------------------------------OPTIONAL-------------------------------------

[lambda] must be a numeric value. see the details in reference
