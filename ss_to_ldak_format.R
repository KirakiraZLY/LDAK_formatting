
library(tidyverse)
library(dplyr)
library(data.table) 

# About Input Argument
args <- commandArgs(trailingOnly = TRUE)

inputFile <- NULL
outputFile <- NULL
bimFile <- NULL
N <- NULL
# arg_out <- NULL

for (i in seq_along(args)) {
  if (args[i] == "--inputFile" && i < length(args)) {
    inputFile <- args[i + 1]
  } else if (args[i] == "--outputFile" && i < length(args)) {
    outputFile <- args[i + 1]
  } else if (args[i] == "--bfile" && i < length(args)) {
    bimFile <- args[i + 1]
  } else if (args[i] == "--N" && i < length(args)) {
    N <- args[i + 1]
  }
}

cat("================= Input =================", "\n")

#inputFile <- "/data/leyzha/proj1_prs_spa_psa/data/ukbb/categorical-20002-both_sexes-1477.tsv"
#outputFile <- "/data/leyzha/proj1_prs_spa_psa/data/ukbb/categorical-20002-both_sexes-1477.tsv.geno_pseudo"
#bimFile <- "/data/leyzha/proj1_prs_spa_psa/data/geno/pseudo/pseudo_1br_1ws"
#N <- 130000
bimFile <- paste(bimFile, ".bim", sep = "")

# Input File
#df_other1 <- read.table(file = "/data/leyzha/proj1_prs_spa_psa/data/finngen/sumstats/finngen_R12_M13_PSORIARTH", header = FALSE)
df_other1 <- read.table(file = inputFile, sep = "\t", header = TRUE)
print(head(df_other1))

# colnames(df_other1) <- unlist(df_other1[1,])
# df_other1 <- df_other1[-1,]


cat("Number of snp in SS: ", nrow(df_other1), "\n")
cat("================= Reading =================", "\n")

# df_other1 <- read_table2(file = inputFile, col_names = TRUE)

#if ("#chrom" %in% names(df_other1)) {
#  # Rename column
#  df_other1 <- df_other1 %>% rename(Chr = `#chrom`)
#}
#print(head(df_other1$Chr))


cols_to_remove <- colSums(is.na(df_other1)) == nrow(df_other1)
df_other1 <- df_other1[, !cols_to_remove]


df_other1 <- na.omit(df_other1)

df_other <- df_other1

print("###############################################################")
print(head(df_other))

#beta <- as.numeric(df_other$beta)
#df_other$or <- exp(beta)

# Column names could be
# Ref == A2, Alter = A1
names_to_replace <- c("Predictor",
                      "rsids",
                      "SNP", 
                      "#chrom",
                      "chrom",
                      "`#chrom`",
					  "#CHROM",
                      "pos",
                      "BP",
                      "Basepair", 
                      "GENPOS", 
                      "ref",
                      "alt",
                      "zscore",
                      "Direction",
                      "Stat",
                      "pval",
                      "beta",
                      "sebeta",
                      "or",
                      
                      ## Bolt
                      "SNP",
                      "CHR",
                      "BP",
                      "GENPOS",
                      "ALLELE1",
                      "ALLELE0",
                      "A1FREQ",
                      "F_MISS",
                      "BETA",
                      "SE",
                      "P_BOLT_LMM_INF",
                      "P_BOLT_LMM",
                      
                      ## LDAK
                      "Chromosome",
                      "Predictor",
                      "Basepair",
                      "A1",
                      "A2",
                      "Wald_Stat",
                      "Wald_P",
                      "Effect",
                      "SD",
                      "Effect_Liability",
                      "SD_Liability",
                      "A1_Mean",
                      "MAF",
					  
					  ## PGC_schiz
					  "CHROM",
					  "ID",
					  "POS",
					  "A1",
					  "A2",
					  "FCAS",
					  "FCON",
					  "IMPINFO",
					  "BETA",
					  "SE",
					  "PVAL",
					  "NCAS",
					  "NCON",
					  "NEFF",
					  
					  ## PGC_MDD
					  "MarkerName",
					  "A1",
					  "A2",
					  "Freq",
					  "LogOR",
					  "StdErrLogOR",
					  "P",
					  
					  ## PGC_ALZ2
					  "chromosome",
					  "base_pair_location",
					  "effect_allele",
					  "other_allele",
					  "beta",
					  "standard_error",
					  "p_value",
					  "N",
					  
					  ## PGC_ADHD
					  "Nca",
					  "Nco",
					  
            ## mvp dbGAP
					  "Allele1",
					  "Allele2",
					  
					  ## mvp AAA
					  "SNP_ID",
					  "SampleSize",
					  "EA",
					  "NEA",
					  "CHR:BP_hg19",
					  "Pvalue",
					  "EffectEstimate",
					  "Position",
					  "PValue",
					  
					  ## Regenie
					  "LOG10P",
					  "RSID",
					  "pval(-log10)",
					  "a1",
					  "a2",
					  "chr",
					  "rsid",
					  "A1_Mean",
					  "FRQ_A",
					  "FRQ_A_40463",
					  "af_alt",
					  "EAF",
					  
					  "beta_EUR",
					  "se_EUR",
					  "neglog10_pval_EUR",
					  "af_cases_EUR"
					  
)

replacement_rules <- c("Predictor" = "Predictor",
                       "rsids" = "rsids",
                       "#chrom" = "Chr",
                       "chrom" = "Chr",
                      "`#chrom`" = "Chr",
					  "#CHROM" = "Chr",
                       "pos" = "Pos",
                       "Basepair" = "Pos",
                       "alt" = "A1",
                       "ref" = "A2",
                       "zscore" = "Z",
                       "Direction" = "Direction",
                       "Stat" = "Stat",
                       "pval" = "P",
                       "beta" = "Beta",
                       "sebeta" = "SEbeta",
                       "or" = "OR",
                       
                       ## Bolt-inf
                       "SNP" = "rsids",
                       "CHR" = "Chr",
                       "BP" = "Pos",
                       "GENPOS" = "Pos",
                       "ALLELE1" = "A1",
                       "ALLELE0" = "A2",
                       "A1FREQ" = "freq",
                       "F_MISS" = "F_MISS",
                       "BETA" = "Beta",
                       "SE" = "SEbeta",
                       "P_BOLT_LMM_INF" = "P_BOLT_LMM_INF",
                       "P_BOLT_LMM" = "P",
                       
                       ## LDAK
                       "Chromosome" = "Chr",
                       "Predictor" = "Predictor",
                       "Basepair" = "Pos",
                       "A1" = "A1",
                       "A2" = "A2",
                       "Wald_Stat" = "Stat",
                       "Wald_P" = "P",
                       "Effect" = "Beta",
                       "SD" = "SD",
                       "Effect_Liability" = "Effect_Liability",
                       "SD_Liability" = "SD_Liability",
                       "A1_Mean" = "freq",
                       "MAF" = "MAF",
					   
					   ## PGC_schiz
					   "CHROM" = "Chr",
					   "ID" = "rsids",
					   "POS" = "Pos",
					   "A1" = "A1",
					   "A2" = "A2",
					   "FCAS" = "freq",
					   "FCON" = "FCON",
					   "IMPINFO" = "IMPINFO",
					   "BETA" = "Beta",
					   "SE" = "SEbeta",
					   "PVAL" = "P",
					   "NCAS" = "Ncas",
					   "NCON" = "Ncon",
				  	   "NEFF" = "n",
					  
					  ## PGC_MDD
					  "MarkerName" = "rsids",
					  "A1" = "A1",
					  "A2" = "A2",
					  "Freq" = "freq",
					  "LogOR" = "Beta",
					  "StdErrLogOR" = "SEbeta",
					  "P" = "P",
					  
					  ## PGC_ALZ2
					  "chromosome" = "Chr",
					  "base_pair_location" = "Pos",
					  "effect_allele" = "A1",
					  "other_allele" = "A2",
					  "beta" = "Beta",
					  "standard_error" = "SEbeta",
					  "p_value" = "P",
					  "N" = "n",
					  
					  ## PGC_ADHD
					  "Nca" = "Ncas",
					  "Nco" = "Ncon",
					  
					  ## mvp dbGAP
					  "Allele1" = "A1",
					  "Allele2" = "A2",
					  
					  ## mvp AAA
					  "SNP_ID" = "rsids",
					  "SampleSize" = "n",
					  "EA" = "A1",
					  "NEA" = "A2",
					  "CHR:BP_hg19" = "Predictor",
					  "Pvalue" = "P",
					  "EffectEstimate" = "Beta",
					  "Position" = "Pos",
					  "PValue" = "P",
					  
					  ## Regenie
					  "LOG10P" = "Log10P",
					  "RSID" = "rsids",
					  "pval(-log10)" = "Log10P",
					  "a1" = "A1",
					  "a2" = "A2",
					  "chr" = "Chr",
					  "rsid" = "rsids",
					  "A1_Mean" = "freq",
					  "FRQ_A" = "freq",
					  "FRQ_A_40463" = "freq" ,
					  "af_alt" = "freq",
					  "EAF" = "freq",
					  
					  "beta_EUR" = "Beta",
					  "se_EUR"= "SEbeta",
					  "neglog10_pval_EUR" = "Log10P" ,
					  "af_cases_EUR" = "freq"
                       
)


cat("================= Initializing =================", "\n")

df_other_alt <- df_other %>%
  rename_all(~ ifelse(. %in% names_to_replace, replacement_rules[.], .))
# df_other_alt <- mutate_at(df_other_alt, vars(c("Z", "Stat", "P", "Beta", "SEbeta", "OR")), as.numeric)

print(head(df_other))
print("###############################################################")
print(head(df_other_alt))
## I use rsids instead of predictor, for geno3 format
case_predictor <- c("rsids", "Chr", "Pos")
# cat("IF::::::::::::::::::::::::: ", all(case_predictor %in% colnames(df_other_alt)), "\n")

##################################### N
N_checker_case3 <- c("Ncas", "Ncon")
if(!is.null(N)){
	n <- N
} else if ("n" %in% colnames(df_other_alt)) { 
	n <- df_other_alt[1,]$n
} else if (sum(colnames(df_other_alt) %in% N_checker_case3) == 2) { 
	n <- df_other_alt[1,]$Ncas + df_other_alt[1,]$Ncon 
} else {
	stop("Check n existing in inputFile, or input --N")
}

cat("n = ", n, "\n")


###################################### Using Chr:Pos
if(FALSE)
{
	if (all(case_predictor %in% colnames(df_other_alt))) {
	  df_other_alt$Predictor <- paste(df_other_alt$Chr, ":", df_other_alt$Pos, sep = "")
	  df_other_alt$Predictor <- gsub(" ", "", df_other_alt$Predictor)
	  print("IF1")
	} else if (any(grepl("rsids", colnames(df_other_alt)) & !grepl("Predictor", colnames(df_other_alt)))){ #########################################################################################################################################
		df_other_alt$Predictor <- df_other_alt$rsids
		print("IF2")
	} else if(any(grepl("Predictor", colnames(df_other_alt)))) {
		df_other_alt$Predictor <- df_other_alt$Predictor
		print("IF3")
	}
}
#################################### Or using Predictor
cat("================= rsids ?= Predictor =================", "\n")
df_other_alt$Predictor <- df_other_alt$rsids
# df_other_alt$rsids <- df_other_alt$Predictor

df_other_alt$n <- n

print(head(df_other_alt))

cat("================= N added successfully =================", "\n")

cat("================= About P Value =================", "\n")
if ("Log10P" %in% colnames(df_other_alt) & !("P" %in% colnames(df_other_alt))) {
  df_other_alt$Log10P <- as.numeric(df_other_alt$Log10P)
  df_other_alt$P <- 10^(-df_other_alt$Log10P)
  print("Log10P converted to P successfully")
  print(head(df_other_alt))
}


cat("================= Z score =================", "\n")

case1 <- c("Z")
case2 <- c("Beta", "P")
case3 <- c("OR", "P")
case4 <- c("Beta", "SEbeta")
case5 <- c("OR", "SEbeta")
case6 <- c("Direction", "Stat")
# case7 <- c("Direction", "P")

if (sum(colnames(df_other_alt) %in% case1) == 1) {
  df_other_alt <- mutate_at(df_other_alt, vars(c("Z")), as.numeric)
  df_other_alt$Z <- df_other_alt$Z
  print("Case 1: directly by Z")
} else if (sum(colnames(df_other_alt) %in% case2) > 1){
  df_other_alt <- mutate_at(df_other_alt, vars(c("P", "Beta")), as.numeric)
  beta <- df_other_alt$Beta
  pval <- df_other_alt$P
  c <- -qnorm(pval/2)
  zscore <- ifelse(beta>0, c, -c)
  df_other_alt$Z <- zscore
  df_other_alt$SEbeta <- sqrt((df_other_alt$Beta ^ 2) / qchisq(df_other_alt$P, 1, lower.tail=F))
  cat("=== Check if there is any SEbeta == 0 ===", "\n")
  print(df_other_alt[df_other_alt$SEbeta == 0, ])
  print("Case 2: calculating Z by Beta and P")
} else if (sum(colnames(df_other_alt) %in% case3) > 1){
  df_other_alt <- mutate_at(df_other_alt, vars(c("P","OR")), as.numeric)
  OR <- df_other_alt$OR
  pval <- df_other_alt$P
  c <- qnorm(1-pval/2)
  zscore <- ifelse(OR>1, c, -c)
  df_other_alt$Z <- zscore
  df_other_alt$Beta <- log(df_other_alt$OR)
  print(head(df_other_alt))
  print("Case 3: calculating Z by OR and P")
} else if (sum(colnames(df_other_alt) %in% case4) > 1){
  df_other_alt <- mutate_at(df_other_alt, vars(c("Beta", "SEbeta")), as.numeric)
  beta <- df_other_alt$Beta
  se <- df_other_alt$SEbeta
  zscore <- beta/se
  df_other_alt$Z <- zscore
  print("Case 4: calculating Z by Beta and SEbeta")
} else if (sum(colnames(df_other_alt) %in% case5) > 1){
  df_other_alt <- mutate_at(df_other_alt, vars(c("SEbeta", "OR")), as.numeric)
  or <- df_other_alt$OR
  se <- df_other$SEbeta
  zscore <- log(or)/se
  df_other_alt$Z <- zscore
  print("Case 5: calculating Z by OR and SEbeta")
} else if (sum(colnames(df_other_alt) %in% case6) > 1){
  df_other_alt$Stat <- as.numeric(df_other_alt$Stat)
  df_other_alt$Direction <- as.numeric(df_other_alt$Direction)
  stat <- df_other_alt$Stat
  direc <- df_other_alt$Direction
  zscore <- ifelse(direc > 0, sqrt(stat), -sqrt(stat))
  df_other_alt$Z <- zscore
  print("Case 6: calculating Z by Direction and Stat")
} else{
  print("None is matched")
}

cat("================= Data Cleaning =================", "\n")

## Delete len(A1) > 1 and len(A2) > 1
# df_other_alt$Predictor <- df_other_alt$rsids
df_other_alt <- df_other_alt[nchar(df_other_alt$A1) <= 1, ]
df_other_alt <- df_other_alt[nchar(df_other_alt$A2) <= 1, ]
## Delete duplicate rsid
cat("================= Data Cleaning: duplicated rsids =================", "\n")
df_other_alt <- df_other_alt[!duplicated(df_other_alt$Predictor), ]

print(head(df_other_alt))
#df_other_alt$Beta <- log(df_other_alt$OR)

cat("================= Combining with bim file =================", "\n")
if(FALSE)
{
	cat("================= LDpred =================", "\n")
	head(df_other_alt)
	# Output File
	df_other_alt_ldpred <- df_other_alt[, c("rsids", "Chr", "Pos", "A1", "A2", "n","Z", "Beta", "SEbeta", "P")]
	colnames(df_other_alt_ldpred)[1] <- "Predictor"
	print(head(df_other_alt_ldpred))
	cat("================= SBayesRC =================", "\n")
	df_other_alt_cojo <- df_other_alt[, c("rsids", "A1", "A2", "Beta", "SEbeta", "P", "n")]
	df_other_alt$freq <- as.numeric(df_other_alt$freq)
	df_other_alt_cojo$freq <- df_other_alt$freq / 2
	df_other_alt_cojo <- df_other_alt_cojo[,c(1,2,3,8,4:7)]
	print(head(df_other_alt_cojo))
	colnames(df_other_alt_cojo)[1] <- "Predictor"
	df_other_alt <- df_other_alt[, c("Predictor", "A1", "A2", "n","Z")]

	cat("================= LDAK =================", "\n")
	print(head(df_other_alt))
}
cat(" ====== main summary statistics ====== ", "\n")
head(df_other_alt)

# Read bim
bim <- read.table(bimFile)
head(bim)
colnames(bim) <- c("Chr", "Predictor", "cM", "Pos", "A1", "A2")
bim_pred <- data.frame(bim$Predictor)
colnames(bim_pred) <- c("Predictor")

cat("N_bim_snp: ", nrow(bim_pred), " N_ss_snp: ", nrow(df_other_alt), "\n")
cat("class of Predictor in bim: ", class(bim_pred$Predictor), "  in ss: ", class(df_other_alt$Predictor), "\n")

# print(head(df_other_alt))
df_other_alt <- df_other_alt[grepl("^rs", df_other_alt$Predictor), ]
head(df_other_alt %>%
       inner_join(bim_pred, by = "Predictor"))
	   

df_other_alt_1 <- merge(df_other_alt, bim_pred, by = "Predictor", all.x = TRUE)
# df_other_alt_ldpred <- merge(df_other_alt_ldpred, bim_pred, by = "Predictor", all.x = TRUE)
#df_other_alt_ldpred <- df_other_alt_ldpred[grepl("^rs", df_other_alt_ldpred$Predictor), ]
# df_other_alt_cojo <- merge(df_other_alt_cojo, bim_pred, by = "Predictor", all.x = TRUE)
# colnames(df_other_alt_cojo) <- c("SNP", "A1","A2","freq", "b", "se","p","N")


# print(head(df_other_alt_1))
cat("============== Bim: ==============", "\n")
bim_pred <- bim
colnames(bim_pred) <- c("Chr","Predictor","cM","Pos","A1","A2")
df_other_alt_1$A1 <- toupper(df_other_alt_1$A1)
df_other_alt_1$A2 <- toupper(df_other_alt_1$A2)
# df_other_alt_ldpred$A1 <- toupper(df_other_alt_ldpred$A1)
# df_other_alt_ldpred$A2 <- toupper(df_other_alt_ldpred$A2)
# df_other_alt_cojo$A1 <- toupper(df_other_alt_cojo$A1)
# df_other_alt_cojo$A2 <- toupper(df_other_alt_cojo$A2)
bim_pred$A1 <- toupper(bim_pred$A1)
bim_pred$A2 <- toupper(bim_pred$A2)

print(head(bim_pred))

cat("============== After Formatting: ==============", "\n")
cat("============== LDAK: ==============", "\n")
df_other_alt_2 <- df_other_alt_1 %>%
  inner_join(bim_pred, by = "Predictor") %>% 
  filter((A1.x == A1.y & A2.x == A2.y) | (A1.x == A2.y & A2.x == A1.y)) %>%
  select(c("Predictor", "A1.x","A2.x","n","Z"))
colnames(df_other_alt_2) <- c("Predictor", "A1", "A2", "n","Z")
head(df_other_alt_2)

# cat("============== LDpred2: ==============", "\n")
# cat("class of A1 A2 in LDpred2: ", class(df_other_alt_ldpred$A1),class(df_other_alt_ldpred$A2), "  in bim: ", class(bim_pred$A1),class(bim_pred$A2), "\n")
# print(head(df_other_alt_ldpred))
# #df_other_alt_ldpred$A1 <- as.character(df_other_alt_ldpred$A1)
# #df_other_alt_ldpred$A2 <- as.character(df_other_alt_ldpred$A2)
# print(head(df_other_alt_ldpred %>%
# inner_join(bim_pred, by = "Predictor"))
# )
# print(head(df_other_alt_ldpred %>%
# inner_join(bim_pred, by = "Predictor")) %>% 
# filter((A1.x == A1.y & A2.x == A2.y) | (A1.x == A2.y & A2.x == A1.y))
# )
# df_other_alt_ldpred <- df_other_alt_ldpred %>%
# inner_join(bim_pred, by = "Predictor") %>% 
# filter((A1.x == A1.y & A2.x == A2.y) | (A1.x == A2.y & A2.x == A1.y)) %>%
# select(c("Predictor", "Chr.x", "Pos.x", "A1.x", "A2.x", "n", "Z", "Beta", "SEbeta", "P"))
# colnames(df_other_alt_ldpred) <- c("Predictor", "Chr", "Pos", "A1", "A2", "n", "Z", "Beta", "SEbeta", "P")
# #print(df_other_alt_ldpred[df_other_alt_ldpred$SEbeta == 0,])
# df_other_alt_ldpred <- df_other_alt_ldpred[df_other_alt_ldpred$SEbeta != 0,]
# head(df_other_alt_ldpred)

# cat("============== SBayesRC: ==============", "\n")
# cat("class of A1 A2 in SBayesRC: ", class(df_other_alt_cojo$A1),class(df_other_alt_cojo$A2), "  in bim: ", class(bim_pred$A1),class(bim_pred$A2), "\n")
# print(head(df_other_alt_cojo))
# #df_other_alt_cojo$A1 <- as.character(df_other_alt_cojo$A1)
# #df_other_alt_cojo$A2 <- as.character(df_other_alt_cojo$A2)
# colnames(df_other_alt_cojo)[1] <- "Predictor"
# head(df_other_alt_cojo %>%
# inner_join(bim_pred, by = "Predictor")
# )
# head(df_other_alt_cojo %>%
# inner_join(bim_pred, by = "Predictor")%>% 
# filter((A1.x == A1.y & A2.x == A2.y) | (A1.x == A2.y & A2.x == A1.y))
# )
# df_other_alt_cojo <- df_other_alt_cojo %>%
# inner_join(bim_pred, by = "Predictor") %>% 
# filter((A1.x == A1.y & A2.x == A2.y) | (A1.x == A2.y & A2.x == A1.y)) %>%
# select(c("Predictor", "A1.x", "A2.x", "freq", "b", "se", "p", "N"))
# colnames(df_other_alt_cojo) <- c("SNP", "A1", "A2", "freq", "b", "se", "p", "N")
# cat("============== SBayesRC After filtering: ==============", "\n")
# head(df_other_alt_cojo)
# #print(df_other_alt_cojo <- df_other_alt_cojo[df_other_alt_cojo$se == 0, ])
# df_other_alt_cojo <- df_other_alt_cojo[df_other_alt_cojo$se != 0, ]

outputFile_ldak <- paste(outputFile, ".ldak", sep="")
# outputFile_ldpred <- paste(outputFile, ".ldpred.ss", sep="")
# outputFile_cojo <- paste(outputFile, ".sbayesrc.cojo", sep="")
# cat("============== Write file to: ==============", "\n")
print(outputFile_ldak)
# print(outputFile_ldpred)
# print(outputFile_cojo)

write.table(df_other_alt_2, outputFile_ldak, row.names = FALSE, col.names = TRUE, quote = FALSE)
# write.table(df_other_alt_ldpred, outputFile_ldpred, row.names = FALSE, col.names = TRUE, quote = FALSE)
# write.table(df_other_alt_cojo, outputFile_cojo, row.names = FALSE, col.names = TRUE, quote = FALSE)
cat("Output Predictor: ", head(df_other_alt_2$Predictor), "\n")
cat("============== Mission Completed ==============", "\n")
cat("\n","\n","\n","\n","\n","\n")
