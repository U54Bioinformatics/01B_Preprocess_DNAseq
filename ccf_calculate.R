# functions from publications:
# 1. Gopal, Priyanka, et al. "Clonal selection confers distinct evolutionary trajectories in BRAF-driven cancers." Nature communications 10.1 (2019): 1-14.
# 2. Zapata, Luis, et al. "Signatures of positive selection reveal a universal role of chromatin modifiers as cancer driver genes." Scientific reports 7.1 (2017): 1-15.


CCF <- function(sample.mutations, VAF = NULL, ploidy = NULL, CCF_CNV = NULL, purity = NULL, correct=TRUE){
  if (is.atomic(sample.mutations)) {
    sample.mutations <- data.frame(x = sample.mutations)
  } 
  
  if (!is.null(VAF)){
    sample.mutations <- assign.columns(sample.mutations, VAF, "VAF")
  }
  if (!is.null(ploidy)){
    sample.mutations <- assign.columns(sample.mutations, ploidy, "ploidy")
  }
  if (!is.null(CCF_CNV)){
    sample.mutations <- assign.columns(sample.mutations, CCF_CNV, "CCF_CNV")
  }
  if (!is.null(purity)){
    sample.mutations <- assign.columns(sample.mutations, purity, "purity")
  }
  
  # make it not sensitive to lower/upper case in column names
  original.col.names <- colnames(sample.mutations)
  num.col <- ncol(sample.mutations)
  colnames(sample.mutations) <-  tolower(colnames(sample.mutations))
  
  # check if BAF column is there
  if ( 'vaf' %in% colnames(sample.mutations) ){
    if  (!is.numeric(sample.mutations$vaf)){
      stop("VAF column is not numeric!")
    }            
  } else {
    stop("There is no mandatory VAF column!")
  }
  
  if ( 'ploidy' %in% colnames(sample.mutations) ){
    if  (!is.numeric(sample.mutations$ploidy)){
      stop("Ploidy column is not numeric!")
    }   
    if ( 'ccf_cnv' %in% colnames(sample.mutations) ){
      if  (!is.numeric(sample.mutations$ccf_cnv)){
        stop("CCF_CNV column is not numeric!")
      }
      if ('purity' %in% colnames(sample.mutations) ) {
        # calculate CCF as ploidy is 2 
        sample.mutations$CCF <- ccfPloidy(sample.mutations$vaf, sample.mutations$ploidy, sample.mutations$ccf_cnv, purity=sample.mutations$purity)
      } else {
        # calculate CCF! there is baf, ploidy and ccf of cnv
        sample.mutations$CCF <- ccfPloidy(sample.mutations$vaf, sample.mutations$ploidy, sample.mutations$ccf_cnv)
      }
    } else {
      if ('purity' %in% colnames(sample.mutations) ) {
        # calculate CCF as ploidy is 2 
        sample.mutations$CCF <- ccfPloidy(sample.mutations$vaf, sample.mutations$ploidy,  purity=sample.mutations$purity)
      } else {
        # calculate CCF! there is baf, ploidy and ccf of cnv
        sample.mutations$CCF <- ccfPloidy(sample.mutations$vaf, sample.mutations$ploidy)  
      }
    }           
  } else {
    if ('purity' %in% colnames(sample.mutations) ) {
      # calculate CCF as ploidy is 2 
      sample.mutations$CCF <- ccfPloidy(sample.mutations$vaf, purity=sample.mutations$purity)
    } else {
      # calculate CCF as ploidy is 2 
      sample.mutations$CCF <- ccfPloidy(sample.mutations$vaf)
    }
    
  }
  
  if (correct){
    sample.mutations <- ccfCorrection(sample.mutations)
  }
  
  colnames(sample.mutations)[1:num.col] <- original.col.names
  
  sample.mutations
  
}


ccfPloidy <- function (vaf, ploidy = 2, ccf_cnv = 1, purity = 1) {
  if (sum(is.na(ploidy))){
    ploidy[is.na(ploidy)] <- 2
  }
  if (sum(is.na(ccf_cnv))){
    ccf_cnv[is.na(ccf_cnv)] <- 1
  }  
  if (sum(is.na(purity))){
    purity[is.na(purity)] <- 1
  }  
  ccf <- ((2 + (ploidy-2)*ccf_cnv)*vaf)/purity    
  return(ccf)
}

ccfCorrection <- function(sample.mutations){  
  if  (!'purity' %in% colnames(sample.mutations)){
    
    # correct BAF between 0.5 and 0.6 and diploid
    if ( 'ploidy' %in% colnames(sample.mutations) ){
      condition <- (sample.mutations$vaf > 0.5 & sample.mutations$vaf <=0.6) & sample.mutations$ploidy == 2   
    } else {
      condition <- (sample.mutations$vaf > 0.5 & sample.mutations$vaf <=0.6 ) 
    }
    if (sum(condition,  na.rm = T)) {
      condition[is.na(condition)] <- FALSE
      sample.mutations[condition, ]$CCF <- 1        
    }
    
    # correct BAF between 0.6 and 1 and diploid
    if ( 'ploidy' %in% colnames(sample.mutations) ){    
      condition <- sample.mutations$vaf > 0.6  &  (sample.mutations$ploidy == 2 | is.na(sample.mutations$ploidy ))
    } else {
      condition <- sample.mutations$vaf > 0.6 
    }
    if (sum(condition,  na.rm = T)) {
      condition[is.na(condition)] <- FALSE
      sample.mutations[condition,]$CCF <- ccfPloidy(sample.mutations[condition ,]$vaf, ploidy=1)
    }
    
    # correct ploidy != 2 and ccf >1
    if ( 'ploidy' %in% colnames(sample.mutations) ){   
      condition <- sample.mutations$CCF > 1  & (sample.mutations$ploidy != 2   | is.na(sample.mutations$ploidy ))
    } else {
      condition <- sample.mutations$CCF > 1      
    }
    if (sum(condition,  na.rm = T)) {
      condition[is.na(condition)] <- FALSE
      sample.mutations[condition, ]$CCF <- 1
    }
    
  } else {
    if (sum(is.na(sample.mutations$purity))){
      sample.mutations[is.na(sample.mutations$purity),'purity'] <- 1
    } 
    
    # correct BAF between 0.5 and 0.6 and diploid
    if ( 'ploidy' %in% colnames(sample.mutations) ){
      condition <- (sample.mutations$vaf > 0.5 & sample.mutations$vaf <=0.6) & (sample.mutations$ploidy == 2 | is.na(sample.mutations$ploidy ))  
    } else {
      condition <- (sample.mutations$vaf > 0.5 & sample.mutations$vaf <=0.6 ) 
    }
    if (sum(condition,  na.rm = T)) {
      condition[is.na(condition)] <- FALSE
      sample.mutations[condition, ]$CCF <-  min( (sample.mutations[condition, ]$vaf*2  / sample.mutations[condition, ]$purity ) , 1)    
    }
    
    # correct BAF between 0.6 and 1 and diploid
    if ( 'ploidy' %in% colnames(sample.mutations) ){    
      condition <- sample.mutations$CCF > 1.2  &  (sample.mutations$ploidy == 2  | is.na(sample.mutations$ploidy ))
    } else {
      condition <- sample.mutations$CCF > 1.2 
    }
    if (sum(condition,  na.rm = T)) {
      condition[is.na(condition)] <- FALSE
      sample.mutations[condition,]$CCF <- ccfPloidy(sample.mutations[condition ,]$vaf, ploidy=1, purity=sample.mutations[condition ,]$purity)
    }
    
    # correct ploidy != 2 and ccf >1
    if ( 'ploidy' %in% colnames(sample.mutations) ){   
      condition <- sample.mutations$CCF > 1  #& sample.mutations$ploidy != 2   
    } else {
      condition <- sample.mutations$CCF > 1      
    }
    if (sum(condition,  na.rm = T)) {
      condition[is.na(condition)] <- FALSE
      sample.mutations[condition, ]$CCF <- 1
    }
    
  }
  
  
  sample.mutations
}

get_ccf <- function(sample_to_use, vaf, cna, purity){
  
  vaf_cols <- c(colnames(vaf)[1:8], sample_to_use)
  vaf_dt <- vaf[, ..vaf_cols]
  colnames(vaf_dt)[9] <- 'VAF'
  
  
  cna_cols <- c('Gene Symbol', 'chrom', 'chromStart', 'chromEnd', sample_to_use)
  cna_dt <- cna[, ..cna_cols]
  colnames(cna_dt)[length(colnames(cna_dt))] <- 'CN'
  
  pur <- purity[Sample == sample_to_use, Purity]
  
  # match vaf with cna by coordinate, since there could be different CNs within one gene
  # matching by gene name will remove 1/2 of mutations
  
  #vaf_dt[, pos_key:=paste0(Chrom, ":", floor(Pos/100000))]
  #cna_dt[, pos_key:=paste0(chrom, ":", floor(chromStart/100000))]
  #ccf_input <- merge(vaf_dt, cna_dt, by.x = 'pos_key', by.y = 'pos_key')
  ccf_input <- merge(vaf_dt, cna_dt, by.x = 'AnnovarGene.refGene', by.y = 'Gene Symbol')
  #check A2ML1 
  
  ccf_input <- ccf_input[Chrom == chrom & Pos<=chromEnd & Pos>=chromStart, ]
  
  
  ccf_input[, purity:=pur]     
  ccf_input[, seg_mean:=log2(CN)]
  ccf_input[, ploidy:=(2*(2^seg_mean))]
  
  ccf_input[, CCF_CNV:=VAF*(((purity * ploidy) + 2*(1-purity))/purity)]  
  ccf_input[ploidy==2, CCF_CNV:=1]
  return (data.table(CCF(data.frame(ccf_input, stringsAsFactors = F))))
  
}



# examples of use


source('calculate_CCF.R')
## preprocess your vaf file by Mutation X (Annotation+Sample) data.frame or data.table
## For example, mine is 
#   Chrom    Pos Ref             Alt AnnovarFunc.refGene AnnovarGene.refGene AnnovarGeneDetail.refGene
#     1  15118   A               G      ncRNA_intronic              WASH7P                          
#     1  16257   G               C      ncRNA_intronic              WASH7P                          
#     1 139382   G GCCCCTCCAGGCCCA                                                                  
   AnnovarExonicFunc.refGene BRST004_A BRST004_AAA
#                               0.123       0.149
#                               0.086       0.242
#                               0.001       0.002


vaf <- fread('cDriver/150x//BRST004_unfiltered_mutations.table.txt')
select_cols <- c(colnames(vaf)[c(1:8)], colnames(vaf)[grepl("Coverage", colnames(vaf))])
vaf <- vaf[, ..select_cols]
colnames(vaf) <- gsub("Coverage|______", "", colnames(vaf))

vaf <- melt.data.table(vaf, id.vars = colnames(vaf)[1:8], value.factor = F, variable.factor = F)
vaf <- vaf[variable %in% samples, ]
vaf[, value:=tstrsplit(value, "/", keep = 3)]
vaf[, value:=as.numeric(value)]
vaf <- vaf[!is.na(value), ]
vaf <- dcast(vaf, Chrom+Pos+Ref+Alt+AnnovarFunc.refGene+AnnovarGene.refGene+AnnovarGeneDetail.refGene+AnnovarExonicFunc.refGene~variable)

##copy number data
cna <- fread('cDriver/150x/cna-BRST004.150x.WES/copy_number.by_gene.txt')

#purity data
purity <- fread('cDriver/150x/cna-BRST004.150x.WES/tumor_purity.txt')

#cancer cell frequency ccf
ccf_dt <- data.table()

for (i in 1:length(samples)){
  ccf <- get_ccf(samples[i], vaf, cna, purity)
  sample <- paste0("SNV_CCF_", samples[i])
  if (i == 1){
    ccf_dt <- ccf
    colnames(ccf_dt)[ncol(ccf_dt)] <- sample
  }else{
    
    ccf_dt[, (sample):=ccf[, CCF]]
  }
}


#ccf_dt contains the final output

