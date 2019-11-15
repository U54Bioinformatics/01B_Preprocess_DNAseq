
# Do copy number analysis.
betsy_run.py --network_png net.pdf --num_cores 20 \
  --input GenericCNVResults --input_file cnv \
  --input GenericCNVModelSelectionFile --input_file model.txt \
  --input ReferenceGenome --input_file Broad.hg19 \
  --input GTFGeneModel --input_file gtf.txt \
  --output CopyNumberAnalysis —output_file cnv.out \
  —dattr FACETSModelSelectionFile.model_selection=adhoc \
  —mattr facets_gbuild=hg19 \
  —mattr cn_header=tcn.em \
  —mattr total_cn_header=tcn.em \
  —mattr minor_cn_header=lcn.em \
  —mattr cn_header2=CNt \
  —mattr total_cn_header2=CNt \
  —mattr minor_cn_header2=B \
  —mattr discard_chrom_with_prefix=GL,JH,KB,KE,NC_,MT


# Do variant calls from FASTQ files with freebayes.

A=Mills_and_1000G_gold_standard
betsy_run.py --num_cores 8 --network_png call02.pdf --receipt call03.txt \
  --input FastqFolder --input_file proc21 \
  --input SampleGroupFile --input_file samp41.txt \
  --input ReferenceGenome \
  --input_file genomes/Broad.hg19/Homo_sapiens_assembly19.fasta \
  --output SimpleVariantMatrix --output_file call01.txt \
  --also_save_highest ManyCallerVCFFolders,call05 \
  --dattr BamFolderChunked.base_quality_recalibrated=yes \
  --dattr VCFFolder.caller=freebayes \
  --dattr VCFFolder.vartype=snp \
  --dattr SimpleVariantMatrix.caller_suite=single \
  --mattr wgs_or_wes=wgs \
  --mattr realign_known_sites1=v/$A.indels.b37.vcf.gz \
  --mattr realign_known_sites2=v/1000G_phase1.indels.b37.vcf.gz \
  --mattr recal_known_sites1=v/$A.indels.b37.vcf.gz \
  --mattr recal_known_sites2=v/1000G_phase1.indels.b37.vcf.gz \
  --mattr recal_known_sites3=v/dbsnp_138.b37.vcf.gz 

The input files you need:
- “proc21” is a directory that contains all your Fastq files
- “samp41.txt” is a file that you create that tells which fastq file goes with which sample.  It can be a tab-delimited text file, or an excel file.  I have attached an example file.  You can use the same format, but use your data instead.
- “Homo_sapiens_assembly19.fasta” is the reference genome.
- Mills_and_1000G_gold_standard.indels.b37.vcf.gz
 1000G_phase1.indels.b37.vcf.gz
 dbsnp_138.b37.vcf.gz
 These are files that come with the reference genome.  They tell you where the SNPs and indels are.  It is used for indel realignment.

There’s a line:
--mattr wgs_or_wes=wgs
If you are running with exome sequencing, please change to:
--mattr wgs_or_wes=wes


# Do variant calls from BAM files with freebayes.
A=Mills_and_1000G_gold_standard
betsy_run.py --num_cores 8 --network_png call02.pdf --receipt call03.txt \
  --input BamFolder --input_file proc21 \
  --dattr BamFolder.aligner=bwa_mem \
  --input SampleGroupFile --input_file samp41.txt \
  --input ReferenceGenome \
  --input_file genomes/Broad.hg19/Homo_sapiens_assembly19.fasta \
  --output SimpleVariantMatrix --output_file call01.txt \
  --also_save_highest ManyCallerVCFFolders,call05 \
  --dattr BamFolderChunked.base_quality_recalibrated=yes \
  --dattr VCFFolder.caller=freebayes \
  --dattr VCFFolder.vartype=snp \
  --dattr SimpleVariantMatrix.caller_suite=single \
  --mattr wgs_or_wes=wgs \
  --mattr realign_known_sites1=v/$A.indels.b37.vcf.gz \
  --mattr realign_known_sites2=v/1000G_phase1.indels.b37.vcf.gz \
  --mattr recal_known_sites1=v/$A.indels.b37.vcf.gz \
  --mattr recal_known_sites2=v/1000G_phase1.indels.b37.vcf.gz \
  --mattr recal_known_sites3=v/dbsnp_138.b37.vcf.gz 



# Do variant calls from FASTQ files with freebayes, no indel
# realignment.
betsy_run.py --num_cores 8 --network_png call02.pdf --receipt call03.txt \
 --input FastqFolder --input_file proc21 \
 --input SampleGroupFile --input_file samp41.txt \
 --input ReferenceGenome \
 --input_file genomes/Broad.hg19/Homo_sapiens_assembly19.fasta \
 --output SimpleVariantMatrix --output_file call01.txt \
 --also_save_highest ManyCallerVCFFolders,call05 \
 --dattr BamFolderChunked.base_quality_recalibrated=no \
 --dattr BamFolderChunked.indel_realigned=no \
 --dattr VCFFolder.caller=freebayes \
 --dattr VCFFolder.vartype=snp \
 --dattr SimpleVariantMatrix.caller_suite=single \
 --mattr wgs_or_wes=wgs


# Do variant calls from BAM files with freebayes, no indel
# realignment.
betsy_run.py --num_cores 8 --network_png call02.pdf --receipt call03.txt \
--input BamFolder --input_file proc21 \
--dattr BamFolder.aligner=bwa_mem \
--input SampleGroupFile --input_file samp41.txt \
--input ReferenceGenome \
--input_file genomes/Broad.hg19/Homo_sapiens_assembly19.fasta \
--output SimpleVariantMatrix --output_file call01.txt \
--also_save_highest ManyCallerVCFFolders,call05 \
--dattr BamFolderChunked.base_quality_recalibrated=no \
--dattr BamFolderChunked.indel_realigned=no \
--dattr VCFFolder.caller=freebayes \
--dattr VCFFolder.vartype=snp \
--dattr SimpleVariantMatrix.caller_suite=single \
--mattr wgs_or_wes=wgs



# Filter mutations for PyClone.
betsy_run.py --network_png pc02.pdf --num_cores 40 \
   --input SimpleVariantMatrix --input_file mut07.txt \
   --dattr SimpleVariantMatrix.with_coverage=yes \
   --input SequenzaResults --input_file cn01 \
   --input SequenzaModelSelectionFile --input_file mod02.xls \
   --output PyCloneMutationsFolder --output_file pc01 \
   --mattr cn_header=CNt \
   --mattr total_cn_header=CNt \
   --mattr minor_cn_header=B \
   --mattr filter_by_min_total_reads=20 \
   --mattr filter_by_min_alt_reads=5 \
   --mattr filter_by_min_vaf=0.05 \
   --mattr max_copynum_for_pyclone=8 \
   --mattr use_only_consistent_cn=yes

# Run PyClone
GTF=genomes/Broad.hg19.RefSeq.NM_only.no_isoforms.170703.gtf
betsy_run.py --network_png pc02.pdf --num_cores 40 \
  --input SimpleVariantMatrix --input_file mut07.txt \
  --dattr SimpleVariantMatrix.with_coverage=yes \
  --input SequenzaResults --input_file cn01 \
  --input SequenzaModelSelectionFile --input_file mod02.xls \
  --input GTFGeneModel --input_file $GTF \
  --output PyCloneAnalysis --output_file pc11 \
  --mattr cn_header=CNt \
  --mattr total_cn_header=CNt \
  --mattr minor_cn_header=B \
  --mattr filter_by_min_total_reads=20 \
  --mattr filter_by_min_alt_reads=5 \
  --mattr filter_by_min_vaf=0.05 \
  --mattr max_copynum_for_pyclone=8 \
  --mattr use_only_consistent_cn=yes
  
  
  # Run Sequenza.
betsy_run.py --network_png cn12.pdf --receipt cnv13.txt --num_cores 20 \
   --input BamFolder --input_file bam11 \
   --dattr BamFolder.has_read_groups=yes \
   --dattr BamFolder.sorted=coordinate \
   --input ReferenceGenome --input_file genomes/Broad.hg19 \
   --input NormalCancerFile --input_file nc01.xls \
   --output SequenzaResults --output_file cn11 \
   --mattr sequenza_assembly=hg19 \
   --mattr discard_chrom_with_prefix=GL,JH,KB,KE,NC_,MT
   

# Somatic variant calls from BAM files.
COSMIC=cosmic.v79.grch37.mutation_data.txt.gz
betsy_run.py --num_cores 20 --network_png call02.pdf --receipt call03.txt \
  --input BamFolder --input_file bam01 \
  --dattr BamFolder.sorted=coordinate \
  --dattr BamFolder.indexed=yes \
  --dattr BamFolder.has_read_groups=yes \
  --dattr BamFolder.duplicates_marked=yes \
  --dattr BamFolder.indel_realigned=yes \
  --dattr BamFolder.base_quality_recalibrated=yes \
  --dattr BamFolder.aligner=bwa_mem \
  --input ReferenceGenome \
  --input_file genomes/Broad.hg19/Homo_sapiens_assembly19.fasta \
  --input SampleGroupFile --input_file samp02.xlsx \
  --input NormalCancerFile --input_file nc01.xlsx \
  --output SimpleVariantMatrix --output_file call01.txt \
  --also_save_highest ManyCallerVCFFolders,call05 \
  --dattr SimpleVariantMatrix.duplicates_marked=yes \
  --dattr SimpleVariantMatrix.caller_suite=lance3 \
  --mattr wgs_or_wes=wgs \
  --mattr mutect_dbsnp_vcf=MuTect/dbsnp_132_b37.leftAligned.vcf \
  --mattr mutect_cosmic_vcf=MuTect/b37_cosmic_v54_120711.vcf \
  --dattr SimpleVariantMatrix.filtered_calls=yes \
  --mattr filter_by_min_total_reads=20 \
  --dattr SimpleVariantMatrix.annotated_with_annovar=yes \
  --mattr annovar_buildver=hg19 \
  --dattr SimpleVariantMatrix.annotated_with_snpeff=yes \
  --mattr snpeff_genome=GRCh37.75 \
  --dattr SimpleVariantMatrix.variant_annotations=cancer_cosmic \
  --mattr cancer_genes_file="008 Cancer Genes/cancer_genes.txt" \
  --mattr cosmic_variants_file="008 Cancer Genes/${COSMIC}" \
  --dattr SimpleVariantMatrix.with_coverage=yes 

