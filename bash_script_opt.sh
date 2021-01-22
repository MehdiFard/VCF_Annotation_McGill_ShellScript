#!/bin/bash
START=`date +%s`
echo ""

# Making directories 
mkdir -p intermediate_sh
mkdir -p log_sh
mkdir -p results_sh

# Defining variables
INPUT="data/info_files/metadata.csv"
VCF_DIR="data/vcfs/"
DIR="intermediate_sh"
LOG="log_sh"
RES="results_sh"

# Parsing Info file
echo "# Parsing Info file ..."
cut -d"," -f1 $INPUT | sed '1d' | while read line; do echo $VCF_DIR${line}; done > $DIR/vcfnames.list # VCFs names
cut -d"," -f1 --complement $INPUT | awk 'NR==1 || NR==2' > $DIR/qcthr.csv # Thresholds

# A while loop to assign each threshold value to its name (assigns value of each field to field name)
nr=1
while IFS=, read -r -a ary; do
    if (( nr == 1 )); then
        col_names=("${ary[@]}")
    else
        for (( i = 0; i < ${#ary[@]}; i++ )); do
            printf -v "${col_names[i]}" "${ary[i]}"
        done
    fi
    (( nr++ ))
done < $DIR/qcthr.csv

# Merging VCF files
echo "# Merging VCF files ..."
PicardCommandLine MergeVcfs I=$DIR/vcfnames.list O=$DIR/merged.vcf.gz 2> $LOG/merge.log
# rm $DIR/vcfnames.list $DIR/qcthr.csv

# Quality control
echo "# Quality control ..."
echo "Thresholds:"
echo "Minor Allel Frequency: $MAF"
echo "Hardy-Weinberg test: $HWE"
echo "Missing genotype per site: $Missing_per_Mrks"
echo "Missing genotype per individual: $Missing_per_Inds"
vcftools --gzvcf $DIR/merged.vcf.gz --missing-indv --out $DIR/miss_data_per_inds 2> $LOG/miss_per_inds.log
awk '$5 > $Missing_per_Inds {print $1}' $DIR/miss_data_per_inds.imiss | sed '1d' > $DIR/high_miss_inds
vcftools --gzvcf $DIR/merged.vcf.gz --remove $DIR/high_miss_inds --max-missing $Missing_per_Mrks --maf $MAF --hwe $HWE --recode --recode-INFO-all --out $DIR/merged_qcok.vcf --stdout | bgzip -c > $DIR/merged_qcok.vcf.gz 2> $LOG/qc.log
sed -n 19p $DIR/merged_qcok.vcf.log; sed -n 21p $DIR/merged_qcok.vcf.log
# rm $DIR/merged.vcf.gz $DIR/merged.vcf.gz.tbi

# Identifing unrelated Individuals
echo "# Identifing unrelated Individuals, Cutoff threshold for kinship coefficient: < 0.0442 ..."
plink2 --vcf $DIR/merged_qcok.vcf.gz --king-cutoff 0.0442 --out $DIR/filter  &>  $LOG/plink_kin.filter.log #/dev/null
grep -v "IID" $DIR/filter.king.cutoff.in.id > $RES/list_of_unrelated_individuals.txt
echo "List of unrelated individuals created in $RES directory."
echo "# Check remained individuals for kinship ..."
plink2 --vcf $DIR/merged_qcok.vcf.gz --remove $DIR/filter.king.cutoff.out.id --make-king-table --out $DIR/merged_qc_unrel_kin &> $LOG/plink_kin.table.log
Rscript "src/kinship_check _opt_bash.R"
echo "# Removing related Individulas from merged VCF file"
bcftools view -S "$RES/list_of_unrelated_individuals.txt" "$DIR/merged_qcok.vcf.gz" | bcftools annotate -x INFO | bgzip -c > $DIR/merged_clean.vcf.gz 2> $LOG/filter_related.log
# rm $DIR/merged_qcok.vcf.gz

# Annotating merged vcf file
Rscript "src/info_field_cal_opt _bash.R"

# Filter based on AF
echo "# Filter based on AF (keeping 0<AF<1)"
bcftools view -e "AF=0 | AF=1"  $DIR/merged_clean_annotated.vcf.bgz | bgzip -c > $DIR/merged_clean_annotated_filtered.vcf.gz 2> $LOG/filter_AF.log
# rm $DIR/merged_clean.vcf.gz $DIR/merged_clean_annotated.vcf.bgz $DIR/merged_clean_annotated.vcf.bgz.tbi

# Split and reheader merged VCF file
echo "# Split and reheader merged VCF file into chromosomes"
tabix -f -p vcf $DIR/merged_clean_annotated_filtered.vcf.gz # Index VCF file
for i in {19..22};do
bcftools view $DIR/merged_clean_annotated_filtered.vcf.gz --regions chr$i | sed -e "/ID=chr${i}/p" -e '/contig/d' -e '/##bcftools/d' | bgzip -c > $RES/$i.out.vcf.gz 2> $LOG/chr${i}_split_reheader.log
done
# rm $DIR/merged_clean_annotated_filtered.vcf.gz $DIR/merged_clean_annotated_filtered.vcf.gz.tbi

# Run Time
END=`date +%s`
RUN_TIME=$((END-START))
echo "All Done! Run Time: $RUN_TIME sec"

### END
