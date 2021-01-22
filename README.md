# Answer to skills_test
A bash script for annotation and filtering VCF files.

## List of used software tools:
* R =4.0.3
* java =11.0.9.1
* picard =2.23.3
* vcftools =0.1.16
* bcftools =1.10.2
* plink2 =2.00a2.3LM
* tabix =1.10.2-3
* VariantAnnotation r package =1.36.0

## Operating system requirements
* moreutils =0.63-1

## Compilation instructions
* No compilation needed.

## To run, execute:
     ./bash_script_opt
### Expected output
* A directory called `results` contains output VCFs and a list of unrelated individuals.
* A directory called `intermediate` contains intermediate files.
* A directory called `LOG` contains log files.
### Workflow steps
![alt text](https://github.com/MehdiFard/Skill_test_McGill_Snakemake/blob/main/DAG.svg?raw=true)
