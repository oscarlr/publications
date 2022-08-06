# README
## Data
1. [SRA accessions for raw HiFi capture data](https://github.com/oscarlr/publications/blob/main/Targeted%20long-read%20sequencing%20facilitates%20phased%20diploid%20assembly%20and%20genotyping%20of%20the%20human%20T%20cell%20receptor%20alpha%2C%20delta%20and%20beta%20loci/metadata-11518572-processed-ok.tsv)
2. [TRA/D and TRB assemblies](https://github.com/oscarlr/publications/raw/main/Targeted%20long-read%20sequencing%20facilitates%20phased%20diploid%20assembly%20and%20genotyping%20of%20the%20human%20T%20cell%20receptor%20alpha%2C%20delta%20and%20beta%20loci/assemblies.tar.gz)
3. [Novel VDJ alleles](https://raw.githubusercontent.com/oscarlr/publications/main/Targeted%20long-read%20sequencing%20facilitates%20phased%20diploid%20assembly%20and%20genotyping%20of%20the%20human%20T%20cell%20receptor%20alpha%2C%20delta%20and%20beta%20loci/vdj_alleles.fasta)
4. [Novel constant alleles](https://raw.githubusercontent.com/oscarlr/publications/main/Targeted%20long-read%20sequencing%20facilitates%20phased%20diploid%20assembly%20and%20genotyping%20of%20the%20human%20T%20cell%20receptor%20alpha%2C%20delta%20and%20beta%20loci/constant_alleles.fasta)

## Scripts
1. Testing the accuracy of alleles by checking the mendelian inheritance of the alleles in the proband using the alleles in the parents
```
python scripts/test_allele_acc/trio_validate_alleles.py \
  NA18506_assembly_alleles.bed \
  NA18507_assembly_alleles.bed \
  NA18508_assembly_alleles.bed
  
python scripts/test_allele_acc/trio_validate_alleles.py \
  HG02059_assembly_alleles.bed \
  HG02060_assembly_alleles.bed \
  HG02061_assembly_alleles.bed
```
