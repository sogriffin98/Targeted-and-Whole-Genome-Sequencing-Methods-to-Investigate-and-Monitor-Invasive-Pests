# Yellow Legged Hornets Combined MiSeq Runs Analysis
Each sequencing run was analysed separately first but then using all the bam files which had been aligned to the published reference genome (NCBI Accession: PRJEB46979) I analysed them as one group and created different subsets depending on sample location and I used a 1 per nest dataset for population genetics analysis. These were the instructions used to create the VCF file:

## 1. Merge the bam files from multiple sequencing runs
Some of the samples were repeated and therefore I combined their bam files using the merge function of samtools:
```
conda activate samtools
samtools merge <>_combined.bam <1>.bam <2>.bam
```

## 2. Sort and index the bam files
```
conda activate samtools
chmod +x bam_sort_index.sh
./bam_sort_index.sh 
```

## 3. Filter the bam files
```
conda activate bamtools
chmod +x filter_bam.sh
./filter_bam.sh
```

## 4. Create VCF file using the filtered bam files and reference genome
```
conda activate bcftools
bcftools mpileup -f ./reference/GCF_912470025.1_iVesVel2.1_genomic.fna --annotate AD,DP ./bam/*primaryalignment.bam | bcftools call -m -v -f GQ --skip-variants indels -o ./vcf/AH_400_inds.vcf
```

## 5. Filter and Quality Check the VCF
```
chmod +x vcf_checks.sh
