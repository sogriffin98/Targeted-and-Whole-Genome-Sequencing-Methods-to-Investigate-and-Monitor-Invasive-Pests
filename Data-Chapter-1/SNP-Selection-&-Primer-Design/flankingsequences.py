from __future__ import print_function
import pysam


# open  text file
f = open("SNP_CDS.txt", "w")
# open vcf file
vcf = pysam.VariantFile("AH_subset_7inds_minDP6GQ18_biallelic_014missing_mac2_maxDP_nomissing_thin10k_04_06freq_exhetero_LDpruned.recode.vcf")
# open fasta file
genome = pysam.FastaFile("GCF_912470025.1_iVesVel2.1_genomic.fna")
L=[]
# define by how many bases the variant should be flanked
flank = 100

# iterate over each variant
for record in vcf:
    # extract sequence
    #
    # The start position is calculated by subtract the number of bases
    # given by 'flank' from the variant position. The position in the vcf file
    # is 1-based. pysam's fetch() expected 0-base coordinate. That's why we
    # need to subtract on more base.
    #
    # The end position is calculated by adding the number of bases
    # given by 'flank' to the variant position. We also need to add the length
    # of the REF value and subtract again 1 due to the 0-based/1-based thing.
    #
    # Now we have the complete sequence like this:
    # [number of bases given by flank]+REF+[number of bases given by flank]
    seq = genome.fetch(record.chrom, record.pos-1-flank,record.pos-1+len(record.ref)+flank)
    t=(seq, record.chrom, str(record.pos), record.id, record.ref, record.alts[0])
    L.append(t)

    # print out tab seperated columns:
    # CRHOM, POS, REF, ALT, flanking sequencing with variant given in the format '[REF/ALT]'
for i in L:

    seq,chr,pos,id,ref,alt=i

    sequence='{}[{}/{}]{}'.format(seq[:flank],ref,alt, seq[flank+len(ref):])

    line='{}\t{}\t{}\t{}'.format(chr, pos, id, sequence)

    f.write(line+'\n')



f.close()
