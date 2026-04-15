#https://isoseq.how/classification/workflow #pigeon
#Begin with the bulk workflow which ends at isoseq cluster2,
#then continue to pigeon workflow for transcript mapping, collapse, and classification.

# 0. Generate segmented reads
skera hifi_reads.bam mas16_primers.fasta segmented.bam 

# 1. Input (reads should be segmented) otherwise run skera prior to these steps
ls DATA-EXAMPLE/human_80k_Sreads.segmented.bam

# 2. Primer removal and demux
lima DATA-EXAMPLE/human_80k_Sreads.segmented.bam REF-primers/IsoSeq_v2_primers_12.fasta movieX.fl.bam --isoseq --peek-guess

# 3a. Refine (trimming polyA & concatemer)
isoseq refine movieX.fl.IsoSeqX_bc11_5p--IsoSeqX_3p.bam REF-primers/IsoSeq_v2_primers_12.fasta movieX.flnc.bam --require-polya

# 3b. Merge SMART cells (if necessary)
#`ls movie*.flnc.bam movie*.flnc.bam movie*.flnc.bam > flnc.fofn`

# 4. Cluster isoforms (CPU intense ? akin to a de novo assembly)
isoseq cluster2 -j 1 movieX.flnc.bam transcripts.bam

# 5. Check sequences (generate a .fastq)
samtools fastq transcripts.bam >transcripts.fastq

# 6. ppmm2: map reads to human genome
pbmm2 align -j 1 --preset ISOSEQ --sort transcripts.bam ../reference/Homo_sapiens/Homo_sapiensChr.GRCh38.dna.primary_assembly.fa mapped.bam

# 7. Collapse into single isoforms
isoseq collapse --do-not-collapse-extra-5exons mapped.bam movieX.flnc.bam collapsed.gff

# 8. Prepare reference files for pigeon
#pigeon prepare ../reference/gencode/gencode.v49.annotation.gtf ../reference/Homo_sapiens/Homo_sapiens.GRCh38.dna.primary_assembly.fa
pigeon prepare ../reference/Human_hg38_Gencode_v39//gencode.v39.annotation.gtf  ../reference/Human_hg38_Gencode_v39/human_GRCh38_no_alt_analysis_set.fasta

# 9.Prepare .gff (sort)
pigeon prepare collapsed.gff

# 10. Classify isoforms
#pigeon classify collapsed.sorted.gff ../reference/gencode/gencode.v49.annotation.sorted.gtf ../reference/Homo_sapiens/Homo_sapiens.GRCh38.dna.primary_assembly.fa 
pigeon classify collapsed.sorted.gff ../reference/Human_hg38_Gencode_v39/gencode.v39.annotation.sorted.gtf ../reference/Human_hg38_Gencode_v39/human_GRCh38_no_alt_analysis_set.fasta --fl collapsed.flnc_count.txt #add count data outputted from isoseq collapse
pigeon filter classification.txt --isoforms collapsed.sorted.gff # filter .gff

# 11. Gene saturation to check if you sequenced enough
pigeon report --exclude-singletons classification.filtered_lite_classification.txt saturation.txt



