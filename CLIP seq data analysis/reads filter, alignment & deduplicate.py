fastp --thread 8 \
      --detect_adapter_for_pe \
      -i SRR16936854_1.fastq.gz -I SRR16936854_2.fastq.gz \
      -o trimmed_SRR16936854_1.fastq.gz -O trimmed_SRR16936854_2.fastq.gz  \
      --adapter_sequence=AGATCGGAAGAGCACACGTC --adapter_sequence_r2=AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
      -q 20 -u 50 -n 5 --dedup -U --umi_loc=read1 --umi_len=8
STAR --runThreadN 8 \
     --runMode genomeGenerate 
     --genomeDir hs38_STAR 
     --genomeFastaFiles Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa  \
     --sjdbGTFfile Homo_sapiens.GRCh38.105.gtf \
     --limitGenomeGenerateRAM 31000000000 \
     --genomeSAsparseD 2 
     --sjdbOverhang 100
STAR --genomeDir ~/hs38/hs38_STAR \
     --readFilesCommand zcat  \
     --readFilesIn trimmed_SRR16936854_1.fastq.gz trimmed_SRR16936854_2.fastq.gz  \
     --outFileNamePrefix SRR16936854_  \
     --runThreadN 8 \
     --outBAMsortingThreadN 4 \
     --outSAMstrandField intronMotif \
     --outFilterIntronMotifs RemoveNoncanonicalUnannotated  \
     --outSAMattributes All  \
     --outSAMtype BAM SortedByCoordinate  \
     --outFilterScoreMinOverLread 0.2  \
     --outFilterMatchNminOverLread 0.2  \
     --outFilterMismatchNmax 2
samtools index SRR16936854_Aligned.sortedByCoord.out.bam
umi_tools dedup --stdin=SRR16936854_Aligned.sortedByCoord.out.bam   \
     --log=UMItoolsdedup_logfile > SRR16936854_deduplicated.bam  \
     --read-length --umi-separator=':'

