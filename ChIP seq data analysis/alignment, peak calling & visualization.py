# alignment
bowtie2-build -f Mus_musculus.GRCm39.dna_sm.primary_assembly.fa mouseChIPindex
gunzip filtered_SRR3020039.fastq.gz
bowtie2 -p 2 -q --local \   
      -x ~/GSE76078_ChIPseq/GCRm39_index/mouseChIPindex \
      -U filtered_SRR3020039.fastq   \
      -S filtered_SRR3020039.sam
samtools view -h -S -b \
      -o filtered_SRR3020039.bam \
      filtered_SRR3020039.sam
sambamba sort -t 8 \
      -o filtered_SRR3020039_sorted.bam \
      filtered_SRR3020039.bam
sambamba view -h -t 2 -f bam \
      -F "[XS] == null and not unmapped and not duplicate" \
      filtered_SRR3020039_sorted.bam > sambamba_SRR3020039.bam 

# peak calling
macs2 callpeak -t sambamba_SRR3020039.bam \
	    -c /input/sambamba_SRR3020076.bam \
 	    -f BAM \
      --gsize mm \            
	    -n SRR3020039 \
	    --outdir SRR3020039_polII_tumor

# visualizing in IGV
bamCoverage -b SRR3020039.bam \
      -o bigWig/SRR3020039.bw \
      --binSize 20 \
      --normalizeUsing BPM \
      --smoothLength 60 \
      --extendReads 150 \
      --centerReads \
      -p 4 2> bigWig/SRR3020039.log
computeMatrix reference-point --referencePoint TSS \
     -b 1000 -a 1000 \
     -R SRR3020039_polII_tumor/SRR3020039_summits.bed \
     -S bigWig/SRR3020039.bw \
     --skipZeros \
     -o bigWig/SRR3020039_matrix_TSS.gz \
     -p 4 \
     --outFileSortedRegions bigWig/SRR3020039_regions_TSS.bed
plotHeatmap -m bigWig/SRR3020039_matrix.gz \
     -out bigWig/SRR3020039_matrix.pdf \
     --colorMap 'Blues' \
     --whatToShow 'heatmap and colorbar' \
     --zMin 0 --zMax 1


