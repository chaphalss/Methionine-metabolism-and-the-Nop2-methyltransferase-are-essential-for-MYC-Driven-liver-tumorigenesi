fastqc -t 8 DMF_2_1.fq.gz DMF_2_2.fq.gz
fastp --thread 8 --detect_adapter_for_pe -i DMF_2_1.fq.gz -I DMF_2_2.fq.gz -o trimmed_DMF_2_1.fq.gz -O trimmed_DMF_2_2.fq.gz -q 20 -u 50 -n 5
wget http://ftp.ensembl.org/pub/release-104/fasta/mus_musculus/dna/Mus_musculus.GRCm39.dna_sm.primary_assembly.fa.gz
gunzip Mus_musculus.GRCm39.dna_sm.primary_assembly.fa.gz
wget http://ftp.ensembl.org/pub/release-104/gtf/mus_musculus/Mus_musculus.GRCm39.104.gtf.gz
gunzip Mus_musculus.GRCm39.104.gtf.gz
STAR --runThreadN 8 --runMode genomeGenerate --genomeDir EC4 --genomeFastaFiles Mus_musculus.GRCm39.dna_sm.primary_assembly.fa --sjdbGTFfile 
Mus_musculus.GRCm39.104.gtf --limitGenomeGenerateRAM 31000000000 --genomeSAsparseD 2 --sjdbOverhang 100
STAR --genomeDir ~/EC4/Mus_musculus.GRCm39 --readFilesCommand gunzip -c --readFilesIn trimmed_DMF_2_1.fq.gz trimmed_DMF_2_2.fq.gz --runThreadN 8 outFilterIntronMotifs RemoveNoncanonical --outSAMattributes All --outSAMtype BAM SortedByCoordinate
samtools sort -n Aligned.sortedByCoord.out.bam -o Namesorted_Aligned.sortedByCoord.out.bam
samtools index Aligned.sortedByCoord.out.bam
gtfToGenePred ~/EC4/Mus_musculus.GRCm39.104.gtf ~/EC4/Mus_musculus.GRCm39.genePred
genePredToBed ~/EC4/Mus_musculus.GRCm39.genePred ~/EC4/Mus_musculus.GRCm39.bed
bam_stat.py -i Namesorted_Aligned.sortedByCoord.out.bam
junction_saturation.py -i Namesorted_Aligned.sortedByCoord.out.bam -r ~/EC4/Mus_musculus.GRCm39.bed -o output
read_distribution.py -i Namesorted_Aligned.sortedByCoord.out.bam -r ~/EC4/Mus_musculus.GRCm39.bed
geneBody_coverage.py -r ~/EC4/Mus_musculus.GRCm39.bed -i Namesorted_Aligned.sortedByCoord.out.bam -o output
inner_distance.py -i Namesorted_Aligned.sortedByCoord.out.bam -o output -r ~/EC4/Mus_musculus.GRCm39.bed
htseq-count --stranded=reverse --order=name --idattr=gene_id --type=exon --mode=union --format=bam --additional-attr=gene_name Namesorted_Aligned.sortedByCoord.out.bam ~/EC4/Mus_musculus.GRCm39.104.gtf > DMF_2_counts.txt
