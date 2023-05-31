#!/bin/sh
threads=6

REF="../references/CP015498.fasta"
TRIMMED="../outputs/trimmed"
SPADES="../outputs/spades"
MEGAHIT="../outputs/megahit"
RAGTAG="../outputs/ragtag"
SCAFFOLDS="ragtag.scaffold.fasta"
BLAST="../outputs/blast"
CONTIGS="final.contigs.fa"
GEPARD="../outputs/gepard"
BWA="../outputs/bwa"
ERR="ERR204044"
SRR15="SRR15131330"
SRR18="SRR18214264"
QUAST2="../outputs/quast2"
JAR="../../gepard/dist/Gepard-1.40.jar"
Q="../outputs/quast/all"
Q1="../outputs/quast/ERR204044"
Q2="../outputs/quast/SRR15131330"
Q3="../outputs/quast/SRR18214264"

#//////////////////Getting data//////////////////////////

#Atsisiunciam duomenis
#wget -P ../inputs "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR204/ERR204044/ERR204044_1.fastq.gz"
#wget -P ../inputs "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR204/ERR204044/ERR204044_2.fastq.gz"
#wget -P ../inputs "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR151/030/SRR15131330/SRR15131330_1.fastq.gz"
#wget -P ../inputs "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR151/030/SRR15131330/SRR15131330_2.fastq.gz"
#wget -P ../inputs "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR182/064/SRR18214264/SRR18214264_1.fastq.gz"
#wget -P ../inputs "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR182/064/SRR18214264/SRR18214264_2.fastq.gz"
#efetch -db nucleotide -id CP015498 -format fasta > ../references/CP015498.fasta

#///////////////////Data QA/QC///////////////////////////

# Data quality evaluation
#fastqc -t $threads ../inputs/* -o ../outputs/raw_data/

#ERR204044 kokybė gera.
#SRR15131330 kokybė labai gera.
#SRR18214264 kokybė prasčiausia iš visų, suprastėjo ties 3' galu.

# Trimming process

#for i in ../inputs/*1.fastq.gz
#do
#R1=$i
#R2="../inputs/$(basename $R1 1.fastq.gz)2.fastq.gz"
#trim_galore -paired $R1 $R2 --fastqc -o ../outputs/trimmed/ --length 20 -q 30
#done

# KOMENTARAS
#ERR204044 kokybė gera, šiek tiek apkarpytos sekos.
#SRR15131330 kokybė labai gera, minimaliai apkarpytas 3' galas.
#SRR18214264 kokybė ties 3' galu pakilo, tačiau dar matome prastesnę galo kokybę.

#multiqc -p ../outputs/raw_data ../outputs/trimmed -o ../outputs/multiqc

#/////////////////Genome assembly//////////////////////

# genome assembly with spades

#for i in $TRIMMED/*_1_val_1.fq.gz
#do
#R1=$i
#R2=$TRIMMED/"$(basename $R1 _1_val_1.fq.gz)_2_val_2.fq.gz"
#/usr/lib/spades/bin/spades.py -o $SPADES/$(basename $R1 _1_val_1.fq.gz) -1 $R1 -2 $R2
#done

# genome assembly with megahit

#for i in $TRIMMED/*_1_val_1.fq.gz
#do
#R1=$i
#R2=$TRIMMED/"$(basename $R1 _1_val_1.fq.gz)_2_val_2.fq.gz"
#megahit -1 $R1 -2 $R2 -o $MEGAHIT/$(basename $R1 _1_val_1.fq.gz)
#done

# quast evaluation

#python3 ~/quast/quast/quast.py -o $Q1 -r $REF --gene-finding $SPADES/$ERR/contigs.fasta $MEGAHIT/$ERR/final.contigs.fa
#python3 ~/quast/quast/quast.py -o $Q2 -r $REF --gene-finding $SPADES/$SRR15/contigs.fasta $MEGAHIT/$SRR15/final.contigs.fa
#python3 ~/quast/quast/quast.py -o $Q3 -r $REF --gene-finding $SPADES/$SRR18/contigs.fasta $MEGAHIT/$SRR18/final.contigs.fa
#python3 ~/quast/quast/quast.py -o $Q -r $REF --gene-finding $SPADES/$ERR/contigs.fasta $MEGAHIT/$ERR/final.contigs.fa $SPADES/$SRR15/contigs.fasta $MEGAHIT/$SRR15/final.contigs.fa $SPADES/$SRR18/contigs.fasta $MEGAHIT/$SRR18/final.contigs.fa

# ERR204044 genome fraction šiek tiek geresnė su megahit, turi didesni bendra ilgi ir contig'ų skaičių, tik N50 didesnis su spades.
# SRR15131330 Tiek genome fraction, total length, number of contigs didesni su megahit, bet N50 - su spades, tačiau skirtumas nedidelis. 
# SRR18214264 genome fraction, total length, number of contigs yra didesni su megahit, spades turi didesnį contig'ų sk ir mažiau misassemblies.

# scaffolds with ragtag

#for i in $SPADES/*
#do
#R1=$i/contigs.fasta
#ragtag.py correct $REF $R1 -o $RAGTAG/spades/$(basename $i)
#done

#for i in $MEGAHIT/*
#do
#R1=$i/final.contigs.fa
#ragtag.py correct $REF $R1 -o $RAGTAG/megahit/$(basename $i)
#done

#for i in $RAGTAG/spades/*
#do
#R1=$i/ragtag.correct.fasta
#ragtag.py scaffold $REF $R1 -o spades_$i
#done

#for i in $RAGTAG/megahit/*
#do
#R1=$i/ragtag.correct.fasta
#ragtag.py scaffold $REF $R1 -o megahit_$i
#done

#python3 ~/quast/quast/quast.py $RAGTAG/megahit/$SRR18/$SCAFFOLDS $RAGTAG/megahit/$SRR15/$SCAFFOLDS $RAGTAG/megahit/$ERR/$SCAFFOLDS $RAGTAG/spades/$SRR18/$SCAFFOLDS $RAGTAG/spades/$SRR15/$SCAFFOLDS $RAGTAG/spades/$ERR/$SCAFFOLDS

#python3 ~/quast/quast/quast.py  -r $REF -o $QUAST2/$SRR18 --gene-finding $RAGTAG/megahit/$SRR18/$SCAFFOLDS $RAGTAG/spades/$SRR18/$SCAFFOLDS
#python3 ~/quast/quast/quast.py  -r $REF -o $QUAST2/$SRR15 --gene-finding $RAGTAG/megahit/$SRR15/$SCAFFOLDS $RAGTAG/spades/$SRR15/$SCAFFOLDS
#python3 ~/quast/quast/quast.py  -r $REF -o $QUAST2/$ERR --gene-finding $RAGTAG/megahit/$ERR/$SCAFFOLDS $RAGTAG/spades/$ERR/$SCAFFOLDS

# KOMENTARAS
#Renkuosi assembly po megahit. Genome fraction, largest alignment, N50, NG50, total lenght yra didesni nei spades (geresni įverčiai, nors ir nežymiai) 

# mapping

#for i in $RAGTAG/megahit/*
#do
#R1=$i/ragtag.scaffold.fasta
#bwa index $R1
#done

#for i in $TRIMMED/*_1_val_1.fq.gz
#do
#R1=$i
#R2=$TRIMMED/"$(basename $R1 _1_val_1.fq.gz)_2_val_2.fq.gz"
#R3="$(basename $R1 _1_val_1.fq.gz)"
#echo "$R3"
#echo "$R1 $R2"
#bwa index $RAGTAG/megahit/"$R3"/$SCAFFOLDS
#bwa mem $RAGTAG/megahit/"$R3"/$SCAFFOLDS $R1 $R2 > $BWA/$(basename $R1 _1_val_1.fq.gz).sam
#done

#for i in $BWA/*
#do
#echo $(basename $i)
#samtools view -bS $i -@ 6 | samtools sort -@ 6  -o $BWA/$(basename $i .sam)_sorted.bam
#done

#for i in $BWA/*.bam
#do
#echo $i
#samtools index -b $i
#samtools flagstat $i
#done

#samtools faidx $REF

#for i in $BWA/*.bam
#do
#bedtools bamtobed -i $BWA/"$(basename $i .bam).bam" | sort -k1,1 -k2,2n > $BWA/"$(basename $i .bam).bed"
#bedtools genomecov -ibam $BWA/"$(basename $i .bam).bam" -g $../references/CP015498.fasta.fai > $BWA/"genome_coverage_$(basename $i .bam)_assembly.txt"
#done

# gepard

#java -cp $JAR org.gepard.client.cmdline.CommandLine -seq1 $RAGTAG/megahit/$ERR/$SCAFFOLDS -seq2 $RAGTAG/megahit/$SRR15/$SCAFFOLDS -matrix ../../gepard/resources/matrices/edna.mat  -outfile $GEPARD/E_S15_output.png
#java -cp $JAR org.gepard.client.cmdline.CommandLine -seq1 $RAGTAG/megahit/$SRR15/$SCAFFOLDS -seq2 $RAGTAG/megahit/$SRR18/$SCAFFOLDS -matrix ../../gepard/resources/matrices/edna.mat -outfile $GEPARD/S15_S18_output.png
#java -cp $JAR org.gepard.client.cmdline.CommandLine -seq1 $RAGTAG/megahit/$ERR/$SCAFFOLDS -seq2 $RAGTAG/megahit/$SRR18/$SCAFFOLDS -matrix ../../gepard/resources/matrices/edna.mat -outfile $GEPARD/E_S18_output.png

# KOMENTARAS
# Gepard results: Panašiausi yra ERR204044 ir SRR18214264, nes įstrižainė yra tiesiausia. Mažiausias panašumas yra tarp SRR18214264 ir SRR15131330.

#makeblastdb -in $REF -dbtype nucl

#makeblastdb -in $MEGAHIT/$ERR/$CONTIGS -dbtype nucl
#makeblastdb -in $MEGAHIT/$SRR15/$CONTIGS -dbtype nucl
#makeblastdb -in $MEGAHIT/$SRR18/$CONTIGS -dbtype nucl

# KOMENTARAS
# Busco results: 
#ERR204044 turime complete and single copy genu - 99%, fragmented and missing - po 0.5%.
#SRR15131330 turime complete and single copy genu - 98.5%, fragmented - 0.5%, missing - 1%.
#SRR18214264 complete and single copy genu - 99%, fragmented and missing - po 0.5%.

#blastn -query ../references/genes.txt -db $MEGAHIT/$ERR/$CONTIGS > $BLAST/$ERR/blast_gene_info.txt
#blastn -query ../references/genes.txt -db $MEGAHIT/$SRR15/$CONTIGS > $BLAST/$SRR15/blast_gene_info.txt
#blastn -query ../references/genes.txt -db $MEGAHIT/$SRR18/$CONTIGS > $BLAST/$SRR18/blast_gene_info.txt

#tblastn -query ../references/prot.txt -db $MEGAHIT/$ERR/$CONTIGS > $BLAST/$ERR/blast_prot_info.txt;
#tblastn -query ../references/prot.txt -db $MEGAHIT/$SRR15/$CONTIGS > $BLAST/$SRR15/blast_prot_info.txt;
#tblastn -query ../references/prot.txt -db $MEGAHIT/$SRR18/$CONTIGS > $BLAST/$SRR18/blast_prot_info.txt;

#blastn -query ../references/genes.txt -db  $MEGAHIT/$ERR/$CONTIGS -out $BLAST/$ERR/gene.fasta.blastn -outfmt 6 -evalue 1e-5
#blastn -query ../references/genes.txt -db  $MEGAHIT/$SRR15/$CONTIGS -out $BLAST/$SRR15/gene.fasta.blastn -outfmt 6 -evalue 1e-5
#blastn -query ../references/genes.txt -db  $MEGAHIT/$SRR18/$CONTIGS -out $BLAST/$SRR18/gene.fasta.blastn -outfmt 6 -evalue 1e-5

#tblastn -query ../references/prot.txt -db $MEGAHIT/$ERR/$CONTIGS -out $BLAST/$ERR/prot.fasta.blastn -outfmt 6 -evalue 1e-5
#tblastn -query ../references/prot.txt -db $MEGAHIT/$SRR15/$CONTIGS -out $BLAST/$SRR15/prot.fasta.blastn -outfmt 6 -evalue 1e-5
#tblastn -query ../references/prot.txt -db $MEGAHIT/$SRR18/$CONTIGS -out $BLAST/$SRR18/prot.fasta.blastn -outfmt 6 -evalue 1e-5

# KOMENTARAS
#Su RAST: ERR204044 - rasta 2435 genų, SRR15131330 - 2605, SRR18214264 - 2446.
#Su GeneMarkS2: ERR204044 - 2290 genų, SRR15131330 - 2517, SRR18214264 - 2309.
#
