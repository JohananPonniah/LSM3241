echo "------------------------------------------"
echo "Code for LSM3241 CA2"
echo "By: Johanan Ponniah (A0155676R)"
echo "------------------------------------------"
echo "Common Step to both methods: examine read quality using FastQC"
gunzip A0155676R.tgz"
tar -xvf A0155676R.tar" #extract sequenced data results
mkdir -p fastqc-results #make a directory for fastqc results
cd Raw-Data-and-Instruction
fastqc A0155676R*.fq -o ../fastqc-results #generate fastqc files
ls ../fastqc-results
cd ../fastqc-results
for filename in *.zip; do
unzip $filename
done #unzip both fastqc files
cat */summary.txt > fastqc-summary.txt #combine summary report for both files
less fastqc-summary.txt
echo "------------------------------------------"

echo "Method A"
echo "Step 1: create and index a combined yeast-transposon reference sequence"
cat Raw-Data-and-Instruction/ty5_6p.fa Raw-Data-and-Instruction/sacCer3.fa > Raw-Data-and-Instruction/YandT.fa #combine transposon and yeast genome
bowtie2-build Raw-Data-and-Instruction/YandT.fa Raw-Data-and-Instruction/YandT #index the combined reference file
echo "---------------------------------------------"
echo "Step 2: align sequence data to yeast-transposon reference"
bowtie2 -x Raw-Data-and-Instruction/YandT \
-1 Raw-Data-and-Instruction/A0155676R_1.fq \
-2 Raw-Data-and-Instruction/A0155676R_2.fq \
-S align-results/YeastTrans-align.sam #Align sequenced data with reference transposon and yeast sequence
echo "---------------------------------------------"
echo "Step 3: verify sequencing quality of reads"
samtools view -S -b align-results/YeastTrans-align.sam > align-results/YeastTrans-align.bam #convert SAM file to BAM file
samtools view -q 10 -b align-results/YeastTrans-align.bam > align-results/YeastTrans-align-filt.bam #filter only reads that have mapping quality of at least 10
samtools fastq align-results/YeastTrans-align-filt.bam > align-results/YeastTrans-align-filt.fastq #convert BAM file to fastq file
fastqc align-results/YeastTrans-align-filt.fastq -o fastqc-results #generate fastqc results
echo "---------------------------------------------"
echo "Step 4: prepare file for viewing in IGV (filtering of mapping quality done within IGV)"
samtools sort align-results/YeastTrans-align.bam -o align-results/YeastTrans-align-sort.bam #Sort BAM file by coordinates
samtools flagstat align-results/YeastTrans-align-sort.bam #view information about sorted BAM file
samtools index align-results/YeastTrans-align-sort.bam #index the sorted BAM file
samtools view -q 5 align-results/YeastTrans-align.sam | awk '$7!="="' | wc -l #verify results of flagstat (number of mates mapped to different chromosome)
echo "---------------------------------------------"
echo "Step 5: preview potential transposon insertion sites"
samtools view -q 10 align-results/YeastTrans-align.sam | awk '$7!="="' | cut -f 3,4 | sort -n #use more stringent mapping quality threshold (-q 10) and view positions of alignment
echo "---------------------------------------------"

echo "Method B"
echo "Step 1: index the reference genomes"
mkdir align-results
bowtie2-build Raw-Data-and-Instruction/sacCer3.fa Raw-Data-and-Instruction/sacCERs288c #Index yeast reference genome
bowtie2-build Raw-Data-and-Instruction/ty5_6p.fa Raw-Data-and-Instruction/ty5-6p #Index transposon reference genome
echo "------------------------------------------"
echo "Step 2: align sequenced data with transposon reference sequence"
bowtie2 --fr -x Raw-Data-and-Instruction/ty5-6p \
-1 Raw-Data-and-Instruction/A0155676R_1.fq \
-2 Raw-Data-and-Instruction/A0155676R_2.fq \
-S align-results/transposon-align.sam #Align sequenced data with reference transposon sequence

samtools view -S -b align-results/transposon-align.sam > align-results/transposon-align.bam #covert SAM file to BAM format
samtools sort align-results/transposon-align.bam -o align-results/transposon-align-sorted.bam #Sort BAM file by coordinates
samtools flagstat align-results/transposon-align-sorted.bam #view information about sorted BAM file

samtools index align-results/transposon-align-sorted.bam #index the sorted BAM file
samtools tview align-results/transposon-align-sorted.bam Raw-Data-and-Instruction/ty5_6p.fa #view the alignment in the terminal
echo "-------------------------------------------"
echo "Step 3: identify mapped reads with unmapped mates"
samtools view -F 4 -f 8 -b align-results/transposon-align-sorted.bam >align-results/transposon-yeast.bam #extract reads that map to transposon but have unmapped mates
samtools index align-results/transposon-yeast.bam #index the bam file
echo "upload indexed 'transposon-yeast.bam' to IGV" #view alignment in IGV
echo "in IGV: filter out low mapping quality reads (threshold 10)" #only keep alignments with good mapping quality
echo "26 reads remaining -> exported as 'align-results/transposon-yeast-mapq.sam'" #export filtered alignments
echo "-------------------------------------------"
echo "Step 4: Quality control and filtering of the 26 reads"
samtools view -S -b align-results/transposon-yeast-mapq.sam > align-results/transposon-yeast-mapq.bam #convert 26-read SAM file to BAM file
samtools fastq align-results/transposon-yeast-mapq.bam > align-results/transposon-yeast-mapq.fastq #convert 26-read BAM file to fastq file
fastqc align-results/transposon-yeast-mapq.fastq -o fastqc-results #generate fastqc results

samtools index align-results/transposon-yeast-mapq.bam #index the bam file
samtools view -b align-results/transposon-yeast-mapq.bam TY5:1-1000 TY5:4000-5388 > align-results/filt22.bam #remove reads aligning to middle of transposon
samtools fastq align-results/filt22.bam > align-results/filt22.fastq #convert 22-read BAM file to fastq file
sed -n /@chr/p align-results/filt22.fastq > align-results/filt22.txt #extract all 22 read IDs from fastq file
echo "-------------------------------------------"
echo "Step 5: map mates of 22 reads to yeast reference genome"
mkdir yeast-mates #make directory for mates who did not map to transposon
grep "/1" align-results/filt22.txt | sed 's/1$/2/g' > yeast-mates/mates2.txt #identify unmapped mates in file 2
grep "/2" align-results/filt22.txt | sed 's/2$/1/g' > yeast-mates/mates1.txt #identify unmapped mates in file 1

seqtk subseq Raw-Data-and-Instruction/A0155676R_1.fq yeast-mates/mates1.txt > yeast-mates/mates1.fq #filter out unmapped mates from file 1
seqtk subseq Raw-Data-and-Instruction/A0155676R_2.fq yeast-mates/mates2.txt > yeast-mates/mates2.fq #filter out unmapped mates from file 2

bowtie2 --local --score-min G,10,8 -x Raw-Data-and-Instruction/sacCERs288c \
-U yeast-mates/mates1.fq,yeast-mates/mates2.fq \
-S align-results/yeast-align.sam #Align unmapped mates to reference yeast genome using local alignment

samtools view -S -b align-results/yeast-align.sam > align-results/yeast-align.bam #convert SAM file to BAM file
samtools sort align-results/yeast-align.bam -o align-results/yeast-align-sorted.bam #sort the BAM file
samtools flagstat align-results/yeast-align-sorted.bam #view details of BAM file
samtools index align-results/yeast-align-sorted.bam #index the bam file
echo "---------------------------------------------"
