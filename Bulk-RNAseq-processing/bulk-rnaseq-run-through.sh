#!/usr/bin/env bash

# set up working directory
# workdir='/path/to/workdir'
workdir=`pwd`

# download SRA toolkit
wget -c https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/current/sratoolkit.current-centos_linux64.tar.gz
tar zxvf ./sratoolkit.current-centos_linux64.tar.gz

# download htseq if needed
case :$PATH:
  in *:htseq-count:*) ;;
     *) pip install HTSeq;;
esac


# Refer to the GEO and SRA accession page. https://www.ncbi.nlm.nih.gov/Traces/study/?acc=PRJNA893220
# download SRA file and convert to fastq
bulkf4=("SRR22013626" "SRR22013628" "SRR22013630" "SRR22013632" "SRR22013634" "SRR22013636" "SRR22013627" "SRR22013629" "SRR22013631" "SRR22013633" "SRR22013635" "SRR22013637")
./sratoolkit.3.1.0-centos_linux64/bin/prefetch `( IFS=$' '; echo "${bulkf4[*]}" )`
for j in ${!bulkf4[@]}; do
  ./sratoolkit.3.1.0-centos_linux64/bin/fastq-dump ${j}.sra --split-files --gzip -O ${j}/
done

# check md5 sums
md5sum SRR22013*/*.fastq.gz
# d9e6453a0f3b2d8cb60759da876a3180  SRR22013626/SRR22013626_1.fastq.gz
# 78feb9e7e557d4c8bb9f1e9d290ea4a3  SRR22013626/SRR22013626_2.fastq.gz
# 906270b9673b20ca756b48c921ba0ba0  SRR22013627/SRR22013627_1.fastq.gz
# bebc153c73cadc4eecbf1a5cc293733c  SRR22013627/SRR22013627_2.fastq.gz
# cb9224b78157c9b015e25b002d6482cb  SRR22013628/SRR22013628_1.fastq.gz
# 974f767ffd9b5c7eb85e565da5f3fe5a  SRR22013628/SRR22013628_2.fastq.gz
# a430c00897a259bf21b2d2238e248f73  SRR22013629/SRR22013629_1.fastq.gz
# 11db85bd6c6286b5778defdf2c3f370b  SRR22013629/SRR22013629_2.fastq.gz
# 567a746346fda2d6477d58d16baa927a  SRR22013630/SRR22013630_1.fastq.gz
# 84ee3e1a1abedf16bd2a020cdb1ee2b7  SRR22013630/SRR22013630_2.fastq.gz
# 33b89e926664f8fea503c5c89a230b60  SRR22013631/SRR22013631_1.fastq.gz
# 5daa9f13223462e41c0e08892e6d37b7  SRR22013631/SRR22013631_2.fastq.gz
# 24057a43496a8b4b242a82ecd55a90a3  SRR22013632/SRR22013632_1.fastq.gz
# df23a49d562189cdc42385a9b3bbcc7f  SRR22013632/SRR22013632_2.fastq.gz
# 891939257846f0b9f3a530e32bdc1c92  SRR22013633/SRR22013633_1.fastq.gz
# bd2de1eac68adf565d62c5f4e2e239fc  SRR22013633/SRR22013633_2.fastq.gz
# 830e05b6128453998ef4daa2a0192798  SRR22013634/SRR22013634_1.fastq.gz
# 7b4ab6c332c6b0f183aa87837dcff622  SRR22013634/SRR22013634_2.fastq.gz
# 33c8fe953dc9ec9545610136f43401ac  SRR22013635/SRR22013635_1.fastq.gz
# bf21f33555bf5eba4f9efde9f1350b97  SRR22013635/SRR22013635_2.fastq.gz
# d1c83fe8a44d184d1b3a0e9423998b84  SRR22013636/SRR22013636_1.fastq.gz
# 7931c6568cf92e489c2edd213673bdcc  SRR22013636/SRR22013636_2.fastq.gz
# ec8464424287c9779840c9b848060efa  SRR22013637/SRR22013637_1.fastq.gz
# 283a0084bda3196586c7a7b89b8d75ac  SRR22013637/SRR22013637_2.fastq.gz

# run alignment script
# set up reference genome directory
# refstar='/path/to/starindex'
refstar='./references/Homo_sapiens/UCSC/hg19/Sequence/STARIndex'
# reffasta='/path/to/fasta'
reffasta='./references/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta'

curl "https://awscli.amazonaws.com/awscli-exe-linux-x86_64.zip" -o "awscliv2.zip"
unzip awscliv2.zip
sudo ./aws/install --bin-dir ${PWD}/bin --install-dir ${PWD}/aws-cli
./bin/aws s3 --no-sign-request --region eu-west-1 sync s3://ngi-igenomes/igenomes/Homo_sapiens/UCSC/hg19/Sequence/STARIndex/ ${refstar}/
./bin/aws s3 --no-sign-request --region eu-west-1 sync s3://ngi-igenomes/igenomes/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/ ${reffasta}/
wget -c https://github.com/alexdobin/STAR/releases/download/2.7.11b/STAR_2.7.11b.zip
unzip STAR_2.7.11b.zip
wget -c https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_38/GRCh37_mapping/gencode.v38lift37.annotation.gtf.gz
gunzip gencode.v38lift37.annotation.gtf.gz

./STAR_2.7.11b/Linux_x86_64/STAR \
  --runMode alignReads \
  --runThreadN 16 \
  --genomeDir ${refstar} \
  --limitBAMsortRAM 10000000000 \
  --outSAMstrandField intronMotif \
  --outSAMattrRGline ID:SRR22013626SRR22013627 PL:ILLUMINA SM:SRR22013626SRR22013627 \
  --readFilesIn \
    SRR22013626/SRR22013626_1.fastq.gz,SRR22013627/SRR22013627_1.fastq.gz \
    SRR22013626/SRR22013626_2.fastq.gz,SRR22013627/SRR22013627_2.fastq.gz \
  --readFilesCommand zcat \
  --outFileNamePrefix SRR22013626SRR22013627/ \
  --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
  --outSAMattributes All \
  --outSAMunmapped Within \
  --outSAMtype BAM Unsorted \
  --outBAMcompression -1 \
  --outFilterMultimapNmax 1 \
  --outFilterMismatchNoverLmax 0.06 \
  --outFilterMatchNmin 30 \
  --outSJfilterOverhangMin 30 10 10 10 \
  --seedSearchStartLmax 30 \
  --alignIntronMin 20 \
  --alignIntronMax 20000 \
  --alignEndsType Local
./STAR_2.7.11b/Linux_x86_64/STAR \
  --runMode alignReads \
  --runThreadN 16 \
  --genomeDir ${refstar} \
  --limitBAMsortRAM 10000000000 \
  --outSAMstrandField intronMotif \
  --outSAMattrRGline ID:SRR22013628SRR22013629 PL:ILLUMINA SM:SRR22013628SRR22013629 \
  --readFilesIn \
    SRR22013628/SRR22013628_1.fastq.gz,SRR22013629/SRR22013629_1.fastq.gz \
    SRR22013628/SRR22013628_2.fastq.gz,SRR22013629/SRR22013629_2.fastq.gz \
  --readFilesCommand zcat \
  --outFileNamePrefix SRR22013628SRR22013629/ \
  --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
  --outSAMattributes All \
  --outSAMunmapped Within \
  --outSAMtype BAM Unsorted \
  --outBAMcompression -1 \
  --outFilterMultimapNmax 1 \
  --outFilterMismatchNoverLmax 0.06 \
  --outFilterMatchNmin 30 \
  --outSJfilterOverhangMin 30 10 10 10 \
  --seedSearchStartLmax 30 \
  --alignIntronMin 20 \
  --alignIntronMax 20000 \
  --alignEndsType Local
./STAR_2.7.11b/Linux_x86_64/STAR \
  --runMode alignReads \
  --runThreadN 16 \
  --genomeDir ${refstar} \
  --limitBAMsortRAM 10000000000 \
  --outSAMstrandField intronMotif \
  --outSAMattrRGline ID:SRR22013630SRR22013631 PL:ILLUMINA SM:SRR22013630SRR22013631 \
  --readFilesIn \
    SRR22013630/SRR22013630_1.fastq.gz,SRR22013631/SRR22013631_1.fastq.gz \
    SRR22013630/SRR22013630_2.fastq.gz,SRR22013631/SRR22013631_2.fastq.gz \
  --readFilesCommand zcat \
  --outFileNamePrefix SRR22013630SRR22013631/ \
  --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
  --outSAMattributes All \
  --outSAMunmapped Within \
  --outSAMtype BAM Unsorted \
  --outBAMcompression -1 \
  --outFilterMultimapNmax 1 \
  --outFilterMismatchNoverLmax 0.06 \
  --outFilterMatchNmin 30 \
  --outSJfilterOverhangMin 30 10 10 10 \
  --seedSearchStartLmax 30 \
  --alignIntronMin 20 \
  --alignIntronMax 20000 \
  --alignEndsType Local
./STAR_2.7.11b/Linux_x86_64/STAR \
  --runMode alignReads \
  --runThreadN 16 \
  --genomeDir ${refstar} \
  --limitBAMsortRAM 10000000000 \
  --outSAMstrandField intronMotif \
  --outSAMattrRGline ID:SRR22013632SRR22013633 PL:ILLUMINA SM:SRR22013632SRR22013633 \
  --readFilesIn \
    SRR22013632/SRR22013632_1.fastq.gz,SRR22013633/SRR22013633_1.fastq.gz \
    SRR22013632/SRR22013632_2.fastq.gz,SRR22013633/SRR22013633_2.fastq.gz \
  --readFilesCommand zcat \
  --outFileNamePrefix SRR22013632SRR22013633/ \
  --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
  --outSAMattributes All \
  --outSAMunmapped Within \
  --outSAMtype BAM Unsorted \
  --outBAMcompression -1 \
  --outFilterMultimapNmax 1 \
  --outFilterMismatchNoverLmax 0.06 \
  --outFilterMatchNmin 30 \
  --outSJfilterOverhangMin 30 10 10 10 \
  --seedSearchStartLmax 30 \
  --alignIntronMin 20 \
  --alignIntronMax 20000 \
  --alignEndsType Local
./STAR_2.7.11b/Linux_x86_64/STAR \
  --runMode alignReads \
  --runThreadN 16 \
  --genomeDir ${refstar} \
  --limitBAMsortRAM 10000000000 \
  --outSAMstrandField intronMotif \
  --outSAMattrRGline ID:SRR22013634SRR22013635 PL:ILLUMINA SM:SRR22013634SRR22013635 \
  --readFilesIn \
    SRR22013634/SRR22013634_1.fastq.gz,SRR22013635/SRR22013635_1.fastq.gz \
    SRR22013634/SRR22013634_2.fastq.gz,SRR22013635/SRR22013635_2.fastq.gz \
  --readFilesCommand zcat \
  --outFileNamePrefix SRR22013634SRR22013635/ \
  --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
  --outSAMattributes All \
  --outSAMunmapped Within \
  --outSAMtype BAM Unsorted \
  --outBAMcompression -1 \
  --outFilterMultimapNmax 1 \
  --outFilterMismatchNoverLmax 0.06 \
  --outFilterMatchNmin 30 \
  --outSJfilterOverhangMin 30 10 10 10 \
  --seedSearchStartLmax 30 \
  --alignIntronMin 20 \
  --alignIntronMax 20000 \
  --alignEndsType Local
./STAR_2.7.11b/Linux_x86_64/STAR \
  --runMode alignReads \
  --runThreadN 16 \
  --genomeDir ${refstar} \
  --limitBAMsortRAM 10000000000 \
  --outSAMstrandField intronMotif \
  --outSAMattrRGline ID:SRR22013636SRR22013637 PL:ILLUMINA SM:SRR22013636SRR22013637 \
  --readFilesIn \
    SRR22013636/SRR22013636_1.fastq.gz,SRR22013637/SRR22013637_1.fastq.gz \
    SRR22013636/SRR22013636_2.fastq.gz,SRR22013637/SRR22013637_2.fastq.gz \
  --readFilesCommand zcat \
  --outFileNamePrefix SRR22013636SRR22013637/ \
  --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
  --outSAMattributes All \
  --outSAMunmapped Within \
  --outSAMtype BAM Unsorted \
  --outBAMcompression -1 \
  --outFilterMultimapNmax 1 \
  --outFilterMismatchNoverLmax 0.06 \
  --outFilterMatchNmin 30 \
  --outSJfilterOverhangMin 30 10 10 10 \
  --seedSearchStartLmax 30 \
  --alignIntronMin 20 \
  --alignIntronMax 20000 \
  --alignEndsType Local

# filter sjdb.
cat SRR22013*/*SJ.out.tab | awk '($1 != "chrM" && $5 > 0 && $6 == 0 && $7 > 2)' | cut -f1-6 | sort | uniq > SJ.filtered.tab

# genome generate
./STAR_2.7.11b/Linux_x86_64/STAR \
  --runMode genomeGenerate \
  --runThreadN 16 \
  --genomeDir `pwd` \
  --genomeFastaFiles ${reffasta} \
  --sjdbGTFfile gencode.v38lift37.annotation.gtf \
  --sjdbOverhang 100 \
  --sjdbFileChrStartEnd SJ.filtered.tab

# second alignment
./STAR_2.7.11b/Linux_x86_64/STAR \
  --runMode alignReads \
  --runThreadN 16 \
  --genomeDir `pwd` \
  --limitBAMsortRAM 10000000000 \
  --outSAMstrandField intronMotif \
  --outSAMattrRGline ID:SRR22013626SRR22013627 PL:ILLUMINA SM:SRR22013626SRR22013627 \
  --readFilesIn \
    SRR22013626/SRR22013626_1.fastq.gz,SRR22013627/SRR22013627_1.fastq.gz \
    SRR22013626/SRR22013626_2.fastq.gz,SRR22013627/SRR22013627_2.fastq.gz \
  --readFilesCommand zcat \
  --outFileNamePrefix SRR22013626SRR22013627/Second. \
  --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
  --outSAMattributes All \
  --outSAMunmapped Within \
  --outSAMtype BAM SortedByCoordinate --outBAMcompression -1 \
  --outFilterMultimapNmax 1 \
  --outFilterMismatchNoverLmax 0.06 \
  --outFilterMatchNmin 30 \
  --outSJfilterOverhangMin 30 10 10 10 \
  --seedSearchStartLmax 30 \
  --alignIntronMin 20 \
  --alignIntronMax 20000 \
  --alignEndsType Local
./STAR_2.7.11b/Linux_x86_64/STAR \
  --runMode alignReads \
  --runThreadN 16 \
  --genomeDir `pwd` \
  --limitBAMsortRAM 10000000000 \
  --outSAMstrandField intronMotif \
  --outSAMattrRGline ID:SRR22013628SRR22013629 PL:ILLUMINA SM:SRR22013628SRR22013629 \
  --readFilesIn \
    SRR22013628/SRR22013628_1.fastq.gz,SRR22013629/SRR22013629_1.fastq.gz \
    SRR22013628/SRR22013628_2.fastq.gz,SRR22013629/SRR22013629_2.fastq.gz \
  --readFilesCommand zcat \
  --outFileNamePrefix SRR22013628SRR22013629/Second. \
  --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
  --outSAMattributes All \
  --outSAMunmapped Within \
  --outSAMtype BAM SortedByCoordinate --outBAMcompression -1 \
  --outFilterMultimapNmax 1 \
  --outFilterMismatchNoverLmax 0.06 \
  --outFilterMatchNmin 30 \
  --outSJfilterOverhangMin 30 10 10 10 \
  --seedSearchStartLmax 30 \
  --alignIntronMin 20 \
  --alignIntronMax 20000 \
  --alignEndsType Local
./STAR_2.7.11b/Linux_x86_64/STAR \
  --runMode alignReads \
  --runThreadN 16 \
  --genomeDir `pwd` \
  --limitBAMsortRAM 10000000000 \
  --outSAMstrandField intronMotif \
  --outSAMattrRGline ID:SRR22013630SRR22013631 PL:ILLUMINA SM:SRR22013630SRR22013631 \
  --readFilesIn \
    SRR22013630/SRR22013630_1.fastq.gz,SRR22013631/SRR22013631_1.fastq.gz \
    SRR22013630/SRR22013630_2.fastq.gz,SRR22013631/SRR22013631_2.fastq.gz \
  --readFilesCommand zcat \
  --outFileNamePrefix SRR22013630SRR22013631/Second. \
  --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
  --outSAMattributes All \
  --outSAMunmapped Within \
  --outSAMtype BAM SortedByCoordinate --outBAMcompression -1 \
  --outFilterMultimapNmax 1 \
  --outFilterMismatchNoverLmax 0.06 \
  --outFilterMatchNmin 30 \
  --outSJfilterOverhangMin 30 10 10 10 \
  --seedSearchStartLmax 30 \
  --alignIntronMin 20 \
  --alignIntronMax 20000 \
  --alignEndsType Local
./STAR_2.7.11b/Linux_x86_64/STAR \
  --runMode alignReads \
  --runThreadN 16 \
  --genomeDir `pwd` \
  --limitBAMsortRAM 10000000000 \
  --outSAMstrandField intronMotif \
  --outSAMattrRGline ID:SRR22013632SRR22013633 PL:ILLUMINA SM:SRR22013632SRR22013633 \
  --readFilesIn \
    SRR22013632/SRR22013632_1.fastq.gz,SRR22013633/SRR22013633_1.fastq.gz \
    SRR22013632/SRR22013632_2.fastq.gz,SRR22013633/SRR22013633_2.fastq.gz \
  --readFilesCommand zcat \
  --outFileNamePrefix SRR22013632SRR22013633/Second. \
  --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
  --outSAMattributes All \
  --outSAMunmapped Within \
  --outSAMtype BAM SortedByCoordinate --outBAMcompression -1 \
  --outFilterMultimapNmax 1 \
  --outFilterMismatchNoverLmax 0.06 \
  --outFilterMatchNmin 30 \
  --outSJfilterOverhangMin 30 10 10 10 \
  --seedSearchStartLmax 30 \
  --alignIntronMin 20 \
  --alignIntronMax 20000 \
  --alignEndsType Local
./STAR_2.7.11b/Linux_x86_64/STAR \
  --runMode alignReads \
  --runThreadN 16 \
  --genomeDir `pwd` \
  --limitBAMsortRAM 10000000000 \
  --outSAMstrandField intronMotif \
  --outSAMattrRGline ID:SRR22013634SRR22013635 PL:ILLUMINA SM:SRR22013634SRR22013635 \
  --readFilesIn \
    SRR22013634/SRR22013634_1.fastq.gz,SRR22013635/SRR22013635_1.fastq.gz \
    SRR22013634/SRR22013634_2.fastq.gz,SRR22013635/SRR22013635_2.fastq.gz \
  --readFilesCommand zcat \
  --outFileNamePrefix SRR22013634SRR22013635/Second. \
  --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
  --outSAMattributes All \
  --outSAMunmapped Within \
  --outSAMtype BAM SortedByCoordinate --outBAMcompression -1 \
  --outFilterMultimapNmax 1 \
  --outFilterMismatchNoverLmax 0.06 \
  --outFilterMatchNmin 30 \
  --outSJfilterOverhangMin 30 10 10 10 \
  --seedSearchStartLmax 30 \
  --alignIntronMin 20 \
  --alignIntronMax 20000 \
  --alignEndsType Local
./STAR_2.7.11b/Linux_x86_64/STAR \
  --runMode alignReads \
  --runThreadN 16 \
  --genomeDir `pwd` \
  --limitBAMsortRAM 10000000000 \
  --outSAMstrandField intronMotif \
  --outSAMattrRGline ID:SRR22013636SRR22013637 PL:ILLUMINA SM:SRR22013636SRR22013637 \
  --readFilesIn \
    SRR22013636/SRR22013636_1.fastq.gz,SRR22013637/SRR22013637_1.fastq.gz \
    SRR22013636/SRR22013636_2.fastq.gz,SRR22013637/SRR22013637_2.fastq.gz \
  --readFilesCommand zcat \
  --outFileNamePrefix SRR22013636SRR22013637/Second. \
  --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
  --outSAMattributes All \
  --outSAMunmapped Within \
  --outSAMtype BAM SortedByCoordinate --outBAMcompression -1 \
  --outFilterMultimapNmax 1 \
  --outFilterMismatchNoverLmax 0.06 \
  --outFilterMatchNmin 30 \
  --outSJfilterOverhangMin 30 10 10 10 \
  --seedSearchStartLmax 30 \
  --alignIntronMin 20 \
  --alignIntronMax 20000 \
  --alignEndsType Local


htseq-count -f bam -r pos -s no -t CDS -i gene_name -m union SRR22013626SRR22013627/Second.bam gencode.v38lift37.annotation.gtf > SRR22013626SRR22013627/SRR22013626SRR22013627-htseq-count-matrix.txt
htseq-count -f bam -r pos -s no -t CDS -i gene_name -m union SRR22013628SRR22013629/Second.bam gencode.v38lift37.annotation.gtf > SRR22013628SRR22013629/SRR22013628SRR22013629-htseq-count-matrix.txt
htseq-count -f bam -r pos -s no -t CDS -i gene_name -m union SRR22013630SRR22013631/Second.bam gencode.v38lift37.annotation.gtf > SRR22013630SRR22013631/SRR22013630SRR22013631-htseq-count-matrix.txt
htseq-count -f bam -r pos -s no -t CDS -i gene_name -m union SRR22013632SRR22013633/Second.bam gencode.v38lift37.annotation.gtf > SRR22013632SRR22013633/SRR22013632SRR22013633-htseq-count-matrix.txt
htseq-count -f bam -r pos -s no -t CDS -i gene_name -m union SRR22013634SRR22013635/Second.bam gencode.v38lift37.annotation.gtf > SRR22013634SRR22013635/SRR22013634SRR22013635-htseq-count-matrix.txt
htseq-count -f bam -r pos -s no -t CDS -i gene_name -m union SRR22013636SRR22013637/Second.bam gencode.v38lift37.annotation.gtf > SRR22013636SRR22013637/SRR22013636SRR22013637-htseq-count-matrix.txt


Rscript bulk-rnaseq-run-through.R





