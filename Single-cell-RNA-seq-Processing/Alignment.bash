#Select fastq files and use script similar to this. Fastq files are located in the box
cellranger count --id=Sample \
--refdata-gex-GRCh38-2020-A/ \
--fastqs=Fastq Files
--expect-cells=8000 \
--localmem=80 \
--localcores=16
