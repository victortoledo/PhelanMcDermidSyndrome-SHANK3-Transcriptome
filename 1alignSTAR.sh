STAR/STAR-2.7.10b/bin/Linux_x86_64_static/STAR --genomeDir STAR/STAR_hg38_index \
--runThreadN 12 \
--readFilesIn ${INPUT_FOLDER}/${i}_R1_paired.fq.gz ${INPUT_FOLDER}/${i}_R2_paired.fq.gz \
--readFilesCommand zcat --outFileNamePrefix ${OUTPUT_FOLDER}/${i}_STAR_GRCh38 \
--outSAMunmapped Within \
--outFilterType BySJout \
--outSAMattributes NH HI AS NM MD \
--outFilterMultimapNmax 20 \
--outFilterMismatchNmax 999 \
--outFilterMismatchNoverLmax 0.04 \
--alignIntronMin 20 \
--alignIntronMax 1000000 \
--alignMatesGapMax 1000000 \
--alignSJoverhangMin 8 \
--alignSJDBoverhangMin 1 \
--sjdbScore 1 \
--genomeLoad NoSharedMemory \
--outSAMtype BAM Unsorted \
--quantMode TranscriptomeSAM \
--outSAMheaderHD @HD VN:1.4 SO:unsorted
