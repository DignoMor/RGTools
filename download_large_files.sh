
# K562 GROcap tracks
wget \
    -O large_files/ENCFF775FNU.mn.bw \
    https://www.encodeproject.org/files/ENCFF775FNU/@@download/ENCFF775FNU.bigWig

wget \
    -O large_files/ENCFF565BWR.pl.bw \
    https://www.encodeproject.org/files/ENCFF565BWR/@@download/ENCFF565BWR.bigWig

# K562 GROcap peaks
wget \
    -O large_files/ENCFF156JSS.bedTRE.gz \
    https://www.encodeproject.org/files/ENCFF156JSS/@@download/ENCFF156JSS.bed.gz
gunzip large_files/ENCFF156JSS.bedTRE.gz

# hg38 reference genome
wget \
    -O large_files/hg38.fa.gz \
    https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz
gunzip large_files/hg38.fa.gz

wget \
    -O large_files/hg38.chrom.sizes \
    https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.chrom.sizes

