# Extract missense mutations from .vcf files for ALL file in input_files director -> scratch directory
for f in input_files/vcf_files/*.vcf
do
  egrep '##|#CHROM|MISSENSE' $f > 'scratch/'${f#input_files/vcf_files/}
done

# gzip all files in scratch directory
for f in scratch/*.vcf
do
  gzip $f
done

for f in scratch/*.gz
do
  mv -f $f $f'.Z'
done

# Convert VCF to genotype format
for f in scratch/*.Z
do
  zcat $f | vcf-to-tab > ${f%.vcf.gz.Z}'.tab'
done

rm scratch/*.gz.Z

# Create labels: 'Somatic', 'Germline', 'LOH', 'OnlyNormal' for each mutation 
# NOTE: Each of the normal and tumor file names need to be individually added to the bash script 
# perl scripts/tumnorm_groups.pl 2 scratch/<normal-filename>.tab scratch/<tumor-filename>.tab scratch/<normal-filename>-Ann.tab scratch/<tumor-filename>-Ann.tab

perl scripts/tumnorm_groups.pl 2 scratch/Normal-snpeffect_Cleaveland038_normal_recal.raw.tab scratch/Tumor-snpeffect_Cleaveland038_tumor_recal.raw.tab scratch/Normal-snpeffect_Cleaveland038_normal_recal.raw-Ann.tab scratch/Tumor-snpeffect_Cleaveland038_tumor_recal.raw-Ann.tab

# Convert to PPH2, SIFT, Xvar batch input format
for f in scratch/*-Ann.tab
do
  _out_file=${f#scratch/}
  perl scripts/convert_algorithms.pl 0 $f > 'input_files/PPH2_input/'${_out_file%-Ann.tab}'-PPH2-in-hg19.txt'
  perl scripts/convert_algorithms.pl 1 $f > 'input_files/SIFT_input/'${_out_file%-Ann.tab}'-SIFT-in-hg19.txt'
  perl scripts/convert_algorithms.pl 2 $f > 'input_files/Xvar_input/'${_out_file%-Ann.tab}'-Xvar-in-hg19.txt'
done

find scratch -type f -not -name '*-Ann.tab' | xargs rm