rm -r input_files/postMUT_input/*.txt

# Extract the annotated file in 'scratch/*.tab' and create genotype file
for f in scratch/*.tab
do
  _temp_file=${f#scratch/}	
  awk -F"\t" '{if (NR!=1) {print $1 "\t" $2 "\t" $3 "\t" substr($4,3,1) "\t" $4}}' $f > 'input_files/postMUT_input/'${_temp_file%.tab}'_geno.txt'
done

# Using new genotype file, collapse all functional predictions from SIFT, PPH2 and Xvar into one file as input in postMUT.  Files are placed in 'input_files/postMUT_input/*-postMUT-in.txt'
for f in input_files/postMUT_input/*_geno.txt
do 
  _temp_file=${f#input_files/postMUT_input/}
  _temp=${_temp_file%_geno.txt}
  perl scripts/collapse.pl notum $f 'output_files/SIFT_output/'${_temp}'-SIFT-out.txt' 'output_files/Xvar_output/'${_temp}'-Xvar-out.txt' 'output_files/PPH2_output/'${_temp}'-PPH2-out.txt' > 'input_files/postMUT_input/'${_temp}'-postMUT-in.txt'
done