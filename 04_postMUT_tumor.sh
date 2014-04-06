# Run postMUT.R; saves results in output_files/postMUT_output
for f in input_files/postMUT_input/*-postMUT-in.txt
do 
  _temp_file=${f#input_files/postMUT_input/}
  _temp=${_temp_file%-postMUT-in.txt}
  ./scripts/postMUT_tumor.R $_temp 
done

