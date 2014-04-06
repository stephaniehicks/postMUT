# Run script to obtain all functional predictions from SIFT (~ 10min / file)
for f in input_files/SIFT_input/*.txt
do 
  _temp_file=${f#input_files/SIFT_input/}
  _temp=${_temp_file%-SIFT-in-hg19.txt}
  perl scripts/SIFT_predictions.pl $f 'output_files/SIFT_output/'${_temp}'-SIFT-out.txt'
done

# Run script to obtain all functional predictions from PPH2 (~ 30min / file)
for f in input_files/PPH2_input/*.txt
do 
  _temp_file=${f#input_files/PPH2_input/}
  _temp=${_temp_file%-PPH2-in-hg19.txt}
  perl scripts/pph2_predictions.pl $f 'output_files/PPH2_output/'${_temp}'-PPH2-out.txt'
done

# Run script to obtain all functional predictions from Xvar (~ 6-8min / file)
for f in input_files/Xvar_input/*.txt
do 
  _temp_file=${f#input_files/Xvar_input/}
  _temp=${_temp_file%-Xvar-in-hg19.txt}
  perl scripts/xvar_predictions.pl $f 'output_files/Xvar_output/'${_temp}'-Xvar-out.txt'
done
