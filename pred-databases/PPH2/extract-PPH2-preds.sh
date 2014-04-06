touch polyphen-2.2.2-whess-2011_12/00_PPH2.tab

# PPH2 database. This script combines the features.tab and scores.tab files and extracts only the necessary information into the 00_PPH2.tab
for f in polyphen-2.2.2-whess-2011_12/*.features.tab
do 
 _name=${f%.features.tab}
 paste $f ${_name}'.scores.tab' | cut -d$'\t' -f1-10,16-20,56,58,62,64 >> polyphen-2.2.2-whess-2011_12/00_PPH2.tab
done

# This splits the 00_PPH2.tab into chromosomal files
for ((i = 1; i <= 22; i++))
do 
  grep -a -m 1 '#chr_pos' polyphen-2.2.2-whess-2011_12/00_PPH2.tab | sed 's/ //g' >> 'PPH2.chr'$i'.tab'
  grep -a 'chr'$i':' polyphen-2.2.2-whess-2011_12/00_PPH2.tab | sed 's/ //g' >> 'PPH2.chr'$i'.tab'
done

grep -a -m 1 '#chr_pos' polyphen-2.2.2-whess-2011_12/00_PPH2.tab | sed 's/ //g' >> 'PPH2.chrX.tab'
grep -a 'chrX' polyphen-2.2.2-whess-2011_12/00_PPH2.tab | sed 's/ //g' >> 'PPH2.chrX.tab'

grep -a -m 1 '#chr_pos' polyphen-2.2.2-whess-2011_12/00_PPH2.tab | sed 's/ //g' >> 'PPH2.chrY.tab'
grep -a 'chrY' polyphen-2.2.2-whess-2011_12/00_PPH2.tab | sed 's/ //g' >> 'PPH2.chrY.tab'
