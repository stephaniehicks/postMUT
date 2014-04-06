#!/usr/bin/perl -w
use warnings; 
use strict; 

# This script takes in SEARCHLIST which contains the CHROM	POS	REF_AA	MUT_AA and converts it to a PPH2, SIFT, Xvar batch input format
# 0 = PPH2 Batch format
# 1 = SIFT Batch format
# 2 = Xvar Batch format
# perl convert_algorithms.pl 0 data/FCP416.tab
# perl convert_algorithms.pl 1 data/FCP416.tab
# perl convert_algorithms.pl 2 data/FCP416.tab
# input.tab is in the form: 	CHROM	POS	REF	GENOTYPE

my $algorithm = shift @ARGV; 
my $source = shift @ARGV; # pops top element in command line into $source
my $pattern1 = ">";
my $pattern2 = "#";
my (@finalarray, @subgeno);  

# Open the search list and output the information 
open SEARCHLIST, $source; 
 	while (my $searchl = <SEARCHLIST>){
		next if($searchl =~ /$pattern2/);

		chomp $searchl; 
		@finalarray = split(/\t/, $searchl);
		@subgeno = split(/\//, $finalarray[3]); 
 
		if($algorithm == 0){
	 		print join("\t", join(":", join("", "chr", $finalarray[0]), $finalarray[1]), join("\/", $finalarray[2], $subgeno[1])), "\n"; 
			# print join("\t", join(":", join("", $finalarray[0]), $finalarray[1]), join("\/", $finalarray[2], $subgeno[1])), "\n"; 
		}
		
		if($algorithm == 1){ 
			# my $siftfinalarray = substr($finalarray[0], 3);
	 		print join(",", $finalarray[0], $finalarray[1], 1, join("\/", $finalarray[2], $subgeno[1])), "\n"; 
		}

		if($algorithm == 2){ 
			# my $xvarfinalarray = substr($finalarray[0], 3);
	 		print join(",","hg19", $finalarray[0], $finalarray[1], $finalarray[2], $subgeno[1]), "\n"; 
		}		
	} 
close(SEARCHLIST); 

