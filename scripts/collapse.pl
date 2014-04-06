#!/usr/bin/perl
use strict; 

# This script takes in multiple text files: geno (complete list of variants), and SIFT, PPH2, Xvar predictions. 
# to run script: 
# perl scripts/collapse.pl tum input_files/postMUT_input/<filename>_geno.txt output_files/SIFT_output/<filename>-SIFT-out.txt output_files/Xvar_output/<filename>-Xvar-out.txt output_files/PPH2_output/<filename>-PPH2HD-out.txt > input_files/postMUT_input/<filename>-postMUT-in.txt

my $ftum = shift @ARGV;
my $ffull = shift @ARGV;
my $fsift= shift @ARGV;
my $fxvar= shift @ARGV;
my $fpph2HD= shift @ARGV;
my (%sift_pred, %sift_score, %xvar_pred, %xvar_score, %pph2hd_pred, %pph2hd_score);  

if($ftum eq "tum"){
	print join("\t", "CHROM", "POS", "REF_AA", "MUT_AA", "GENOTYPE", "LABEL", "SIFT_pred", "SIFT_score", "Xvar_pred", "Xvar_score", "PPH2_HD_pred", "PPH2_HD_score"), "\n";
} 

if($ftum eq "notum"){
	print join("\t", "CHROM", "POS", "REF_AA", "MUT_AA", "GENOTYPE", "SIFT_pred", "SIFT_score", "Xvar_pred", "Xvar_score", "PPH2_HD_pred", "PPH2_HD_score"), "\n";
}

# slurp the sift predictions into a hash
open FSIFT, $fsift;
	while (<FSIFT>){
		chomp $_;
 		my @siftID = split /\t/, $_;

		if(defined $siftID[1]){
			my @sift_coord = split(",", $siftID[0]);
			my @sift_geno = split("\/", $sift_coord[3]); 
			my $sift_in = join(",", $sift_coord[0], $sift_coord[1], $sift_geno[0], $sift_geno[1]); 
			chomp $sift_in;

			$sift_score{$sift_in} = $siftID[1];
			if($siftID[1] < 0.05){ 
				$sift_pred{$sift_in} = "DAMAGING";
			} else { 
				$sift_pred{$sift_in} = "TOLERATED"; 
			}
		} 
  	}
close(FSIFT); 

# slurp the Xvar predictions into a hash
open FXVAR, $fxvar;
	while (<FXVAR>){
		# next if($_ =~ /issue/);
		
		chomp $_;
 		my @xvarID = split /\t/, $_;
		# shift @xvarID; 

		if(defined $xvarID[6]){
			my @xvar_coord = split(",", $xvarID[0]);
			shift @xvar_coord;  
			my $xvar_in = join(",", @xvar_coord); 
			chomp $xvar_in; 
			$xvar_pred{$xvar_in} = $xvarID[6];
			$xvar_score{$xvar_in} = $xvarID[7];
		} 
  	}
close(FXVAR); 

# slurp the PPH2HD predictions into a hash
open FPPH2HD, $fpph2HD;
	while (<FPPH2HD>){
		
		chomp $_;
		my @pph2HDID = split /\t/, $_;
 		my @org = split /\t/, $_;

		my @header = split(/:/, $pph2HDID[0]);
		my @aa = split(/\//, $pph2HDID[1]); 
		
		my $chrom = substr($header[0],3); 
		my $pph2hd_in = join(",", $chrom, $header[1], $aa[0], $aa[1]);
		chomp $pph2hd_in;

		$pph2hd_pred{$pph2hd_in} = $pph2HDID[15];
		$pph2hd_score{$pph2hd_in} = $pph2HDID[16];
	}
close(FPPH2HD); 

# Open Complete list of variants file
open FFULL, $ffull; 
	while (<FFULL>){
		my @finalID = split /\t/, $_;
		my $temp = join(",", $finalID[0], $finalID[1], $finalID[2], $finalID[3]); 
		chomp $_;
			chomp @finalID;  
			print join("\t", join("\t", @finalID), $sift_pred{$temp}, $sift_score{$temp}, $xvar_pred{$temp}, $xvar_score{$temp}, $pph2hd_pred{$temp}, $pph2hd_score{$temp}), "\n";
	} 
close FFULL;

