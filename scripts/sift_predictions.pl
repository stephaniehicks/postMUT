#!/usr/bin/perl -w
use warnings; 
use strict; 
use File::Slurp; 

# This script takes in the input file for SIFT and returns the functional predictions. 
# to run script: 
# perl scripts/sift_predictions.pl input_files/SIFT_input/<filename>-SIFT-in-hg19.txt output_files/SIFT_output/<filename>-SIFT-out.txt

my $finput = shift @ARGV;  
my $foutput = shift @ARGV;
my @temp;
my %counter;
my @values = qw(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y); 
my %final; 

my @order = read_file($finput); 

# Slurp the SIFT-in-hg19.txt file into a hash
my %siftID;
open FINPUT, $finput;
while(<FINPUT>){
	chomp $_; 
	$siftID{$_} = $_; 
}
close(FINPUT); 

# Open destination folder
open FOUTPUT, ">$foutput" or die "Can't read file $foutput: $!\n";

foreach my $vals (@values){
	open FILETEMP, "/Users/stephaniehicks/Documents/statgen/postMUT/pred-databases/SIFT/SIFT.chr$vals.tab"; 

 	while (my $outerfile = <FILETEMP>){
		next if($outerfile =~ /CHROM/); 
		chomp $outerfile;
		my @outer = split(/\t/, $outerfile); 
		my $ID = join(",",substr($outer[0],3),$outer[1],1,join("\/",$outer[2],$outer[3]));
		if( defined($siftID{$ID}) ){
			$final{$ID} = join("\t", $ID, $outer[4]);
			delete $siftID{$ID}; 
		}
 	}
	close(FILETEMP);
}

while (my ($key, $value) = each(%siftID)) {
	$final{$key} = $siftID{$key};       
}

foreach my $ord (@order){ 
	chomp $ord; 
	print FOUTPUT $final{$ord}, "\n";
}
