#!/usr/bin/perl -w
use warnings; 
use strict; 
use File::Slurp; 

# This script takes in the input file for Xvar and returns the functional predictions. 
# to run script: 
# perl scripts/xvar_predictions.pl input_files/Xvar_input/<filename>-Xvar-in-hg19.txt output_files/Xvar_output/<filename>-Xvar-out.txt

my $finput = shift @ARGV;  
my $foutput = shift @ARGV;
my @temp;
my %counter;
my @values = qw(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y); 
my %final; 

my @order = read_file($finput); 

# Slurp the Xvar-in-hg19.txt file into a hash
my %xvarID;
open FINPUT, $finput; 
while(<FINPUT>){
	chomp $_; 
	$xvarID{$_} = $_; 
}
close(FINPUT); 

# Open destination folder
open FOUTPUT, ">$foutput" or die "Can't read file $foutput: $!\n";

foreach my $vals (@values){
 	open FILETEMP, "/Users/stephaniehicks/Documents/statgen/postMUT/pred-databases/Xvar/MA.hg19/MA.chr$vals.txt"; 

 	while (my $outerfile = <FILETEMP>){
		next if($outerfile =~ /@/); 
		chomp $outerfile;
		my @outer = split(/\t/, $outerfile);
		if( defined($xvarID{$outer[0]}) ){
			$final{$outer[0]} = $outerfile; 
			delete $xvarID{$outer[0]}; 
		}
 	}
	close(FILETEMP);
}

while (my ($key, $value) = each(%xvarID)) {
	$final{$key} = $xvarID{$key};       
}

foreach my $ord (@order){ 
	chomp $ord; 
	print FOUTPUT $final{$ord}, "\n";
}
