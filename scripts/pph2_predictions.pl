#!/usr/bin/perl -w
use warnings; 
use strict; 
use File::Slurp; 

# This script takes in the input file for PPH2 and returns the functional predictions. 
# to run script: 
# perl scripts/pph2_predictions.pl input_files/PPH2_input/<filename>-PPH2-in-hg19.txt output_files/PPH2_output/<filename>-PPH2-out.txt

my $finput = shift @ARGV;  
my $foutput = shift @ARGV;
my @temp;
my %counter;
my @values = qw(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y); 
my %final; 

my @order = read_file($finput); 

# Slurp the PPH2-in-hg19.txt file into a hash
my %pph2ID;
open FINPUT, $finput; 
while(<FINPUT>){
	chomp $_; 
	$pph2ID{$_} = $_; 
}
close(FINPUT); 

# Open destination folder
open FOUTPUT, ">$foutput" or die "Can't read file $foutput: $!\n";

foreach my $vals (@values){
 	open FILETEMP, "/Users/stephaniehicks/Documents/statgen/postMUT/pred-databases/PPH2/PPH2.chr$vals.tab"; 

 	while (my $outerfile = <FILETEMP>){
		next if($outerfile =~ /#/); 
		chomp $outerfile;
		my @outer = split(/\t/, $outerfile);
		my $ID = join("\t",$outer[0],$outer[1]);
		if( defined($pph2ID{$ID}) ){
			$final{$ID} = $outerfile; 
			delete $pph2ID{$ID}; 
		}
 	}
	close(FILETEMP);
}

while (my ($key, $value) = each(%pph2ID)) {
	$final{$key} = $pph2ID{$key};       
}

foreach my $ord (@order){ 
	chomp $ord; 
	print FOUTPUT $final{$ord}, "\n";
}
