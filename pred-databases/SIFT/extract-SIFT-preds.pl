#!/usr/bin/perl -w
use warnings; 
use strict; 
use DBI;
use DBD::SQLite; 

# To run script: 
# perl scripts/extract-SIFT-preds.pl

my @values = qw(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y); 

foreach my $vals (@values){

	# Open destination folder
	open(FOUTPUT,'>', "SIFT.chr$vals.tab");
	print FOUTPUT join("\t", "CHROM", "POS", "REF", "ALT", "SIFT_score"), "\n";

	# Connect to SQLite database
	my $db_chr = DBI->connect( "dbi:SQLite:dbname=Human_CHR$vals.sqlite","", "",{ RaiseError => 1, AutoCommit => 1 });

	# Prepare SQLite statement to query database
	my $sth_db_name = $db_chr -> prepare("SELECT name FROM sqlite_master WHERE type='table' ORDER BY name");
	$sth_db_name -> execute();
	
	# For each table in the SQLite database, print chromosome, coordinate position, REF & ALT nucleotide, and SIFT score
	while (my $rows = $sth_db_name -> fetchrow()){
		foreach my $table_name ( $rows ){
			
			my $sth_db_chr = $db_chr -> prepare("SELECT CHR,COORD2,NT1,NT2,SCORE FROM $table_name");
			$sth_db_chr -> execute();
			
			my @rows;
			my @query_result;
			while (@rows = $sth_db_chr -> fetchrow_array() ){
				push @query_result, join("\t", @rows);
			}

			# Print output
			foreach my $row (@query_result){
				chomp $row;
				my @arrayID = split(/\t/, $row); 
				if( defined $arrayID[4]){ 
					print FOUTPUT $row, "\n";
				}
			}
		}
	}
	close(FOUTPUT); 
}

