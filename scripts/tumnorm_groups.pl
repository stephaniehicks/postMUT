#!/usr/bin/perl -w
use warnings; 
use strict; 
use File::Slurp; 
use Array::Utils qw(:all);

# This script takes in a Normal and Tumor .tab file and produced three files: intersection and two unique sets of mutations to the Normal and Tumor
# perl scripts/tumnorm_groups.pl 2 scratch/010-HCC1187N.txt scratch/011-HCC1187T.txt scratch/010-HCC1187N-Annotated.txt scratch/010-HCC1187T-Annotated.txt

my $numSAMPLES = shift @ARGV; # pops top element in command line into $source

if($numSAMPLES == 2){
	my $sourceN = shift @ARGV; 
	my $sourceT = shift @ARGV; 
	my $destination02 = shift @ARGV;
	my $destination03 = shift @ARGV;

	my @des00 = read_file($sourceN); 
	my @des01 = read_file($sourceT);

	my (%Norm_HomREF, %Norm_HetREF, %Norm_HomALT, %Norm_HetCOMP, %NormFile); 
	my (%Tum_HomREF, %Tum_HetREF, %Tum_HomALT, %Tum_HetCOMP, %TumFile);

	open SOURCEN, $sourceN; 
	while (<SOURCEN>){ 
		next if ($_ =~ /#/); 
		chomp $_; 

		my @normID = split /\t/, $_;
		my @geno = split /\//, $normID[3];
		my $headID = join("\t", $normID[0], $normID[1], $normID[2]); 

		if($geno[0] eq $geno[1]){ 
			if($geno[0] ne $normID[2]){ 
				$Norm_HomALT{$headID} = $_; 
			} else { 
				$Norm_HomREF{$headID} = $_; 
			}
		} else { 
			if( ($normID[2] eq $geno[0]) || ($normID[2] eq $geno[1]) ){
				$Norm_HetREF{$headID} = $_; 
			} else { 
				$Norm_HetCOMP{$headID} = $_; 
			}
		}
	}
	close SOURCEN; 

	open SOURCET, $sourceT; 
	while (<SOURCET>){ 
		next if ($_ =~ /#/); 
		chomp $_; 

		my @tumID = split /\t/, $_;
		my @geno = split /\//, $tumID[3]; 
		my $headID = join("\t", $tumID[0], $tumID[1], $tumID[2]); 

		if($geno[0] eq $geno[1]){ 
			if($geno[0] ne $tumID[2]){ 
				$Tum_HomALT{$headID} = $_; 
			} else { 
				$Tum_HomREF{$headID} = $_; 
			}
		} else { 
			if( ($tumID[2] eq $geno[0]) || ($tumID[2] eq $geno[1]) ){
				$Tum_HetREF{$headID} = $_; 
			} else { 
				$Tum_HetCOMP{$headID} = $_; 
			}
		}
	}
	close SOURCET; 


	foreach (keys %Norm_HomREF){ 
		$NormFile{$_} = 'Ref'; 
	}

	foreach (keys %Norm_HetREF){ 
		$NormFile{$_} = 'OnlyNormal' if (not exists $Tum_HetREF{$_} || not exists $Tum_HomALT{$_} || not exists $Tum_HetCOMP{$_} ); 
		$NormFile{$_} = 'Germline' if (exists $Tum_HetREF{$_});
		$NormFile{$_} = 'LOH' if (exists $Tum_HomALT{$_} || exists $Tum_HetCOMP{$_}); 
	}

	foreach (keys %Norm_HomALT){
		$NormFile{$_} = 'OnlyNormal' if (not exists $Tum_HetREF{$_} || not exists $Tum_HetCOMP{$_} || not exists $Tum_HomALT{$_}); 
		$NormFile{$_} = 'rLOH' if (exists $Tum_HetREF{$_} || exists $Tum_HetCOMP{$_}); 
		$NormFile{$_} = 'Germline' if (exists $Tum_HomALT{$_}); 	
	}

	foreach (keys %Norm_HetCOMP){
		$NormFile{$_} = 'Unknown' if ( exists $Tum_HetREF{$_} || exists $Tum_HomALT{$_} || not exists $Tum_HetCOMP{$_}); 
		$NormFile{$_} = 'Germline' if exists $Tum_HetCOMP{$_}; 
	}



	foreach (keys %Tum_HomREF){ 
		$TumFile{$_} = 'Ref'; 
	}

	foreach (keys %Tum_HetREF){ 
		$TumFile{$_} = 'Somatic' if (not exists $Norm_HetREF{$_} || not exists $Norm_HomALT{$_} || not exists $Norm_HetCOMP{$_} ); 
		$TumFile{$_} = 'Germline' if (exists $Norm_HetREF{$_});
		$TumFile{$_} = 'rLOH' if (exists $Norm_HomALT{$_}); 
		$TumFile{$_} = 'Unknown' if (exists $Norm_HetCOMP{$_}); 
	}

	foreach (keys %Tum_HomALT){
		$TumFile{$_} = 'Somatic' if (not exists $Norm_HetREF{$_} || not exists $Norm_HomALT{$_} || not exists $Norm_HetCOMP{$_}); 
		$TumFile{$_} = 'LOH' if (exists $Norm_HetREF{$_}); 
		$TumFile{$_} = 'Germline' if (exists $Norm_HomALT{$_}); 	
		$TumFile{$_} = 'Unknown' if (exists $Norm_HetCOMP{$_}); 
	}

	foreach (keys %Tum_HetCOMP){
		$TumFile{$_} = 'Somatic' if (not exists $Norm_HetREF{$_} || not exists $Norm_HomALT{$_} || not exists $Norm_HetCOMP{$_} ); 
		$TumFile{$_} = 'LOH' if (exists $Norm_HetREF{$_}); 
		$TumFile{$_} = 'rLOH' if (exists $Norm_HomALT{$_}); 
		$TumFile{$_} = 'Germline' if (exists $Norm_HetCOMP{$_}); 
	}


	open OUT02, ">$destination02" or die "Can't read file $destination02: $!\n";
	open OUT03, ">$destination03" or die "Can't read file $destination03: $!\n";

	open SOURCEN, $sourceN; 
	while (<SOURCEN>){
		if ($_ =~ /#/){
			print OUT02 $_;
			next;
		}
		chomp $_; 

		my @normID = split /\t/, $_;
		my @geno = split /\//, $normID[3];
		my $headID = join("\t", $normID[0], $normID[1], $normID[2]); 

		push (@normID, $NormFile{$headID}); 
		print OUT02 join("\t", @normID), "\n"; 
	} 
	close SOURCEN; 

	open SOURCET, $sourceT;
	while (<SOURCET>){
		if ($_ =~ /#/){
			print OUT03 $_;
			next;
		}
		chomp $_; 

		my @tumID = split /\t/, $_;
		my @geno = split /\//, $tumID[3];
		my $headID = join("\t", $tumID[0], $tumID[1], $tumID[2]); 

		push @tumID, $TumFile{$headID}; 
		print OUT03 join("\t", @tumID), "\n"; 
	} 
	close SOURCET; 

}