#!c:/perl/bin/perl.exe


# script to parse the results of reciprocal blast
# proj 3983
# 16 xii 2017

#v. 2.1: solving a problem with reading from the current dir
# perl integrate_homologues_counts_v0.2.pl –samples 3983_sample_fastq.txt –dir_annotated . –dir_counts . –outfile ./tst3

# AgataSMBP:counts agatasmialowska$ perl integrate_homologues_counts_v0.2.1.pl –samples 3983_sample_fastq.txt –dir_annotated . –dir_counts . –outfile ./tst3
# parameters: –samples 3983_sample_fastq.txt –dir_annotated . –dir_counts . –outfile ./tst3
# Use of uninitialized value $path2outfile in concatenation (.) or string at
# 	integrate_homologues_counts_v0.2.1.pl line 79 (#1)
#     (W uninitialized) An undefined value was used as if it were already
#     defined.  It was interpreted as a "" or a 0, but maybe it was a mistake.
#     To suppress this warning assign a defined value to your variables.
    
#     To help you figure out what was undefined, perl will try to tell you
#     the name of the variable (if any) that was undefined.  In some cases
#     it cannot do this, so it also tells you what operation you used the
#     undefined value in.  Note, however, that perl optimizes your program
#     anid the operation displayed in the warning may not necessarily appear
#     literally in your program.  For example, "that $foo" is usually
#     optimized into "that " . $foo, and the warning will refer to the
#     concatenation (.) operator, even though there is no . in
#     your program.
    
# Use of uninitialized value $path2outfile in concatenation (.) or string at
# 	integrate_homologues_counts_v0.2.1.pl line 80 (#1)
# Use of uninitialized value $path2outfile in concatenation (.) or string at
# 	integrate_homologues_counts_v0.2.1.pl line 81 (#1)
# Uncaught exception from user code:
# 	invalid top directory at /System/Library/Perl/5.18/File/Find.pm line 472.
# 	File::Find::_find_opt('HASH(0x7f9ed4803ff0)', undef) called at /System/Library/Perl/5.18/File/Find.pm line 1079
# 	File::Find::find('HASH(0x7f9ed4803ff0)', undef) called at integrate_homologues_counts_v0.2.1.pl line 94


#change in line 
#my @files = glob "$inputDir*";
#http://www.perlmonks.org/?node_id=926958

use warnings;
use strict;
use diagnostics;
use File::Basename;
use Getopt::Long;
use List::MoreUtils qw(any uniq);
use File::Find;


my $script_name="integrate_homologues_counts.pl";



if ($ARGV[0]eq qw '-h'){
	print "please provide arguments for: \n perl $script_name [arguments]\n";
	print "arguments:\n";
	print "--samples: /path/to/file/with/sample_description: sample_number - tab - sample_id\n";
	print "--dir_annotated: /path/to/dir/with/results_annotated_with_arp7_gene_id\n";
	print "--dir_counts: /path/to/dir/with/count_tables\n";
	print "--outfile: /path/to/outfile\n\n";
}

else{
	my $parameters=join(' ', @ARGV);
	print "parameters: $parameters\n";
	
	#commandline parsing for parameters
	GetOptions(
		'samples=s'		=>	\(my $samples_file),
		'dir_annotated=s'		=>	\(my $dir_annotated),
		'dir_counts=s'		=>	\(my $dir_counts),
		'outfile=s'	=>	\(my $path2outfile)
	) or die "Error in command line arguments";


	#outfiles
	my $path2outfile1="$path2outfile\.counts"; #just a counts table
	my $path2outfile2="$path2outfile\.homologues"; #file with homologues info
	my $path2outfile3="$path2outfile\.homologues_multiple"; #file with homologues info, multiple genes assigned to the same ARP7 consensus gene group


	#############################################
	#get all count table files
	my @counts_files;
	sub find_files_counts {
    	my $F = $File::Find::name;

    	if ($F =~ /.+corset-counts.filtered.header.txt/ ) { #Pfa_corset-counts.filtered.header.txt
    		push @counts_files, $F;
    	}
	}
	find({ wanted => \&find_files_counts, no_chdir=>1}, $dir_counts);


	#############################################
	#get all files with annotation
	my @annot_files;
	sub find_files_annot {
    	my $F = $File::Find::name;

    	if ($F =~ /.+\.arp7.reblast.clust.annot.txt/ ) { #psa.arp7.reblast.clust.annot.txt
    		push @annot_files, $F;
    	}
	}
	find({ wanted => \&find_files_annot, no_chdir=>1}, $dir_annotated);


	#############################################
	# get sample info
	my %samples;
	open (INFILE_SAMPLES, "<$samples_file") or die "Cannot open input file $samples_file $!"; 
	while(<INFILE_SAMPLES>){
		chomp $_;
		if ($_ =~m/(\d+)\t(\w+)/){
			my ($sample_id,$sample_name)=($1,$2);
			$samples{$sample_id}=$sample_name;
		}
	}
	close(INFILE_SAMPLES);




	#############################################
	#outfiles


	open (OUTFILE2, ">$path2outfile2") or die "Cannot open output file $path2outfile2 $!"; 

	my $header2="ARP7_gene_id\tARP7_description\tNxa_denovo_geneID\tNxa_clusterID\tPca_denovo_geneID\tPca_clusterID\tPfa_denovo_geneID\tPfa_clusterID\tPsa_denovo_geneID\tPsa_clusterID";
	print OUTFILE2 "$header2\n";

	open (OUTFILE3, ">$path2outfile3") or die "Cannot open output file $path2outfile3 $!"; 

	my $header3="ARP7_gene_id\tARP7_description\tNxa_n\tPca_n\tPfa_n\tPsa_n\tNxa_denovo_geneID::Nxa_clusterID\tPca_denovo_geneID::Pca_clusterID\tPfa_denovo_geneID::Pfa_clusterID\tPsa_denovo_geneID::Psa_clusterID";
	print OUTFILE3 "$header3\n";



	#############################################
	#process annoatated results

	my %arp7_descriptions;
	my %annotations;

	foreach my $path (@annot_files){
		my($infile_name, $dirs, $suffix) = fileparse($path);
		print "processing file $infile_name\n";
		open(FILE,"<$path") or die "Cannot open input file $path $!";

		#parse file name to get the species
		my $species;
		if ($infile_name =~m/nxa/i){
			$species="Nxa";
		}elsif ($infile_name =~m/pca/i){
			$species="Pca";
		}elsif ($infile_name =~m/psa/i){
			$species="Psa";
		}elsif ($infile_name =~m/pfa/i){
			$species="Pfa";
		}

#trinity_gene_id	corset_cluster_id	reblast_confirmed_arp7	arp7_consensus_group	arp7_description
#TRINITY_DN27327_c0_g1	Cluster-553.0	PogchaEGm007091t1	ARP7f_G1952	Leucine-rich repeat serine/threonine-protein kinase

		while(<FILE>){
			chomp $_;
			if ($_ !~m/trinity_gene_id/){#header line
				my @line=split/\t/,$_;
				my $trinity_gene_id=$line[0];
				my $corset_cluster_id=$line[1];
				my $reblast_hit=$line[2];
				my $arp7_consensus_group=$line[3];
				my $arp7_description=$line[4];

				$arp7_descriptions{$arp7_consensus_group}=$arp7_description;

				#my $annot="$trinity_gene_id::$corset_cluster_id::$reblast_hit";
				my $annot="$trinity_gene_id::$corset_cluster_id"; 

				#$annotations{$arp7_consensus_group}{$species}=$annot;#this does not work; multiple transcripts / clusters with the same arp7 consensus group exist (e.g. ARP7f_G0)
				push @{$annotations{$arp7_consensus_group}{$species}}, $annot;
				
			}
		}
		close(FILE);
	}



	#############################################
	#process annoatated results -> output file *.homologues

	
	my %clusterID_arp7cg; #for easier access downstream, should only contain checked annotations

	foreach my $arp7_consensus_group (sort keys %annotations){
		
		my $arp7_description=$arp7_descriptions{$arp7_consensus_group};

		my $nxa_denovo;
		my $nxa_clust;
		my $pca_denovo;
		my $pca_clust;
		my $pfa_denovo;
		my $pfa_clust;
		my $psa_denovo;
		my $psa_clust;
		my %group_species_multiple;
		foreach my $species (keys $annotations{$arp7_consensus_group}){
			
			my @hits=@{$annotations{$arp7_consensus_group}{$species}};

			#identify (and collect info) if the $arp7_consensus_group has multiple trinity genes assigned for any species
			if ( scalar (@hits)>1 ){ #multiple transcripts annotated with one ARP7 consensus group
				#my $number=scalar(@hits);
				my $genes_clusters=join',',@hits;
				my $species_info="$genes_clusters";
				$group_species_multiple{$species}=$species_info;
			}
			elsif ( scalar (@hits)==1 ){ #unique pairing transcript - arp7 consensus group
				my $info=$hits[0];
				my ($trinity_gene_id,$corset_cluster_id)=split/::/,$info;

				if ($species eq qw /Nxa/){
				 	$nxa_denovo=$trinity_gene_id;
				 	$nxa_clust=$corset_cluster_id;
				 }
				if ($species eq qw /Pca/){
				 	$pca_denovo=$trinity_gene_id;
				 	$pca_clust=$corset_cluster_id;
				 }
				if ($species eq qw /Pfa/){
				 	$pfa_denovo=$trinity_gene_id;
				 	$pfa_clust=$corset_cluster_id;
				 }
				if ($species eq qw /Psa/){
				 	$psa_denovo=$trinity_gene_id;
				 	$psa_clust=$corset_cluster_id;
				 }
			}
		}

		#check if the $arp7_consensus_group is present in each species
		unless (defined $pca_denovo){
			$pca_denovo="NA";
			$pca_clust="NA";
		}
		unless (defined $psa_denovo){
			$psa_denovo="NA";
			$psa_clust="NA";
		}
		unless (defined $pfa_denovo){
			$pfa_denovo="NA";
			$pfa_clust="NA";
		}
		unless (defined $nxa_denovo){
			$nxa_denovo="NA";
			$nxa_clust="NA";
		}

		#check if $arp7_consensus_group has multiple genes assigned in any of the species
		my @seen; #this is only for checks if multiple pairings exist in an easy way
		if ( (exists $group_species_multiple{"Nxa"}) | ( exists $group_species_multiple{"Pca"}) | (exists $group_species_multiple{"Psa"}) | ( exists $group_species_multiple{"Pfa"}) ){
			#print "some species_infiles have multiples: $arp7_consensus_group\n"; #wc -l 460
			push @seen,$arp7_consensus_group; 

			my $n_nxa;
			if (exists $annotations{$arp7_consensus_group}{"Nxa"}){
				$n_nxa=scalar( @{$annotations{$arp7_consensus_group}{"Nxa"}} );
			}else{
				$n_nxa="na";
			}

			my $n_pca;
			if (exists $annotations{$arp7_consensus_group}{"Pca"}){
				$n_pca=scalar( @{$annotations{$arp7_consensus_group}{"Pca"}} );
			}else{
				$n_pca="na";
			}

			my $n_pfa;
			if (exists $annotations{$arp7_consensus_group}{"Pfa"}){
				$n_pfa=scalar( @{$annotations{$arp7_consensus_group}{"Pfa"}} );
			}else{
				$n_pfa="na";
			}

			my $n_psa;
			if (exists $annotations{$arp7_consensus_group}{"Psa"}){
				$n_psa=scalar( @{$annotations{$arp7_consensus_group}{"Psa"}} );
			}else{
				$n_psa="na";
			}

			my $line_nxa;
			 if($n_nxa eq qw /1/){
			 	$line_nxa="$nxa_denovo::$nxa_clust";
			 }elsif($n_nxa eq qw /na/){
			 	$line_nxa="na";
			 }else {#($n_nxa>1){
			 	$line_nxa=$group_species_multiple{"Nxa"};
			}
			 
			 my $line_pca;
			 if($n_pca eq qw /1/){
			 	$line_pca="$pca_denovo::$pca_clust";
			 }elsif($n_pca eq qw /na/){
			 	$line_pca="na";
			 }else {#($n_nxa>1){
			 	$line_pca=$group_species_multiple{"Pca"};
			}


			 my $line_pfa;
			 if($n_pfa eq qw /1/){
			 	$line_pfa="$pfa_denovo::$pfa_clust";
			 }elsif($n_pfa eq qw /na/){
			 	$line_pfa="na";
			 }else {#($n_nxa>1){
			 	$line_pfa=$group_species_multiple{"Pfa"};
			}


			 my $line_psa;
			 if($n_psa eq qw /1/){
			 	$line_psa="$psa_denovo::$psa_clust";
			 }elsif($n_psa eq qw /na/){
			 	$line_psa="na";
			 }else {#($n_nxa>1){
			 	$line_psa=$group_species_multiple{"Psa"};
			}

			my $line="$arp7_consensus_group\t$arp7_description\t$n_nxa\t$n_pca\t$n_pfa\t$n_psa\t$line_nxa\t$line_pca\t$line_pfa\t$line_psa";
			print OUTFILE3 "$line\n";

		}

		if (scalar (@seen)==0){
			#print "none has multiples: $arp7_consensus_group\n"; # wc -l 7039

			$clusterID_arp7cg{"Nxa"}{$nxa_clust}=$arp7_consensus_group;
			$clusterID_arp7cg{"Pca"}{$pca_clust}=$arp7_consensus_group;
			$clusterID_arp7cg{"Pfa"}{$pfa_clust}=$arp7_consensus_group;
			$clusterID_arp7cg{"Psa"}{$psa_clust}=$arp7_consensus_group;

			my $line="$arp7_consensus_group\t$arp7_description\t$nxa_denovo\t$nxa_clust\t$pca_denovo\t$pca_clust\t$pfa_denovo\t$pfa_clust\t$psa_denovo\t$psa_clust";
			print OUTFILE2 "$line\n";
		}

	}
	close(OUTFILE2);
	close(OUTFILE3);


	#############################################
	#process counts
	
	my %sample_ids;
	my %counts;

	#check how many tab-delimited fields a line in the counts table contains
	#add NA when some arp7 consensus groups are missing in individual counts tables
	#otherwise the values are all shifted to the left when merging counts from individual species
	my $sample_count_nxa;
	my $sample_count_pca;
	my $sample_count_pfa;
	my $sample_count_psa;

	foreach my $path (@counts_files){
	 	my($infile_name, $dirs, $suffix) = fileparse($path);
	 	print "processing file $infile_name\n";
	 	open(FILE,"<$path") or die "Cannot open input file $path $!";

	 	#parse file name to get the species
	 	my $species;
	 	if ($infile_name =~m/nxa/i){
	 		$species="Nxa";
	 	}elsif ($infile_name =~m/pca/i){
	 		$species="Pca";
	 	}elsif ($infile_name =~m/psa/i){
	 		$species="Psa";
	 	}elsif ($infile_name =~m/pfa/i){
	 		$species="Pfa";
	 	}

	 	##
		#############################################
		## general strategy:
		## build hash of hashes of arrays
		## %superhash{gene = arp7 consensus group}{species}counts
		## where counts holds the line from the counts table
		## for arp7 groups with data in all species, no multiple assignments
		## additional hash with lines with header (sample order) (sample ids)

		while(<FILE>){
			chomp $_;

			if ($_ !~m/Cluster/){#header line, does not contain the first field, R-style
				$sample_ids{$species}=$_;

				my @header_line=split/\t/,$_;
				my $sample_count=scalar(@header_line);

				$sample_count_nxa=$sample_count if ($species eq qw /Nxa/);
				$sample_count_pca=$sample_count if ($species eq qw /Pca/);
				$sample_count_pfa=$sample_count if ($species eq qw /Pfa/);
				$sample_count_psa=$sample_count if ($species eq qw /Psa/);

				print "$sample_count\n";

			}
			if ($_ =~m/(Cluster\S+)\t(.+)$/){#line with counts
				my($clusterID,$counts)=($1,$2);
				
				#get arp7 consensus group
				if (exists $clusterID_arp7cg{$species}{$clusterID}){
					my $arp7_consensus_group=$clusterID_arp7cg{$species}{$clusterID};

					$counts{$arp7_consensus_group}{$species}=$counts;
				}
				

			}
		}
		close(FILE);
	}
	
	#process the hashes to print the counts table
	open (OUTFILE1, ">$path2outfile1") or die "Cannot open output file $path2outfile1 $!"; 

	#header (sample names)
	my @header_ids;
	foreach my $species (sort keys %sample_ids){
		my $ids=$sample_ids{$species};
		push @header_ids,$ids;
	}

	my @new_headerIDs;
	my @header_ids_numbers;

	foreach my $id_line (@header_ids){
		my @ids_single=split/\t/,$id_line;

		push @header_ids_numbers,@ids_single;
	}

	foreach my $id (@header_ids_numbers){
		my $smpl_name=$samples{$id};
		push @new_headerIDs,$smpl_name;
	}
	my $header_ids=join "\t",@new_headerIDs;
	my $header="ARP7_gene_id\tDescription\t$header_ids";
	print OUTFILE1 "$header\n";


	#counts
	foreach my $arp7_consensus_group (keys %counts){
		my @line_counts;
		
		if (exists $counts{$arp7_consensus_group}{"Nxa"}){
			my $line_counts_species=$counts{$arp7_consensus_group}{"Nxa"};
			push @line_counts,$line_counts_species;
		}else{
			for (1..$sample_count_nxa) { push @line_counts, "na"; }
		}
		if (exists $counts{$arp7_consensus_group}{"Pca"}){
			my $line_counts_species=$counts{$arp7_consensus_group}{"Pca"};
			push @line_counts,$line_counts_species;
		}else{
			for (1..$sample_count_pca) { push @line_counts, "na"; }
		}

		if (exists $counts{$arp7_consensus_group}{"Pfa"}){
			my $line_counts_species=$counts{$arp7_consensus_group}{"Pfa"};
			push @line_counts,$line_counts_species;
		}else{
			for (1..$sample_count_pfa) { push @line_counts, "na"; }
		}

		if (exists $counts{$arp7_consensus_group}{"Psa"}){
			my $line_counts_species=$counts{$arp7_consensus_group}{"Psa"};
			push @line_counts,$line_counts_species;
		}else{
			for (1..$sample_count_psa) { push @line_counts, "na"; }
		}
		
		my $description=$arp7_descriptions{$arp7_consensus_group};
		my $line_cnts=join "\t",@line_counts;
		my $line="$arp7_consensus_group\t$description\t$line_cnts";
		print OUTFILE1 "$line\n";
	}
	close(OUTFILE1);

}
exit;

