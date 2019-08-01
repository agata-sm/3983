#!c:/perl/bin/perl.exe


# script to parse the results of reciprocal blast
# proj 3983
# 16 xii 2017


use warnings;
use strict;
use diagnostics;
use File::Basename;
use Getopt::Long;
use List::Util qw(min max);
use List::MoreUtils qw(any uniq);
use File::Find;


my $script_name="parse_reblast.pl";



if ($ARGV[0]eq qw '-h'){
	print "please provide arguments for: \n perl $script_name [arguments]\n";
	print "arguments:\n";
	print "--dir_blast: /path/to/dir/with/blast\n";
	print "--dir_reblast: /path/to/dir/with/reblast\n";
	print "--outfile: /path/to/outfile\n\n";
}

else{
	my $parameters=join(' ', @ARGV);
	print "parameters: $parameters\n";
	
	#commandline parsing for parameters
	GetOptions(
		'dir_blast=s'		=>	\(my $dir_blast),
		'dir_reblast=s'		=>	\(my $dir_reblast),
		'outfile=s'	=>	\(my $path2outfile)
	) or die "Error in command line arguments";


	#outfiles
	my $path2outfile1="$path2outfile\.all_transcripts_matches";
	my $path2outfile2="$path2outfile\.genes_best_matches_confirmed_reblast_best";
	my $path2outfile3="$path2outfile\.genes_best_matches_confirmed_reblast_any";


	#############################################
	#get all blastp files
	my @blastp_files;
	sub find_files_blastp {
    	my $F = $File::Find::name;

    	if ($F =~ /.+\.arp7.blastp.outfmt6/ ) {
    		push @blastp_files, $F;
    	}
	}
	find({ wanted => \&find_files_blastp, no_chdir=>1}, $dir_blast);


	#############################################
	#get all reblast files
	my @reblast_files;
	sub find_files_reblast {
    	my $F = $File::Find::name;

    	if ($F =~ /.+\.arp7.15xii2017.blastp.outfmt6/ ) {
    		push @reblast_files, $F;
    	}
	}
	find({ wanted => \&find_files_reblast, no_chdir=>1}, $dir_reblast);



	#############################################
	#outfiles
	open (OUTFILE1, ">$path2outfile1") or die "Cannot open output file $path2outfile1 $!"; 

	my $header="gene_id_trinity\tisoform_id_trinity\tbest_hit_blastp\tmatch_length_blastp\tperc_identity_blastp\te_val_blastp\tbit_score_blastp\tfound_by_reblast_query\tmatch_length_reblast\tperc_identity_reblast\te_val_reblast\tbit_score_reblast";
	print OUTFILE1 "$header\n";

	#############################################
	#process results from reblast
	my %reblast_hits;	#trinity isoform id as key
	my %reblast_best_hits; #trinity gene id as key

	foreach my $path (@reblast_files){
		my($infile_name, $dirs, $suffix) = fileparse($path);
		print "processing file $infile_name\n";
		open(FILE,"<$path") or die "Cannot open input file $path $!";
		my $last_seen="0::0::0";
		my $tophit="yes";

		while(<FILE>){
			chomp $_;
			my @line=split/\t/,$_;
			my $query_id=$line[0];
			my $trinotate_id=$line[1];
			my @trinotateid=split/::/,$trinotate_id;
			my $isoform_id=$trinotateid[1];
			my $gene_id=$trinotateid[0];
			my $identity=$line[2];
			my $match_length=$line[3];
			my $bitscore=$line[11];
			my $evalue=$line[10];

			my $isoform_id_info="$query_id::$identity::$match_length::$bitscore::$evalue";

			push @{$reblast_hits{$isoform_id}},$isoform_id_info;

			#choose the best hit for each query based on highest bitscore
			#the hits are sorted based on e value / bitscore, take the first entry (line) per query
			#if there are more than one best hits, choose all
			my ($prev_query_id,$prev_bitscore,$prev_gene_id)=split/::/,$last_seen;

			#a way first derived from taking the first line of matches reported for the same query
			unless ($query_id eq $prev_query_id){
				#unless (exists $reblast_best_hits{$query_id}){#redundant
					#$reblast_best_hits{$query_id}=$gene_id;
					push @{$reblast_best_hits{$query_id}},$gene_id; #because there are more than one "best hits" of equal e value for some genes
					$tophit="yes";
				#}
			}else{ #if query id eq prev query id
				if ( ($bitscore eq $prev_bitscore) && ($gene_id ne $prev_gene_id) ) {
					if ($tophit eq qw /yes/){
						push @{$reblast_best_hits{$query_id}},$gene_id;
					}	
				}
				$tophit="no";
			}

			#or an alternative for choosing best hits in reblast (esp when multiple best hits of equal bitscore exist)



			$last_seen="$query_id::$bitscore::$gene_id";

		}
		close(FILE);
	}

	#############################################
	#collect results from blast
	my %blastp_hits;

	foreach my $path (@blastp_files){
		my($infile_name, $dirs, $suffix) = fileparse($path);
		print "processing file $infile_name\n";
		open(FILE,"<$path") or die "Cannot open input file $path $!";
		while(<FILE>){
			chomp $_;
			my @line=split/\t/,$_;
			my $query_id_trinotate=$line[0]; #TRINITY_DN27975_c0_g1::TRINITY_DN27975_c0_g1_i1::g.13::m.13
			my $hit_id_blastp=$line[1];
			my @queryid_trinotate=split/::/,$query_id_trinotate;
			my $isoform_id=$queryid_trinotate[1];
			my $gene_id=$queryid_trinotate[0];			
			my $identity=$line[2];
			my $match_length=$line[3];
			my $bitscore=$line[11];
			my $evalue=$line[10];
			my $hit_id_blastp_bitscore="$hit_id_blastp::$bitscore";


			$blastp_hits{$isoform_id}=$hit_id_blastp_bitscore;

			#print info to file with all info
			if (exists $reblast_hits{$isoform_id}){ #@{$reblast_hits{$isoform_id}} 
				my @hits_reblast=@{$reblast_hits{$isoform_id}};

				foreach my $reblast_hit (@hits_reblast){
					my ($query_id_reblast,$identity_rb,$match_length_rb,$bitscore_rb,$evalue_rb)=split/::/,$reblast_hit;
					my $line="$gene_id\t$isoform_id\t$hit_id_blastp\t$match_length\t$identity\t$evalue\t$bitscore\t$query_id_reblast\t$match_length_rb\t$identity_rb\t$evalue_rb\t$bitscore_rb";
					print OUTFILE1 "$line\n";

				}
			}else{
				my $line="$gene_id\t$isoform_id\t$hit_id_blastp\t$match_length\t$identity\t$evalue\t$bitscore\t.\t.\t.\t.\t.";
				print OUTFILE1 "$line\n";
			}

		}
		close(FILE);
	}
	close(OUTFILE1);


	#get best matches from blastp on gene level (in case if one of the isoforms gets best match to another sequence; probably a homolog also, but may lead to false negative)
	my %blastp_all_matches;
	my %blastp_best_matches;

	foreach my $isoform_id (keys %blastp_hits){
		my $hit_id_blastp_bts=$blastp_hits{$isoform_id};
		$isoform_id=~m/(\w+_c\d+_g\d+)_i\d+/;
		my $gene_id=$1; #trinity gene id
		push @{$blastp_all_matches{$gene_id}},$hit_id_blastp_bts;
	}


	foreach my $gene_id (keys %blastp_all_matches){
		my @blastp_matches=@{$blastp_all_matches{$gene_id}};
					
		my @bitscores;
		my @prot_ids;

		foreach my $match (@blastp_matches){
			my ($hit_id_blastp,$bitscore)=split/::/,$match;
			push @bitscores,$bitscore;
			push @prot_ids,$hit_id_blastp;
		}

		my $best_bitscore=max(@bitscores);
		my @idxs_bitscore= List::MoreUtils::indexes {$_ eq $best_bitscore} @bitscores;

		#get the prot ids with highest bitscore
		my @ids_best_match_gene_lvl;
		foreach my $idx (@idxs_bitscore){
			push @ids_best_match_gene_lvl, $prot_ids[$idx];
		}

		my @uniq_hit_ids=uniq(@ids_best_match_gene_lvl); #what if there are two best hits of equal bitscore? and there are, see below!
		#my $number=scalar(@uniq_hit_ids);
		#print "$gene_id\t$number\t@uniq_hit_ids\n";# TRINITY_DN32739_c0_g2	2	pogcha:PogchaEGm005043t1 tribcas:XP_973009
		
		#treat both hits equally, without further filtering
		foreach my $id (@uniq_hit_ids){
			push @{$blastp_best_matches{$gene_id}},$id;

		}
		

		undef(@bitscores);
		undef(@prot_ids);
	}




	#get best matches on gene level which are confirmed by reblast
	#if there are >1 best matches per gene, accept any if found by reblast
	# if more than one match of equal bitscore was found by reblast, list all ( from @{$reblast_best_hits{$query_id}},$gene_id;) query id - db id; gene id - trinity gene id
	my $header2="gene_id_trinity\tbest_hit_blastp_reblast";

	open (OUTFILE2, ">$path2outfile2") or die "Cannot open output file $path2outfile2 $!"; #best in reblast
	print OUTFILE2 "$header2\n";



	my %reblast_hits_confirmed_any;
	foreach my $gene_id (keys %blastp_best_matches){ # $gene_id - trinity gene id
		
		#on a gene level - if different isoforms blastp to two db_ids equally well
		my @blastp_matches=@{$blastp_best_matches{$gene_id}}; #@{$blastp_best_matches{$gene_id}},$id; gene_id - trinity gene id ; id - hit id in blastp = db id

		#print "geneid: $gene_id\tblastp_matches: @blastp_matches\n";
		foreach my $db_id (@blastp_matches){ #db ids
			
			###
			# reblast best matches
			if (exists $reblast_best_hits{$db_id}){ #@{$reblast_best_hits{$query_id}},$gene_id; query id - db id
				#print "blastp: $gene_id\treblast: $_\n"; #this does not print best hits! prints reblast id if it was found in the hash, without checking the value

				#print "$db_id\n";
				my @reblast_best_matches=@{$reblast_best_hits{$db_id}};  #@{$reblast_best_hits{$query_id}},$gene_id; query id - db id


				#print "$gene_id\tdbid $db_id\t@reblast_best_matches\n"; #prints hits confirmed by reblast as best hits in @
				
				if (any { $_ eq $gene_id} (@reblast_best_matches) )	{
					#print "success!\n";
					#print "reblast: $gene_id\t$_\n";
					my $line="$gene_id\t$db_id";
					print OUTFILE2 "$line\n";
				}
			}

			# reblast_matches (all in top 5)
			# @{$reblast_hits{$isoform_id}},$isoform_id_info; $isoform_id_info="$query_id::$identity::$match_length::$bitscore::$evalue";
			foreach my $trinity_trx_id (keys %reblast_hits){
				my @reblast_matches=@{$reblast_hits{$trinity_trx_id}};

				$trinity_trx_id=~m/(TRINITY_\w+_\w+_g\d+)_i\d+/;
				my $trinity_gene_id_reblast=$1;

				if ($trinity_gene_id_reblast eq $gene_id){
					foreach $_ (@reblast_matches){
						my ($query_id,$identity,$match_length,$bitscore,$evalue)=split/::/,$_;
						
						if ($query_id eq $db_id){ #this prints repetitions - lines are duplicated (because of isoforms! - fix this)
							$reblast_hits_confirmed_any{$gene_id}=$query_id;
						}
					}

				}
			}
		}
			
	open (OUTFILE3, ">$path2outfile3") or die "Cannot open output file $path2outfile3 $!"; #any in reblast
	print OUTFILE3 "$header2\n";


	while (my ($denovo_gene_id,$db_id) = each %reblast_hits_confirmed_any){
		#print "$denovo_gene_id\n";
		my $line="$denovo_gene_id\t$db_id";
		print OUTFILE3 "$line\n";

	}

#geneid: TRINITY_DN32549_c0_g1	blastp_matches: pogcha:PogchaEGm016353t1
#TRINITY_DN32549_c0_g1	dbid pogcha:PogchaEGm016353t1	TRINITY_DN32549_c0_g1 TRINITY_DN32549_c2_g1


		# #old
		# OBS! this prints the pair if the db id - trinity gene id was found as a best match in reblast in ANY db id - trinity id pair, not necessarily the one being examined
		# my @blastp_matches=@{$blastp_best_matches{$gene_id}};

		# while (my ($query_id_reblast,$gene_id_reblast) = each %reblast_best_hits){ #unless ($query_id eq $last_seen); $reblast_best_hits{$query_id}=$gene_id;

		# 	if ( any { $_ eq $query_id_reblast} (@blastp_matches) ) {
		# 		#reblast confirms match found by blast
		# 		#print "success!\n";
		# 		my $line="$gene_id\t$query_id_reblast";
		# 		print OUTFILE2 "$line\n";

		# 	}
		# } 
 


	}

	close(OUTFILE2);
	close(OUTFILE3);

}
exit;


