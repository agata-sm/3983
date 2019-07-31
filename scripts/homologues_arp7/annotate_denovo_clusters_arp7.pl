#!c:/perl/bin/perl.exe


# script to annotate results of reciprocal blast and trx clustering with consensus gene groups from ARP7
# proj 3983
# 10 ii 2018

#warning still not fixed
#Useless use of private variable in void context at annotate_denovo_clusters_arp7.pl line 94 (#1)




use warnings;
use strict;
use diagnostics;
use File::Basename;
use Getopt::Long;
#use List::Util qw(min max);
use List::MoreUtils qw(any uniq);
#use File::Find;


my $script_name="annotate_denovo_clusters_arp7.pl";



if ($ARGV[0]eq qw '-h'){
	print "please provide arguments for: \n perl $script_name [arguments]\n";
	print "arguments:\n";
	print "--clusters: /path/to/file/with/clusters/*.corset-clusters.txt\n";
	print "--reblast: /path/to/dir/with/reblast/*.genes_best_matches_confirmed_reblast_best\n";
	print "--outfile: /path/to/outfile\n\n";
}

else{
	my $parameters=join(' ', @ARGV);
	print "parameters: $parameters\n";
	
	#commandline parsing for parameters
	GetOptions(
		'clusters=s'		=>	\(my $file_clusters),
		'reblast=s'		=>	\(my $file_reblast),
		'outfile=s'	=>	\(my $path2outfile)
	) or die "Error in command line arguments";

	#this one is hardcoded: arp7s10f_genes.ugp.txt
	my $file_arp7_cons="arp7s10f_genes.ugp.txt";

	open (OUTFILE, ">$path2outfile") or die "Cannot open output file $path2outfile $!"; 
	my $header="trinity_gene_id\tcorset_cluster_id\treblast_confirmed_arp7\tarp7_consensus_group\tarp7_cg_description";
	print OUTFILE "$header\n";




	################################################################################
	#gene names info
	
	my %consensus_ids_arp7;
	
	my @chunks;
	open (INFILE_CLUST, "<$file_arp7_cons") or die "Cannot open output file $file_arp7_cons $!"; 

	undef $/ ;
	@chunks = split (/GeneSummary id="euGenes:/, <INFILE_CLUST>);
	close INFILE_CLUST;

	#remove the first empty element
	my $foo=shift @chunks;
	#print " this is foo $foo\n";

	my %homologues;
	my %definitions;
	foreach my $record (@chunks){
		#print "$record\n\n";
		my @lines=split/\n/, $record;
		#print "lines for chunk \n";
		#print "@lines\n";



		#my $gene_id=$lines[2];
		#$gene_id=s/"://,$gene_id;
		#print "$gene_id\n";

		#print "\n\n";
		#print "\n*****\n NOWY CHUNK\n";
		
		#if ($lines[0] =~m/ARP7/){
			

		
		my $gene_id=$lines[0];
		$gene_id=~s/"://,$gene_id;
		#print "$gene_id\t";

		my $gene_def=$lines[12];
		$gene_def=~s/\s+description:\s+//g,$gene_def;
		$gene_def=~s/;//,$gene_def;
		#print "$gene_def\n";

		#my @similar_genes;
		foreach my $element (@lines){
			chomp $element;
			if ($element=~m/similarity/){
				if ( ($element=~m/acc: (\w+);/) | ($element=~m/acc: (\w+\.\d+);/) ){
					my $homolog=$1;
					#print "$homolog\n";
					push @{$homologues{$gene_id}}, $homolog;
				}
			}
		}
		$definitions{$gene_id}=$gene_def;
	}

	#just for sanity checks later
	foreach my $key (keys %definitions){
		#print "$key\t$definitions{$key}\t";
		my @array=@{$homologues{$key}};
		#print "@array\n";
	}


	################################################################################
	$/="\n" ;


	################################################################################
	#reblast

	my %reblast;
	open (INFILE_REBLAST, "<$file_reblast") or die "Cannot open output file $file_reblast $!"; 
	while(<INFILE_REBLAST>){
		chomp $_;
		if ($_ =~m/(TRINITY\w+)\t(\S+)/){
			my ($trinity_gene_id,$reblast_hit)=($1,$2);
			$reblast_hit=~m/(\w+):(\S+)/;
			my ($org,$accession)=($1,$2);
			#print "$trinity_gene_id\t$reblast_hit\n";
			$reblast{$trinity_gene_id}=$accession;
			#print "$trinity_gene_id\n";
		}
	}
	close(INFILE_REBLAST);

	#foreach my $key  (keys %reblast){
	 #	print "$key\t$reblast{$key}\n";
	 #}

	################################################################################
	#trinity_id - clusters



	my %clusters;
	open (INFILE_CLUST, "<$file_clusters") or die "Cannot open output file $file_clusters $!"; 
	while(<INFILE_CLUST>){
		chomp $_;
		if ($_ =~m/(TRINITY\w+)\t(\S+)/){
			my ($trinity_trx_id,$cluster_id)=($1,$2);
			$trinity_trx_id=~m/(\w+)_i\d+$/;
			my $trinity_gene_id=$1;
			push @{$clusters{$trinity_gene_id}},$cluster_id;
			#print "$trinity_gene_id\n";
		}
	}
	close(INFILE_CLUST);

	################################################################################

	#superclusters that are assigned more than 1 to 1 trinity gene id
	# i.e. TRINITY_DN659_c0_g1	Cluster-1291.23814 Cluster-1291.29106
	# 3234 genes, 5 superclusters - in superclusters.multiple
	# deal with this later
	my $foo2_file="tmp.trinity_geneids.multiple_clusters";
	open (OUTFILE_FOO2, ">$foo2_file") or die "Cannot open output file $foo2_file $!"; 

	my %supercluster;
	foreach my $trinity_gene_id (keys %clusters){
		my @clusters_=@{$clusters{$trinity_gene_id}};
		my @clusters=uniq(@clusters_);
		my $number=scalar(@clusters);
		if ($number>1){
			print OUTFILE_FOO2 "$trinity_gene_id\t@clusters\n";
			foreach my $cluster_id (@clusters){ #Cluster-1291.18482 or Cluster-0.0
				#print "$cluster_id\n";

				if ($cluster_id=~m/(Cluster\S\d+)\..+/){
					#print "$cluster_id\n";
					my $supercluster_id=$1;
					$supercluster{$supercluster_id}=$cluster_id;
				}
				elsif($cluster_id=~m/Cluster\S(\d+)$/){
					my $supercluster_id=$1;
					$supercluster{$supercluster_id}=$cluster_id;
				}
				
			}
		}
	}
	close(OUTFILE_FOO2);

	my $foo_file="tmp.superclusters.multiple";
	open (OUTFILE_FOO, ">$foo_file") or die "Cannot open output file $foo_file $!"; 

	foreach my $supercluster (keys %supercluster){
		print OUTFILE_FOO "$supercluster\n";
	}
	close (OUTFILE_FOO);


	################################################################################
	################################################################################
	# connect all info together

	foreach my $trinity_gene_id (sort keys %clusters){
		my @clusters_=@{$clusters{$trinity_gene_id}};
		my @clusters=uniq(@clusters_);
		my $number=scalar(@clusters);
		if ($number>1){ #trinity_gene_id with ambiguous cluster assignment (several clusters per gene id)
			#print OUTFILE "$trinity_gene_id\t@clusters\n";
		}if($number==1){
			my $cluster_id=$clusters[0];
			
			#print "$trinity_gene_id\t";

			if (exists ( $reblast{$trinity_gene_id}) ) {
				my $reblast_hit=$reblast{$trinity_gene_id};
				#print "$reblast_hit\n";
				#print "$trinity_gene_id\t$reblast_hit\n";

				#$definitions{$gene_id}=$gene_def;
				#push @{$homologues{$gene_id}}, $homolog;

				#my @reblast_matches=@{$reblast_hits{$trinity_trx_id}};

				foreach my $arp7_gene_id (keys %homologues){
					#rint "$arp7_gene_id\n";
					my @homologues=@{$homologues{$arp7_gene_id}};
									#my @reblast_matches=@{$reblast_hits{$trinity_trx_id}};

					#print "@homologues\n";
					#if ($reblast_hit eq $arp7_gene_id){
					if (any { $_ eq $reblast_hit} (@homologues) )	{
						#print "$trinity_gene_id\t$reblast_hit\t$arp7_gene_id\n";

						my $gene_desc=$definitions{$arp7_gene_id};
						my $line="$trinity_gene_id\t$cluster_id\t$reblast_hit\t$arp7_gene_id\t$gene_desc";
						print OUTFILE "$line\n";
					}

				}

				#if (any { $_ eq $gene_id} (@reblast_best_matches) )	{

				#my $line="$trinity_gene_id\t$cluster_id\tcorset_super_cluster_id\t$reblast_hit\tarp7_gene_id\tgene_desc";
				#		print OUTFILE "$line\n";
			}
			
			



		}
	}



#if (any { $_ eq $gene_id} (@reblast_best_matches) )	{


}
exit;