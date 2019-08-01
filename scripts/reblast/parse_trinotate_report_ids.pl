#!c:/perl/bin/perl.exe


# script to parse the report from trinotate to generate a list of IDs for reciprocal blast


use warnings;
use strict;
use diagnostics;
use File::Basename;
use Getopt::Long;


my $script_name="parse_trinotate_report_ids.pl";



if ($ARGV[0]eq qw '-h'){
	print "please provide arguments for: \n perl $script_name [arguments]\n";
	print "arguments:\n";
	print "--infile: /path/to/infile\n";
	print "--out_prefix: output prefix\n\n";
}

else{
	my $parameters=join(' ', @ARGV);
	print "parameters: $parameters\n";
	
	#commandline parsing for parameters
	GetOptions(
		'infile=s'		=>	\(my $path2infile),
		'out_prefix=s'	=>	\(my $out_prefix)
	) or die "Error in command line arguments";


	my($infilename, $dirs, $suffix) = fileparse($path2infile);
	print "processing file $infilename\n";


	my $outfile_sprot="$out_prefix\.$infilename\.ids_sprot.txt";
	my $outfile_trembl="$out_prefix\.$infilename\.ids_trembl.txt";
	my $outfile_arp7="$out_prefix\.$infilename\.ids_arp7.txt";


	### code

	open (INFILE, "<$path2infile") or die "Cannot open input file $path2infile $!";
	open (OUT_S, ">$outfile_sprot") or die "Cannot open output file $outfile_sprot $!"; #sprot IDs (consensus from blastp and blastx)
	open (OUT_T, ">$outfile_trembl") or die "Cannot open output file $outfile_trembl$!"; #trembl IDs for genes with no sprot ID
	open (OUT_A, ">$outfile_arp7") or die "Cannot open output file $outfile_arp7 $!"; #arp7 IDs (all)


	my %sprot_p;
	my %sprot_x;
	my %trembl;
	my %arp7;


	while (<INFILE>){
	chomp $_;

	unless($_ =~m/^\#/ ){ 
		my @entry=split/\t/,$_;
		
		my $gene_id=$entry[0];
		my $transcript_id=$entry[1];
		my $prot_id=$entry[4];

		my $sprot_BLASTX=$entry[2];
		my $sprot_BLASTP=$entry[6];

		my $ARP7_BLASTP=$entry[9];
		my $TREMBL_BLASTP=$entry[10];

		my $pfam=$entry[11];
		my $signalP=$entry[12];
		my $tmm=$entry[13];

		my @hit_sprot_blastx=split/\^/,$sprot_BLASTX;
		my $sprot_blastxID=$hit_sprot_blastx[0];

		my @hit_sprot_blastp=split/\^/,$sprot_BLASTP;
		my $sprot_blastpID=$hit_sprot_blastp[0];

		my @hit_arp7_blastp=split/\^/,$ARP7_BLASTP;
		my $arp7_blastpID=$hit_arp7_blastp[0];

		my @hit_trembl_blastp=split/\^/,$TREMBL_BLASTP;
		my $trembl_blastpID=$hit_trembl_blastp[0];


		#to retain only unique ids
		if ($sprot_blastpID ne qw /./) {
			$sprot_p{$sprot_blastpID}=$gene_id;
		}
		if ($sprot_blastxID ne qw /./) {
			$sprot_x{$sprot_blastxID}=$gene_id;
		}
		if ($arp7_blastpID ne qw /./) {
			$arp7{$arp7_blastpID}=$gene_id;
		}
		if ($trembl_blastpID ne qw /./) {
			$trembl{$trembl_blastpID}=$gene_id;
		}

	}
	}


	#printing to files only unique ids
	for my $id (keys %sprot_p){
	print OUT_S "$id\n";
	}
	for my $id (keys %sprot_x){
	unless(exists($sprot_p{$id})){
		print OUT_S "$id\n";
	}
	}
	for my $id (keys %trembl){
	unless( (exists($sprot_p{$id})) | (exists($sprot_x{$id})) ){
		print OUT_T "$id\n";
	}
	}
	for my $id (keys %arp7){
	print OUT_A "$id\n";
	}


	close (INFILE);
	close (OUT_A);
	close (OUT_S);
	close (OUT_T);
}
exit;


