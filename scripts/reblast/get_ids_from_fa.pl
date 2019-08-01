#!c:/perl/bin/perl.exe


# script to get full IDs from fa file based on a list of partial IDs


use warnings;
use strict;
use diagnostics;
use File::Basename;
use Getopt::Long;


my $script_name="get_ids_from_fa.pl";



if ($ARGV[0]eq qw '-h'){
	print "please provide arguments for: \n perl $script_name [arguments]\n";
	print "arguments:\n";
	print "--list: /path/to/infile.txt\n";
	print "--fa: /path/to/infile.fa\n";
	print "--outfile: /path/to/outfile\n\n";
}

else{
	my $parameters=join(' ', @ARGV);
	print "parameters: $parameters\n";
	
	#commandline parsing for parameters
	GetOptions(
		'list=s'		=>	\(my $path2infile_lst),
		'fa=s'		=>	\(my $path2infile_fa),
		'outfile=s'	=>	\(my $path2outfile)
	) or die "Error in command line arguments";


	my($infilename, $dirs, $suffix) = fileparse($path2infile_lst);
	print "processing file $infilename\n";


	### code

	open (INFILE_LIST, "<$path2infile_lst") or die "Cannot open input file $path2infile_lst $!";
	open (INFILE_FA, "<$path2infile_fa") or die "Cannot open input file $path2infile_fa $!";
	open (OUTFILE, ">$path2outfile") or die "Cannot open output file $path2outfile $!"; 

	my %ids;



	while (<INFILE_LIST>){
		chomp $_;
		$ids{$_}=1;
	}
	close(INFILE_LIST);

	while (<INFILE_FA>) {
		if ($_=~m/^>/){
			$_=~ s/>//;
			my @header=split/\s/,$_;
			my $name=$header[0];
			if (exists $ids{$name} ){
				print OUTFILE "$_";
			}

		}
	}

	close (INFILE_FA);
	close (OUTFILE);

}
exit;


