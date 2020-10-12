#!/usr/bin/perl -w
use strict;
use Getopt::Std; 
use File::Spec;
my %opts;
my $scriptdir=$ENV{PWD};
getopts('f:s:h',\%opts);

if ($opts{'h'}) {
print "STD  output
      -f   fasta file.
      -s   secondary structure file (optional)(format like example/1ddy_A_dot.txt).
      -h   help.\n";
}

exit() unless ($opts{'f'} || $opts{'s'});

if($opts{'f'}) {
	my ($file1, $out_dir, $fasta_name) = File::Spec->splitpath($opts{'f'});
	$fasta_name=~s/\.fasta//g;
	my $name=$fasta_name;

	`/var/www/rnaligands/ViennaRNA/bin/RNAfold --noPS -i $out_dir/$name.fasta >$out_dir/$name\_dot.txt`;#This should be replaced with your RNAfold directory
	`$scriptdir/dot2motif.pl $out_dir/$name`;
	`$scriptdir/dot2motif2.pl $out_dir/$name`;
	`$scriptdir/process_motif.pl $out_dir/$name`;
	`$scriptdir/process_motif2.pl $out_dir/$name`;
	`$scriptdir/search_motif.pl $out_dir/$name`;
	`$scriptdir/search_R-BIND.pl $out_dir/$name`;
	#`$scriptdir/search_miRBase.pl $out_dir/$name`;
}

if($opts{'s'}) {
	my ($file, $out_dir, $ss_name) = File::Spec->splitpath($opts{'s'});
	$ss_name=~s/\.(.+)//g;
	my $name=$ss_name;my $format=$1;

	`cp "$out_dir/$name.$format" "$out_dir/$name\_dot.txt"`;
	`$scriptdir/dot2motif.pl $out_dir/$name`;
        `$scriptdir/dot2motif2.pl $out_dir/$name`;
        `$scriptdir/process_motif.pl $out_dir/$name`;
        `$scriptdir/process_motif2.pl $out_dir/$name`;
        `$scriptdir/search_motif.pl $out_dir/$name`;
        `$scriptdir/search_R-BIND.pl $out_dir/$name`;
        `$scriptdir/search_miRBase.pl $out_dir/$name`;
}

