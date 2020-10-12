#!/usr/bin/perl -w
use strict;
use File::Spec;

#######################################
#This code is used to trim motifs decoded from a RNA chain.

#Saisai Sun
#Date:10-24-2019

#######################################


my $infile=$ARGV[0];
#my $outfile=$ARGV[1];

my ($file, $output_dir, $query_id) = File::Spec->splitpath($infile);

my @motif=`cat $output_dir/$query_id\_motif2.out`;chomp(@motif);
open(OUT,">$output_dir/$query_id\_motif2_final.out");

my @multi_branch=();
my @multi_branch2=();
my @internal=();
my @hairpin=();
my @bugle=();

my $multi_num=0;
my $multi_num2=0;
my $hairpin_num=0;
my $internal_num=0;
my $bugle_num=0;
for(my $i=0;$i<@motif;$i++)
{
	#internal or hairpin or multi-branch
	if($motif[$i] =~ /\((\d+\,\d+)\)(\s+.+)\((\d+\,\d+)\)(\s+.+)\((\d+\,\d+)\)\s+$/)#maybe \s+.+ can be replaced by \s+.*
	{
		my $pair1=$1;
		my $dot1=$2;
		my $pair2=$3;
		my $dot2=$4;
		my $pair3=$5;
		#print "$pair1##$dot1##$pair2##$dot2##$pair3\n";
		#my $pair2_re=reverse($pair2);
		#if($pair1 =~ /\)/)
		if($motif[$i] =~ s/!//g)
		{
			$motif[$i] =~ /\((\d+\,\d+)\)(\s+.+)\((\d+\,\d+)\)(\s+.+)\((\d+\,\d+)\)/;
			$pair1=$1;$dot1=$2;$pair2=$3;$dot2=$4;
			push(@multi_branch,"($pair1)$dot1($pair2)$dot2");
			$multi_num++;
			#print "Multi-branch loop: ($pair1)$dot1($pair2)$dot2\n";
		} 
		#elsif($dot1 eq $dot2 && $pair1 eq $pair2_re && length($dot1)>2)
		elsif($dot1 =~ s/#//g)
		{
			#push(@hairpin,"($pair1)$dot1");
			$hairpin_num++;
			print OUT "Hairpin loop: ($pair1)$dot1\n";
		}
		else{
			$internal_num++;
			#next if($internal_num % 2 eq 0);
			#push(@internal,"($pair1)$dot1($pair2)$dot2");
			print OUT "Internal loop: ($pair1)$dot1($pair2)$dot2\n" if($hairpin_num%2==0);
			#print OUT "Internal loop: ($pair1)$dot1($pair2)$dot2\n";
		}
		#last;
	}
	#bugle
	elsif($motif[$i] =~ /\((\d+\,\d+)\)(\s+.+)\((\d+\,\d+)\)\s+\((\d+\,\d+)\)\s+$/)
	{
		my $pair1=$1;
		my $dot1=$2;
		my $pair2=$3;
		my $pair3=$4;
		$bugle_num++;
		#print $pair1,"\t",$dot1,"\t",$pair2,"\t",$pair3,"\n";
		#push(@bugle,"($pair1)$dot1($pair2)-");
		print OUT "Bulge loop: ($pair1)$dot1($pair2) \n";
	}
	elsif($motif[$i] =~ /\((\d+\,\d+)\)\s+\((\d+\,\d+)\)(\s+.+)\((\d+\,\d+)\)\s+$/)
	{
		my $pair1=$1;
                my $pair2=$2;
                my $dot1=$3;
                my $pair3=$4;
		$bugle_num++;
		#print $pair1,"\t",$dot1,"\t",$pair2,"\t",$pair3,"\n";
		#push(@bugle,"($pair1)-($pair2)$dot1");
		print OUT "Bulge loop: ($pair1) ($pair2)$dot1\n";
	}
	#Multi-branch2
	elsif($motif[$i] =~ /\((\d+\,\d+)\)(\s+.*)\((\d+\,\d+)\)(\s+.*)/)
        {
		my $pair1=$1;
		my $dot1=$2;
		my $pair2=$3;
		my $dot2=$4;
		$dot2 =~ s/[\d+\s+]+$//;
		$multi_num2++;
		push(@multi_branch2,"($pair1)$dot1($pair2)$dot2 ");
		#print OUT "Multi-branch loop2: ($pair1)$dot1($pair2)\n";
	}
=pod
	elsif($motif[$i] =~ /^\((\d+\,\d+)\)([\s+\d+\s+]+)$/)
	{
		my $pair1=$1;
		my $dot1=$2;
		print OUT "Terminal ends:$dot1\n";
	}
	elsif($motif[$i] =~ /^\d+.+/)
	{
		print OUT "Terminal ends: $motif[$i]\n";
	}
=cut
}
#for(my $i=0;$i<scalar(@internal)/2;$i++)
#{
#	print OUT "Internal loop: $internal[$i]\n";
#}
$internal_num=int(($internal_num+1)/2);
for(my $i=0;$i<scalar(@multi_branch);$i++)
{
	my $string=$multi_branch[$i];
	if(my $count= ($string =~ s/\)\s+\d+/#/g))
	{
		#print $count,"\n";
		$multi_num=$multi_num-($count-1);
		splice(@multi_branch,$i+1,$count-1);
	}
	print OUT "Multi-branch loop: $multi_branch[$i]\n";
}
for(my $i=0;$i<scalar(@multi_branch2);$i++)
{
        my $string=$multi_branch2[$i];
        if(my $count= ($string =~ s/\)\s+\d+/#/g))
        {
		if($count>2)
                {$multi_num2=$multi_num2-($count-2);}
		else{
		$multi_num2=1;}
                splice(@multi_branch2,$i+1,$count-2);
        }
        print OUT "Multi-branch loop: $multi_branch2[$i]\n";
}
#print OUT "Hairpin loop: $hairpin_num, Internal loop: $internal_num, Bulge number: $bugle_num, Multi-branch loop: $multi_num, Multi-branch loop2: $multi_num2\n"; 
close(OUT);
