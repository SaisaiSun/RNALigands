#!/usr/bin/perl -w
use strict;
use File::Spec;

#######################################
#This code is used to trim motifs decoded from a RNA chain.

#Saisai Sun
#Date:10-18-2019

#######################################


my $infile=$ARGV[0];
#my $outfile=$ARGV[1];

my ($file, $output_dir, $query_id) = File::Spec->splitpath($infile);

my @motif=`cat $output_dir/$query_id\_motif.out`;chomp(@motif);
open(OUT,">$output_dir/$query_id\_motif_final.out");
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
	if($motif[$i] =~ /\((\S+)\)(\S+)\((\S+)\)(\S+)\((\S+)\)$/)
	{
		my $pair1=$1;
		my $dot1=$2;
		my $pair2=$3;
		my $dot2=$4;
		my $pair3=$5;
		#print OUT $pair1,"\t",$dot1,"\t",$pair2,"\t",$dot2,"\t",$pair3,"\n";
		my $pair2_re=reverse($pair2);
		#if($pair1 =~ /\)/)
		if($motif[$i] =~ s/!//g)
		{
			$motif[$i] =~ /\((\S+)\)(\S+)\((\S+)\)(\S+)\((\S+)\)/;
			$pair1=$1;$dot1=$2;$pair2=$3;$dot2=$4;
			push(@multi_branch,"($pair1)$dot1($pair2)$dot2");
			$multi_num++;
			#print OUT "Multi-branch loop: ($pair1)$dot1($pair2)$dot2\n";
		} 
		#elsif($dot1 eq $dot2 && $pair1 eq $pair2_re && length($dot1)>2)
		elsif($dot1 =~ s/#//g)
		{
			#push(@hairpin,"($pair1)$dot1");
			$hairpin_num++;
			#my @pair=split(//,$pair1);
			print OUT "Hairpin loop: ($pair1) $dot1\n";
		}
		else{
			$internal_num++;
			#next if($internal_num % 2 eq 0);
			#push(@internal,"($pair1)$dot1($pair2)$dot2");
			#my @pair1=split(//,$pair1);
			#my @pair2=split(//,$pair2);
			print OUT "Internal loop: ($pair1) $dot1 ($pair2) $dot2\n" if($hairpin_num%2==0);
		}
		#last;
	}
	#bugle
	elsif($motif[$i] =~ /\((\S+)\)(\S+)\((\S+)\)\((\S+)\)$/)
	{
		my $pair1=$1;
		my $dot1=$2;
		my $pair2=$3;
		my $pair3=$4;
		$bugle_num++;
		#print OUT $pair1,"\t",$dot1,"\t",$pair2,"\t",$pair3,"\n";
		#push(@bugle,"($pair1)$dot1($pair2)-");
		#my @pair1=split(//,$pair1);
		#my @pair2=split(//,$pair2);
		print OUT "Bulge loop: ($pair1) $dot1 ($pair2)\n";
	}
	elsif($motif[$i] =~ /\((\S+)\)\((\S+)\)(\S+)\((\S+)\)$/)
	{
		my $pair1=$1;
                my $pair2=$2;
                my $dot1=$3;
                my $pair3=$4;
		$bugle_num++;
		#push(@bugle,"($pair1)-($pair2)$dot1");
		my @pair1=split(/,/,$pair1);
		my @pair2=split(/,/,$pair2);
		print OUT "Bulge loop: ($pair2[1],$pair2[0]) $dot1 ($pair1[1],$pair1[0])\n";
	}
	#Multi-branch2
	elsif($motif[$i] =~ /\((\S+)\)(.+)\((\S+)\)/)
	{
		my $pair1=$1;
                my $dot1=$2;
                my $pair2=$3;
                $multi_num2++;
		push(@multi_branch2,"($pair1)$dot1($pair2)");
                #print OUT "Multi-branch loop2: ($pair1)$dot1($pair2)\n";
	}
}
#for(my $i=0;$i<scalar(@internal)/2;$i++)
#{
#	print OUT "Internal loop: $internal[$i]\n";
#}
$internal_num=int(($internal_num+1)/2);
for(my $i=0;$i<scalar(@multi_branch);$i++)
{
	my $string=$multi_branch[$i];
	if(my $count= ($string =~ s/\)[a-zA-Z]+/#/g))
	{
		#print OUT $count,"\n";
		$multi_num=$multi_num-($count-1);
		splice(@multi_branch,$i+1,$count-1);
	}
	$multi_branch[$i]=~s/\)/\) /g;
	$multi_branch[$i]=~s/\(/ \(/g;
	print OUT "Multi-branch loop: $multi_branch[$i]\n";
}
for(my $i=0;$i<scalar(@multi_branch2);$i++)
{
        my $string=$multi_branch2[$i];
        if(my $count= ($string =~ s/\)[a-zA-Z]+/#/g))
        {
                $multi_num2=$multi_num2-($count-1);
                splice(@multi_branch2,$i+1,$count-1);
        }
	$multi_branch2[$i]=~s/\)/\) /g;
	$multi_branch2[$i]=~s/\(/ \(/g;
        print OUT "Multi-branch loop: $multi_branch2[$i]\n";
}
#print OUT "Hairpin loop: $hairpin_num, Internal loop: $internal_num, Bulge number: $bugle_num, Multi-branch loop: $multi_num, Multi-branch loop2: $multi_num2\n";
close(OUT);
