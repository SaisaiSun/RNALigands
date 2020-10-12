#!/usr/bin/perl -w
use strict;

#This program is to generate .bpseq file from ss file.

#Saisai Sun
#Date:6-5-2020

###################################


my $query_id=$ARGV[0];
my $output_dir="/var/www/rnaligands";

	open(OUT,">$output_dir/$query_id.bpseq");
	my @ss=`cat $output_dir/$query_id\_dot.txt`;chomp(@ss);
	$ss[2]=~/(\S+)\s+/;my $SS=$1;
	#my $SS=$ss[2];
	my @seq=split(//,$ss[1]);
	my @arr=split(//,$SS);
	my @pair1=();
        my @pair2=();
        my @pair3=();
	my @dot=();
	my @matrix=();
        for(my $i=0;$i<@arr;$i++)
        {
                for(my $j=0;$j<@arr;$j++)
                {
                        $matrix[$i][$j]=0;
                }
	}
	for(my $i=0;$i<@arr;$i++)
	{
		if($arr[$i] eq '(')
		{
			push(@pair1,$i+1);
		}
                elsif($arr[$i] eq '[')
                {
                        push(@pair2,$i+1);
                }
                elsif($arr[$i] eq '{')
                {
                        push(@pair3,$i+1);
                }
		elsif($arr[$i] eq '.')
		{
			push(@dot,$i);
		}
		elsif($arr[$i] eq ')')
                {
			my $positon=pop(@pair1);
			#print OUT $positon,"\t",$i+1,"\n";
			$matrix[$positon-1][$i]=1;
		}
                elsif($arr[$i] eq ']')
                {
                        my $positon=pop(@pair2);
                        #print OUT $positon,"\t",$i+1,"\n";
			$matrix[$positon-1][$i]=1;
                }
                elsif($arr[$i] eq '}')
                {
                        my $positon=pop(@pair3);
                        #print OUT $positon,"\t",$i+1,"\n";
			$matrix[$positon-1][$i]=1;
                }
	}
	for(my $i=0;$i<@arr;$i++)
        {
		for(my $j=0;$j<$i;$j++)
                {
			$matrix[$i][$j]=$matrix[$j][$i];
		}
	}
	OUTSIDE: for(my $i=0;$i<@arr;$i++)
        {
                for(my $j=0;$j<@arr;$j++)
                {
			if($matrix[$i][$j] eq '1')
			{
				print OUT $i+1,"\t",$seq[$i],"\t",$j+1,"\n";
				next OUTSIDE;
			}
                        #print OUT $i+1,"\t",$j+1,"\t",$matrix[$i][$j],"\n";
                }
		if($matrix[$i][0] eq '0')
		{
			print OUT $i+1,"\t",$seq[$i],"\t",0,"\n";
		}
        }
	close(OUT);
