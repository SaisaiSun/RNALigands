#!/usr/bin/perl -w
use strict;
use File::Spec;

############################
#This code is to extract RNA structural motif from dot and bracket file.

#Saisai Sun
#Date:10-24-2019

###########################

#my $dir=$ENV{PWD};
my $infile=$ARGV[0];
#my $outfile=$ARGV[1];

my ($file, $output_dir, $query_id) = File::Spec->splitpath($infile);

my @file=`cat $output_dir/$query_id\_dot.txt`;chomp(@file);
$file[0]=~s/\r//g;
$file[1]=~s/\r//g;
$file[2]=~s/\r//g;
my $length=length($file[1]);
my $SS=substr($file[2],0,$length);
my @seq=split(//,$file[1]);
my @ss=split(//,$SS);
my @pair1=();
my @pair2=();
my @pair3=();
my @dot=();
my %seq_ss=();

for(my $i=0;$i<@ss;$i++)
{
	if($ss[$i] eq '(')
	{
		push(@pair1,$i);
	}
	elsif($ss[$i] eq '<')
	{
		push(@pair2,$i);
	}
	elsif($ss[$i] eq '{')
        {
                push(@pair3,$i);
        }
	elsif($ss[$i] eq '.')
	{
		push(@dot,$i);
		$seq_ss{$i}=$seq[$i];
	}
	elsif($ss[$i] eq ')')
	{
		my $position=pop(@pair1);
		$seq_ss{$position}="($seq[$position]$seq[$i])-$i";
		$seq_ss{$i}="($seq[$i]$seq[$position])-$position";
	}
        elsif($ss[$i] eq '>')
        {
                my $position=pop(@pair2);
                #$seq_ss[$position][$i]="($seq[$position]$seq[$i])";
        }
        elsif($ss[$i] eq '}')
        {
                my $position=pop(@pair3);
                #$seq_ss[$position][$i]="($seq[$position]$seq[$i])";
        }
}
=pod
for(my $i=0;$i<@ss;$i++)
{
	print $i,"\t",$seq_ss{$i},"\n";
}
=cut

#print scalar(@dot),"\n";
if(scalar(@dot)==scalar(@seq))
{
	print "No paired nucleotide!\n" ;
	exit;
}

open(OUT,">$output_dir/$query_id\_motif2.out");
my $pair_num1=@seq-1;
my $pair_num2=0;
for(my $i=0;$i<@dot;$i++)
{
	#print OUT $dot[$i]+1," ";
#=pod
	if($dot[$i] == 0)
	{
		#print OUT $dot[$i]+1," ";
		until($dot[$i+1]!=$dot[$i]+1)
		{
			print OUT $dot[$i]+1," ";
			$i++;
		#	break;
		}
		print OUT $dot[$i]+1," ";
		print OUT "\n";
	}
	else{
	if($seq_ss{$dot[$i]-1} =~ /(.+)\-(\d+)/ )
	{
		my $pair=$1;
		$pair_num1=$2;
		#print OUT $pair;
		#print $seq_ss{$dot[$i]-1};
		print OUT "(",$dot[$i]-1+1,",",$pair_num1+1,") ";
	}
	#print $seq_ss{$dot[$i]};
	print OUT $dot[$i]+1," ";
	if($i<@dot-1)
	{
		next if($dot[$i+1]==$dot[$i]+1);
		if($seq_ss{$dot[$i]+1} =~ /(.+)\-(\d+)/ )
		{
			print OUT "# " if($dot[$i]+1 eq $pair_num1);
			my $pair=$1;
			$pair_num2=$2;
			#print OUT $pair;
			#print $seq_ss{$dot[$i]+1};
			print OUT "(",$dot[$i]+1+1,",",$pair_num2+1,") ";
			my $first_num=$pair_num1;
			my $last_num=$pair_num2;
			my $j;
			for($j=$pair_num2+1;$j<@seq;$j++)
			{
				print OUT $j+1," " if($j ~~ @dot);
				#print $seq_ss{$j} if($j ~~ @dot);
				if($seq_ss{$j} =~ /(.+)\-(\d+)/ )
				{
					$first_num=$j;
					$last_num=$2;
					my $pair=$1;
					#print OUT $pair;
					#print $seq_ss{$first_num};
					print OUT "(",$first_num+1,",",$last_num+1,") ";
					#$j=$last_num;
					last;
				}
			}
			for(my $e=0;;$e++)
			{
				last if($j==scalar(@seq));
				last if($first_num == $pair_num1);
				for($j=$last_num+1;$j<@seq;$j++)
				{
					print OUT $j+1," " if($j ~~ @dot);
					#print $seq_ss{$j} if($j ~~ @dot);
					if($seq_ss{$j} =~ /(.+)\-(\d+)/ )
					{
						print OUT "! ";
						$first_num=$j;
						$last_num=$2;
						my $pair=$1;
						#print OUT $pair;
						#print $seq_ss{$first_num};
						print OUT "(",$first_num+1,",",$last_num+1,") ";
						last;
					}
				}
				#last if($j==scalar(@seq));
                                #last if($first_num == $pair_num1);
			}
			print OUT "\n";
		}
	}
	else
	{
		if($dot[$i]+1 < @seq && $seq_ss{$dot[$i]+1} =~ /(.+)\-(\d+)/ )
		{
			my $pair=$1;
			my $last_num=$2;
			#print OUT $pair;
			#print $seq_ss{$dot[$i]+1};
			print OUT "# " if($dot[$i]+1 eq $pair_num1);
			print OUT "(",$dot[$i]+1+1,",",$2+1,") ";
                        for(my $e=$last_num+1;$e<@seq;$e++)
                        {
                                if($seq_ss{$e}=~ /(.+)\-(\d+)/ )
                                {
					print OUT "(",$e+1,",",$2+1,") ";
                                        last;
                                }
				print OUT $e+1," ";
                        }
			
		}
	}
	}
#=cut
	#last;
	#print $seq_ss{$dot[$i]+1},"\n" if($seq_ss{$dot[$i]+1}=~/[A-Z][A-Z]/);
}
print OUT "\n";
close(OUT);
