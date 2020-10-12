#!/usr/bin/perl -w
use strict;

######################################
#Usage: ./Motif-align.pl 1.motif 2.motif

######################################
#my $dir=$ENV{'PWD'};
#my $motif1=$ARGV[0];
#my $motif2=$ARGV[1];
#my $ligand=$ARGV[2];

sub motif_align{
	my ($dir,$motif1,$motif2,$ligand,$id) =@_;

my @length_cost=();
my @length_cost_file=`cat $dir/Loop_length_cost.matrix`;chomp(@length_cost_file);
for(my $i=2;$i<@length_cost_file;$i++)
{
	my @arr=split(/\s+/,$length_cost_file[$i]);
	for(my $j=1;$j<@arr;$j++)
	{
		$length_cost[$i-2][$j-1]=$arr[$j];
		#print $length_cost[$i-2][$j-1],"\t";
	}
	#print "\n";
}
my %base_pair_sub=();
my @base_pair_sub_file=`cat $dir/Base_pair_substitution.matrix`;chomp(@base_pair_sub_file);
my @base1=split(/\s+/,$base_pair_sub_file[1]);
my @base2=split(/\s+/,$base_pair_sub_file[2]);
my @base=();
for(my $i=0;$i<@base1;$i++)
{
	$base[$i]="$base1[$i]$base2[$i]";
}
for(my $i=3;$i<@base_pair_sub_file;$i++)
{
	my @bp_score=split(/\s+/,$base_pair_sub_file[$i]);
	for(my $j=0;$j<@bp_score-2;$j++)
	{
		$base_pair_sub{$base[$i-3]}{$base[$j]}=$bp_score[$j];
	}
}

my $ss_sub=0;
my $length_sub=0;
my $bp_sub=0;
my $score=0;
if($motif1=~/^(\S+)\s+loop:\s+(.+)/)
{
	my $type1=substr($1,0,1);
	my $seq1=$2;
	my $seq1_1=$2;
	if($motif2=~/^([A-Z])\d+\s+(.+)/)
	{
		my $type2=$1;
		my $seq2=$2;
		my $seq2_2=$2;
		exit if($type1 ne $type2);
		if($type1 eq $type2 && $type1 eq 'E')
		{
			my $length1=length($seq1);
			my $length2=length($seq2);
			exit if($length1 <5 || $length2 <5);
			$ss_sub=&NWalign($dir,$seq1,$seq2);
			$score=$ss_sub;
			#print $motif1,"\t",$seq2_2,"\t",$ligand,"\t",$score,"\n" if($score >0);
			if($score >0){
				my $result="$type1: $seq1_1;$type1: $seq2_2;$ligand;$id;$score\n";
				return $result;
			}
		}
		elsif($type1 eq $type2 && $type1 eq 'X')
                {
                        my $length1=length($seq1);
                        my $length2=length($seq2);
                        exit if($length1 <5 || $length2 <5);
                        $ss_sub=&NWalign($dir,$seq1,$seq2);
                        $score=$ss_sub;
			#print $motif1,"\t",$seq2_2,"\t",$ligand,"\t",$score,"\n" if($score >0);
			if($score >0){
				my $result="$type1: $seq1_1;$type1: $seq2_2;$ligand;$id;$score\n";
                                return $result;
                        }
                }
		elsif($type1 eq $type2 && $type1 eq 'H')
		{
			$seq1 =~ s/\(//g;$seq1=~s/\)//g;$seq1=~s/\,//g;$seq1=~s/\s+//g;
			$seq2 =~ s/\(//g;$seq2=~s/\)//g;$seq2=~s/\,//g;$seq2=~s/\s+//g;
			my $length1=length($seq1)-2;
			my $length2=length($seq2)-2;
			my $basepair1=substr($seq1,0,2);
			my $basepair2=substr($seq2,0,2);
			$seq1=substr($seq1,2);
			$seq2=substr($seq2,2);
			$ss_sub=&NWalign($dir,$seq1,$seq2);
			$length_sub=-abs($length_cost[$length1-1][1]-$length_cost[$length2-1][1]);
			$bp_sub=$base_pair_sub{$basepair1}{$basepair2};
			$score=$ss_sub+$length_sub+$bp_sub;
			#print $motif1,"\t",$seq2_2,"\t",$ligand,"\t",$score,"\n" if($score >0);
			if($score >0){
				my $result="$type1: $seq1_1;$type1: $seq2_2;$ligand;$id;$score\n";
                                return $result;
                        }
			
		}
		elsif($type1 eq $type2 && $type1 eq 'I')
		{
			my $basepair1_1="";
			my $basepair1_2="";
			my $basepair2_1="";
			my $basepair2_2="";
			if($seq1=~/\((.+)\)\s+(.+)\s+\((.+)\)\s+(.+)/)
			{
				$basepair1_1=$1;
				$basepair1_2=$3;
				$seq1="$2$4";
				$basepair1_1=~s/\,//g;$basepair1_2=~s/\,//g;
				#print $seq1,"\t",$basepair1_1,"\t",$basepair1_2,"\n";
			}
			if($seq2=~/\((.+)\)\s+(.+)\s+\((.+)\)\s+(.+)/)
                        {
                                $basepair2_1=$1;
                                $basepair2_2=$3;
                                $seq2="$2$4";
				$basepair2_1=~s/\,//g;$basepair2_2=~s/\,//g;
				#print $seq2,"\t",$basepair2_1,"\t",$basepair2_2,"\n";
                        }
                        my $length1=length($seq1);
                        my $length2=length($seq2);
                        $ss_sub=&NWalign($dir,$seq1,$seq2);
                        $length_sub=-abs($length_cost[$length1-1][3]-$length_cost[$length2-1][3]);
                        $bp_sub=$base_pair_sub{$basepair1_1}{$basepair2_1}+$base_pair_sub{$basepair1_2}{$basepair2_2};
                        $score=$ss_sub+$length_sub+$bp_sub;
			#print $motif1,"\t",$seq2_2,"\t",$ligand,"\t",$score,"\n" if($score >0);
			if($score >0){
				my $result="$type1: $seq1_1;$type1: $seq2_2;$ligand;$id;$score\n";
                                return $result;
                        }
		}
		elsif($type1 eq $type2 && $type1 eq 'B')
                {
                        my $basepair1_1="";
                        my $basepair1_2="";
                        my $basepair2_1="";
                        my $basepair2_2="";
                        if($seq1=~/\((.+)\)\s+(.+)\s+\((.+)\)/)
                        {
                                $basepair1_1=$1;
                                $basepair1_2=$3;
                                $seq1=$2;
                                $basepair1_1=~s/\,//g;$basepair1_2=~s/\,//g;
                                #print $seq1,"\t",$basepair1_1,"\t",$basepair1_2,"\n";
                        }
                        if($seq2=~/\((.+)\)\s+(.+)\s+\((.+)\)/)
                        {
                                $basepair2_1=$1;
                                $basepair2_2=$3;
                                $seq2=$2;
                                $basepair2_1=~s/\,//g;$basepair2_2=~s/\,//g;
                                #print $seq2,"\t",$basepair2_1,"\t",$basepair2_2,"\n";
                        }
                        my $length1=length($seq1);
                        my $length2=length($seq2);
                        $ss_sub=&NWalign($dir,$seq1,$seq2);
                        $length_sub=-abs($length_cost[$length1-1][2]-$length_cost[$length2-1][2]);
                        $bp_sub=$base_pair_sub{$basepair1_1}{$basepair2_1}+$base_pair_sub{$basepair1_2}{$basepair2_2};
                        $score=$ss_sub+$length_sub+$bp_sub;
			#print $motif1,"\t",$seq2_2,"\t",$ligand,"\t",$score,"\n" if($score >0);
			if($score >0){
				my $result="$type1: $seq1_1;$type1: $seq2_2;$ligand;$id;$score\n";
                                return $result;
                        }
                }
		elsif($type1 eq $type2 && $type1 eq 'M')
		{
			my @basepairs1=$seq1=~/\(([A-Z]\,[A-Z])\)/g;
			my @basepairs2=$seq2=~/\(([A-Z]\,[A-Z])\)/g;
			my @loop1=$seq1=~/\)\s+([A-Z]+)\s+\(/g;
			my @loop2=$seq2=~/\)\s+([A-Z]+)\s+\(/g;
			$seq1="";$seq2="";
			foreach my $lo1(@loop1)
			{
				$seq1="$seq1$lo1";
			}
			foreach my $lo2(@loop2)
			{
				$seq2="$seq2$lo2";
			}
			#print $seq1,"\n",$seq2,"\n",scalar(@basepairs1),"\t",scalar(@basepairs2),"\n";
			my $length1=length($seq1);
                        my $length2=length($seq2);
			if(($length1 == 0 && $length2 !=0) || ($length1 != 0 && $length2 ==0))
			{
				$ss_sub=-55;
			}
			elsif($length1 == 0 && $length2 ==0)
			{
				$ss_sub=0;
			}
			else
			{
				$ss_sub=&NWalign($dir,$seq1,$seq2);
			}
			my $num_bp=scalar(@basepairs1);
			$num_bp=scalar(@basepairs2) if(scalar(@basepairs1)>scalar(@basepairs2));
			for(my $i=0;$i<$num_bp;$i++)
			{
				$basepairs1[$i]=~s/\,//g;
				$basepairs2[$i]=~s/\,//g;
				#print $basepairs1[$i],"\t",$basepairs2[$i],"\n";
				$bp_sub+=$base_pair_sub{$basepairs1[$i]}{$basepairs2[$i]};
			}
			$score=$ss_sub+$bp_sub;
			#print $motif1,"\t",$seq2_2,"\t",$ligand,"\t",$score,"\n" if($score >0);
			if($score >0){
                                my $result="$type1: $seq1_1;$type1: $seq2_2;$ligand;$id;$score\n";
                                return $result;
                        }
		}
		
	}
}
}

sub NWalign{

my ($dir,$seq1,$seq2)=@_;
my $infinite=-1e+06;
my $gap_open=-110;
my $gap_extension=-55;

my %mutation=();
my @blosum=`cat $dir/Single_strand_substitution.matrix`;chomp(@blosum);shift(@blosum);
my @seq=split(/\s+/,$blosum[0]);
for(my $b=1;$b<@blosum;$b++)
{
        my @column=split(/\s+/,$blosum[$b]);
        for(my $c=1;$c<@column;$c++)
        {
                $mutation{$seq[$b],$seq[$c]}=$column[$c];
        }
}

my @M=();
my @m=();
my @Ix=();
my @ix=();
my @Iy=();
my @iy=();

my @seq1=split(//,$seq1);
my @seq2=split(//,$seq2);

my $i=0;
my $j=0;
$M[$i][$j]=0;
$Ix[$i][$j]=$gap_open;
$Iy[$i][$j]=$gap_open;
$m[$i][$j]='M';
$ix[$i][$j]='Ix';
$iy[$i][$j]='Iy';

for($j=1;$j<@seq2+1;$j++)
{
	$M[$i][$j]=$infinite;
	$Ix[$i][$j]=$infinite;
	$Iy[$i][$j]=$gap_open+$gap_extension*$j;
	$m[$i][$j]='M';
	$ix[$i][$j]='Ix';
	$iy[$i][$j]='Iy';
}
$j=0;
for($i=1;$i<@seq1+1;$i++)
{
	$M[$i][$j]=$infinite;
	$Ix[$i][$j]=$gap_open+$gap_extension*$i;
	$Iy[$i][$j]=$infinite;
        $m[$i][$j]='M';
        $ix[$i][$j]='Ix';
        $iy[$i][$j]='Iy';
}


my $score0=0;
my $score1=0;
my $score2=0;
for($i=1;$i<@seq1+1;$i++)
{
	for($j=1;$j<@seq2+1;$j++)
	{
		if($seq1[$i-1] ~~ @seq)
		#if(grep /^$seq1[$i-1]$/, @seq)
		{
			#print $seq1[$i-1],"\t";
		}
		else{
			$seq1[$i-1]='*';
			#print $seq1[$i-1],"\t";
		}
		if($seq2[$j-1] ~~ @seq)
		#if(grep /^$seq2[$j-1]$/, @seq)
                {
			#print $seq2[$j-1],"\n";
		}
		else{
                        $seq2[$j-1]='*';
			#print $seq2[$j-1],"\n";
                }
		$score2=$M[$i-1][$j-1]+$mutation{$seq1[$i-1], $seq2[$j-1]};
                $score1=$Ix[$i-1][$j-1]+$mutation{$seq1[$i-1], $seq2[$j-1]};
                $score0=$Iy[$i-1][$j-1]+$mutation{$seq1[$i-1], $seq2[$j-1]};
		($M[$i][$j],my $m)=&max($score0,$score1,$score2);
		if($m==0){
			$m[$i][$j]='Iy'
		}
		elsif($m==1){
			$m[$i][$j]='Ix'
		}
                elsif($m==2){
                        $m[$i][$j]='M'
                }
		#print $M[$i][$j],"\t";
		
		$score1=$M[$i-1][$j]+$gap_open;
		$score0=$Ix[$i-1][$j]+$gap_extension;
		($Ix[$i][$j],my $ix)=&max($score0,$score1);
                if($ix==0){
                        $ix[$i][$j]='Ix'
                }
                elsif($ix==1){
                        $ix[$i][$j]='M'
                }
		#print $ix[$i][$j],"\t";

		$score1=$M[$i][$j-1]+$gap_open;
                $score0=$Iy[$i][$j-1]+$gap_extension;
                ($Iy[$i][$j],my $iy)=&max($score0,$score1);
                if($iy==0){
                        $iy[$i][$j]='Iy'
                }
                elsif($iy==1){
                        $iy[$i][$j]='M'
                }
		#print $iy[$i][$j],"\t";

	}
	#print "\n";
}


my $max_score=$Iy[@seq1][@seq2];
$max_score=$Ix[@seq1][@seq2] if($Ix[@seq1][@seq2]>$max_score);
$max_score=$M[@seq1][@seq2] if($M[@seq1][@seq2]>$max_score);

return ($max_score);
#print "The alignment score is $max_score\n";

}

sub max{
	my @score=@_;
	my $max=$score[0];
	my $index=0;
	for(my $i=1;$i<@score;$i++)
	{
		if($score[$i]>$max)
		{
			$max=$score[$i];
			$index=$i;
		}	
	}
	return ($max,$index);
}
1;
