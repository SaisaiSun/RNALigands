#!/usr/bin/perl -w
#use strict;

################################
#This program is to generate position weight matrix from .pwm file for each motif.

#Saisai Sun
#Date: 22-06-2020

################################

my $dir="$ENV{PWD}";
my $datadir="$dir/../dataset";
my @list=`cat $dir/../RNA_HET_bp_site_cdhit_list2`;chomp(@list);
open(OUT1,">$dir/All_H_PWM.txt");
open(OUT2,">$dir/All_B_PWM.txt");
open(OUT3,">$dir/All_I_PWM.txt");
open(OUT4,">$dir/All_M_PWM.txt");
open(OUT5,">$dir/All_E_PWM.txt");
open(OUT6,">$dir/All_X_PWM.txt");
open(OUT7,">$dir/All_PK_PWM.txt");

foreach my $l(@list)
{
	#$l="1arj_N";
	#print $l,"\n";
#=pod
	my @pwm=`cat $datadir/$l/$l.pwm`;chomp(@pwm);
	my @motif=`cat $datadir/$l/$l.motif`;chomp(@motif);
	my @motif_num=`cat $datadir/$l/$l\_num.motif`;chomp(@motif_num);
	my @motif_ligand=`cat $datadir/$l/$l\_motif_ligand`;chomp(@motif_ligand);
	for(my $i=0;$i<@motif_num;$i++)
	{
		if($motif[$i] ~~ @motif_ligand)
		{
		if($motif_num[$i]=~/^H\d+\s+\((\d+)\,(\d+)\)/)
		{
			my $position1=$1;
			my $position2=$2;
			my @position=();
			for(my $i=$position1;$i<$position2+1;$i++)
			{
				push(@position,$i);
			}
			my @result=&match_pwm($motif[$i],\@pwm,\@position);
			print OUT1 $motif[$i],"\n";
			foreach my $line(@result)
			{
				print OUT1 $line,"\n";
			}
		}
		elsif($motif_num[$i]=~/^B\d+\s+\((\d+)\,(\d+)\).+\((\d+)\,(\d+)\)/)
		{
			my $position1=$1;
			my $position2=$3;
			my $position3=$4;
			my $position4=$2;
			my @position=();
			for(my $i=$position1;$i<$position2+1;$i++)
                        {
                                push(@position,$i);
                        }
			push(@position,$position3);push(@position,$position4);
			my @result=&match_pwm($motif[$i],\@pwm,\@position);
			print OUT2 $motif[$i],"\n";
                        foreach my $line(@result)
                        {
                                print OUT2 $line,"\n";
                        }
		}
		elsif($motif_num[$i]=~/^I\d+\s+\((\d+)\,(\d+)\).+\((\d+)\,(\d+)\).+/)
		{
			my $position1=$1;
			my $position2=$4;
			my $position3=$3;
			my $position4=$2;
			my @position=();
			for(my $i=$position1;$i<$position2+1;$i++)
                        {
                                push(@position,$i);
                        }
			for(my $i=$position3;$i<$position4+1;$i++)
                        {
                                push(@position,$i);
                        }
			my @result=&match_pwm($motif[$i],\@pwm,\@position);
			print OUT3 $motif[$i],"\n";
                        foreach my $line(@result)
                        {
                                print OUT3 $line,"\n";
                        }
		}
		elsif($motif_num[$i]=~/^M/)
		{
			my @position1=$motif_num[$i]=~/(\d+)\,/g;
			my @position2=$motif_num[$i]=~/\,(\d+)/g;
			my @position=();
			for(my $i=$position1[0];$i<$position1[1]+1;$i++)
                        {
                                push(@position,$i);
                        }			
			for(my $j=1;$j<@position2-1;$j++)
			{
				for(my $i=$position2[$j];$i<$position1[$j+1]+1;$i++)
				{
					push(@position,$i);
				}
			}
			my @result=&match_pwm($motif[$i],\@pwm,\@position);
			print OUT4 $motif[$i],"\n";
                        foreach my $line(@result)
                        {
                                print OUT4 $line,"\n";
                        }
		}
		elsif($motif_num[$i]=~/^E\d+\s+(\d+)\.\.(\d+)/)
		{
			my $position1=$1;
			my $position2=$2;
			my @position=();
                        for(my $i=$position1;$i<$position2+1;$i++)
                        {
                                push(@position,$i);
                        }
			my @result=&match_pwm($motif[$i],\@pwm,\@position);
			print OUT5 $motif[$i],"\n";
                        foreach my $line(@result)
                        {
                                print OUT5 $line,"\n";
                        }
		}
		elsif($motif_num[$i]=~/^X\d+\s+(\d+)\.\.(\d+)/)
                {       
                        my $position1=$1;
                        my $position2=$2;
                        my @position=();
                        for(my $i=$position1;$i<$position2+1;$i++)
                        {       
                                push(@position,$i);
                        }
                        my @result=&match_pwm($motif[$i],\@pwm,\@position);
			print OUT6 $motif[$i],"\n";
                        foreach my $line(@result)
                        {
                                print OUT6 $line,"\n";
                        }
                }
		elsif($motif_num[$i]=~/^PK/)
                {
                        my @position1=$motif_num[$i]=~/(\d+)\,/g;
                        my @position2=$motif_num[$i]=~/\,(\d+)/g;
                        my @position=();
                        for(my $i=$position1[0];$i<$position1[-1]+1;$i++)
                        {
                                push(@position,$i);
                        }
                        for(my $i=$position2[0];$i<$position2[-1]+1;$i++)
                        {
                                push(@position,$i);
                        }
			my @result=&match_pwm($motif[$i],\@pwm,\@position);
                        print OUT7 $motif[$i],"\n";
                        foreach my $line(@result)
                        {
                                print OUT7 $line,"\n";
                        }
		}
		}

	}
	#close(OUT);
#=cut
	#last;
}
close(OUT1);
close(OUT2);
close(OUT3);
close(OUT4);
close(OUT5);
close(OUT6);
close(OUT7);

sub match_pwm
{
	(my $motif,my $pwm,my $position)=@_;
	my @PWM=@$pwm;
	my @result=();
	#print $motif,"\n";
	for(my $i=0;$i<@PWM;$i++)
	{
		#print $PWM[$i],"\n" if($i+1 ~~ @$position);
		push(@result,$PWM[$i]) if($i+1 ~~ @$position);
	}
	return(@result);	
}
