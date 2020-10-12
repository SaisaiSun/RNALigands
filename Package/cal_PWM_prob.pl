#!/usr/bin/perl -w
use strict;

###################################
#This code is used to calculate the probability for a query motif from PWM.

#Saisai Sun
#Date:22-6-2020

###################################


my $query_id=$ARGV[0];
my $output_dir=$ARGV[1];
my $datadir="/home/zhanglab3/sssun/RNA_motif-ligand";
my @id_motif=`cat $output_dir/$query_id.motif`;chomp(@id_motif);

	foreach my $line(@id_motif)
	{
		if($line=~/^(M)\d+\s+(.+)/)
		{
			my $type=$1;
			my $seq=$2;
			$seq =~ s/\(//g;$seq=~s/\)//g;$seq=~s/\,//g;$seq=~s/\s+//g;
			my $length=length($seq);
			my $first=substr($seq,0,1); my $second=substr($seq,1,1); my $mid=substr($seq,2,$length-4);
			$seq="";$seq=$first.$mid.$second;
			&cal_weight($line,$seq,$type);
			#print $seq,"\n";
		}
		elsif($line=~/^(PK)\d+\s+(.+)/)
                {
                        my $type=$1;
                        my $seq=$2;
			my @base1=$seq=~/([A-Z])\,/g;
                        my @base2=$seq=~/\,([A-Z])/g;
			my $seq_new="";
			foreach my $b1(@base1){	$seq_new.=$b1;	}
			foreach my $b2(@base2){ $seq_new.=$b2;  }
			my $length=length($seq_new);
			&cal_weight($line,$seq_new,$type);
		}
		elsif($line=~/^(E)\d+\s+(.+)/)
		{
			my $type=$1;
                        my $seq=$2;
			next if(length($seq)<3);
			&cal_weight($line,$seq,$type);
		}
		elsif($line=~/^(X)\d+\s+(.+)/)
                {
                        my $type=$1;
                        my $seq=$2;
			next if(length($seq)<3);
                        &cal_weight($line,$seq,$type);
                }
		elsif($line=~/^([A-Z])\d+\s+(.+)/)
		{
			my $type=$1;
                        my $seq=$2;
			next if(length($seq)<3);
                        $seq =~ s/\(//g;$seq=~s/\)//g;$seq=~s/\,//g;$seq=~s/\s+//g;
                        my $first=substr($seq,0,1); my $second=substr($seq,1,1); my $mid=substr($seq,2);
                        $seq="";$seq=$first.$mid.$second;
			&cal_weight($line,$seq,$type);
			#print $seq,"\n";
		}
	}

sub cal_weight{

	my $datadir="/home/zhanglab3/sssun/RNA_motif-ligand/package";
	my ($line,$seq,$type)=@_;
	my $query_len=length($seq);
	my @seq=split(//,$seq);
	my @Weight=();
	my @Motif=();
	#if($type=~/^H/)
	#{
		#print "Query motif: $line\t";
		my @database=`cat $datadir/All_$type\_PWM.txt`;chomp(@database);
		my $flag=0;
		my $weight=0;
		my $j=0;
		for(my $i=0;$i<@database;$i++)
		{
			if($database[$i]=~/^$type\d+\s+(.+)/)
                        {
				my $motif=$1;
				$motif =~ s/\(//g;$motif=~s/\)//g;$motif=~s/\,//g;$motif=~s/\s+//g;
				my $target_len=length($motif);
				next if($query_len != $target_len);
				$flag=1;
				push(@Motif,$database[$i]);
				#print $database[$i];
                        }
			elsif($flag == 1)
			{
				my @arr=split(/\s+/,$database[$i]);
				if($seq[$j] eq "A")
				{
					$weight=$weight+$arr[1];
				}
				elsif($seq[$j] eq "U")
				{
					$weight=$weight+$arr[2];
				}
				elsif($seq[$j] eq "G")
                                {
                                       	$weight=$weight+$arr[3];
                                }
				elsif($seq[$j] eq "C")
                                {
                                       	$weight=$weight+$arr[4];
                                }
				#print $weight,"\t",$j,"\n";
				$j++;
				if(($i<@database-1 && $database[$i+1]=~/^$type/) || $i==@database-1)
				{
					#print $weight,"\n";
					push(@Weight,$weight);
					$j=0;
					$flag=0;
					$weight=0;
				}
			}
		}
	#}
	for(my $w=0;$w<@Weight;$w++)
	{
		if($Weight[$w]>0)
		{
			print "Query motif: $line\t";
			print "Target motif $Motif[$w]\t";
			print $Weight[$w],"\n";
		}
	}
	#return @Weight;
}
