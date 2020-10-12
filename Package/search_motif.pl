#!/usr/bin/perl -w
use strict;
use File::Spec;
require "./Motif-align.pl";

##########################################
#This code is used to search the same motif of query sequence from the motif database and return the ligand.

#SaisaiSun
#Date:24-04-2020

##########################################

my $infile=$ARGV[0];
my ($file, $output_dir, $query_id) = File::Spec->splitpath($infile);

my $script_dir=$ENV{"PWD"};
my @database=`cat $script_dir/All_site2motif2_new.txt`;chomp(@database);
my @id_motif=`cat $output_dir/$query_id\_motif_final.out`;chomp(@id_motif);

my $num=0;
foreach my $line(@id_motif)
{
	if($line=~/^(\S+)\s+loop:\s+(.+)/)
	{
		my $type=substr($1,0,1);
		my $motif=$2;
		my $result_str="";
		$num++;
		open(OUT2,">$output_dir/$query_id\_$num.out");
		foreach my $data(@database)
		{
			if($data=~/^($type\d+\s+.+)\s+(\S+)/)
			{
				my $motif_db=$1;
				my $ligand=$2;
				$motif_db=~s/ +$//;
				#print $line,"\t",$motif_db,"\n";
				#if($motif eq $motif_db)
				#{
				#	print $data,"\n";
				#}
				#$result_str .= `$script_dir/Motif-align.pl '$line' '$motif_db' '$ligand'`;
				$result_str .=&motif_align($script_dir,$line,$motif_db,$ligand,'111')
			}
		}
		print OUT2 $result_str;
		close(OUT2);
	}
}


open(OUT1,">$output_dir/$query_id\_final.out");
for(my $i=1;$i<$num+1;$i++)
{
	my $result_file="$output_dir/$query_id\_$i.out";
	my $max_score_motif=&max_score($result_file);
	print OUT1 $max_score_motif,"\n";
	#print OUT1 "success\n";
}
close(OUT1);

sub max_score{

	(my $file)=@_;
	my @result=`cat $file`;chomp(@result);
	my $max=0;
	my @index=();
	for(my $i=0;$i<@result;$i++)
	{
		my @arr=split(/;/,$result[$i]);
		if($arr[-1]>$max)
		{
			$max=$arr[-1];
			push(@index,$i);
		}
	}
	#$result[$index[-1]]=~/^[A-Z]:\s+\S+\s+[A-Z]:\s+\S+\s+(\S+)\s+(\d+)$/;
	return $result[$index[-1]];
}
