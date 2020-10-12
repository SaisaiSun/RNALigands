#!/usr/bin/perl -w
use strict;

##########################################
#This code is used to search the same motif of query sequence from the motif database and return the ligand.

#SaisaiSun
#Date:24-04-2020

##########################################


my $query_id=$ARGV[0];
my $output_dir=$ARGV[1];

my $datadir="/home/zhanglab3/sssun/RNA_motif-ligand";
my @database=`cat $datadir/All_site2motif2_new.txt`;chomp(@database);
my @id_motif=`cat $output_dir/$query_id.motif`;chomp(@id_motif);

foreach my $line(@id_motif)
{
	if($line=~/^([A-Z])\d+\s+(.+)/)
	{
		my $type=$1;
		my $motif=$2;
		foreach my $data(@database)
		{
			if($data=~/^$type\d+\s+(.+)\s+(\S+)/)
			{
				my $motif_db=$1;
				my $ligand=$2;
				if($motif eq $motif_db)
				{
					print $data,"\n";
				}
			}
		}
	}
}

