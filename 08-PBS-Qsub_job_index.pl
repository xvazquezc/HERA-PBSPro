#!/usr/bin/perl
#Author: huilong du & YiLei
#Date: 2021/05/31
#Note: Based on 08-qsub_job_index.pl

use warnings;
use strict;

my $infile=shift;
my $queue=shift;
my $genome=shift;
my ${Bwa}=shift;

open IN,"<$infile" or die $!;

my $count=1;

while(<IN>){
    chomp;
    my $line=$_;
    open OUT,">$count.pbs" or die $!;
    print OUT "#!/bin/bash
#PBS -N $genome-INDEX-$count
#PBS -o $count.log
#PBS -e $count.err
#PBS -l nodes=1:ppn=1
#PBS -q $queue

${Bwa}/bwa index $line
echo Make index for every part of the pacbios and non-scaffolded contigs finished at `date`
";
     close(OUT);
     system("qsub $count.pbs");
     $count++;
}
close(IN);
