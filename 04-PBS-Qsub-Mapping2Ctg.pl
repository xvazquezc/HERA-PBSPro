#!/usr/bin/perl
#Author: YiLei
#Date: 2021/05/31
#Note: Based on 04-Qsub-Mapping2Ctg.pl
#Note: align the corrected pacbios and non-scaffolded contigs to scaffolded contigs

use warnings;
use strict;

my $infile=shift;
my $ref=shift;
my $output=shift;
my $script=shift;
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
#PBS -N $genome-Map-$count
#PBS -o $count.log
#PBS -e $count.err
#PBS -l nodes=1:ppn=1
#PBS -q $queue

${Bwa}/bwa mem -a -t 1 $ref $line >$output/Part_Alignment_$count.sam
echo bwa aligment finished at `date`

perl $script/sam2blasr.pl $output/Part_Alignment_$count.sam $output/Part_Alignment_$count.txt
echo sam2blasr finished at `date`

rm -f $output/Part_Alignment_$count.sam
";
     close(OUT);
     system("qsub $count.pbs");
     $count++;
}
close(IN);
