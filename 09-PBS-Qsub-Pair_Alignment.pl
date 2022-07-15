#!/usr/bin/perl
#Author: huilong du & YiLei
#Date: 2021/05/31
#Note: Based on 04-Qsub-Mapping2Ctg.pl
#Note: align the corrected pacbios and non-scaffolded contigs to scaffolded contigs

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
#PBS -l select=1:ncpus=1
#PBS -q $queue

${Bwa}/bwa index $line
echo Make index for every part of the pacbios and non-scaffolded contigs finished at `date`
";
     close(OUT);
     system("qsub $count.pbs");
     $count++;
}
close(IN);

(base) [stu_wuyilei@mn01 HERA]$ vim 08-PBS-Qsub_job_index.pl
(base) [stu_wuyilei@mn01 HERA]$ clear
(base) [stu_wuyilei@mn01 HERA]$ cat 09-PBS-Qsub-Pair_Alignment.pl
#!/usr/bin/perl
#Author: huilong du
#Note: align the corrected pacbios and non-scaffolded contigs selfly
use warnings;
use strict;

my $infile=shift;
my $script=shift;
my $queue=shift;
my $genome=shift;
my $outfolder=shift;
my ${Bwa}=shift;

open IN,"<$infile" or die $!;


my $count=1;
my @part=();
while(<IN>){
        chomp;
        my $line=$_;
        push(@part,$line);
}
close(IN);
for(my $i=0;$i<@part;$i++){
    for(my $j=$i;$j<@part;$j++){
        open OUT,">$count.pbs" or die $!;
        if($count<350){
                print OUT "#!/bin/bash
#PBS -N $genome-Pair-$i-$j
#PBS -o $count.log
#PBS -e $count.err
#PBS -l select=1:ncpus=1
#PBS -q $queue
";
        }
        elsif($count>=350){
                print OUT "#!/bin/bash
#PBS -N $genome-Pair-$i-$j
#PBS -o $count.log
#PBS -e $count.err
#PBS -l select=1:ncpus=1
#PBS -q $queue
";
        }
        if($i==$j){
           print OUT "

${Bwa}/bwa mem -a -e -t 1 $part[$i] $part[$j] >${outfolder}/03-Pacbio-SelfAlignment/Part_SelfAlignment_$i-$j.sam
echo Bwa SelfAlignment finished at `date`

perl $script/sam2blasr.pl ${outfolder}/03-Pacbio-SelfAlignment/Part_SelfAlignment_$i-$j.sam ${outfolder}/03-Pacbio-SelfAlignment/Part_SelfAlignment_$i-$j.txt
echo Sam2blasr finished at `date`

rm -f ${outfolder}/03-Pacbio-SelfAlignment/Part_SelfAlignment_$i-$j.sam
";
        }
        else{
           print OUT "
${Bwa}/bwa mem -a -t 1 $part[$i] $part[$j] >${outfolder}/03-Pacbio-SelfAlignment/Part_SelfAlignment_$i-$j.sam
echo Bwa SelfAlignment finished at `date`

perl $script/sam2blasr.pl ${outfolder}/03-Pacbio-SelfAlignment/Part_SelfAlignment_$i-$j.sam ${outfolder}/03-Pacbio-SelfAlignment/Part_SelfAlignment_$i-$j.txt
echo Sam2blasr finished at `date`

rm -f ${outfolder}/03-Pacbio-SelfAlignment/Part_SelfAlignment_$i-$j.sam
";
        }
        close(OUT);
        system("qsub $count.pbs");
        $count++;
    }
}
