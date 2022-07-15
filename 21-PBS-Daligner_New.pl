#!/usr/bin/perl
#######################################################################################################
#   This script is used to align the path-contigs to the corresponding gap with Daligner
#                      Author: huilong du & Yilei
#                      Date: 2021/05/31 
#   Usage: perl $0 Scaffold2Ctg_Gap.txt Prosudo_ScaffoldNonEnzyme2Contig.fasta PathContigRename.fasta
#######################################################################################################
use warnings;
use strict;
my $infile1=shift;            #Scaffold2Ctg_Gap.txt
my $infile2=shift;            #Prosudo_ScaffoldNonEnzyme2Contig.fasta
my $infile3=shift;            #PathContigRename.fasta
my $queue=shift;
my $script=shift;
my $genome_name=shift;
my $DAZZ_DB=shift;            #DAZZ_DB
my $DALIGNER=shift;           #DALIGNER
my ${outfolder}=shift;

open IN1,"<$infile1" or die $!;
open IN2,"<$infile2" or die $!;
open IN3,"<$infile3" or die $!;

my %ScaffoldContig=();
my $sign="";
#Super-Scaffold_3.1
while(<IN2>){
        chomp;
        my $line=$_;
        if($line=~/^>(\S+)/){
                $sign=$1;
        }
        else{
                $ScaffoldContig{$sign}=$line;
        }
}
close(IN2);

my %PathContig=();
#Super-Scaffold_262-5-6.1
my $gap="";
while(<IN3>){
        chomp;
        my $line=$_;
        if($line=~/^>(\S+)\.(\d+)/){
                $gap=$1;
                $sign=$2;
        }
        else{
                $PathContig{$gap}{$sign}=$line;
        }
}
close(IN3);

while(<IN1>){
        chomp;
        my $line=$_;
        my @content=split /\s+/,$line;
        $content[0]=~/^(\S+)\.(\d+)/;
        my $scaffold=$1;
        my $First=$2;
        my $Second=$First+1;
        my $Gap_Info=$scaffold."-".$First."-".$Second;
        print "$Gap_Info\n";
        my $commond1=`mkdir $Gap_Info`;
#       my $commond1=`mkdir $Gap_Info`;
#       my $Abs_path="${outfolder}/06-Daligner"
        open OUT,">$Gap_Info/$Gap_Info.fasta" or die $!;
        print OUT ">$content[0]\n";
        print OUT "$ScaffoldContig{$content[0]}\n";
        print OUT ">$content[1]\n";
        print OUT "$ScaffoldContig{$content[1]}\n";
        my $temp=$PathContig{$Gap_Info};
        foreach my $key (keys %$temp){
                next if(length($PathContig{$Gap_Info}{$key})>800000);
                print OUT ">$Gap_Info.$key\n$PathContig{$Gap_Info}{$key}\n";
        }
        close(OUT);
        open OUT,">$Gap_Info.pbs" or die $!;
        print OUT "#!/bin/bash
#PBS -N $genome_name-DALIGNER-$Gap_Info
#PBS -o $Gap_Info.log
#PBS -e $Gap_Info.err
#PBS -l select=1:ncpus=1
#PBS -q $queue

cd ${outfolder}/06-Daligner/$Gap_Info

perl $script/lines_to_split.pl ${outfolder}/06-Daligner/$Gap_Info/$Gap_Info.fasta \
    ${outfolder}/06-Daligner/$Gap_Info/$Gap_Info-formated.fasta

#perl $script/SplitRef.pl $Gap_Info.fasta $Gap_Info-formated.fasta

$DAZZ_DB/fasta2DAM ${outfolder}/06-Daligner/$Gap_Info/$Gap_Info ${outfolder}/06-Daligner/$Gap_Info/$Gap_Info-formated.fasta

$DAZZ_DB/DBdust ${outfolder}/06-Daligner/$Gap_Info/$Gap_Info.dam

$DAZZ_DB/DBsplit -x1000 -s50 ${outfolder}/06-Daligner/$Gap_Info/$Gap_Info.dam

$DALIGNER/HPC.daligner ${outfolder}/06-Daligner/$Gap_Info/$Gap_Info.dam > ${outfolder}/06-Daligner/$Gap_Info/$Gap_Info.sh

time sh ${outfolder}/06-Daligner/$Gap_Info/$Gap_Info.sh

rm -f ${outfolder}/06-Daligner/$Gap_Info/$Gap_Info.*.$Gap_Info.*.?*.las

$DAZZ_DB/DBdump -rh ${outfolder}/06-Daligner/$Gap_Info/$Gap_Info.dam | perl $script/ParseDAZZDB.pl >${outfolder}/06-Daligner/$Gap_Info/ParseDAZZDB.txt

cat ${outfolder}/06-Daligner/$Gap_Info/$Gap_Info*.las >${outfolder}/06-Daligner/$Gap_Info/All.las

$DALIGNER/LAdump -cd ${outfolder}/06-Daligner/$Gap_Info/$Gap_Info.dam ${outfolder}/06-Daligner/$Gap_Info/All.las | perl $script/ParseLA.pl > ${outfolder}/06-Daligner/$Gap_Info/$Gap_Info-Final.txt
perl $script/Daligner_Reformate.pl ${outfolder}/06-Daligner/$Gap_Info/$Gap_Info-Final.txt ${outfolder}/06-Daligner/$Gap_Info/$Gap_Info-Final_Reformated.txt
";
                close(OUT);
                sleep(3);
                #my $commond=`qsub $Gap_Info.pbs`;
                system("qsub $Gap_Info.pbs");
}
close(IN1);
