#!usr/bin/perl -w
=head
1. clean reads data;
2. mapping;
3. ANGsd in all;
Update: 2015-12-31 15:24:30
Version 3
This will just do the mapping.
=cut
use strict;
use Getopt::Long;
#die "perl mapping.pl mapping.list; the format of mappling.list:\$dir/\$readsName	\$newName\n" unless @ARGV==1;

my $reference;
my $suffix = "";

GetOptions(
	"indica" => \$reference
);

my $cpu=8;
my $maxDep=150;
my $fa="/home/hongru/work/database/mouse/genomes/osativa/japonica/nipponbare7/np7.fasta";

if ($reference){
	$fa="/mnt/mouse/genomes/osativa/indica/minghui63/minghui63.fasta";
	$suffix = "i";
}
my $phred="phred33";
my $samtools="samtools";
my $samtools_old="/usr/biobin/samtools-0.1.19/samtools";
my $bwa="bwa";
my $gatk="/usr/biobin/GenomeAnalysisTK-v3.5.jar";
my $angsd = "angsd0614";
my $trimmomatic="/usr/biobin/trimmomatic-0.36.jar";
my $name;

while (<>){
	chomp;
	if(/^#/){next;}
	my @tmp=split;
	my $reads=$tmp[0];
	$name=$tmp[1];
	$name = $name.$suffix;

    system "java -jar $trimmomatic PE -threads $cpu -$phred $reads\_1.fq.gz $reads\_2.fq.gz  $name.1p.fq $name.1u.fq $name.2p.fq $name.2u.fq SLIDINGWINDOW:4:20 LEADING:20 TRAILING:20 MINLEN:30 TOPHRED33 && echo trimmomatic_ok";
	## check file before and after filtering;
=head	
	my @all; my @fq1; my @fq2; my $add;
	my $out = `du -b $reads\_1.fq.gz`; my $fq1 =  (split /\s/, $out)[0];	
	$out = `du -b $reads\_1.fq.gz`;    my $fq2 =  (split /\s/, $out)[0];
	push @all, $fq1+$fq2; push @fq1, $fq1; push @fq2, $fq2;	
	$out = `du -b $name.1p.fq`; $fq1 = (split /\s/, $out)[0]; push @fq1, $fq1; $add = $fq1;
    $out = `du -b $name.2p.fq`; $fq2 = (split /\s/, $out)[0]; push @fq2, $fq2; $add+=$fq2;
	$out = `du -b $name.1u.fq`; $fq1 = (split /\s/, $out)[0]; push @fq1, $fq1; $add+=$fq1;
    $out = `du -b $name.2u.fq`; $fq2 = (split /\s/, $out)[0]; push @fq2, $fq2; $add+=$fq2;
	push @all, $add;
	print "@all\n@fq1\n@fq2\n";
=cut
	## Start mapping here.
    system "$bwa  mem -M -t $cpu -v 1 -R '\@RG\tID:$name\tSM:$name\tLB:$name'  $fa $name.1p.fq $name.2p.fq |$samtools view -bS -  > $name.p.bam && echo bwa_pair_ok && rm $name.1p.fq $name.2p.fq";
    system "$bwa  mem -M -t $cpu -v 1 -R '\@RG\tID:$name\tSM:$name\tLB:$name'  $fa $name.1u.fq |$samtools view -bS -  > $name.1u.bam && echo bwa_1u_ok && rm $name.1u.fq";
    system "$bwa  mem -M -t $cpu -v 1 -R '\@RG\tID:$name\tSM:$name\tLB:$name'  $fa $name.2u.fq |$samtools view -bS -  > $name.2u.bam && echo bwa_2u_ok && rm $name.2u.fq";
    system "$samtools merge -h $name.p.bam $name.bam $name.p.bam $name.1u.bam $name.2u.bam && echo merge_ok && rm $name.p.bam $name.1u.bam $name.2u.bam";
    system "$samtools fixmate  $name.bam $name\_fixmate.bam  && echo fixmate_ok && rm $name.bam";
    system "$samtools sort -@ $cpu -O bam -o $name.bam $name\_fixmate.bam && echo fixmate_ok_sort_ok && rm $name\_fixmate.bam";
	system "$samtools rmdup $name.bam $name.rmdup.bam && echo $name\_rmdup_ok && mv $name.rmdup.bam $name.bam";
	system "$samtools index $name.bam && echo $name\_index_OK";
	## I will use GATK to fix gap below
	system "java -Xmx8g -jar $gatk -T RealignerTargetCreator -R $fa -I $name.bam -o $name.intervals && echo $name\_intervals_ok";
	system "java -Xmx8g -jar $gatk -T IndelRealigner  -R $fa -I $name.bam -targetIntervals $name.intervals -o $name.gatk.bam && echo $name\_gatk_ok && mv $name.gatk.bam $name.bam";
	system "$samtools index $name.bam && echo $name\_index_OK";
}
