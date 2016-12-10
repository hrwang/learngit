## Usage: perl snpIndex.pl list
## The format of the list: "$poolfile\t$reffile"
## This is different from the old one in the way that it also produce vcf for reference. and then cop with filterVCF.new.pl
use File::Basename;
use strict;
use warnings;
my $cpu=8;
my $maxDep=150;
my $fa="/home/hongru/work/database/mouse/genomes/osativa/japonica/nipponbare7/nipponbare7";
my $samtools="samtools";
my $angsd = "angsd0614";

while (<>){
  chomp;
  my @tmp=split;
  my $poolFile = $tmp[0];
  my $poolName = basename ($tmp[0]);
  my $refFile = $tmp[1];
  my $refName = basename ($tmp[1]);
	open (FH, ">$poolName.bamlist");
	# my $b= `pwd`;
	# chomp $b;
	my $c=$poolFile.".bam";
	my $d=$refFile.".bam";
	print FH "$c\n";
	print FH "$d\n";
  for (my $i=1;$i<=12;$i++){
    system "$angsd -nThreads $cpu -bam $poolName.bamlist -uniqueOnly 1 -remove_bads 1 -minMapQ 20 -minQ 20 -baq 1 -C 50 -ref $fa -r Chr$i:1- -out $poolName.Chr$i -doCounts 1 -GL 1 -doMajorMinor 2 -doMaf 2 -SNP_pval 0.01";
		if ($i==1) {system "zcat $poolName.Chr$i.mafs.gz | sed \'1d\' | cut -f1-2 > $poolName.pos";
		} else {system "zcat $poolName.Chr$i.mafs.gz | sed \'1d\' | cut -f1-2 >> $poolName.pos";}
  }
  system "rm $poolName.Chr*.arg";
  system "/usr/biobin/samtools-0.1.19/samtools mpileup -C50 -q20 -Q 20  -l $poolName.pos  -ugf ~/work/database/mouse/genomes/osativa/japonica/nipponbare7/nipponbare7 $refFile.bam | /usr/biobin/bcftools-0.1.19/bcftools  view -cg - > $refName-$poolName.vcf";
  system "/usr/biobin/samtools-0.1.19/samtools mpileup -C50 -q20 -Q 20  -l $poolName.pos  -ugf ~/work/database/mouse/genomes/osativa/japonica/nipponbare7/nipponbare7 $poolFile.bam | /usr/biobin/bcftools-0.1.19/bcftools  view -cg - > $poolName.vcf";
}


=head
G

while (<>){
	chomp;
	my @tmp=split;
	my $poolFile = $tmp[0];
	my $poolName = basename ($tmp[0]);
	my $refFile = $tmp[1];
	my $refName = basename ($tmp[1]);
	for (my $chr=1;$chr<=12;$chr++){
		## Below I will call informative sites from both ref and pool.
    system "zcat $poolFile.Chr$chr.mafs.gz | awk '{if (\$6>0.2){print \$0}}'| grep -v chrom > $poolName.Chr$chr.filter.mafs";
    system "zcat $poolFile.Chr$chr.mafs.gz | awk '{if (\$3 != \$5 && \$6<=0.2){print \$0}}'| grep -v chrom >> $poolName.Chr$chr.filter.mafs";
    system "sort -k2,2n $poolName.Chr$chr.filter.mafs > $poolName.Chr$chr.tmp";
    system "mv $poolName.Chr$chr.tmp $poolName.Chr$chr.filter.mafs";

    system "zcat $refFile.Chr$chr.mafs.gz | awk '{if (\$6>0.2){print \$0}}'| grep -v chrom > $refName.Chr$chr.filter.mafs";
    system "zcat $refFile.Chr$chr.mafs.gz | awk '{if (\$3 != \$5 && \$6<=0.2){print \$0}}'| grep -v chrom >> $refName.Chr$chr.filter.mafs";
    system "sort -k2,2n $refName.Chr$chr.filter.mafs > $refName.Chr$chr.tmp";
    system "mv $refName.Chr$chr.tmp $refName.Chr$chr.filter.mafs";
		### Done
		## Below I will subtract ref SNPs from the pool SNPs.
    my %hash=();
		my @sites; # To store all the sites for writing into filehandle.
    open (LS,"$refName.Chr$chr.filter.mafs") or die "give me a bank!~~: $!";
    while(<LS>){
        chomp;
        my @tmp=split;
        $hash{$tmp[1]}=$_;
    }
    close LS;

    open (FILE,"$poolName.Chr$chr.filter.mafs") or die "kidding me?! give me something to filter!: $!";
    while(<FILE>){
        chomp;
        my @tmp=split /\t/;
        if (exists $hash{$tmp[1]}){next;}
        push @sites, "$tmp[0]\t$tmp[1]\n";
    }
    close FILE;
		### Done
		# Below I will write viable sites into a list to supply into samtools mpileup;
		if ($chr==1){open (FH, ">$poolName.pos");} else {open (FH, ">>$poolName.pos");}
        select FH;
				foreach (@sites){
        	print "$_";
				}
        close FH;
		# Below I will do some cleaning;
		system "rm $refName.Chr$chr.filter.mafs";
		system "rm $poolName.Chr$chr.filter.mafs";
		### Cleaning done.
	}
	## Below I will do the samtools running;
	system "/usr/biobin/samtools-0.1.19/samtools mpileup -C50 -q20 -Q 20  -l $poolName.pos  -ugf ~/work/database/mouse/genomes/osativa/japonica/nipponbare7/nipponbare7 $refFile.bam | /usr/biobin/bcftools-0.1.19/bcftools  view -cg - > $refName-$poolName.vcf";
	system "/usr/biobin/samtools-0.1.19/samtools mpileup -C50 -q20 -Q 20  -l $poolName.pos  -ugf ~/work/database/mouse/genomes/osativa/japonica/nipponbare7/nipponbare7 $poolFile.bam | /usr/biobin/bcftools-0.1.19/bcftools  view -cg - > $poolName.vcf";
}
