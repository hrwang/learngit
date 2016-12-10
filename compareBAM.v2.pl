=head1 Command-line Option
  --threads <num>       maximum number of threads used for this program.   
  --list <str>          the list containing reads directory of case and control. 
  --ref <str>           which reference to use for this mapping? should be one of (np7|mh63|r498|tair10).
  --help                output help information to screen  
=cut

use strict;
use warnings;
use File::Basename;
use Getopt::Long;
use Parallel::ForkManager;
my $maxDep=150;
my $angsd = "angsd0614";

my ($Threads, $List, $Reference, $Help);
GetOptions(
	"threads:n"=>\$Threads,
	"list:s"=>\$List,
	"ref:s"=>\$Reference,
	"help"=>\$Help
);

die `pod2text $0` if ($Help);
die `pod2text $0` unless (defined ($Threads) && defined ($List) && defined ($Reference));
$Reference = "\$".$Reference;
$Reference = `echo $Reference`;
chomp ($Reference);

my @list;
open (LS, $List) or die "list should be file that can be opened.\n";
while (<LS>){
	if(/^#/){next;}
	chomp;
	push @list, $_;
}
close LS;

### Find SNP between case and control.
foreach my $list1 (@list) {
	my @tmp=split /\s+/, $list1;
	my $poolFile = $tmp[0];
	my $poolName = basename ($tmp[0]);
	my $refFile = $tmp[1];
	my $refName = basename ($tmp[1]);
	open (FH, ">$poolName.bamlist");
	my $c=$poolFile.".bam";
	my $d=$refFile.".bam";
	print FH "$c\n";
	print FH "$d\n";
	close FH;
	
	my $pm=new Parallel::ForkManager($Threads);
	for (my $i=1;$i<=12;$i++){
		$pm->start and next;
		system "$angsd -nThreads 1 -bam $poolName.bamlist -uniqueOnly 1 -remove_bads 1 -minMapQ 20 -minQ 20 -baq 1 -C 50 -ref $Reference -r Chr$i:1- -out $poolName.Chr$i -doCounts 1 -GL 1 -doMajorMinor 2 -doMaf 2 -SNP_pval 0.01";
		if ($i==1) {
			system "zcat $poolName.Chr$i.mafs.gz | sed \'1d\' | cut -f1-2 > $poolName.pos";
		} else {
			system "zcat $poolName.Chr$i.mafs.gz | sed \'1d\' | cut -f1-2 >> $poolName.pos";
		}
		$pm->finish;
	}
	$pm->wait_all_children;
	system "rm $poolName.Chr*.arg";
}
##########


### Call VCF for case and control at SNP sites.
my $pmm=new Parallel::ForkManager($Threads);
foreach my $list2 (@list) {
	$pmm->start and next;
	my @tmp=split /\s+/, $list2;
	my $poolFile = $tmp[0];
	my $poolName = basename ($tmp[0]);
	my $refFile = $tmp[1];
	my $refName = basename ($tmp[1]);
	system "/usr/biobin/samtools-0.1.19/samtools mpileup -C50 -q20 -Q 20  -l $poolName.pos  -ugf $Reference $refFile.bam | /usr/biobin/bcftools-0.1.19/bcftools  view -cg - > $refName-$poolName.vcf";
	system "/usr/biobin/samtools-0.1.19/samtools mpileup -C50 -q20 -Q 20  -l $poolName.pos  -ugf $Reference $poolFile.bam | /usr/biobin/bcftools-0.1.19/bcftools  view -cg - > $poolName.vcf";
	$pmm->finish;
}
$pmm->wait_all_children;



