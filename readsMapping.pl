=head1 Command-line Option
  --threads <num>       maximum number of threads used for this program.   
  --list <str>          the list containing reads directory and output name. 
  --ref <str>           which reference to use for this mapping? should be one of (np7|mh63|r498|tair10).
  --help                output help information to screen  
=cut

use strict;
use Getopt::Long;
use Parallel::ForkManager;

my $samtools="samtools";
my ($Threads, $List, $Reference, $Help);
GetOptions(
	"threads:n"=>\$Threads,
	"list:s"=>\$List,
	"ref:s"=>\$Reference,
	"help"=>\$Help
);

die `pod2text $0` if ($Help);
die `pod2text $0` unless (defined ($Threads) && defined ($List) && defined ($Reference));
my $tag = $Reference;
$Reference = "\$".$Reference;
$Reference = `echo $Reference`;
chomp ($Reference);

my @list;
open (LS, $List) or die "list should be file that can be opened.\n";
while (<LS>){
	if(/^#/){next;}
	chomp;
	my ($dir, $name) = split;
	$name = $name.".".$tag;
	push @list, "$dir\t$name";
}
close LS;
#print "Reference:$Reference\nLists:@list\nThreads:$Threads\n";


foreach my $list (@list){
	my ($dir, $name) = split /\s+/, $list;
	my @reads = glob("$dir*");
	
	&trimFastq(FASTQ=>[@reads], NAME=>$name, THREADS=>$Threads);
	
	my @bwa_p_reads = ("$name.1p.fq", "$name.2p.fq");
	my @bwa_1u_reads = ("$name.1u.fq");
	my @bwa_2u_reads = ("$name.2u.fq");
	&fastq2bwa2bam(FASTQ=>[@bwa_p_reads], NAME=>"$name.p", REF=>$Reference, THREADS=>$Threads);
	&fastq2bwa2bam(FASTQ=>[@bwa_1u_reads], NAME=>"$name.1u", REF=>$Reference, THREADS=>$Threads);
	&fastq2bwa2bam(FASTQ=>[@bwa_2u_reads], NAME=>"$name.2u", REF=>$Reference, THREADS=>$Threads);
}

my $pm=new Parallel::ForkManager($Threads);
foreach my $list2 (@list) {
	$pm->start and next;
	my ($dir, $name) = split /\s+/, $list2;

    system "$samtools merge -h $name.p.bam $name.bam $name.p.bam $name.1u.bam $name.2u.bam && echo merge_ok && rm $name.p.bam $name.1u.bam $name.2u.bam";
    system "$samtools fixmate  $name.bam $name\_fixmate.bam  && echo fixmate_ok && rm $name.bam";
    system "$samtools sort -@ $cpu -O bam -o $name.bam $name\_fixmate.bam && echo fixmate_ok_sort_ok && rm $name\_fixmate.bam";
	system "$samtools rmdup $name.bam $name.rmdup.bam && echo $name\_rmdup_ok && mv $name.rmdup.bam $name.bam";
	system "$samtools index $name.bam && echo $name\_index_OK";
	
	&fixBamGap(BAM=>"$name.bam", REF=>$Reference, KEEP_INTERVALS=>"NO");
	
	system "$samtools index $name.bam && echo $name\_index_OK";
	$pm->finish;
}
$pm->wait_all_children;


##### Subroutines below.

sub trimFastq {
	my %args = (
		PAIR    => "PE",
		THREADS => 1,
		PHRED   => "phred33",
		FASTQ   => undef,  ## Fastq names should be stored in an array.
		NAME    => undef,
		SLIDINGWINDOW => "4:20",
		LEADING => 20, 
		TRAILING => 20, 
		MINLEN  => 30, 
		OTHER => "TOPHRED33",
		@_,         # actual args override defaults
	);
	if (!defined($args{FASTQ}) || !defined ($args{NAME})){
		die "ERROR in function trimFastq: Fastq array and name are required!\n";
	}
	my $outReads;
	if ($args{PAIR} eq "PE"){$outReads = "$args{NAME}.1p.fq $args{NAME}.1u.fq $args{NAME}.2p.fq $args{NAME}.2u.fq";}
	if ($args{PAIR} eq "SE"){$outReads = "$args{NAME}.fq";}
	my $trimmomatic="/usr/biobin/trimmomatic-0.36.jar";
	system "java -jar $trimmomatic $args{PAIR} -threads $args{THREADS} -$args{PHRED} @{$args{FASTQ}}  $outReads SLIDINGWINDOW:$args{SLIDINGWINDOW} LEADING:$args{LEADING} TRAILING:$args{TRAILING} MINLEN:$args{MINLEN} $args{OTHER}";
	print "$args{NAME}_fastq_trimming_ok\n";
}

sub fastq2bwa2bam {
	my %args = (
		THREADS => 1,
		FASTQ   => undef, ## Fastq names should be stored in an array.
		REF    => undef,
		NAME => undef,
		OTHER => "-M -v 1",
		@_,         # actual args override defaults
	);
	
	if (!defined($args{FASTQ}) || !defined ($args{NAME}) || !defined ($args{REF})) {
		die "ERROR in function fastq2bwa2bam: reference, fastq array and name are required!\n";
	}
	
	my $bwa="bwa";
	my $samtools="samtools";
	system "$bwa  mem -t $args{THREADS}  $args{OTHER} -R '\@RG\tID:$args{NAME}\tSM:$args{NAME}\tLB:$args{NAME}'  $args{REF} @{$args{FASTQ}} | $samtools view -bS -  > $args{NAME}.bam";
	#print "\n";
	print "$args{NAME}_fastq_bwa2bam_ok\n";
}

sub fixBamGap {
	my %args = (
		MEMORY => "-Xmx8g",
		BAM   => undef, 
		REF    => undef,
		KEEP_INTERVALS => "NO",
		@_,         # actual args override defaults
	);
	
	if (!defined($args{BAM}) || !defined ($args{REF})) {
		die "ERROR in function fixBamGap: reference, bam are required!\n";
	}
	
	my $gatk="/usr/biobin/GenomeAnalysisTK-v3.5.jar";

	system "java $args{MEMORY} -jar $gatk -T RealignerTargetCreator -R $args{REF} -I $args{BAM} -o $args{BAM}.intervals";
	return "$args{BAM}.intervals_ok";
	system "java $args{MEMORY} -jar $gatk -T IndelRealigner  -R $args{REF} -I $args{BAM} -targetIntervals $args{BAM}.intervals -o $args{BAM}.gatk.bam";
	return "$args{BAM}_gatk_ok";
	system "mv $args{BAM}.gatk.bam $args{BAM}";
	if ($args{KEEP_INTERVALS} =~ /NO/i) {
		system "rm $args{BAM}.intervals";
	}
}
#     



