use strict;
use warnings;
while (<>){
	if (/^#/) {next;}
	chomp;
	my $out = &vcfParser(VCF=>$_, QUAL=>30, ALT1HOM =>"YES", MAPQ=>20, MINDP=>2, MAXDP=>50);
	if ($out =~ /DP/) {
		print "$_\n";
	}
}


sub vcfParser {

	my %args = (
		VCF    => undef,
		QUAL => 20,
		MAPQ   => 20,
		MINDP  => undef,
		MAXDP  => undef,
		DP4AF => undef,
		ALT1HOM => undef,
		RMINDEL => undef, ## rm indel
		RMTRI => undef,  ## rm more than one alt.
		@_,         # actual args override defaults
	);
	
	if (!defined($args{VCF})) {
		die "ERROR in function vcfFilter: VCF line is required!\n";
	}
	
	my ($chr, $pos, $id, $ref, $alt, $qual, $filter, $info, $format, @genotypes) = split /\t/, $args{VCF};
	if ($ref eq "N"){return 0; last;}
	if (defined $args{RMINDEL}) {
		if ($args{VCF} =~ /INDEL/) {return 0; last;}
	}
	
	if (defined $args{RMTRI}) {
		if ($alt =~ /,/) {return 0; last;}
	}
	
    if ($qual<$args{QUAL}) {return 0; last;}
    if ($args{VCF} =~ /DP=(\d+)/){
        if (defined $args{MINDP}) {if ($1<$args{MINDP}) {return 0; last;}}
        if (defined $args{MAXDP}) {if ($1>$args{MAXDP}) {return 0; last;}}
    }
	if (defined $args{ALT1HOM}) {if (!($genotypes[0] =~ /1\/1/)) {return 0; last;}}

	if (defined $args{DP4AF}) {	
		if ($args{VCF} =~ /DP4=(\d+),(\d+),(\d+),(\d+)/){
			if ($3+$4+$1+$2>=4 && $3>0 && $4>0 ){
				my $af=1-($1+$2)/($3+$4+$1+$2);
				if ($af < $args{DP4AF}) { 
						return 0; last;
				}
			}
		}
	}
	return $args{VCF};
}