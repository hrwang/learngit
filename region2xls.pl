while (<>){
	chomp;
	my @tmp = split;
	if (/HIGH|MODERATE/){
		print $_;
		my @annts = split /,/, $tmp[8];
		foreach my $annt (@annts){
			#print "$annt\n";
			if ($annt=~/HIGH|MODERATE/){
				#print "$annt\n";
				if ($annt =~ /(Os[0-9][0-9]g\d+)/){
					#print "\t$1\n";
					my $info = `grep $1 /mnt/mouse/genomes/osativa/japonica/gff/locus_msu7.gene.gff`;
					#my $info = `grep $1 /mnt/mouse/genomes/osativa/japonica/gff/locus_rapdb_20140625.gff`;
					chomp ($info);
					#$info =~ s/20/ /gi;
					#$info =~ s/2C/ /gi;
					my @info = split /\t/, $info;
					print "\t$info[1]\t";
					$info[8] =~ s/%20/ /gi;
					$info[8] =~ s/%2C/ /gi;
					$info[8] =~ s/%2F/ /gi;
					$info[8] =~ s/\s+/ /gi;
					@anntinfo = split /;/, $info[8];
					foreach my $anntinfo (@anntinfo){
						if ($anntinfo =~ /ID=/){next;}
						if ($anntinfo =~ /Name=/){next;}
						if ($anntinfo =~ /Transcript variants=/){next;}
						print "$anntinfo ";
					}
					print "\n";
					last;
				}
			}
		}
	} else {
		print "$_\tNA\tNA\n"
	}

}
