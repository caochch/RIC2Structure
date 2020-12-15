#!/usr/bin/perl
die "perl $0 all_reads resolution\n" if(@ARGV != 2);
my $all_pairs=shift;
my $resolution=shift;

my $window=5000;

my $factor=$window/$resolution;

open(AP,$all_pairs) || die;
while(my $line=<AP>){
	chomp $line;
	my @sub=split/\s+/,$line;
	
	$sub[3]=$sub[3]*$factor;
	$sub[7]=$sub[7]*$factor;
	my $new=join"\t",@sub;
	print $new,"\n";
}

open(HS,">hs.resolution$resolution.size") || die;
print HS "hs1\t";
print HS 29903*$factor,"\n";
close HS;
