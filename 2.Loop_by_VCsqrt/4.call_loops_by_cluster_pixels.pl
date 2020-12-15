#!/usr/bin/perl
die "perl $0 enriched_pixels_5000.resolution5.bedpe resolution\n" if(@ARGV != 2);
my $all_pixels_bedpe=shift;
my $resolution=shift;

my $slop=int(2*$resolution);

`bedtools pairtopair -a $all_pixels_bedpe -b $all_pixels_bedpe -is -rdn -slop $slop > tmp.pixels_to_pixels.list`;
`perl clustering_bridge_SaveMemory_fast.pl tmp.pixels_to_pixels.list tmp.pixels.cluster tmp.ID_for_pixels_in_cluster.list`;


print "#chr1\tx1\tx2\tchr2\ty1\ty2\tname\tscore\tstrand1\tstrand2\n";

my $manualClusterID;
open(TC,"tmp.pixels.cluster") || die;
while(my $line=<TC>){
	chomp $line;
	my @sub=split/\s+/,$line;
	my $new=join"\t",@sub[0..5];
	$manualClusterID++;
	print $new,"\tManualLoops_$manualClusterID\t255\t+\t+\n";
}

my %pixels_in_manual;
open(TI,"tmp.ID_for_pixels_in_cluster.list") || die;
while(my $line=<TI>){
	chomp $line;
	$pixels_in_manual{$line}=1;
}

my $lonelyID;
open(APB,$all_pixels_bedpe) || die;
<APB>;
while(my $line=<APB>){
	chomp $line;
	my @sub=split/\s+/,$line;
	if($pixels_in_manual{$sub[6]}){
		next;
	}
	else{
		print $line,"\n";
	}
}

`rm -rf tmp.pixels_to_pixels.list tmp.pixels.cluster tmp.ID_for_pixels_in_cluster.list`;
