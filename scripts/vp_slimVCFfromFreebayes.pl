#!/usr/bin/perl
# author: biotang
# version: 1.5
# 2019.10.10
use strict;
use warnings;
use 5.010;
use File::Basename;

# usage
my $script = basename $0;
die "Usage: perl $script raw.fb.vcf slim.fb.vcf\n" if(not @ARGV);

# command line
my $file_vcf = shift @ARGV;
my $file_out = shift @ARGV;

#
open FILE, $file_vcf or die "";
open OUT, ">$file_out" or die "";
while(<FILE>){
	chomp ;
	if($_=~m/^##/){
		print OUT "$_\n" if($_!~m/^##INFO/ and $_!~m/^##FORMAT/);
		print OUT "$_\n" if($_=~m/^##INFO=<ID=DP/);
		print OUT "$_\n" if($_=~m/^##INFO=<ID=RO/);
		print OUT "$_\n" if($_=~m/^##INFO=<ID=AO/);
		print OUT "$_\n" if($_=~m/^##INFO=<ID=CIGAR/);
		print OUT "$_\n" if($_=~m/^##FORMAT=<ID=GT/);
		print OUT "$_\n" if($_=~m/^##FORMAT=<ID=DP/);
		print OUT "$_\n" if($_=~m/^##FORMAT=<ID=RO/);
		print OUT "$_\n" if($_=~m/^##FORMAT=<ID=AO/);
		next;
	}
	if($_=~m/^#/){
		print OUT "$_\n";
		next;
	}
	#
	my ($chr,$loc,$id,$ref,$alt,$qual,$flt,$inf,$tag,@samples) = (split /\t/, $_);
	my $new_inf=&slim_info(\$inf);
	my $new_tag="GT:DP:RO:AO";
	my $new_samples=&slim_tag(\$tag,\@samples);
	print OUT "$chr\t$loc\t$id\t$ref\t$alt\t$qual\t$flt\t$new_inf\t".
	"$new_tag\t".$new_samples."\n";
}
close FILE;
close OUT;

sub slim_info{
	my $inf = shift @_;
	# per-base coverage depth
	my $dp = $1 if($$inf=~m/\bDP=(\S+?)[;\t]/);
	# reference allele obversation count
	my $ro = $1 if($$inf=~m/\bRO=(\S+?)[;\t]/);
	# alternative alleles obversation counts
	my $ao = $1 if($$inf=~m/\bAO=(\S+?)[;\t]/);
	# CIGAR
	my $cigar = $1 if($$inf=~m/\bCIGAR=(\S+?)[;\t]/);
	# type
	#my $type = $1 if($$inf=~m/\bTYPE=(\S+?)[;\t]/);
	#
	my $new="DP=$dp;RO=$ro;AO=$ao;CIGAR=$cigar;";	
	return $new;
}

sub slim_tag{
	my $tag = shift @_;
	my $spl = shift @_;

	my @new = ();
	my @key = split /:/, $$tag;
	foreach my $str (@$spl){
		if($str eq "."){
			push @new, ".";
		}
		else{
			my @value = split /:/, $str;
			my %tag = ();
			for(my $i=0;$i<scalar @key;$i++){
				$tag{$key[$i]}=$value[$i];
			}
			push @new, "$tag{'GT'}:$tag{'DP'}:$tag{'RO'}:$tag{'AO'}";
		}
	}
	my $new_samples=join "\t",@new;
	return $new_samples;
}

