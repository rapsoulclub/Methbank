#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;

my $Program_name = `basename $0`;
chomp $Program_name;
my $usage=<<EOF;
Usage: $Program_name -n name

Options: 
	-n    name (0-8)
	-h    help

Example:
	$Program_name -n name
EOF

my ($name,$help);

GetOptions('name=i'	=> \$name,
			'help'	=> \$help);

die $usage if($help);
die $usage if(!(defined $name));

#chr length
open(F,"/leofs/zhangz_group/sunshx/snp-call-wangl/mouse/chromInfo.txt");
my %chr_length;#chr1-19,chrX,chrY
while (my $content=<F>) {
	my @temp=split(/\t/,$content);
	$chr_length{$temp[0]}=$temp[1];
}
close(F);

#methy level
my @methylevel_name=qw(2cell 4cell e13_5f e13_5m e6_5 e7_5 icm oocyte sperm);#(1 2 3 4 5 6 7 8 9)
my $methylevel={};
for (my $i=0; $i<@methylevel_name; $i++) {
	next unless ($i==$name);
	open(F,"/leofs/zhangz_group/sunshx/snp-call-wangl/mouse/$methylevel_name[$i]"."_mc.bed");
	my $chrom;
	my $k=$i+1;
	while (my $content=<F>) {
		chomp $content;
		my @word=split(/\t/,$content);
		next if (($word[4]+$word[5])<5);
		$methylevel->{$word[0]}->{$k}->{$word[2]}=$word[3];#$methylevel{chr}{methylevel_name}{site}=level
	}
	close(F);
}

#reflat, by gene, sort and store
open(F,"/leofs/zhangz_group/sunshx/snp-call-wangl/mouse/refFlat_new.txt");
open(R1,">/leofs/zhangz_group/sunshx/snp-call-wangl/mouse/split/gene_methylevel_$name");
while (my $content=<F>) {
	my $promoter_start=0;
	my $promoter_end=0;
	my @methy_gene_data=qw(0 0 0 0 0 0 0 0 0 0);
	my @methy_promoter_data=qw(0 0 0 0 0 0 0 0 0 0);
	my $NA=1;
	chomp $content;
	my @word=split(/\t/,$content);
	if (exists $chr_length{$word[2]}) {
		if ($word[3]=~m,\+,) {
			if ($word[4]>2000) {
				$promoter_start=$word[4]-2000;
			}
			else {
				$promoter_start=1;
			}
			$promoter_end=$word[4]-1;
		}
		if ($word[3]=~m,\-,) {
			if ($chr_length{$word[2]}>($word[5]+2000)) {
				$promoter_end=$word[5]+2000;
			}
			else {
				$promoter_end=$chr_length{$word[2]};
			}
			$promoter_start=$word[5]+1;
		}
	}
	else {
		$NA=0;
	}
	#$methylevel{chr}{methylevel_name}{site}=level
	if (exists $methylevel->{$word[2]}) {
		foreach my $number (sort {$a<=>$b} keys %{$methylevel->{$word[2]}}) {
			my $methylevel_gene_number=0;
			my $methylevel_promoter_number=0;
			foreach my $site (sort {$a<=>$b} keys %{$methylevel->{$word[2]}->{$number}}) {
				if ($site>=$word[4] && $site<=$word[5]) {
					$methylevel_gene_number++;
					my $temp=$methy_gene_data[$number];
					my $length=$methylevel->{$word[2]}->{$number}->{$site};
					$methy_gene_data[$number]=$temp+$length;
				}
				if ($NA==1) {
					if ($site>=$promoter_start && $site<=$promoter_end) {
						$methylevel_promoter_number++;
						my $temp=$methy_promoter_data[$number];
						my $length=$methylevel->{$word[2]}->{$number}->{$site};
						$methy_promoter_data[$number]=$temp+$length;
					}
				}
			}
			if ($methy_gene_data[$number]==0) {
				$methylevel_gene_number=1;
			}
			if ($methy_promoter_data[$number]==0) {
				$methylevel_promoter_number=1;
			}
			my $temp=$methy_gene_data[$number]/$methylevel_gene_number;
			$methy_gene_data[$number]=$temp;
			if ($NA==1) {
				$temp=$methy_promoter_data[$number]/$methylevel_promoter_number;
				$methy_promoter_data[$number]=$temp;
			}
			else {
				$methy_promoter_data[$number]="-";
			}
		}
	}
	else {
		for (my $i=0; $i<@methy_gene_data; $i++) {
			$methy_gene_data[$i]="-";
			$methy_promoter_data[$i]="-";
		}
	}
	my $print_i=$name+1;
	if ($methy_gene_data[$print_i]=~m,\-,) {
		print R1 "$methy_gene_data[$print_i]\t";
	}
	else {
		my $print_gene=sprintf("%.4f",$methy_gene_data[$print_i]);
		print R1 "$print_gene\t";
	}
	if ($methy_promoter_data[$print_i]=~m,\-,) {
		print R1 "$methy_promoter_data[$print_i]\n";
	}
	else {
		my $print_promoter=sprintf("%.4f",$methy_promoter_data[$print_i]);
		print R1 "$print_promoter\n";
	}
}
close(R1);
close(F);


