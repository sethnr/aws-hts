#!/usr/bin/perl

use Getopt::Long;
my $vcf;
my $flank=0;
my $output="FASTA";
my $minQ = 0;
my $onlyVars;
GetOptions(
	   'vcf=s' => \$vcf,
	   'header=s' => \$header,
	   'fasta=s' => \$fasta,
	   'flank=s' => \$flank,
	   'output=s' => \$output,
	   'minQ=s' => \$minQ,     # min quality to add variant to file
	   'only_vars|vars=s' => \$onlyVars
          );

$output = uc($output);

my %onlyVars;
if($onlyVars) {
open(VARS,"$onlyVars");
  while(<VARS>) {  
    chomp;
    $onlyVars{$_} = 1;
#    print STDERR $_ ."->". $onlyVars{$_}."\n"; 
  }
}

open(HEAD,$header);
while(<HEAD>) {
$head_str .= $_;
}
close(HEAD);

my $new_fa = "get_gene_consensus.tmp.fa";
open(NEWFA,">$new_fa");
system("samtools faidx $fasta");

my $new_vcf = "get_gene_consensus.tmp.vcf";
open(NEWVCF,">$new_vcf");
print NEWVCF $head_str;

my %SNP_codes = (
	       "[A/C]", "M",
	       "[A/G]", "R",
	       "[A/T]", "W",
	       "[C/G]", "S",
	       "[C/T]", "Y",
	       "[G/T]", "K",

	       "[C/A]", "M",
	       "[G/A]", "R",
	       "[T/A]", "W",
	       "[G/C]", "S",
	       "[T/C]", "Y",
	       "[T/G]", "K");



while(<>) {
  chomp;
  #print STDERR $_."\n";
  my ($gene, $chr, $start, $end, $set, $priority, $notes) = split("\t",$_);
  
  my $region_vcf = `tabix $vcf -r $chr:$start-$end`;
  my @region_vcf = split("\n",$region_vcf);

  my %uniq_line;
  foreach $line (@region_vcf) {
    #print STDERR $line."\n";
    my @F = split("\t",$line);
    #my $snp_id = $F[0]."_".$F[1]."_".$F[3]."/".$F[4];

    $F[0] = $gene;
    my $global_start = $F[1];
    $F[1] = $F[1] - $start +$flank+1;

    #uniquify:
    $uniq_line{$global_start} = \@F;

  }
#  print NEWFA ">".$gene."\t".$chr.":".$start."-".$end."\tfl:".$flank."\n";
  $start = $start - $flank;
  $end = $end + $flank;
  my $region_fa = `samtools faidx $fasta $chr:$start-$end`;
  $region_fa =~ s/\>/\>$gene\t/gi;
  @region_fa = split("\n",$region_fa);
  $fa_head = shift(@region_fa);
  $region_seq = join("",@region_fa);
  my @region_fa = split(//,$region_seq);
  print NEWFA $fa_head."\n";
  print NEWFA $region_fa."\n";
  

  my @keys = keys %uniq_line;
  my @snp_ids = ();
#  print STDERR $gene."\t".$#keys."\t".$#region_vcf."\n";
  foreach $key (sort {$a <=> $b} @keys) {
#    print STDERR $key."\n";
    my @F = @{$uniq_line{$key}};
    
    next if $F[5]+0 < $minQ;
    $snp_id =  $chr."_".$key."_".$F[3]."/".$F[4];
    
    my $newVar = "!";
    if($onlyVars) {
      if($onlyVars{$snp_id} == 1) {
	$newVar = "[".$F[3]."/".$F[4]."]";
	push (@snp_ids, $snp_id);
      }
      else {
	$newVar = $SNP_codes{"[".$F[3]."/".$F[4]."]"};
      }
    }
    
    my $idx = $F[1] -1;
    print NEWVCF join("\t",@F)."\n";

    if(uc($region_fa[$idx]) eq uc($F[3])) {
      $region_fa[$idx] = $newVar;
    }
    else {
      my $fs = $F[7];
      $fs =~ s/^.*FS/FS/g;
      print STDERR $region_fa[$idx] ." != ". $F[3]."\n".$fs."\n";
      $region_fa[$idx] = "[!!!]";
      
    }
  }

  if($output eq "FASTA") {
    print $fa_head."\n";
    for ($i=0; $i <= $#region_fa; $i++) {
      print $region_fa[$i];
      print "\n" if (($i+1) % 60 == 0);
    }
    print "\n";
  }
  elsif ($output eq "SEQUENOM") {
    next if $#snp_ids < 0;
    $keys = join(", ",@snp_ids);
    print  $gene."\t".
           "\t".
           $keys."\t".
	   $set."\t".
	   $priority."\t".
	   $notes."\t".
	   join("",@region_fa)."\n";
  }

  

}
close(NEWVCF);
close(NEWFA);
print STDERR "bgzip $new_vcf\n";
#system("bgzip $new_vcf");
print STDERR "tabix -pvcf $new_vcf.gz\n";
#system("tabix -pvcf $new_vcf.gz");
print STDERR "cat $new_fa | vcf-consensus $new_vcf.gz > $new_fa.consensus.fa\n";
#system("cat $new_fa | vcf-consensus $new_vcf.gz > $new_fa.consensus.fa");



