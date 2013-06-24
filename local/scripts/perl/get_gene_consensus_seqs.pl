#!/usr/bin/perl

use Getopt::Long;
my $vcf;
my $flank=0;
my $output="FASTA";
my $minQ = 0;
my $onlyVars;
my $snp_gap=100;
my $tag="";
GetOptions(
	   'vcf=s' => \$vcf,
	   'header=s' => \$header,
	   'tag=s' => \$tag,
	   'fasta=s' => \$fasta,
	   'flank=s' => \$flank,
	   'snp_join|snp_gap=s', => \$snp_gap,
	   'output=s' => \$output,
	   'minQ=s' => \$minQ,     # min quality to add variant to file
	   'only_vars|vars=s' => \$onlyVars
          );


$output = uc($output);


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


my %only_vars;
my %var_blocks;
if($onlyVars) {
open(VARS,"$onlyVars");
  # snp_gap = 100 
  my $lastEnd = 0; my $lastChr = "";
  my $block = "";
  my $block_priority = 1000;
  my $block_gene;
  my %block_notes;
  my $block_start;
  my $block_chr;
  while(<VARS>) {  
    next unless (m/[\w\d]/gi);
    next if (m/^VAR.*$/i);
    
    chomp;
    my ($snp_id, undef, undef, $gene, $priority, $notes, undef) = split("\t", $_);
    $only_vars{$snp_id} = 1;
    my ($chr,$start,$var) = split(/_/, $snp_id); 
    $start = $start+0;
    # IF CLOSE TO LAST SNP
    if((($lastEnd + $snp_gap) >= $start) && ($chr eq $block_chr)) {
      $block = $block."::".$snp_id;
      $block_priority = $priority if ($block_priority > $priority);
      $block_notes{$notes}=1;
      $block_count++;
    } 
    #if far from last SNP, print previous and restart new block
    else {
      if($block) {
	my $block_name = $block;
	$block_name = $block_chr."_".$block_start."-".$lastEnd."_".$block_count unless ($lastEnd == $block_start);
	$block_notes = join(" ",keys %block_notes);
	my @block = ($block_name, $block, $block_chr, $block_start, $lastEnd, $block_gene, $block_priority, $block_notes);

#	print STDERR ".".$block_name." -> ".join("\t", @block)." \n"; 
	$var_blocks{$block_name} = \@block;
      }
      $block_gene = $gene;
#      $block_notes = $notes;
      %block_notes = undef;
      $block_notes{$notes} = 1;
      $block_start = $start;
      $block_chr = $chr;
      $block_count=1;
      $block_priority = $priority;
      
      $block=$snp_id;
    } 
    # print STDERR $snp_id."\n"; 
    $lastEnd = $start;
  }
  #ADD LAST ELEMENT TO ARRAY
  my $block_name = $block;
  $block_name = $block_chr."_".$block_start."-".$lastEnd."_".$block_count unless ($lastEnd == $block_start);
  my @block = ($block_name, $block, $block_chr, $block_start, $lastEnd, $block_gene, $block_priority, $block_notes);
#  print STDERR "\t".$block_name." -> ".join("\t", @block)." \n"; 
  $var_blocks{$block_name} = \@block;

}




#while(<>) {
foreach my $block (keys %var_blocks) {

  chomp;
#  my ($gene, $chr, $start, $end, $set, $priority, $notes) = split("\t",$_);
  my ($block_name, $snp_set, $chr, $start, $end, $gene, $priority, $notes) = @{$var_blocks{$block}};
  # print STDERR $block_name.":\n".join("\t",@{$var_blocks{$block}})."\n";
  my $fstart = ($start+0-$flank);
  my $fend = ($end+0+$flank);

  #GET REGION FROM VCF
  my $region_vcf = `tabix $vcf -r $chr:$fstart-$fend`;
  my @region_vcf = split("\n",$region_vcf);

  #FOR EACH LINE IN VCF, ADJUST START POSITION AND RESAVE
  my %uniq_line;
  foreach $line (@region_vcf) {
    #print STDERR $line."\n";
    my @F = split("\t",$line);
    #my $snp_id = $F[0]."_".$F[1]."_".$F[3]."/".$F[4];

#    $F[0] = $gene;
    $F[0] = $block_name;
    my $global_start = $F[1];
    $F[1] = $F[1] - $fstart +1;

    #uniquify:
    $uniq_line{$global_start} = \@F;
   }
#  print NEWFA ">".$gene."\t".$chr.":".$start."-".$end."\tfl:".$flank."\n";

  # GET SEQUENCE FROM FASTA
  $start = $start - $flank;
  $end = $end + $flank;
  # print STDERR "samtools faidx $fasta $chr:$start-$end\n";
  my $region_fa = `samtools faidx $fasta $chr:$start-$end`;
  
  #REPLACE POSITION IN GENOME WITH BLOCK NAME
#  $region_fa =~ s/\>/\>$gene\t/gi;
  $region_fa =~ s/\>/\>$block_name\t/gi;
  @region_fa = split("\n",$region_fa);
  $fa_head = shift(@region_fa);
  $region_seq = join("",@region_fa);
  #MAKE INTO CHAR ARRAY
  my @region_fa = split(//,$region_seq);
  #print NEWFA $fa_head."!\n";
  print NEWFA $region_fa."\n";
  

  # CYCLE THROUGH (UNIQ) LINES IN VCF
  my @keys = keys %uniq_line;
  my @snp_ids = ();
#  print STDERR $gene."\t".$#keys."\t".$#region_vcf."\n";
  foreach $key (sort {$a <=> $b} @keys) {
#    print STDERR $key."\n";
    my @F = @{$uniq_line{$key}};
    
    next if $F[5]+0 < $minQ;
    $snp_id =  $chr."_".$key."_".$F[3]."/".$F[4];
    my $newVar = "!";
    # CHECK SNP IS IN SET TO BE ASSAYED:
    # IF SO REPLACE WITH [A/B] STYLE MRKER\
    # IF NOT REPLACE WITH N/Y/R STYLE
    if($onlyVars) {
      if($only_vars{$snp_id} == 1) {
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
    print  $block_name."\t".
           "\t".
           $keys."\t".
	   $tag.$gene."\t".
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



