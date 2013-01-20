#!/usr/bin/perl
use Getopt::Long;

%chrs = ("X", "NC_004818.2",
"2R", "NT_078266.2",
"3R", "NT_078268.4",
"2L", "NT_078265.2",
"3L", "NT_078267.5",
"UNKN","NC_002084.1");

$ncbi_chrs;
$fixed;
$fixStep=1000;
GetOptions(
	   'ncbi' => \$ncbi_chrs,
	   'fixed' => \$fixed,
	   'step=s' => \$fixStep
          );


my $last_pos; my $last_chr; my $step;

my $format = "%s\t%d\t%d\t%s";
my %tracks;
my @tracks;

sub _parse_track_names {
  @tracks = @_;
  for ($i=0; $i<=$#tracks; $i++) {    
    # fixedStep chrom=chr3 start=400601 step=100

    $tracks{$tracks[$i]} .=  "track type=bedGraph name=".$tracks[$i]."\n";
  }
}

sub _start_track {
    for ($i=0; $i<=$#tracks; $i++) {    
      # shiftfixedStep chrom=chr3 start=400601 step=100
      $tracks{$tracks[$i]} .=  "track type=wiggle_0 name=".$tracks[$i]."\n";
      }
  }
sub _start_block {
    my ($chr, $start, $end) =  @_;
    my $middle = (($end - $start)/2)+$start;
    
    $step = ($end - $start) unless $fixed;
    $step = $fixStep if $fixed;
    $last_chr = $chr;

    $chr = $chrs{$chr} if $ncbi_chrs;
       
    print STDERR "starting block $chr : $start\n";
    for ($i=0; $i<=$#tracks; $i++) {    
      # shiftfixedStep chrom=chr3 start=400601 step=100
	$tracks{$tracks[$i]} .=  "fixedStep chrom=".$chr." start=".($middle - ($step/2))." step=".$step." span=".$step."\n" if $fixed;
	$tracks{$tracks[$i]} .=  "variableStep chrom=".$chr." span=".$step."\n" unless $fixed;
    }
}
  



sub _end_block {
  foreach my $track (@tracks) {
#    open (TRACK,">".$ARGV.".".$track.".wig");
    open (TRACK,">".$track.".wig");
    print TRACK $tracks{$track}."\n";
    close(TRACK);
  }
}


while(<>) {
  my @F = split;
  if (uc($F[0]) eq 'CHR') {@tracks = @F[3..$#F];
		       next;}

  if (!$last_chr) {
    print STDERR "starting new chr $F[0] \t";
    _start_track(); 
    _start_block(@F[0..2]); }
  elsif ($last_chr ne $F[0]) {
    print STDERR "CHRs don't match $last_chr -> $F[0]\t";
    _start_block(@F[0..2]); }
  elsif ($fixed && ($last_pos + $step < $F[1])) {
    print STDERR "fixed, jump too large $last_pos + $step > $F[1] \t";
    _start_block(@F[0..2]); }
  elsif ((!$fixed) && ($last_pos - $F[1] != $step)) {
    print STDERR "var, jump too large\t";
    _start_block(@F[0..2]);}

  $last_chr = $F[0]; $last_pos = $F[1];
  for ($i=0; $i<=$#tracks; $i++) {    
    $val = $F[$i+3];
    $val = 0 unless $val =~ m/\d+/gi;
    $tracks{$tracks[$i]} .=  $val."\n" if $fixed;
    $tracks{$tracks[$i]} .=  $F[1]."\t".$val."\n" unless $fixed;
  }
}

_end_block();
