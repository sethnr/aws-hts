#!/usr/bin/perl

my $bed = shift;
my $flank = shift;
$flank = 0 unless($flank);


open(BED,$bed);

my @locs;
while(<BED>) {

  my ($chr,$start,$end,@undef) = split;
  my @loc = ($chr,($start-$flank),($end+$flank));
  push(@locs,\@loc);
}
close(BED);

my ($chr, $pos, $step, $span, $inblock);
while(<>) {
  if (m/^track/i) { print $_;}
  elsif(m/^(.*step) chrom=(.*) start=(\d+) step=(\d+) span=(\d+)/gi) {
   ($steptype,$chr, $pos, $step, $span) = ($1, $2, $3, $4,$5);
   $inblock = undef;
#   print STDERR $chr.":\t:".$pos.":\t:".$step.":\t:".$span."\n";

  }
  elsif(m/^(.*step) chrom=(.*) span=(\d+)/gi) {
   ($steptype,$chr, $span) = ($1, $2, $3);
   $inblock = undef;
#   print STDERR $chr.":\t:".$pos.":\t:".$step.":\t:".$span."\n";

  }
  else {
    if ($_ =~ m/^(\d+)\s+(\d+)/gi) {
      $pos = $1;
    }

    $start = $pos -($span/2);
    $end = $pos + ($span/2);

    my $found = undef;
    foreach $loc (@locs) {
      my ($bed_chr,$bed_start,$bed_end) = @{$loc};
      if($bed_chr eq $chr) {
	#	print STDERR "$chr";
	if(($start < $bed_end) && ($end > $bed_start)) {
	  unless ($inblock) {
	    print "$steptype chrom=$chr start=$pos step=$step span=$span\n" if ($steptype eq "fixedStep"); 
	    print "$steptype chrom=$chr span=$span\n" if ($steptype eq "variableStep"); 
	  }
#	  print "$bed_chr $bed_start $bed_end\n" unless $inblock;
	  print $_;
	  $inblock = 1;
	  $found = 1;
	  last;
	}
      }
    }
    $inblock = undef unless $found;
    $pos = $pos+$step;
#    print STDERR $pos." "
  }
}
  
