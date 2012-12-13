#!/usr/bin/perl

%chrs = ("X", "NC_004818.2",
"2R", "NT_078266.2",
"3R", "NT_078268.4",
"2L", "NT_078265.2",
"3L", "NT_078267.5");

my $format = "%s\t%d\t%d\t%s";
my %tracks;
my @tracks;

$step = 1000;
$half_step = $step/2;

sub _parse_track_names {
  @tracks = @_;
  for ($i=0; $i<=$#tracks; $i++) {    
    $tracks{$tracks[$i]} .=  "track type=bedGraph name=".$tracks[$i]."\n";
  }
}

if (@ARGV) {
_parse_track_names(@ARGVv);
}


while(<>) {
  my @F = split;
  if ($F[0] eq 'CHR') {_parse_track_names(@F[3..$#F]);
		       next;}
  for ($i=0; $i<=$#tracks; $i++) {    
    $middle = ($F[2] - $F[1])/2+$F[1];
    $tracks{$tracks[$i]} .=  sprintf($format, $chrs{$F[0]},($middle - $half_step),($middle+$half_step),$F[$i+3]);
    $tracks{$tracks[$i]} .=  "\n";
  }
}

foreach my $track (@tracks) {
  print STDOUT $tracks{$track}."\n";
  


}
