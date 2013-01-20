#!/usr/bin/perl

%chrs = ("X", "NC_004818.2",
"2R", "NT_078266.2",
"3R", "NT_078268.4",
"2L", "NT_078265.2",
"3L", "NT_078267.5",
"UNKN","NC_002084.1");

$ncbi_chrs = 0;

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
    
#    $step = ($end - $start +1);
    $step = 1000;
    $last_chr = $chr;
    for ($i=0; $i<=$#tracks; $i++) {    
      # shiftfixedStep chrom=chr3 start=400601 step=100
      $tracks{$tracks[$i]} .=  "fixedStep chrom=".$chrs{$chr}." start=".($middle - ($step/2))." step=".$step." span=".$step."\n" if $ncbi_chrs;
      $tracks{$tracks[$i]} .=  "fixedStep chrom=".$chr." start=".($middle - ($step/2))." step=".$step." span=".$step."\n" unless $ncbi_chrs;
      
      }
  }



sub _end_block {
  foreach my $track (@tracks) {
    open (TRACK,">".$ARGV.".".$track.".wig");
    print TRACK $tracks{$track}."\n";
    close(TRACK);
  }
}


my $last_pos; my $last_chr; my $step;
while(<>) {
  my @F = split;
  if (uc($F[0]) eq 'CHR') {@tracks = @F[3..$#F];
		       next;}

  if (!$last_chr) {_start_track(); 
		   _start_block(@F[0..2]);}
  elsif (($F[1] > $last_pos + $step) || ($last_chr ne $F[0])) {_start_block(@F[0..2]);}
  $last_chr = $F[0]; $last_pos = $F[2];
  for ($i=0; $i<=$#tracks; $i++) {    
    $val = $F[$i+3];
    $val = 0 unless $val =~ m/\d+/gi;
    $tracks{$tracks[$i]} .=  $val."\n";
  }
}

_end_block();
