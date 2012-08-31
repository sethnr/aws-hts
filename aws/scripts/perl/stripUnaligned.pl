while(<>) {
  @F = split; 
  print $_ unless $F[2] eq '*';
}
