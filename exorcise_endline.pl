#!/usr/bin/env perl
use strict;

my @Seq = ();
my $contig;
my $first = 0;
my $head = '';

while(<STDIN>)
{
 if($_ =~ /^>/)
  {
      my $head = substr($_, 0, length($_) - 2);
      my $headN = substr($_, length($_) - 2, 1);
      $head =~ s/\(.+\)//g;
      $head =~ s/>.+_contig/>contig/gi;
      print $head, "\n";
  }
 elsif($_ =~ /^[ACGTN]/)
  {
      my $line = substr($_, 0, length($_) - 2);
      print $line, "\n";
  }

}

exit;
