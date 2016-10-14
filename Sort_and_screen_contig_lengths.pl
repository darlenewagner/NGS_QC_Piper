use strict;

my %Contig = ();

my @len = ();
my @header = ();
my @Seq = ();


my $r = 0;

while(<STDIN>)
{  

  if($_ =~ /^>/)
   {
     chomp;
     push @header, $_;
     $r++;     
   }
  elsif($_ =~ /^(A|T|G|C|U|N)/)
   {
     chomp;
     $Seq[$r - 1] = $Seq[$r - 1].$_;
   }

}

my $sum = 0;
my $i = 0;

foreach(@header)
 {
     $len[$i] = length($Seq[$i]);
     
     my $revised_head = '';
     
     if($header[$i] =~ /^>NODE/)
     {
         my $cut = index($header[$i], '_l');
         $revised_head = substr($header[$i], 0, $cut)."_length_".$len[$i];
         $Contig{$revised_head} = $Seq[$i];
         # print $revised_head, "\t", $len[$i], "\n";
     }
    elsif($header[$i] =~ /^>Pair/)
     {
         $revised_head = $header[$i]."_length_".$len[$i];
         $Contig{$revised_head} = $Seq[$i];
         # print $revised_head, "\t", $len[$i], "\n"; 
     }
   elsif($header[$i] =~ /^>.*_?contig.+/i)
     {
         $header[$i] =~ s/contig_/Contig/g;
         
         my $cut = length($header[$i]);
         
         $header[$i] =~ s/\t/ /g;

         if($header[$i] =~ /\s/)
	 {
           $cut = index($header[$i], ' ');
	 }
                 
         $revised_head = substr($header[$i], 0, $cut)."_length_".$len[$i];
         $Contig{$revised_head} = $Seq[$i];
     }
     
     $i++;
 }

## output contigs sorted by sequence length
## suppress printing of contigs < 500 bp in length

foreach my $name (sort { length($Contig{$b}) <=> length($Contig{$a}) } keys %Contig)
{
  if( length($Contig{$name}) > 499 )
  {
     print $name, "\n"; 
     
     my $offset = 0;
     
     my $temp_Seq = $Contig{$name};
     
     while(($offset + 80)  < length($temp_Seq))
     {
	 print substr($temp_Seq, $offset, 80), "\n";
         $offset = $offset + 80;
     }
    print  substr($temp_Seq, $offset, length($temp_Seq) - $offset), "\n";
  }

}



exit;
