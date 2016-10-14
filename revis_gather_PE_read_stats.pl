#!/usr/bin/perl -w
use strict;
use Tie::File;
use Math::Round;

#INPUT: fas_file = my_contigs.fasta sam_file = my_contigs.sam, dep_file = my_contigs_depth_from_bam.txt
#DEPENDENCIES: 
#OUTPUT: Estimated average coverage, Average mapped read length, Average Insert length, Median Insert length, Insert Quartiles. 

my $fas_file = $ARGV[0];
my $sam_file = $ARGV[1];
my $dep_file = $ARGV[2];
my @SAM = ();
my @Depth = ();

my @Path = split(/\//, $fas_file);

print $Path[0], "\t";

open(FAS, $fas_file) || die ("Cannot find fasta file, $fas_file $!");
#my $check_fas = <FAS>;
#print $check_fas;

open(SAM, $sam_file) || die ("Cannot find .sam file, $sam_file $!");
#tie @SAM, 'Tie::File', $sam_file or die "Cannot tie .sam file, $sam_file $!";
#my $check_sam = $SAM[0];
#print $check_sam, "\n";

open(DEP, $dep_file) || die ("Cannot find samtools depth output file, $dep_file $!");
#tie @Depth, 'Tie::File', $dep_file or die "Cannot tie samtools depth file, $dep_file $!";
#my $check_dep = $Depth[0];
#print $check_dep, "\n";

my $match_sum = 0;
my $count = 0;
my @inserts = ();
my $i = 0;

 print  "\tDepth\tRead(Mapped) Length\tMean(Median) Insert\tContigs\tN50\tLongest Contig\ttotal size\n";

Read_Depth(\*DEP);

Read_Stats(\*SAM);

Contig_Stats(\*FAS);

#untie @SAM;
#untie @Depth;  

exit;

sub Read_Depth
{
   my $data_handle = $_[0];  

   my $avg = 0;
   my $sum = 0;
   my $count = 0;
   my $min = 100;
   my $max = 0;
   
 while(<$data_handle>)
   {
     my @data = split(/\t/, $_);
     $sum = $sum + $data[2];
     $count++;
     if($min > $data[2])
      {
	 $min = $data[2];
      }
  
    if($max < $data[2])
      {
        $max = $data[2]; 
      }

    # print $count, "\n";  
   }
    $avg = $sum / $count;
    my $str = sprintf("%.2f", $avg);
    printf $str, $avg;
    print "\t";

    
}

sub Read_Stats
{
    my $stat_handle = $_[0];
   
    my $avg_match = 0;
    my $avg_insert = 0;
    my @Orig_read = ();
    my @Matched_read = ();
    my @Inserts = ();
    my $code = '';
    my $line = 0;

 while(<$stat_handle>)
  {

  if($_ !~ /^@(HD|SQ|PG)/)
    {      

       my @read = split("\t", $_);
       # print $read[5], "\n";
       
      # print $read[0], "\t"; 
      # print $read[1], "\t"; 
      push @Orig_read, length($read[9]);
      
      


       $code = $read[1];
       $line++;

      ## extract CIGAR string        
       my $match_str = $read[5];
       my @R = split('', $match_str);
       
      ## Process CIGAR string to sum matching positions within reads
   
       if(($code !~ /(\b77\b|\b141\b|\b69\b|\b133\b)/) && ($match_str ne '*'))
         {   
             my $seg = 0;
             my $i = 0;
	     my @match_seq = ();

       ## loop to extract match segments
           foreach(@R)
             {
               if($R[$i] =~ /\d/)
	       {
                  $seg = $seg.$R[$i];
               }
              elsif($R[$i] eq 'M')
               { 
		   push @match_seq, $seg;
                   $seg = '';
               }
             else
	       {
                 $seg = ''; 
	       }
              $i++;  
	    }
         
       ## sum lengths of match segments
	     my $sum = 0;
             my $s = 0;
          foreach(@match_seq) 
	    {
		$match_seq[$s] =~ s/^0//g; # get rid of leading zeros
 	       # print $match_seq[$s], "\t";
               $sum = $sum + $match_seq[$s];
               $s++;
	    }
	   push @Matched_read, $sum;
         }
       else
       {
          push @Matched_read, 0;  
       }

           if($read[8] > 0)
	     {
               push @Inserts, abs($read[8]);
             }
         # else
	 # {
         #        push @Inserts, 0;
	 # }
     }

  }
  
    my $o = 0;
    my $o_sum = 0;
   foreach(@Orig_read)
     {
	 $o_sum = $o_sum + $Orig_read[$o];
         $o++;
     }

     my $m = 0;
    my $m_sum = 0;
    
     foreach(@Matched_read)
      {
	  $m_sum = $m_sum + $Matched_read[$m];
          $m++;
      }

    my $avg_orig = $o_sum / $o;
    my $str = sprintf("%.2f", $avg_orig);
    printf $str, $avg_orig;
    print "\t";
    
    $avg_match = $m_sum / $m;
    $str = sprintf("%.2f", $avg_match);
    printf $str, $avg_match;
    print "\t";

   ## CAUTION: Make certain this is a numeric sort   
   my @sInsert = sort {$a <=> $b} @Inserts;
   my $sum_insert = 0;
   my $half = round( scalar @sInsert / 2 );
   my $quarter = round( scalar @sInsert / 4 );
   my $ii = 0; 
 
  foreach(@sInsert)
   {   
       
       $sInsert[$ii] =~ s/\s+//g;
     ###  print $sInsert[$ii], "\n";
       $sum_insert = $sum_insert + $sInsert[$ii];
       $ii++;
    }
    $avg_insert = $sum_insert / $ii;
    $str = sprintf("%.2f", $avg_insert);
  if($sum_insert > 0)
   {
    printf $str, $avg_insert;
    print "\t";
    print $sInsert[$half], "\t"
   }
   else
    {
     print "No-PE!\n";
    }   
 
}

sub Contig_Stats
{
    my $fileHandle = $_[0];
     
    my $total = 0;
    my @lengths = ();
    my $max = 0;
    my $N50 = 0;
    my $count = 0;    
    
    while(<$fileHandle>)
      {
       if($_ =~ /^>/)
         {
             chomp;
	     my $idx = rindex($_, '_');
             push @lengths, substr($_, $idx + 1, length($_) - $idx);
             $count++;
         }
      }
      
      my  $l = 0;
      foreach(@lengths)
      {
          if($lengths[$l] > $max)         
 	     {
		 $max = $lengths[$l];
	     }
          $total = $total + $lengths[$l];
	  $l++;
      }

    ## The contig that resides at (genomeSize / 2), a.k.a., N50
    
    my $halfGenomeSize = $total/2;
    my $currentGenomeLength = 0;    
    
    $l = 0;
    my $sum = 0;  
    
    for(my $i = 0; $i < $count; $i++)
    {
	$currentGenomeLength += $lengths[$i];
        
        if($currentGenomeLength > $halfGenomeSize)
	 {
           $N50 = $lengths[$i];
           last;  
         }
    }
        
    print  $count, "\t", $N50, "\t", $max, "\t", $total, "\n";
}



# Find the N50 of an assembly.  The N50 is the size N such that 50% of the genome is contained in contigs of size N or greater.
# param: fasta filename
# optional param: $$settings
# return int The N50 in bp.
sub N50($;$){
  my($seqs,$settings)=@_;
  my($seqname,@seqs,$numSeqs,$N50);
  # put the seqs into an array. No defline needed in this subroutine, so it can be discarded
  foreach $seqname (keys %$seqs){
    push(@seqs,$$seqs{$seqname});
  }
  $numSeqs=@seqs;
  # order the contigs by size
  # Largest contigs are towards $seqs[0], and smallest contigs are at the last elements
  @seqs=sort{
    length($a)<=>length($b);
  } @seqs;
  # if the size is not provided, assume the size of the assembly is the genome size
  my $genomeSize = $$settings{expectedGenomeLength}||genomeLength($seqs);
  my $halfGenomeSize=$genomeSize/2;
  my $currentGenomeLength=0;
  # find the contig that resides at (size/2)
  for(my $i=0;$i<$numSeqs;$i++){
    my $seqLength=length($seqs[$i]);
    $currentGenomeLength+=$seqLength;
    if($currentGenomeLength>$halfGenomeSize){
      $N50=$seqLength;
      last;
    }
  }
  # return the length of the contig
  $N50||=0.01;
  return $N50;
}


