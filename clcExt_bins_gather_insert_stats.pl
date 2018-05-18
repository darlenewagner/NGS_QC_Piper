#!/usr/bin/perl -w
use strict;
use Math::Round;

#INPUT: sam_file = my_contigs.sam
#OUTPUT: tab-delimited counts of Insert lengths, suitable for plotting a histogram 

my $sam_file = $ARGV[0];
my @SAM = ();
my @Depth = ();

open(SAM, $sam_file) || die ("Cannot find .sam file, $sam_file $!");

my $match_sum = 0;
my $count = 0;
my @inserts = ();
my $i = 0;

Read_Stats(\*SAM);

exit;


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

     my $m = 0;
    my $m_sum = 0;
    

   my @sInsert = sort {$a <=> $b} @Inserts;  ## force ascending order
   my $sum_insert = 0;
   my $half = round( scalar @sInsert / 2 );
   my $quarter = round( scalar @sInsert / 4 );
   my $total = scalar @sInsert;
   my $min = $sInsert[0];
   my $max = $sInsert[$total - 1];
   my $ii = 0; 
 
  foreach(@sInsert)
   {
       $sInsert[$ii] =~ s/\s+//g;
       $sum_insert = $sum_insert + $sInsert[$ii];
       $ii++;
    }
    $avg_insert = $sum_insert / $ii;
    my $str = sprintf("%.2f", $avg_insert);
  if($sum_insert > 0)
   {
    print "Avg.\t";
    printf $str, $avg_insert;
    print "\n";
    print "Lower Quartile\t", $sInsert[$quarter], "\nMedian\t", $sInsert[$half], "\nMinimum\t", $min, "\nMaximum\t", $max, "\n";
    
   
    if( ($max - $min) > 5 )
      { 
         my $std_dev = Standard_Deviation(\@sInsert, $avg_insert);
        
         my $outlier_limit = 3.5 * $std_dev;  ## only remove 'insane' outliers

        ## subroutine assumes input array, @sInsert, is sorted ascending
        Print_Histogram(\@sInsert, $outlier_limit, $avg_insert);  
      }
    else
    {
      print "Insufficient variation in inserts or low number of inserts!\n"; 
    }
     
   }
   else
    {
     print "No-PE!\t";
    }   
 
}



sub Standard_Deviation {
    my @input_ref = @{$_[0]};
    my $avg = $_[1];
    my $ii = 0;
    my $sum_sq = 0;
    my $std_dev = 0;
    
     foreach(@input_ref)
       {
         $sum_sq = $sum_sq + ($input_ref[$ii] - $avg)*($input_ref[$ii] - $avg);
         $ii++;
       }
    $std_dev = sqrt( $sum_sq / $ii );
    
   # my $str = sprintf("%.2f\n", $std_dev);
   # printf $str, $std_dev;    
    
   return $std_dev;  
  }


sub Print_Histogram {
    my @input_ref = @{$_[0]};
    my $outlier_filter = $_[1];
    my $avg = $_[2];
    
    my @Hist_Bins = (0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130, 140, 150, 160, 170, 180, 190, 200, 210, 220, 230, 240, 250, 260, 270, 280, 290, 300, 310, 320, 330, 340, 350, 360, 370, 380, 390, 400, 410, 420, 430, 440, 450, 460, 470, 480, 490, 500, 510, 520, 530, 540, 550, 560, 570, 580, 590, 600, 610, 620, 630, 640, 650, 660, 670, 680, 690, 700, 710, 720, 730, 740, 750, 760, 770, 780, 790, 800, 810, 820, 830, 840, 850, 860, 870, 880, 890, 900, 910, 920, 930, 940, 950, 960, 970, 980, 990, 1000, 1010, 1020, 1030, 1040, 1050, 1060, 1070, 1080, 1090, 1100, 1110, 1120, 1130, 1140, 1150, 1160, 1170, 1180, 1190, 1200);  


    my $bins = scalar @Hist_Bins - 1; 
   
    my $min = $Hist_Bins[0];
    my $length = scalar @input_ref;
    my $max = $Hist_Bins[$bins - 1];
  # print $min, "\t", $max, "\n";    
    print "Counts per insert size range:\n"; 
  
   # adjust maximum to prevent construting bins
   # based on extreme large inserts
   #  if($max > ($avg + $outlier_filter))
   #   {
   #	$max = $avg + $outlier_filter;
   #   }
      
   #  my $range = $max - $min;
   
   # print $offset, "\n";
    
    my @Hist_Counts = ();
    my $lower = $Hist_Bins[0];
    my $upper = $Hist_Bins[1];
    my $offset = $Hist_Bins[1] - $Hist_Bins[0];

 for(my $h = 0; $h < $bins; $h++)
  {

    #  print $Hist_Bins[$bins], "\t";

   if($h < $bins)
      {
       	$lower = $Hist_Bins[$h];
        $upper = $Hist_Bins[$h + 1];
      }

     $Hist_Counts[$h] = 0;
    # $upper = $upper + $offset;

   for(my $i = 0; $i < $length; $i++)
   {
      if(($lower <= $input_ref[$i]) && ($input_ref[$i] < $upper))
      {
	  $Hist_Counts[$h]++;
      }
    # elsif($Hist_Bins[$bins] <= $input_ref[$i])
    #  {
    #   	  $Hist_Counts[$h]++;
    #      $lower = $Hist_Bins[$bins];
    #  } 
   }   
  
   if($h < ($bins - 1))
   {
      print $lower, "-", $upper - 1, "\t", $Hist_Counts[$h], "\n";   
   }
  elsif($h == ($bins - 1))
  {
    print $lower, "<\t";
    print $Hist_Counts[$h], "\n";   
  }
   
           

  }
   print "\n";
   
}
