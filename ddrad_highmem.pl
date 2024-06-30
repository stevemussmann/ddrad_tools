#! /usr/bin/perl

use warnings;
use strict;
use Getopt::Std;
use Chart::Gnuplot;
use Data::Dumper;

# kill program and print help if no command line arguments were given
if( scalar( @ARGV ) == 0 ){
    &help;
    die "Exiting program because no command line options were used.\n\n";
}

# take command line arguments
my %opts;
getopts( '1:2:e:f:hs:', \%opts );

# if -h flag is used, or if no command line arguments were specified, kill program and print help
if( $opts{h} ){
  &help;
  die "Exiting program because help flag was used.\n\n";
}

# declare hash to hold restriction enzyme cut sites
my %re = (
		'AciI' => 'C^CGC',
		'AgeI' => 'A^CCGGT',
		'AluI' => 'AG^CT',
		#'ApeKI' => 'G^CWGC',
		#'ApoI' => 'R^AATY',
		'AseI' => 'AT^TAAT',
		'BamHI' => 'G^GATCC',
		'BfaI' => 'C^TAG',
		'BgIII' => 'A^GATCT',
		'BspDI' => 'AT^CGAT',
		#'BstYI' => 'R^GATCY',
		'ClaI' => 'AT^CGAT',
		#'DdeI' => 'C^TNAG',
		'DpnII' => '^GATC',
		#'EaeI' => 'Y^GGCCR',
		'EcoRV' => 'GAT^ATC',
		'EcoT22I' => 'ATGCA^T',
		'HindIII' => 'A^AGCTT',
		'KpnI' => 'GGTAC^C',
		'MluCI' => '^AATT',
		'MseI' => 'T^TAA',
		'NdeI' => 'CA^TATG',
		'NheI' => 'G^CTAGC',
		'NlaIII' => 'CATG^',
		'NotI' => 'GC^GGCCGC',
		'NsiI' => 'ATGCA^T',
		'RsaI' => 'GT^AC',
		'SacI' => 'GAGCT^C',
		'Sau3AI' => '^GATC',
		#'SexAI' => 'A^CCWGGT',
		#'SgrAI' => 'CR^CCGGYG',
		'SpeI' => 'A^CTAGT',
		'SphI' => 'GCATG^C',
		'TaqI' => 'T^CGA',
		'XbaI' => 'T^CTAGA',
		'XhoI' => 'C^TCGAG',
		'SbfI' => 'CCTGCA^GG',
		'EcoRI' => 'G^AATTC',
		'PstI' => 'CTGCA^G',
		'MspI' => 'C^CGG',
	);

# if -r flag is used, kill program and print supported enzyme list
if( $opts{r} ){
	&printenzymes( \%re );
	die "Exiting program following request to print supported enzyme list.\n\n";
}

# parse the command line options and return values
my( $first, $second, $error, $file, $size ) = &parsecom( \%opts );


# get the restriction site sequences
my $re1 = "$re{$first}" or die "\nMatch for first restriction enzyme not found in the program.\nCheck your input for flag -1.  Input is case sensitive.\nAlternatively, your enzyme may not be supported at this time.\n\n"; # string to hold sequence of first restriction enzyme
my $re2 = "$re{$second}" or die "\nMatch for second restriction enzyme not found in the program.\nCheck your input for flag -2.  Input is case sensitive.\nAlternatively, your enzyme may not be supported at this time.\n\n"; # string to hold sequence of second restriction enzyme
my $cut1 = "$re{$first}"; # string to hold sequence of RE1 with carot in place of cut site
my $cut2 = "$re{$second}"; # string to hold sequence of RE2 with carot in place of cut site
my $select = 0; # counter to hold the number of RE1 x RE2 fragments that fall within the desired size range
my $discard = 0; # counter to hold the number of RE1 x RE2 fragments that do not fall within the desired size range
my $re1xre1select = 0; # counter to hold the number of RE1 x RE1 fragments that fall within the desired size range
my $re1xre1discard = 0; # counter to hold the number of RE1 x RE1 fragments that do not fall within the desired size range
my $re2xre2select = 0; # counter to hold the number of RE2 x RE2 fragments that fall within the desired size range
my $re2xre2discard = 0; # counter to hold the number of RE2 x RE2 fragments that do not fall within the desired size range
my $re1frags = 0; # counts the total number of RE1 x RE1 fragments
my $re2frags = 0; # counts the total number of RE2 x RE2 fragments
my $ddradfrags = 0; # counts the total number of RE1 x RE2 fragments
my $totalfrags = 0; # counts the total number of fragments produced
my $re1out = "$first.fasta";
my $re2out = "$second.fasta";
my $re1xre2out = "$first.by.$second.fasta";
my $genomelength = 0; # calculate total length of genome
my $adjgenomelength = 0; # adjust genome length for Ns in sequences
my $gc = 0; # count gc bases
my $n = 0; # count the Ns in the genome
my @re1xre1; # array to hold RE1 x RE1 matches
my @re1xre2; # array to hold RE1 x RE2 matches
my @re2xre2; # array to hold RE2 x RE2 matches
my @lengths; # array to hold lengths of RE1 x RE2 matches
my @counts; # array to hold counts of lengths for RE1 x RE2 matches
my %hash; # hash to hold counts of all fragment lengths for RE1 x RE2 matches.  key = fragment length, value = count

# output file name
my $out = $first . '_by_' . $second . '_' . $size . '_' . $error . '.log';

# remove carot from restriction enzyme cut sites for matching
$re1 =~ s/\^//;
$re2 =~ s/\^//;

# convert the restriction enzyme cut sites to symbols
$cut1 =~ s/\^/1\^1/;
$cut2 =~ s/\^/2\^2/;

# parse the fasta file and return strings with header and full sequence
my( $head, $seqs ) = &fastaparse( $file );

# get minimum and maximum fragment lengths for size selection
my $min = ( $size - $error );
my $max = ( $size + $error );

# open files


foreach my $item( @$seqs ){
  # get total length of genome
  $genomelength += length( $item );
  $gc += ( $item =~ s/G/G/g );
  $gc += ( $item =~ s/C/C/g );
  $n += ( $item =~ s/N/N/g );
  #$n += $ncount;
  #$gc += $gcount;
  #$gc += $ccount;

  # find all restriction sites for the two enzymes and replace them with the cuts
  $item =~ s/$re1/$cut1/eg;
  $item =~ s/$re2/$cut2/eg;

  # split the string at the cut sites
  my @fragments = split( /\^/, $item );

  # add to the running total of fragments
  $totalfrags += scalar(@fragments);
  
  foreach my $item( @fragments ){
    if( $item =~ /^1.*1$/ ){
      &sizeselect( \$item, \$re1xre1select, \$re1xre1discard, \$min, \$max );
      &transform( \@re1xre1, \$item, \$re1frags );
    }elsif( $item =~ /^1.*2$|^2.*1$/ ){
      &sizeselect( \$item, \$select, \$discard, \$min, \$max );
      &transform( \@re1xre2, \$item, \$ddradfrags );
      my $len = length( $item );
      $hash{$len}++;
    }elsif( $item =~ /^2.*2$/ ){
      &sizeselect( \$item, \$re2xre2select, \$re2xre2discard, \$min, \$max );
      &transform( \@re2xre2, \$item, \$re2frags );
    }
  }
}

# close files

# open output file
open( OUT, '>', $out ) or die "Can't open $out: $!\n\n";

# print results summary statistics
# print genome length
print OUT "Genome length:\t", $genomelength, "\n\n";

# print GC Content
$adjgenomelength = ( $genomelength - $n );
my $gcpercent = ( $gc/$adjgenomelength );
print OUT "GC Content:\t", $gcpercent, "\n\n";

# print summary of all fragments
print OUT "All fragments (any combination):\n";
print OUT "Total:\t", $totalfrags, "\n\n";

# print summary of RE1 x RE2 fragments
print OUT "$re1 by $re2 fragments:\n";
print OUT "Total:\t", $ddradfrags, "\n";
print OUT $min, "bp to ", $max, "bp:\t", $select, "\n";
print OUT "Others:\t", $discard, "\n\n";

# print summary of RE1 x RE1 fragments
print OUT "$re1 by $re1 fragments:\t", $re1frags, "\n";
print OUT $min, "bp to ", $max, "bp:\t", $re1xre1select, "\n";
print OUT "Others:\t", $re1xre1discard, "\n\n";

# print summary of RE2 x RE2 fragments
print OUT "$re2 by $re2 fragments:\t", $re2frags, "\n";
print OUT $min, "bp to ", $max, "bp:\t", $re2xre2select, "\n";
print OUT "Others:\t", $re2xre2discard, "\n\n";


# calculate the percent of sequences that will be retained.
#my $total = ( $select + $discard );
my $unwanted = ($re1xre1select + $re2xre2select );
my $total = ( $select + $unwanted );
my $percent = ( $select / $total );
$percent = sprintf( '%.2f%%', 100 * $percent );
print OUT $percent, " of fragments between ", $min, " and ", $max, " bp will be selected.\n\n";

close OUT;

foreach my $key( sort {$a <=> $b} keys %hash ){
  $key += 0;
  push( @lengths, $key );
  push( @counts, $hash{$key} );
}
  
#print Dumper(\%hash);

#make a file name to send to gnuplot
my $name = join('.', $first, $second, "histogram", "png");

&gnuhistogram( \@lengths, \@counts, $name );

print "There are $re1frags $re1 by $re1 fragments\n\n";
print "There are $re2frags $re2 by $re2 fragments\n\n";
print "There are $ddradfrags $re1 by $re2 fragments\n\n";

exit;

#####################################################################################################
############################################ Subroutines ############################################
#####################################################################################################

# subroutine to print help
sub help{
  
  print "\nddrad.pl is a perl script developed by Steven Michael Mussmann\n\n";
  print "To report bugs send an email to mussmann\@email.uark.edu\n";
  print "When submitting bugs please include all input files, options used for the program, and all error messages that were printed to the screen\n\n";
  
  print "Program Options:\n";
  print "\t\t[ -1 | -2 | -e | -f | -h | -r | -s ]\n\n";
  print "\t-1:\tUse this to specify the first restriction enzyme.\n";
  print "\t\tProgram will terminate if no enzyme is specified.\n\n";

  print "\t-2:\tUse this to specify the second restriction enzyme.\n";
  print "\t\tProgram will terminate if no enzyme is specified.\n\n";

  print "\t-e:\tUse this to specify the error around the median value of your selected size range.\n";
  print "\t\tFor example, if you have an -s value of 300 and an -e value of 25, the program will calcu    late the number of fragments between 275bp and 325bp.\n";
  print "\t\tProgram will terminate if no value is specified.\n\n";

  print "\t-f:\tUse this to specify the input file name.\n";
  print "\t\tProgram will terminate if no file name is provided.\n\n";

  print "\t-h:\tUse this command to display this help message.\n";
  print "\t\tProgram will terminate after displaying this help message.\n\n";

  print "\t-r:\tUse this command to display a list of supported restriction enzymes.\n";
  print "\t\tProgram will terminate after displaying the enzyme list.\n\n";

  print "\t-s:\tUse this to specify the median value of your selected size range.\n";
  print "\t\tFor example, if you have an -s value of 300 and an -e value of 25, the program will calcu    late the number of fragments between 275bp and 325bp.\n";
  print "\t\tProgram will terminate if no value is specified.\n\n"; 
 
}

#####################################################################################################
# subroutine to print supported enzyme list
sub printenzymes{

	my( $re ) = @_;
	my %enz_seq = %$re;

	print "\nEnzyme", "\t\t", "Target sequence (^ denotes the cut site)", "\n\n";

	#print the list of enzymes
	foreach my $enzyme( sort{ $a cmp $b } keys %enz_seq ){
		print $enzyme, "\t\t", $enz_seq{$enzyme}, "\n";
	}

	print "\n";

}

#####################################################################################################
# subroutine to parse the command line options

sub parsecom{ 
  
  my( $params ) =  @_;
  my %opts = %$params;
  
  # set default values for command line arguments
  my $first = $opts{1} or die "\nFirst restriction enzyme not specified\n\n";
  my $second = $opts{2} or die "\nSecond restriction enzyme not specified\n\n";
  my $error = $opts{e} or die "\nError for size selection not specified (+ or - how many bases?)\n\n";
  my $file = $opts{f} or die "\nNo input file specified\n\n";
  my $size = $opts{s} or die "\nNo size selection value provided (how many bases?)\n\n";
  
  return( $first, $second, $error, $file, $size );

}

#####################################################################################################
# subroutine to parse a fasta file

sub fastaparse{
    my( $file ) = @_;
    
    # declare strings to hold header and sequence
    my @head;
    my @seqs;
    my $seq;
    
    # open the fasta file
    open( FASTA, $file ) or die "Can't open $file: $!\n\n";
    
    # loop through fasta file and extract data
    while( my $line = <FASTA> ){
      chomp( $line );
      if( $line =~/^>/ ){
	push( @head, $line );
	# push the previous sequence off to the @seqs array
	if( $seq ){
	  push( @seqs, uc($seq) );
	  # undefine the $seq string
	  undef( $seq );
	}
	# if the line does not start with > then append the sequence to the $seq string
      }else{
	$seq .= uc($line);
      }
    }
    
    # take care of the last sequence by pushing it to @seqs array
    if( $seq ){
      push( @seqs, uc($seq) );
      undef( $seq );
    }
    
    close FASTA;
    
    return( \@head, \@seqs );
    
}

#####################################################################################################
# subroutine to count fragments of a desired size range

sub sizeselect{
  
  my( $fragment, $select, $discard, $min, $max ) = @_;

  if( length( $$fragment ) ~~ [$$min..$$max] ){
    $$select ++;
  }else{
    $$discard++;
  }
  
}

#####################################################################################################
# subroutine to remove numbers from fragments and push sequences to appropriate arrays

sub transform{ 

  my( $array, $seq, $frags ) = @_;

  $$seq =~ s/[12]//g;
  push( @$array, $$seq );
  $$frags++;

}

#####################################################################################################

sub gnuhistogram{

  my( $data, $counts, $name ) = @_;

  # Create the chart object 
  my $chart = Chart::Gnuplot->new(
				  output => "$name", 
				  xlabel => "Fragment Size",
				  ylabel => "Fragment Count",
				  yrange => [0, '*'],
				  xrange => [0, 500],
				  xtics => {
					    labels => [0, 50, 100, 150, 200, 250, 300, 350, 400, 450, 500],
					   },
				 ); 
  # Data set object 
  my $dataSet = Chart::Gnuplot::DataSet->new( 
					     xdata => \@$data,
					     ydata => \@$counts,
					     title => "ddrad sim",
					     style => "impulses",
					    ); 
  # Plot the graph 
  $chart->plot2d($dataSet);
}
#####################################################################################################
