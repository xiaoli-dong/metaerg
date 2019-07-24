#!/usr/bin/env perl

use strict;
use warnings;
use Time::Piece;
use Scalar::Util qw(openhandle);
use List::Util qw[min max];

#Given two ranges [x1,x2], [y1,y2]
sub is_overlapping{
    my ($x1,$x2,$y1,$y2) = @_;
    return max($x1,$y1) <= min($x2,$y2)
}
sub msg {
  my $t = localtime;
  #my $line = "[".$t->hms."] @_\n";
  my $line = "[".$t->cdate."] @_\n";
  print STDERR $line;
}

#----------------------------------------------------------------------

sub err {
    my($txt) = @_;
    msg($txt);
    exit(2);
}

#----------------------------------------------------------------------

sub delfile {
    for my $file (@_) {
	msg("Deleting unwanted file:", $file);
	unlink $file or warn "Could not unlink $file: $!\n";;
    }
}



#----------------------------------------------------------------------

sub runcmd {
    msg("Running:", @_);
    my $cmd = join("\n", @_);
    system(@_)==0 or err("Could not run command:$cmd, $!\n");
}

sub decode {
  my $str = shift;
  $str =~ tr/+/ /;
  $str =~ s/%([0-9a-fA-F]{2})/chr hex($1)/ge;
  return $str;
}
sub encode {
    my $str = shift;
    $str =~ s/([\t\n\r%&\=;,])/sprintf("%%%X",ord($1))/ge;
    return $str;
}

sub mean {
    @_ == 1 or die ('Sub usage: $average = average(\@array);');
    my ($array_ref) = @_;
    my $sum;
    my $count = scalar @$array_ref;
    foreach (@$array_ref) { $sum += $_; }
    return sprintf("%.2f",$sum / $count);
}

sub median {
    @_ == 1 or die ('Sub usage: $median = median(\@array);');
    my ($array_ref) = @_;
    my $count = scalar @$array_ref;
# Sort a COPY of the array, leaving the original untouched
    my @array = sort { $a <=> $b } @$array_ref;
    if ($count % 2) {
	return $array[int($count/2)];
    } else {
	return sprintf("%.2f", ($array[$count/2] + $array[$count/2 - 1]) / 2);
    }
} 

sub stdev{
    my ($array_ref) = @_;
    my $sqsum = 0;
    my $n = @$array_ref;
    my $mean_value = sum(@$array_ref)/$n;
    
    for (@$array_ref) {
	$sqsum += ( $_ ** 2 );
    } 
    $sqsum /= $n;
    $sqsum -= ( $mean_value ** 2 );
    my $stdev = sqrt($sqsum);
    return sprintf("%.2f", $stdev);
}


sub get_N50{
    my ($arrRef) = @_;
    my @sort = sort {$b <=> $a} @$arrRef;
    my $totalLength = sum(@sort);
    my $n50 = 0;
    my $n50_value = 0;
    foreach my $val(@sort){
	$n50+=$val;
	if($n50 >= $totalLength/2){
	    #print "N50 length is $n50 and N50 value is: $val\n";
	    $n50_value = $val;
	    last;
	}
    }
    return $n50_value;
}

###################################################################
# _histogram_bins - calculates the bins usings Scott's algorithm
#
# Arguements:
#
# $data  (Vector)  - Data values
#
# $nbins (Integer) - Number of bins to create. If $nbins is undef
#                    the number of bins is calculated using Scott's
#                    algorithm
#
###################################################################
sub _histogram_bins {
	my ( $data, $nbins ) = @_;

	if( !defined $data ) { return; }

	my $calcBins = ( defined $nbins )? 0 : 1;
	my $cnt = 0;
	my $mean= 0;
	my $max = my $min = $data->[0];
	foreach (@$data) {
		$mean += $_;
		$min = ( $_ < $min )? $_ : $min;
		$max = ( $_ > $max )? $_ : $max;
		$cnt++;
	}
	$mean /= $cnt if( $cnt > 1 );

	my $sumsq = 0;
	$nbins = 1 if( $calcBins );
	my $s = 0;
	if( $cnt > 1 ) {
		foreach (@$data) {
			$sumsq += ( $_ - $mean )**2;
		}
		$s = sqrt( $sumsq / ($cnt - 1));
		$nbins = 3.49 * $s / $cnt**0.33 if( $s > 0 && $calcBins );
	}

	my $binwidth = ( $max - $min ) / $nbins;

	my $lower = $min;
	my $upper = $lower;

	my $bins;
	my @cutPoints;
	my $cntr = 0;
	while ( $upper <= $max && $cntr < $nbins) {
		$upper = $lower + $binwidth;
		push( @cutPoints, [int($lower), int($upper)] );
		$lower = $upper;
		$cntr++;
	}

	return \@cutPoints;
}
###################################################################
# _histogram_frequency - bins the data
#
#     Lower Boundry <= data value < Upper Boundry
#
# Arguements:
#
# $data  (Vector)  - Data values
#
# $nbins (Integer) - Vector containing the cutpoints to bin the data
#
###################################################################
sub _histogram_frequency {
	my ( $data, $cutPoints ) = @_;

	if( !defined $data || !defined $cutPoints ) { return; }

	my @freqs;
	foreach (@$cutPoints) {
		push( @freqs, 0 );
	}

	foreach (@$data) 
	{
		for( my $i = 0; $i < scalar( @$cutPoints ); $i++ ) 
		{
		if( ($_ >= $cutPoints->[$i]->[0] && $_ < $cutPoints->[$i]->[1])
			||
			($i == (scalar (@$cutPoints) - 1) && $_ >= $cutPoints->[$i]->[1]) ) 
			{	

				$freqs[$i]++;
			}
		}
	}
	return \@freqs;
}



sub calcgc {
    
    my ($seq) = @_;
    my $count = 0;
    my $len   = length($seq);
    for (my $i = 0;$i<$len; $i++) {
	my $base = substr $seq, $i, 1;
	$count++ if $base =~ /[G|C]/i;
    }
    return sprintf("%.2f",($count / $len) * 100);
}
sub run_cmds {

        my ($max_cmds, @cmds) = @_;

        my ($num_children, $pid);

        return unless @cmds; # get out of the function if there is nothing in @cmds

        for($num_children = 1; $num_children < $max_cmds && @cmds; $num_children++){
	    
	    # initialize the number of child processes at 1, and increment it by one
	    #while it is less than $max_cmds

	    my $cmd = shift (@cmds);
	    if($pid = fork) {
		# do nothing if parent
	    }
	    elsif(defined $pid) { # $pid is zero here if defined
		
		system $cmd;
		
		exit;
	    }
	    else {
		#weird fork error
		die "Can't fork: $!\n";
	    }
        }
	
        while(@cmds) {
	    undef $pid;
	  
	  FORK: {
	      my $cmd = shift (@cmds);
	      if($pid = fork) {
		  # parent here
		  $num_children++;
		  wait;
		  $num_children--;
		  next;
		  
	      }
	      elsif(defined $pid) { # $pid is zero here if defined
		
		  system $cmd;
				  
		  exit;
		  
	      }
	      else {
		  #weird fork error
		  die "Can't fork: $!\n";
	      }
	    }
        }
        wait while $num_children--;
}



1;

__END__
    
