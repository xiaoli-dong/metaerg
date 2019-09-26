#!/usr/bin/env perl
@ARGV == 4 or die "Usage: $0 <HMM file> <number of files to split> <output directory><output file prefix>\n",
                  "split HMM model input file into multiple smaller files\n";

$infile = $ARGV[0];
$num_files = $ARGV[1];
$output_dir = $ARGV[2];
$output_file_prefix = $ARGV[3];
chomp ($infile);

open (INFILE, "$infile") or die "Could not open $infile to read, $!\n";

#count sequencs it has in file

$/="//\n";
my $total = 0;
while(<INFILE>){
    $total++;
}
$/="\n";
print "$infile has $total modles\n";

$remainder = $total % $num_files;
$seq_per_file = ($total - $remainder)/$num_files;

open (INFILE, "<$infile");
$count = 0;
$file_count = 1;

$output_file = "$output_dir/$output_file_prefix$file_count.hmm";

open (OUTPUT, ">$output_file") || die "cannot open the file $output_file\n";
#split input file into multiple files
$/="//\n";
while(<INFILE>){
    print OUTPUT "$_";
    $count++;
    
    if($file_count != $num_files && $count == $seq_per_file){
	close(OUTPUT);
	$count = 0;
	
	if($file_count < $num_files){
	    $file_count++;
	    $output_file = "$output_dir/$output_file_prefix$file_count.hmm";
	    open (OUTPUT, ">$output_file") || die "cannot open the file $output_file\n";
	}
    }
}
#print OUTPUT "\n\//\n";
close(INFILE);
close(OUTPUT);
