#!/usr/bin/env perl
use Data::Dumper;

use Getopt::Long;
my ($input, $type, $level, $header);

&GetOptions("i=s" =>\$input,
	    "t=s" =>\$type,
	    "l=s" =>\$level
    );

($input && $type) ||
    die "usage: output json file for tree visulization\n" .
    "$0 OPTIONS where options are:\n" .
    "-i <tab delimited summary file in format of \"a;b;c number\"\n" .
    "-t <data type: e.g, Taxonomy, Subsystem...>\n".
    "-l <level exapanded: default is 3>\n";

$level ||= 2;
my $leading = 1;

open(INPUT, "<$input") or die "Could not open $input to read, $!\n";

my %kid_parent = ();
my %terms = ();
my %level_1 = ();

while (<INPUT>){
    chmop;
    if (/^#/){
	s/\s+$//;
	@titles = split("\t", $_);
	shift @titles;
	next;
    }

    my @cols = split(/\t/, $_);
    my $nodes_str = shift @cols;
    my @l_iterms = ();
   
    my $index = 1;
    $nodes_str =~ s/^\s+|\s+$//g;
    $nodes_str =~ s/\s+/\-/g;
    
    foreach my $item (split(/\";|;\t|;/, $nodes_str)){
	#this is for taxonomic, taxon database is messing, sometimes, one order is assigned to different class.
	next if $item =~ /unclassified|\-\[class\]|Terrabacteria-group|Cyanobacteria\/Melainabacteria-group|-group/i;
	$item =~ s/^\"|\"$//g;
	push(@l_iterms, "MYLEVEL$index\_$item");
	#print join("\t", @l_iterms), "\n";
	$index++;
   }

    my $len = @l_iterms;

    #print STDERR "***********$nodes_str len=$len 0=$terms[0];1=$terms[1];2=$terms[2]\n";

    if(@l_iterms && not exists $level_1{$l_iterms[0]}){
	if($l_iterms[0] eq "root"){
	    shift (@l_iterms);
	}
	$level_1{$l_iterms[0]} = 1;
	#print  "***********$l_iterms[0];\n";
    }

    for(my $i = 0; $i < $len; $i++){
	#each kid can have multiple parents
	#push(@{$kid_parent {trim($l_iterms[$i+1])}}, trim($l_iterms[$i])) if $i < $len -1;
	$kid_parent{trim($l_iterms[$i+1])}->{trim($l_iterms[$i])} = 1 if $i < $len -1;
	for (my $index = 0; $index < @cols; $index++){
	    if($cols[$index] =~ /N\/A/){
		$terms{trim($l_iterms[$i])}->{$titles[$index]} = "N/A";
		next;
	    }
	    $terms{trim($l_iterms[$i])}->{$titles[$index]} += $cols[$index];
	    #print "************$nodes_str\n";
	    #print "************ L72= $terms[$i] title=$titles[$index] data=$cols[$index]\n";
	}

    }

}
close(INPUT);

foreach my $kid (keys %kid_parent){
    foreach my $parent (keys %{$kid_parent{$kid}}){
	#print "kid=$kid, parent=$parent\n";
    }
}

foreach my $term (sort {$terms{$b} <=> $terms{$a}} keys %terms){
    foreach my $title (@titles){
	#print "$term\t", $terms{$term}->{$title}, "\n";
    }
}

my %parent_kid = ();


foreach my $kid (keys %kid_parent){

    #my @parents = keys %{$kid_parent{$kid}};
    #print  "##############kid=$kid parents=@parents\n";
    foreach my $parent (keys %{$kid_parent{$kid}}){
	#print  "##############kid=$kid parents=$parent\n";
	$parent_kid{$parent}->{$kid} = 1;
    }
}

#print (Dumper(%parent_kid));
foreach my $parent (keys %parent_kid){
    
    foreach my $kid (keys %{$parent_kid{$parent}}){
	#print "*****************parent=$parent, kid=$kid\n";
    }
}

print "var mydata = {\n";
print ' ' x 1;
print "\"text\":\"$type\",\n";
print "\"state\": {\"opened\": true},\n";
print ' ' x 1,  "\"children\": [\n";

my $root_term_count = keys %level_1;

my $count = 0;
foreach (keys %level_1){
    $count++;
    
    getKids(1, $_, \%parent_kid, \%terms, 2);
    if($count < $root_term_count){
	print ' ' x 2,  "},\n";
    }
    else{
	print ' ' x 2,  "}\n";
    }
}
print " ]\n";
print "}\n";
sub getKids{

    my ($leading, $current, $parent_kid, $terms, $spaceoffset) = @_;
    #print "*******current=$current\n";

    my $local_leading = $leading;
    
    my $node_str = $current;

    $node_str =~ s/\"/\\"/g;
    $node_str =~ s/^MYLEVEL\d+\_//;
    print ' ' x ($spaceoffset - 1), "{\n";
    print ' ' x $spaceoffset, "\"text\": \"$node_str\",\n";
    print ' ' x $spaceoffset, "\"data\": {";

    my @cols = ();
    foreach my $title (@titles){
	push(@cols, "\"$title\":" . sprintf("%.2f", $terms->{$current}->{$title}));
    }
    print join(",", @cols);

    my @kids = keys %{$parent_kid->{$current}};
    if (!@kids){
	print "}\n";
	return;
    }
    else{
	print "},\n";
	
    }

    #expand tree to class level
    if($leading < $level){
	print "\"state\": {\"opened\": true},\n";
    }


    print ' ' x $spaceoffset, "\"children\": [\n";
    my $num = 0;
    foreach my $kid (@kids){
	$num++;
	#print "current=$current, kid=$kid\n";
	getKids($leading+1, $kid, $parent_kid, $terms, $spaceoffset + 2);
	if($num < @kids){
	    print ' ' x ($spaceoffset + 2 -1), "},\n";
	}
	else{
	    print ' ' x ($spaceoffset + 2 -1),"}\n";
	}
    }


    print ' ' x $spaceoffset, "]\n";

}

#trim leading and trailing space
sub trim($)
{
    my $string = shift;
    $string =~ s/^\s+//;
    $string =~ s/\s+$//;
    $string =~ s/\'//g;
    $string =~ s/uncult ured/uncultured/g;
    $string =~ s/uncul tured/uncultured/g;
    $string =~ s/uncultu red/uncultured/g;
    $string =~ s/unculture d/uncultured/g;
    $string =~ s/;uncultured\s?;/;/g;
    $string =~ s/uncultured bacterium//g;


    return $string;
}
