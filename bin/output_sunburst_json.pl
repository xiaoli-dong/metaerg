#!/usr/bin/env perl

use Getopt::Long;

my ($input, $type, $prefix);

&GetOptions("i=s" =>\$input,
	    "t=s" =>\$type,
	    "p=s" =>\$prefix
    );

($input && $type && $prefix) ||
    die "usage: output json file for sunburst visulization\n" .
    "$0 OPTIONS where options are:\n" .
    "-i <tab delimited summary file in format of \"a;b;c number\"\n" .
    "-p <output file prefix>\n" .
    "-t <data type: e.g, Taxonomy, Subsystem...>\n";

my $leading = 1;

open(INPUT, "$input") or die "Could not open $input to read, $!\n";
open(OUT_C, ">$prefix.count.json") or die "Could not open $prefix.count.json to write, $!\n";


my %kid_parent = ();
my %terms = ();
my %level_1 = ();
my $root_count = 0;

while (<INPUT>){
    chmop;
    next if/^#/;
    next if /^\"?UNKNOWN/i;
    
    if($_ !~ /(\w+.*?)\t+(\d+)/){
	next;
    }
    
    my $str = $1;
    my $count = $2;
    
    my @l_iterms = ();
    my $index = 1;
    #print STDERR "**********$str\n";
    foreach my $item (split(/\";|;|;\t/, $str)){
	
	#this is for taxonomic, taxon database is messing, sometimes, one order is assigned to different class.
	next if $item =~ /unclassified|\-\[class\]|Terrabacteria-group|Cyanobacteria\/Melainabacteria-group|-group/i;
	$item =~ s/^\"//;
	$item =~ s/\"$//;
	push(@l_iterms, "MYLEVEL$index\_$item");
	#print STDERR join("\t", @l_iterms), "\n";
	$index++;
    }
    my $len = @l_iterms;
    
    $root_count += $count;
    
    if(@l_iterms && not exists $level_1{$l_iterms[0]}){
	if($l_iterms[0] eq "root"){
	    shift (@l_iterms);
	}
	$level_1{$l_iterms[0]} = 1;

    }
    
    #for(my $i = 0; $i < $len - 1; $i++){
	
	#$kid_parent {trim($l_iterms[$i+1])} = trim($l_iterms[$i]);
    #}
    for(my $i = 0; $i < $len; $i++){
	#each kid can have multiple parents
	$kid_parent{trim($l_iterms[$i+1])}->{trim($l_iterms[$i])} = 1 if $i < $len -1;
	
    }

    
    $terms{trim($l_iterms[$len -1])}->{count} += $count;
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

foreach my $kid (keys %kid_parent){

    #my @parents = keys %{$kid_parent{$kid}};
    #print  "##############kid=$kid parents=@parents\n";
    foreach my $parent (keys %{$kid_parent{$kid}}){
	#print  "##############kid=$kid parents=$parent\n";
	$parent_kid{$parent}->{$kid} = 1;
    }
}
print OUT_C "data = {\n";
print OUT_C ' ' x 1;
print OUT_C "\"text\":\"$type\",\n";
print OUT_C ' ' x 1,  "\"children\": [\n";

my $root_term_count = keys %level_1;
my $count = 0;
foreach (keys %level_1){
    $count++;
    getKids(1, $_, \%parent_kid, \%terms, 2);
    if($count < $root_term_count){
	print OUT_C ' ' x 2,  "},\n";
	
    }
    else{
	print OUT_C ' ' x 2,  "}\n";
	
    }
}
print OUT_C "]}\n";


sub getKids{

    my ($leading, $current, $parent_kid, $terms, $spaceoffset) = @_;
    my $local_leading = $leading;
    my $node_str = $current;
    
    $node_str =~ s/\"/\\"/g;
    $node_str =~ s/MYLEVEL\d+\_//;
    
    print OUT_C ' ' x ($spaceoffset - 1), "{\n";
    #print OUT_A ' ' x ($spaceoffset - 1), "{\n";
    
    
    
    my $current_count = $terms->{$current}->{count};
    my $current_count_pct = sprintf("%.2f", $current_count*100/$root_count);

    #my @kids = split (/\t/, $parent_kid->{$current});
    my @kids = keys %{$parent_kid->{$current}};
    if(@kids){
	print OUT_C ' ' x $spaceoffset, "\"text\": \"$node_str\",\n";
	
	if(exists $terms{$current}->{count}){
	    print OUT_C ' ' x $spaceoffset,  "\"size\": ", $terms->{$current}->{count}, ",\n";
	   
	}

	print OUT_C ' ' x $spaceoffset, "\"children\":[\n";
	
	
	my $num = 0;
	foreach my $kid (@kids){
	    $num++;

	    getKids($leading+1, $kid, $parent_kid, $terms, $spaceoffset + 2);

	    if($num < @kids){
		print OUT_C ' ' x ($spaceoffset + 2 -1), "},\n";
	
	    }
	    else{
		print OUT_C ' ' x ($spaceoffset + 2 -1), "}\n";
	
	    }
	}

	print OUT_C ' ' x $spaceoffset, "]\n";
	
    }
    else{

	if(exists $terms{$current}->{count}){
	    print OUT_C ' ' x $spaceoffset, "\"text\": \"$node_str\",\n";
	    print OUT_C ' ' x $spaceoffset,  "\"size\": ", $terms->{$current}->{count}, "\n";
	
	}
	else{
	    print OUT_C ' ' x $spaceoffset, "\"text\": \"$node_str\"\n";
	
	}

    }

}

#trim leading and trailing space
sub trim{
    my ($string) = @_;
     $string =~ s/^\s+//;
    $string =~ s/\s+$//;
    $string =~ s/\'//g;
    $string =~ s/(unclassified|uncultured)//g;

    return $string;
}
