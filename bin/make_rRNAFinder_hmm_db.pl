#!/usr/bin/env perl
use File::Fetch;
use Archive::Extract;
use FindBin qw($Bin);
use strict;
use warnings;

my $now_string = localtime;  # e.g., "Thu Oct 13 04:54:34 1994"
print STDERR "build database at: $now_string\n";
my $bin = "$FindBin::Bin";
my $tmpdir = "$bin/../db/tmp";
my $hmm_dir = "$bin/../db/hmm";
print STDERR "hmmdb_dir=$hmm_dir\n";

my $fp = find_exe("hmmbuild");
if(!$fp){
    err("Cannot find hmmbuild tool");
}

my $rfam = "Rfam.seed";
if(! -e "$tmpdir/$rfam.gz"){
    my $ff = File::Fetch->new(uri => "ftp://ftp.ebi.ac.uk/pub/databases/Rfam/CURRENT/$rfam.gz");
    
### fetch the uri to local directory###
    print STDERR "Fetching $rfam.gz\n";
    my $where = $ff->fetch(to => $tmpdir) or die $ff->error;
}
if(! -e "$tmpdir/$rfam"){
    my $ae = Archive::Extract->new(archive =>"$tmpdir/$rfam.gz");
    $ae->extract(to => $tmpdir);
}
my %queries = (
     
    "5S_rRNA" => "$tmpdir/5S_rRNA.rfam.align.fasta",
    "5_8S_rRNA" =>"$tmpdir/5_8S_rRNA.rfam.align.fasta",
    "SSU_rRNA_bacteria" =>"$tmpdir/SSU_rRNA_bacteria.rfam.align.fasta",
    "SSU_rRNA_archaea" =>"$tmpdir/SSU_rRNA_archaea.rfam.align.fasta",
    "SSU_rRNA_eukarya" =>"$tmpdir/SSU_rRNA_eukarya.rfam.align.fasta",
    "LSU_rRNA_archaea" =>"$tmpdir/LSU_rRNA_archaea.rfam.align.fasta",
    "LSU_rRNA_bacteria" =>"$tmpdir/LSU_rRNA_bacteria.rfam.align.fasta",
    "LSU_rRNA_eukarya" =>"$tmpdir/LSU_rRNA_eukarya.rfam.align.fasta"
    );
open(STOCK, "$tmpdir/$rfam") or die "cannot open $rfam for reading: $!";

my $i = 0;
$/ = "\n\//";

foreach my $query (keys %queries){
    
    my $outname = $queries{$query};
    open(OUT, ">$outname") or die "cannot open $outname for writting: $!\n";
    while(<STOCK>){
	chomp;
	if(/\#=GF\s+ID\s+$query/s){
	    my @lines = split(/\n/, $_);
	    foreach my $line (@lines){
		if($line =~ /^\#/ || $line !~ /\S/){
		    next;
		} 	
		else{	
		    #seqid alignment
		    my @a = split(/\s+/, $line);
		    print OUT ">$a[0]\n$a[1]\n";
		}	
		
	    }
	    
	}
    }
    close(OUT);
    seek STOCK, 0, 0;
}
$/ = "\n";
close(STOCK);

fasta2domain("5S_rRNA.rfam.align.fasta", $tmpdir);
fasta2domain("5_8S_rRNA.rfam.align.fasta", $tmpdir);

my %hmms = (
    "arc_5SrRNA" => "$tmpdir/arc_5S_rRNA.rfam.align.fasta",
    "euk_5SrRNA" => "$tmpdir/euk_5S_rRNA.rfam.align.fasta",
    "bac_5SrRNA" => "$tmpdir/bac_5S_rRNA.rfam.align.fasta",
    "euk_5_8SrRNA" => "$tmpdir/euk_5_8S_rRNA.rfam.align.fasta",
    "arc_16SrRNA" => "$tmpdir/SSU_rRNA_archaea.rfam.align.fasta",
    "bac_16SrRNA" => "$tmpdir/SSU_rRNA_bacteria.rfam.align.fasta",
    "euk_18SrRNA" => "$tmpdir/SSU_rRNA_eukarya.rfam.align.fasta",
    "arc_23SrRNA" => "$tmpdir/LSU_rRNA_archaea.rfam.align.fasta",
    "bac_23SrRNA" => "$tmpdir/LSU_rRNA_bacteria.rfam.align.fasta",
    "euk_28SrRNA" => "$tmpdir/LSU_rRNA_eukarya.rfam.align.fasta"
    );

my $cmd = "";
foreach my $name (keys %hmms){
    my $alignFile = $hmms{$name};
    $cmd .= "hmmbuild --rna -n $name $tmpdir/$name.hmm $alignFile;"
}

print STDERR "cmd=$cmd\n";
system ($cmd) >> 8 and  die "Could not execute cmd=$cmd, $!\n";
$cmd = "cat $tmpdir/arc_16SrRNA.hmm $tmpdir/arc_23SrRNA.hmm $tmpdir/arc_5SrRNA.hmm > $hmm_dir/arc.hmm;";
$cmd .= "cat $tmpdir/bac_16SrRNA.hmm $tmpdir/bac_23SrRNA.hmm $tmpdir/bac_5SrRNA.hmm > $hmm_dir/bac.hmm;";
$cmd .= "cat $tmpdir/euk_18SrRNA.hmm $tmpdir/euk_28SrRNA.hmm $tmpdir/euk_5_8SrRNA.hmm $tmpdir/euk_5SrRNA.hmm > $hmm_dir/euk.hmm;";

print STDERR "cmd=$cmd\n";
system ($cmd) >> 8 and  die "Could not execute cmd=$cmd, $!\n";

my @tmp_files = glob "$tmpdir/*";
unlink glob "'$tmpdir/*'";
#rmdir $tmpdir or warn "couldn't rmdir $tmpdir: $!\n";


sub fasta2domain{
    my ($fasta, $tmpdir) = @_;
    use Bio::DB::EUtilities;
    use Bio::SeqIO;
    
    open (FASTA, "$tmpdir/$fasta") or die "Could not open $tmpdir/$fasta to read, $!\n";
    $/ = "\n>";
    open (BAC, ">$tmpdir/bac\_$fasta") or die "Could not open $tmpdir\/bac\_$fasta for writting, $!\n";
    open (ARC, ">$tmpdir\/arc\_$fasta") or die "Could not open $tmpdir\/arc\_$fasta for writting, $!\n";;
    open (EUK, ">$tmpdir\/euk\_$fasta") or die "Could not open $tmpdir\/euk\_$fasta for writting, $!\n";;
    
    while(<FASTA>){
	chomp;
	if(my ($seqid,$other, $align) =  /^>?(\S+?)(\/.*?)\n(.*)/s){
	    my $factory = Bio::DB::EUtilities->new(-eutil   => 'efetch',
						   -db      => 'nucleotide',
						   -rettype => 'gb',
						   -email   => 'xdong@ucalgary.ca',
						   -id      => $seqid);
	    my $file = "$tmpdir/myseqs.gb";
	    
	    # dump HTTP::Response content to a file (not retained in memory)
	    $factory->get_Response(-file => $file);
	    my $seqin = Bio::SeqIO->new(-file   => $file,
					-format => 'genbank');
	    
	    while (my $seq = $seqin->next_seq) {
		my @classification = $seq->species->classification;
		my $lineage = join('\t', @classification);
		if($lineage =~ /Bacteria/i){
		    print  BAC ">$seqid$other\n$align\n";
		}
		elsif($lineage =~ /Archaea/i){
		    print  ARC ">$seqid$other\n$align\n";
		}
		elsif($lineage =~ /Eukaryota/i){
		    print EUK ">$seqid$other\n$align\n";
		}
	    }
	
	}
	
    }
    
    $/ = "\n";
    close(FASTA);
    close(BAC);
    close(ARC);
    close(EUK);
}
sub find_exe {
    my($bin) = shift;
  for my $dir (File::Spec->path) {
    my $exe = File::Spec->catfile($dir, $bin);

    if(-x $exe){

	return $exe;
    }
  }
  return;
}
