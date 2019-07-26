#!/usr/bin/env perl
use strict;
use warnings;
use FindBin;
use IO::Uncompress::Gunzip qw(gunzip $GunzipError) ;
use File::Fetch;
use lib "$FindBin::Bin";
use Archive::Extract;
#use SWISS::Entry;
#use SWISS::KW;
#use XML::Simple;
#use XML::Parser;
use Getopt::Long;
use Time::Piece;
#use File::Basename;

my ($tmp_dir);

&GetOptions(
    "tmpdir=s" =>\$tmp_dir
    );

($tmp_dir) ||
    die "usage: $0 OPTIONS
where options are:\n -tmpdir  <tmp direcoty for intermedidate files>\n";

my $EXE = $FindBin::RealScript;
my $bin = "$FindBin::RealBin";
my $rna_hmm_dir = "$bin/../hmmrna";

msg("construct tmperorary $tmp_dir directories for intermedidate output");
runcmd("mkdir -p \Q$tmp_dir\E") if (! -d $tmp_dir);

build_rRNAFinder_hmmdb();
rmdir $tmp_dir or die "Could not delte $tmp_dir, $!\n";

sub build_rRNAFinder_hmmdb{

    #annotated seed alignments in STOCKHOLM format
    my $rfam = "Rfam.seed";
    if(! -e "$tmp_dir/$rfam.gz"){
	my $ff = File::Fetch->new(uri => "ftp://ftp.ebi.ac.uk/pub/databases/Rfam/CURRENT/$rfam.gz");

### fetch the uri to local directory###
	msg("Fetching $rfam.gz");
	my $where = $ff->fetch(to => $tmp_dir) or die $ff->error;
    }

    if(! -e "$tmp_dir/$rfam"){
	msg("Uncompressing $tmp_dir/$rfam.gz");
	my $ae = Archive::Extract->new(archive =>"$tmp_dir/$rfam.gz");
	$ae->extract(to => $tmp_dir);
    }
    my @rna_types = ("5S_rRNA","5_8S_rRNA","SSU_rRNA_bacteria","SSU_rRNA_archaea","SSU_rRNA_eukarya","LSU_rRNA_archaea","LSU_rRNA_bacteria","LSU_rRNA_eukarya");
    my %rna_align_filehandlers = ();


    my %queries = (
	"5S_rRNA" => "$tmp_dir/5S_rRNA.rfam.align.fasta",
	"5_8S_rRNA" =>"$tmp_dir/5_8S_rRNA.rfam.align.fasta",
	"SSU_rRNA_bacteria" =>"$tmp_dir/SSU_rRNA_bacteria.rfam.align.fasta",
	"SSU_rRNA_archaea" =>"$tmp_dir/SSU_rRNA_archaea.rfam.align.fasta",
	"SSU_rRNA_eukarya" =>"$tmp_dir/SSU_rRNA_eukarya.rfam.align.fasta",
	"LSU_rRNA_archaea" =>"$tmp_dir/LSU_rRNA_archaea.rfam.align.fasta",
	"LSU_rRNA_bacteria" =>"$tmp_dir/LSU_rRNA_bacteria.rfam.align.fasta",
	"LSU_rRNA_eukarya" =>"$tmp_dir/LSU_rRNA_eukarya.rfam.align.fasta"
	);
    foreach (keys %queries){
	open $rna_align_filehandlers{$_}, ">", $queries{$_} or die "cannot open $queries{$_} to write, $!\n";
    }

    open(STOCK, "$tmp_dir/$rfam") or die "cannot open $rfam for reading: $!";

    my $i = 0;
    $/ = "\n\//";
    while(<STOCK>){
	chomp;
	if(/\#=GF\s+ID\s+(\S+)/s){
	    my $id = $1;
	    if(exists $queries{$id}){
		my @lines = split(/\n/, $_);
		foreach my $line (@lines){
		    if($line =~ /^\#/ || $line !~ /\S/){
			next;
		    }
		    else{
			#seqid alignment
			my @a = split(/\s+/, $line);
			#msg("id=$id");
			next if $a[0] eq "AY544260.1/501-381";
			my $handler = $rna_align_filehandlers{$id};
			print $handler ">$a[0]\n$a[1]\n";
		    }

		}
	    }

	}
    }

    $/ = "\n";
    close(STOCK);
    #foreach (keys %rna_align_filehandlers){
#	close($_);
 #   }

    fasta2domain("5S_rRNA.rfam.align.fasta", $tmp_dir);
    fasta2domain("5_8S_rRNA.rfam.align.fasta", $tmp_dir);

    my %hmms = (
	"arc_5SrRNA" => "$tmp_dir/arc_5S_rRNA.rfam.align.fasta",
	"euk_5SrRNA" => "$tmp_dir/euk_5S_rRNA.rfam.align.fasta",
	"bac_5SrRNA" => "$tmp_dir/bac_5S_rRNA.rfam.align.fasta",
	"euk_5_8SrRNA" => "$tmp_dir/euk_5_8S_rRNA.rfam.align.fasta",
	"arc_16SrRNA" => "$tmp_dir/SSU_rRNA_archaea.rfam.align.fasta",
	"bac_16SrRNA" => "$tmp_dir/SSU_rRNA_bacteria.rfam.align.fasta",
	"euk_18SrRNA" => "$tmp_dir/SSU_rRNA_eukarya.rfam.align.fasta",
	"arc_23SrRNA" => "$tmp_dir/LSU_rRNA_archaea.rfam.align.fasta",
	"bac_23SrRNA" => "$tmp_dir/LSU_rRNA_bacteria.rfam.align.fasta",
	"euk_28SrRNA" => "$tmp_dir/LSU_rRNA_eukarya.rfam.align.fasta"
	);

    my $cmd = "";
    foreach my $name (keys %hmms){
	my $alignFile = $hmms{$name};
	$cmd .= "hmmbuild --rna -n $name $tmp_dir/$name.hmm $alignFile;"
    }

    msg("Start running:$cmd");
    runcmd($cmd);

    $cmd = "cat $tmp_dir/arc_16SrRNA.hmm $tmp_dir/arc_23SrRNA.hmm $tmp_dir/arc_5SrRNA.hmm > $rna_hmm_dir/arc.hmm;";
    $cmd .= "cat $tmp_dir/bac_16SrRNA.hmm $tmp_dir/bac_23SrRNA.hmm $tmp_dir/bac_5SrRNA.hmm > $rna_hmm_dir/bac.hmm;";
    $cmd .= "cat $tmp_dir/euk_18SrRNA.hmm $tmp_dir/euk_28SrRNA.hmm $tmp_dir/euk_5_8SrRNA.hmm $tmp_dir/euk_5SrRNA.hmm > $rna_hmm_dir/euk.hmm;";

    msg("Start running:$cmd");
    runcmd($cmd);

    #my @tmp_files = glob "$tmp_dir/*";
    #unlink glob "'$tmp_dir/*'";
#rmdir $tmp_dir or warn "couldn't rmdir $tmp_dir: $!\n";
}

sub fasta2domain{
    my ($fasta, $tmp_dir) = @_;
    use Bio::DB::EUtilities;
    use Bio::SeqIO;
    use Time::HiRes;
    open (FASTA, "$tmp_dir/$fasta") or die "Could not open $tmp_dir/$fasta to read, $!\n";
    $/ = "\n>";
    open (BAC, ">$tmp_dir/bac\_$fasta") or die "Could not open $tmp_dir\/bac\_$fasta for writting, $!\n";
    open (ARC, ">$tmp_dir\/arc\_$fasta") or die "Could not open $tmp_dir\/arc\_$fasta for writting, $!\n";;
    open (EUK, ">$tmp_dir\/euk\_$fasta") or die "Could not open $tmp_dir\/euk\_$fasta for writting, $!\n";;

    while(<FASTA>){
	chomp;
	if(my ($seqid,$other, $align) =  /^>?(\S+?)(\/.*?)\n(.*)/s){
	    
	    
	    
	    Time::HiRes::sleep(0.1);
	    #can get too many request problems. By default, the request from NCBI cannot be >3 times/second.
	    #with the api_key, it can be increased to 10 times/second perl email. 
	    
	    my $factory = Bio::DB::EUtilities->new(-eutil   => 'efetch',
						   -db      => 'nucleotide',
						   -rettype => 'gb',
						   -email   => 'xiaolid@gmail.com',
						   -api_key => 'f0f0c0d0dedaa42991657439b7c577b0ce08',
						   -id      => $seqid);
	    my $file = "$tmp_dir/myseqs.gb";
	    
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


sub err {
    my($txt) = @_;
    msg($txt);
    exit(2);
}

sub msg {
    my $t = localtime;
    #my $line = "[".$t->hms."] @_\n";
    my $line = "[".$t->cdate."] @_\n";
    print STDERR $line;
}

sub runcmd {
    msg("Running:", @_);
    my $cmd = join("\n", @_);
    system(@_)==0 or err("Could not run command:$cmd, $!\n");
}
