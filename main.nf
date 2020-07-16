$HOSTNAME = ""
params.outdir = 'results'  

//pipeline defaults
params.genome_build =  ""  //* @dropdown @options:"human_hg38_gencode_v30, custom"
params.nucleicAcidType = "rna"
params.run_ncRNA_Removal =  "yes"  //* @dropdown @options:"yes, no" @show_settings:"ncRNA_Removal"
params.run_STAR =  "yes"  //* @dropdown @options:"yes, no" @show_settings:"STAR_align_Single_Best_Multimapper"
params.run_riboseq_workflow =  "yes"  //* @dropdown @options:"yes, no" @show_settings:"Densebuilder"


_species = ""
_build = ""
_share = ""
_annotation = ""  
//* autofill
if (params.genome_build == "human_hg38_gencode_v30"){
    _species = "human"
    _build = "hg38"
    _annotation = "gencode_v30"
} 


if ($HOSTNAME == "default"){
    _share = "${params.DOWNDIR}/genome_data"
    $SINGULARITY_IMAGE = "onuryukselen/ribosome-profiling:1.0"
    $DOCKER_IMAGE = "onuryukselen/ribosome-profiling:1.0"
    $DEFAULT_IMAGE = "singularity"
}
//* platform
if ($HOSTNAME == "ghpcc06.umassrc.org"){
    _share = "/share/data/umw_biocore/genome_data/ribosome-profiling"
    $SINGULARITY_IMAGE = "/project/umw_biocore/singularity/UMMS-Biocore-ribosome-profiling-1.0.simg"
    $TIME = 240
    $CPU  = 1
    $MEMORY = 32 
    $QUEUE = "short"
} 
//* platform
if (params.genome_build && $HOSTNAME){
    params.genome       = "${_share}/${_species}/${_build}/main/genome.fa"
    params.genome2bit   = "${_share}/${_species}/${_build}/main/genome.2bit"
    params.single_transcript_gtf = "${_share}/${_species}/${_build}/${_annotation}/gencodeV30_protCode_TermStopCodon_validUTRs.gtf"
    params.single_transcript_mRNA_csv = "${_share}/${_species}/${_build}/${_annotation}/gencodeV30_protCode_TermStopCodon_validUTRs_mRNAseqs.csv"
    params.single_transcript_UTR_csv = "${_share}/${_species}/${_build}/${_annotation}/gencodeV30_protCode_TermStopCodon_validUTRs_UTRs.csv"
    params.single_transcript_stopcodons_csv = "${_share}/${_species}/${_build}/${_annotation}/gencodeV30_protCode_TermStopCodon_validUTRs_stopcodons.csv"
    params.motif_file_root_path = "${_share}/${_species}/${_build}/${_annotation}/gencodeV30_protCode_TermStopCodon_validUTRs_codons/gencodeV30_protCode_TermStopCodon_validUTRs_"
    params.star_index    = "${_share}/${_species}/${_build}/${_annotation}/STARIndex/"
    params.ncRNA_star_index = "${_share}/${_species}/${_build}/${_annotation}/ncRNA_STARIndex/"

}
if ($HOSTNAME){
    params.densebuilder_main_path = "/usr/local/bin/G418_readthrough/riboseq/densebuilder_main.py"
    params.makeavggene_main_path = "/usr/local/bin/G418_readthrough/riboseq/makeavggene_main.py"
    params.riboseq_buildDenseTables_rpkm_path = "/usr/local/bin/G418_readthrough/riboseq/riboseq_buildDenseTables_rpkm.py"
    params.riboseq_buildDenseTables_rpkm_utr3adj_path = "/usr/local/bin/G418_readthrough/riboseq/riboseq_buildDenseTables_rpkm_utr3adj.py"
    params.rphelper_path = "/usr/local/bin/G418_readthrough/riboseq/"
    params.codon_occupancy_main_path = "/usr/local/bin/G418_readthrough/riboseq/codon_occupancy_main.py"
}
//* autofill


if (!params.mate){params.mate = ""} 
if (!params.reads){params.reads = ""} 

Channel.value(params.mate).into{g_1_mate_g_3;g_1_mate_g_4;g_1_mate_g0_3;g_1_mate_g0_11;g_1_mate_g0_16;g_1_mate_g0_18;g_1_mate_g0_19;g_1_mate_g0_20;g_1_mate_g0_21}
Channel
	.fromFilePairs( params.reads , size: (params.mate != "pair") ? 1 : 2 )
	.ifEmpty { error "Cannot find any reads matching: ${params.reads}" }
	.into{g_2_reads_g0_3;g_2_reads_g0_18}


params.run_Adapter_Removal =   "no"   //* @dropdown @options:"yes","no" @show_settings:"Adapter_Removal"
//* @style @multicolumn:{seed_mismatches, palindrome_clip_threshold, simple_clip_threshold} @condition:{Tool_for_Adapter_Removal="trimmomatic", seed_mismatches, palindrome_clip_threshold, simple_clip_threshold}, {Tool_for_Adapter_Removal="fastx_clipper", discard_non_clipped}

//* autofill
//* platform
if ($HOSTNAME == "ghpcc06.umassrc.org"){
    $TIME = 1000
    $CPU  = 1
    $MEMORY = 24
    $QUEUE = "long"
}
//* platform
//* autofill
if (!((params.run_Adapter_Removal && (params.run_Adapter_Removal == "yes")) || !params.run_Adapter_Removal)){
g_2_reads_g0_18.into{g0_18_reads_g0_19}
g0_18_log_file_g0_11 = Channel.empty()
} else {


process Adapter_Trimmer_Quality_Module_Adapter_Removal {

input:
 set val(name), file(reads) from g_2_reads_g0_18
 val mate from g_1_mate_g0_18

output:
 set val(name), file("reads/*.fastq")  into g0_18_reads_g0_19
 file "*.{fastx,trimmomatic}.log"  into g0_18_log_file_g0_11

errorStrategy 'retry'

when:
(params.run_Adapter_Removal && (params.run_Adapter_Removal == "yes")) || !params.run_Adapter_Removal

shell:
phred = params.Adapter_Trimmer_Quality_Module_Adapter_Removal.phred
Tool_for_Adapter_Removal = params.Adapter_Trimmer_Quality_Module_Adapter_Removal.Tool_for_Adapter_Removal
Adapter_Sequence = params.Adapter_Trimmer_Quality_Module_Adapter_Removal.Adapter_Sequence
//trimmomatic_inputs
min_length = params.Adapter_Trimmer_Quality_Module_Adapter_Removal.min_length
seed_mismatches = params.Adapter_Trimmer_Quality_Module_Adapter_Removal.seed_mismatches
palindrome_clip_threshold = params.Adapter_Trimmer_Quality_Module_Adapter_Removal.palindrome_clip_threshold
simple_clip_threshold = params.Adapter_Trimmer_Quality_Module_Adapter_Removal.simple_clip_threshold

//fastx_clipper_inputs
discard_non_clipped = params.Adapter_Trimmer_Quality_Module_Adapter_Removal.discard_non_clipped
    
remove_previous_reads = params.Adapter_Trimmer_Quality_Module_Adapter_Removal.remove_previous_reads
discard_non_clipped_text = ""
if (discard_non_clipped == "yes") {discard_non_clipped_text = "-c"}
nameAll = reads.toString()
nameArray = nameAll.split(' ')
file2 = ""
if (nameAll.contains('.gz')) {
    newName =  nameArray[0] - ~/(\.fastq.gz)?(\.fq.gz)?$/
    file1 =  nameArray[0] - '.gz' 
    if (mate == "pair") {file2 =  nameArray[1] - '.gz'}
    runGzip = "ls *.gz | xargs -i echo gzip -df {} | sh"
} else {
    newName =  nameArray[0] - ~/(\.fastq)?(\.fq)?$/
    file1 =  nameArray[0]
    if (mate == "pair") {file2 =  nameArray[1]}
    runGzip = ''
}
'''
#!/usr/bin/env perl
 use List::Util qw[min max];
 use strict;
 use File::Basename;
 use Getopt::Long;
 use Pod::Usage;
 use Cwd qw();
 
runCmd("mkdir reads adapter unpaired");

open(OUT, ">adapter/adapter.fa");
my @adaps=split(/\n/,"!{Adapter_Sequence}");
my $i=1;
foreach my $adap (@adaps)
{
 print OUT ">adapter$i\\n$adap\\n";
 $i++;
}
close(OUT);

system("!{runGzip}");
my $quality="!{phred}";
print "fastq quality: $quality\\n";
print "tool: !{Tool_for_Adapter_Removal}\\n";

if ("!{mate}" eq "pair") {
    if ("!{Tool_for_Adapter_Removal}" eq "trimmomatic") {
        runCmd("trimmomatic PE -threads 1 -phred${quality} !{file1} !{file2} reads/!{name}.1.fastq unpaired/!{name}.1.fastq.unpaired reads/!{name}.2.fastq unpaired/!{name}.2.fastq.unpaired ILLUMINACLIP:adapter/adapter.fa:!{seed_mismatches}:!{palindrome_clip_threshold}:!{simple_clip_threshold} MINLEN:!{min_length} 2> !{name}.trimmomatic.log");
    } elsif ("!{Tool_for_Adapter_Removal}" eq "fastx_clipper") {
        print "Fastx_clipper is not suitable for paired reads.";
    }
} else {
    if ("!{Tool_for_Adapter_Removal}" eq "trimmomatic") {
        runCmd("trimmomatic SE -threads 1  -phred${quality} !{file1} reads/!{name}.fastq ILLUMINACLIP:adapter/adapter.fa:!{seed_mismatches}:!{palindrome_clip_threshold}:!{simple_clip_threshold} MINLEN:!{min_length} 2> !{name}.trimmomatic.log");
    } elsif ("!{Tool_for_Adapter_Removal}" eq "fastx_clipper") {
        runCmd("fastx_clipper  -Q $quality -a !{Adapter_Sequence} -l !{min_length} !{discard_non_clipped_text} -v -i !{file1} -o reads/!{name}.fastq > !{name}.fastx.log");
    }
}
if ("!{remove_previous_reads}" eq "true") {
    my $currpath = Cwd::cwd();
    my @paths = (split '/', $currpath);
    splice(@paths, -2);
    my $workdir= join '/', @paths;
    splice(@paths, -1);
    my $inputsdir = join '/', @paths;
    $inputsdir .= "/inputs";
    print "INFO: inputs reads will be removed if they are located in the workdir inputsdir\\n";
    my @listOfFiles = `readlink -e !{file1} !{file2}`;
    foreach my $targetFile (@listOfFiles){
        if (index($targetFile, $workdir) != -1 || index($targetFile, $inputsdir) != -1) {
            system("rm -f $targetFile");
            print "INFO: $targetFile deleted.\\n";
        }
    }
}


##Subroutines
sub runCmd {
    my ($com) = @_;
    my $error = system($com);
    if   ($error) { die "Command failed: $error $com\\n"; }
    else          { print "Command successful: $com\\n"; }
}
'''

}
}


params.run_Trimmer =   "no"   //* @dropdown @options:"yes","no" @show_settings:"Trimmer"
//* @style @multicolumn:{trim_length_5prime,trim_length_3prime}, {trim_length_5prime_R1,trim_length_3prime_R1}, {trim_length_5prime_R2,trim_length_3prime_R2} @condition:{single_or_paired_end_reads="single", trim_length_5prime,trim_length_3prime}, {single_or_paired_end_reads="pair", trim_length_5prime_R1,trim_length_3prime_R1,trim_length_5prime_R2,trim_length_3prime_R2}

//* autofill
//* platform
if ($HOSTNAME == "ghpcc06.umassrc.org"){
    $TIME = 500
    $CPU  = 1
    $MEMORY = 8
    $QUEUE = "long"
}
//* platform
//* autofill
if (!((params.run_Trimmer && (params.run_Trimmer == "yes")) || !params.run_Trimmer)){
g0_18_reads_g0_19.into{g0_19_reads_g0_20}
g0_19_log_file_g0_21 = Channel.empty()
} else {


process Adapter_Trimmer_Quality_Module_Trimmer {

input:
 set val(name), file(reads) from g0_18_reads_g0_19
 val mate from g_1_mate_g0_19

output:
 set val(name), file("reads/*q")  into g0_19_reads_g0_20
 file "*.log" optional true  into g0_19_log_file_g0_21

errorStrategy 'retry'

when:
(params.run_Trimmer && (params.run_Trimmer == "yes")) || !params.run_Trimmer

shell:
phred = params.Adapter_Trimmer_Quality_Module_Trimmer.phred
single_or_paired_end_reads = params.Adapter_Trimmer_Quality_Module_Trimmer.single_or_paired_end_reads
trim_length_5prime = params.Adapter_Trimmer_Quality_Module_Trimmer.trim_length_5prime
trim_length_3prime = params.Adapter_Trimmer_Quality_Module_Trimmer.trim_length_3prime
trim_length_5prime_R1 = params.Adapter_Trimmer_Quality_Module_Trimmer.trim_length_5prime_R1
trim_length_3prime_R1 = params.Adapter_Trimmer_Quality_Module_Trimmer.trim_length_3prime_R1
trim_length_5prime_R2 = params.Adapter_Trimmer_Quality_Module_Trimmer.trim_length_5prime_R2
trim_length_3prime_R2 = params.Adapter_Trimmer_Quality_Module_Trimmer.trim_length_3prime_R2
remove_previous_reads = params.Adapter_Trimmer_Quality_Module_Trimmer.remove_previous_reads

nameAll = reads.toString()
nameArray = nameAll.split(' ')
file2 = ""
if (nameAll.contains('.gz')) {
    newName =  nameArray[0] - ~/(\.fastq.gz)?(\.fq.gz)?$/
    file1 =  nameArray[0] - '.gz' 
    if (mate == "pair") {file2 =  nameArray[1] - '.gz'}
    runGzip = "ls *.gz | xargs -i echo gzip -df {} | sh"
} else {
    newName =  nameArray[0] - ~/(\.fastq)?(\.fq)?$/
    file1 =  nameArray[0]
    if (mate == "pair") {file2 =  nameArray[1]}
    runGzip = ''
}
'''
#!/usr/bin/env perl
 use List::Util qw[min max];
 use strict;
 use File::Basename;
 use Getopt::Long;
 use Pod::Usage; 
 use Cwd qw();
 
system("mkdir reads");
system("!{runGzip}");
my $file1 = "";
my $file2 = "";
if ("!{mate}" eq "pair") {
    $file1 = "!{file1}";
    $file2 = "!{file2}";
    my $trim1 = "!{trim_length_5prime_R1}:!{trim_length_3prime_R1}";
    my $trim2 = "!{trim_length_5prime_R2}:!{trim_length_3prime_R2}";
    my $len=getLength($file1);
    print "length of $file1: $len\\n";
    trimFiles($file1, $trim1, $len);
    my $len=getLength($file2);
    print "INFO: length of $file2: $len\\n";
    trimFiles($file2, $trim2, $len);
} else {
    $file1 = "!{file1}";
    my $trim1 = "!{trim_length_5prime}:!{trim_length_3prime}";
    my $len=getLength($file1);
    print "INFO: length of file1: $len\\n";
    trimFiles($file1, $trim1, $len);
}
if ("!{remove_previous_reads}" eq "true") {
    my $currpath = Cwd::cwd();
    my @paths = (split '/', $currpath);
    splice(@paths, -2);
    my $workdir= join '/', @paths;
    splice(@paths, -1);
    my $inputsdir= join '/', @paths;
    $inputsdir .= "/inputs";
    print "INFO: inputs reads will be removed if they are located in the workdir inputsdir\\n";
    my @listOfFiles = `readlink -e !{file1} !{file2}`;
    foreach my $targetFile (@listOfFiles){
        if (index($targetFile, $workdir) != -1 || index($targetFile, $inputsdir) != -1) {
            system("rm -f $targetFile");
            print "INFO: $targetFile deleted.\\n";
        }
    }
}



sub trimFiles
{
  my ($file, $trim, $len)=@_;
    my @nts=split(/[,:\\s\\t]+/,$trim);
    my $inpfile="";
    my $com="";
    my $i=1;
    my $outfile="";
    my $param="";
    my $quality="-Q!{phred}";

    if (scalar(@nts)==2)
    {
      $param = "-f ".($nts[0]+1) if (exists($nts[0]) && $nts[0] >= 0 );
      $param .= " -l ".($len-$nts[1]) if (exists($nts[0]) && $nts[1] > 0 );
      $outfile="reads/$file";  
      $com="fastx_trimmer $quality -v $param -o $outfile -i $file > !{name}.fastx_trimmer.log" if ((exists($nts[0]) && $nts[0] > 0) || (exists($nts[0]) && $nts[1] > 0 ));
      print "INFO: $com\\n";
      if ($com eq ""){
          print "INFO: Trimmer skipped for $file \\n";
          system("mv $file reads/.");
      } else {
          runCmd("$com");
          print "INFO: Trimmer executed for $file \\n";
      }
    }

    
}


sub getLength
{
   my ($filename)=@_;
   open (IN, $filename);
   my $j=1;
   my $len=0;
   while(my $line=<IN>)
   {
     chomp($line);
     if ($j >50) { last;}
     if ($j%4==0)
     {
        $len=length($line);
     }
     $j++;
   }
   close(IN);
   return $len;
}

sub runCmd {
    my ($com) = @_;
    my $error = system($com);
    if   ($error) { die "Command failed: $error $com\\n"; }
    else          { print "Command successful: $com\\n"; }
}

'''

}
}



process Adapter_Trimmer_Quality_Module_Trimmer_Removal_Summary {

input:
 file logfile from g0_19_log_file_g0_21.collect()
 val mate from g_1_mate_g0_21

output:
 file "trimmer_summary.tsv"  into g0_21_outputFileTSV_g_26

shell:
'''
#!/usr/bin/env perl
use List::Util qw[min max];
use strict;
use File::Basename;
use Getopt::Long;
use Pod::Usage;
use Data::Dumper;

my @header;
my %all_files;
my %tsv;
my %tsvDetail;
my %headerHash;
my %headerText;
my %headerTextDetail;

my $i = 0;
chomp( my $contents = `ls *.log` );

my @files = split( /[\\n]+/, $contents );
foreach my $file (@files) {
    $i++;
    my $mapOrder = "1";
    if ($file =~ /(.*)\\.fastx_trimmer\\.log/){
        $file =~ /(.*)\\.fastx_trimmer\\.log/;
        my $mapper   = "fastx_trimmer";
        my $name = $1;    ##sample name
        push( @header, $mapper );
        my $in;
        my $out;
        chomp( $in =`cat $file | grep 'Input:' | awk '{sum+=\\$2} END {print sum}'` );
        chomp( $out =`cat $file | grep 'Output:' | awk '{sum+=\\$2} END {print sum}'` );

        $tsv{$name}{$mapper} = [ $in, $out ];
        $headerHash{$mapOrder} = $mapper;
        $headerText{$mapOrder} = [ "Total Reads", "Reads After Trimmer" ];
    }
}

my @mapOrderArray = ( keys %headerHash );
my @sortedOrderArray = sort { $a <=> $b } @mapOrderArray;

my $summary          = "trimmer_summary.tsv";
writeFile( $summary,          \\%headerText,       \\%tsv );

sub writeFile {
    my $summary    = $_[0];
    my %headerText = %{ $_[1] };
    my %tsv        = %{ $_[2] };
    open( OUT, ">$summary" );
    print OUT "Sample\\t";
    my @headArr = ();
    for my $mapOrder (@sortedOrderArray) {
        push( @headArr, @{ $headerText{$mapOrder} } );
    }
    my $headArrAll = join( "\\t", @headArr );
    print OUT "$headArrAll\\n";

    foreach my $name ( keys %tsv ) {
        my @rowArr = ();
        for my $mapOrder (@sortedOrderArray) {
            push( @rowArr, @{ $tsv{$name}{ $headerHash{$mapOrder} } } );
        }
        my $rowArrAll = join( "\\t", @rowArr );
        print OUT "$name\\t$rowArrAll\\n";
    }
    close(OUT);
}

'''
}

params.run_Quality_Filtering =   "no"   //* @dropdown @options:"yes","no" @show_settings:"Quality_Filtering"
//* @style @multicolumn:{window_size,required_quality}, {leading,trailing,minlen}, {minQuality,minPercent} @condition:{tool="trimmomatic", minlen, trailing, leading, required_quality_for_window_trimming, window_size}, {tool="fastx", minQuality, minPercent}

//* autofill
//* platform
if ($HOSTNAME == "ghpcc06.umassrc.org"){
    $TIME = 240
    $CPU  = 1
    $MEMORY = 8
    $QUEUE = "short"
}
//* platform
//* autofill
if (!((params.run_Quality_Filtering && (params.run_Quality_Filtering == "yes")) || !params.run_Quality_Filtering)){
g0_19_reads_g0_20.into{g0_20_reads_g_3}
g0_20_log_file_g0_16 = Channel.empty()
} else {


process Adapter_Trimmer_Quality_Module_Quality_Filtering {

input:
 set val(name), file(reads) from g0_19_reads_g0_20
 val mate from g_1_mate_g0_20

output:
 set val(name), file("reads/*q")  into g0_20_reads_g_3
 file "*.{fastx,trimmomatic}_quality.log" optional true  into g0_20_log_file_g0_16

errorStrategy 'retry'

when:
(params.run_Quality_Filtering && (params.run_Quality_Filtering == "yes")) || !params.run_Quality_Filtering    

shell:
tool = params.Adapter_Trimmer_Quality_Module_Quality_Filtering.tool
phred = params.Adapter_Trimmer_Quality_Module_Quality_Filtering.phred
window_size = params.Adapter_Trimmer_Quality_Module_Quality_Filtering.window_size
required_quality_for_window_trimming = params.Adapter_Trimmer_Quality_Module_Quality_Filtering.required_quality_for_window_trimming
leading = params.Adapter_Trimmer_Quality_Module_Quality_Filtering.leading
trailing = params.Adapter_Trimmer_Quality_Module_Quality_Filtering.trailing
minlen = params.Adapter_Trimmer_Quality_Module_Quality_Filtering.minlen


// fastx parameters
minQuality = params.Adapter_Trimmer_Quality_Module_Quality_Filtering.minQuality
minPercent = params.Adapter_Trimmer_Quality_Module_Quality_Filtering.minPercent

remove_previous_reads = params.Adapter_Trimmer_Quality_Module_Quality_Filtering.remove_previous_reads

nameAll = reads.toString()
nameArray = nameAll.split(' ')
file2 ="";
if (nameAll.contains('.gz')) {
    newName =  nameArray[0] - ~/(\.fastq.gz)?(\.fq.gz)?$/
    file1 =  nameArray[0] - '.gz' 
    if (mate == "pair") {file2 =  nameArray[1] - '.gz'}
    runGzip = "ls *.gz | xargs -i echo gzip -df {} | sh"
} else {
    newName =  nameArray[0] - ~/(\.fastq)?(\.fq)?$/
    file1 =  nameArray[0]
    if (mate == "pair") {file2 =  nameArray[1]}
    runGzip = ''
}
'''
#!/usr/bin/env perl
 use List::Util qw[min max];
 use strict;
 use File::Basename;
 use Getopt::Long;
 use Pod::Usage; 
 use Cwd qw();
 
system("mkdir reads unpaired");
system("!{runGzip}");
my $param = "SLIDINGWINDOW:"."!{window_size}".":"."!{required_quality_for_window_trimming}";
$param.=" LEADING:"."!{leading}";
$param.=" TRAILING:"."!{trailing}";
$param.=" MINLEN:"."!{minlen}";

my $quality="!{phred}";

print "INFO: fastq quality: $quality\\n";
     
if ("!{tool}" eq "trimmomatic") {
    if ("!{mate}" eq "pair") {
        runCmd("trimmomatic PE -phred${quality} !{file1} !{file2} reads/!{name}.1.fastq unpaired/!{name}.1.fastq.unpaired reads/!{name}.2.fastq unpaired/!{name}.1.fastq.unpaired $param 2> !{name}.trimmomatic_quality.log");
    } else {
        runCmd("trimmomatic SE -phred${quality} !{file1} reads/!{name}.fastq $param 2> !{name}.trimmomatic_quality.log");
    }
} elsif ("!{tool}" eq "fastx") {
    if ("!{mate}" eq "pair") {
        print("WARNING: Fastx option is not suitable for paired reads. This step will be skipped.");
        system("mv !{file1} !{file2} reads/.");
    } else {
        runCmd("fastq_quality_filter  -Q $quality -q !{minQuality} -p !{minPercent} -v -i !{file1} -o reads/!{name}.fastq > !{name}.fastx_quality.log");
    }
}
if ("!{remove_previous_reads}" eq "true") {
    my $currpath = Cwd::cwd();
    my @paths = (split '/', $currpath);
    splice(@paths, -2);
    my $workdir= join '/', @paths;
    splice(@paths, -1);
    my $inputsdir= join '/', @paths;
    $inputsdir .= "/inputs";
    print "INFO: inputs reads will be removed if they are located in the workdir inputsdir\\n";
    my @listOfFiles = `readlink -e !{file1} !{file2}`;
    foreach my $targetFile (@listOfFiles){
        if (index($targetFile, $workdir) != -1 || index($targetFile, $inputsdir) != -1) {
            system("rm -f $targetFile");
            print "INFO: $targetFile deleted.\\n";
        }
    }
}

##Subroutines
sub runCmd {
    my ($com) = @_;
    my $error = system($com);
    if   ($error) { die "Command failed: $error $com\\n"; }
    else          { print "Command successful: $com\\n"; }
}


'''

}
}


params.ncRNA_star_index =  ""  //* @input
if (!((params.run_ncRNA_Removal && (params.run_ncRNA_Removal == "yes")) || !params.run_ncRNA_Removal)){
g0_20_reads_g_3.into{g_3_reads_g_4}
g_3_logOut_g_29 = Channel.empty()
} else {


process ncRNA_Removal {

input:
 set val(name), file(reads) from g0_20_reads_g_3
 val mate from g_1_mate_g_3

output:
 set val(name), file("*_no_ncRNA.fastq")  into g_3_reads_g_4
 set val(name), file("${newName}Log.final.out")  into g_3_logOut_g_29

errorStrategy 'retry'

when:
(params.run_ncRNA_Removal && (params.run_ncRNA_Removal == "yes")) || !params.run_ncRNA_Removal

script:
params_STAR = params.ncRNA_Removal.params_STAR
nameAll = reads.toString()
nameArray = nameAll.split(' ')

if (nameAll.contains('.gz')) {
    newName =  nameArray[0] - ~/(\.fastq.gz)?(\.fq.gz)?$/
    file =  nameAll - '.gz' - '.gz'
    runGzip = "ls *.gz | xargs -i echo gzip -df {} | sh"
} else {
    newName =  nameArray[0] - ~/(\.fastq)?(\.fq)?$/
    file =  nameAll 
    runGzip = ''
}

"""
$runGzip
STAR ${params_STAR}  --genomeDir ${params.ncRNA_star_index} --readFilesIn $file  --outReadsUnmapped Fastx  --outSAMtype BAM SortedByCoordinate     --outFileNamePrefix ${newName}
echo "Alignment completed."
mv ${newName}Unmapped.out.mate1 ${newName}_no_ncRNA.fastq
"""


}
}



process riboseq_ncRNA_Removal_Summary {

input:
 set val(name), file(alignSum) from g_3_logOut_g_29.groupTuple()

output:
 file "*.tsv"  into g_29_outputFile_g_28
 val "ncRNA_removal_sum"  into g_29_name_g_28

shell:
'''
#!/usr/bin/env perl
use List::Util qw[min max];
use strict;
use File::Basename;
use Getopt::Long;
use Pod::Usage; 
use Data::Dumper;

my %tsv;
my @headers = ();
my $name = "!{name}";

alteredAligned();

my @keys = keys %tsv;
my $summary = "$name"."_ncRNA_removal_sum.tsv";
my $header_string = join("\\t", @headers);
`echo "$header_string" > $summary`;
foreach my $key (@keys){
	my $values = join("\\t", @{ $tsv{$key} });
	`echo "$values" >> $summary`;
}


sub alteredAligned
{
	my @files = qw(!{alignSum});
	my $multimappedSum;
	my $alignedSum;
	my $inputCountSum;
	my $afterFiltering;
	push(@headers, "Sample");
    push(@headers, "Total Reads");
	push(@headers, "Reads After ncRNA Removal");
	foreach my $file (@files){
		my $multimapped;
		my $aligned;
		my $inputCount;
		chomp($inputCount = `cat $file | grep 'Number of input reads' | awk '{sum+=\\$6} END {print sum}'`);
		chomp($aligned = `cat $file | grep 'Uniquely mapped reads number' | awk '{sum+=\\$6} END {print sum}'`);
		chomp($multimapped = `cat $file | grep 'Number of reads mapped to multiple loci' | awk '{sum+=\\$9} END {print sum}'`);
		$multimappedSum += int($multimapped);
        $alignedSum += int($aligned);
        $inputCountSum += int($inputCount);
	}
	$afterFiltering = $inputCountSum - $alignedSum - $multimappedSum;
	$tsv{$name} = [$name, $inputCountSum];
	push(@{$tsv{$name}}, $afterFiltering);
}
'''

}


process Merge_TSV_Files_ncRNA_removal {

input:
 file tsv from g_29_outputFile_g_28.collect()
 val outputFileName from g_29_name_g_28.collect()

output:
 file "${name}.tsv"  into g_28_outputFileTSV_g_26

script:
name = outputFileName[0]
"""    
awk 'FNR==1 && NR!=1 {  getline; } 1 {print} ' *.tsv > ${name}.tsv
"""
}

params.star_index =  ""  //* @input

process STAR_align_Single_Best_Multimapper {

input:
 set val(name), file(reads) from g_3_reads_g_4
 val mate from g_1_mate_g_4

output:
 set val(name), file("${newName}.{bam,bam.bai}")  into g_4_mapped_reads_g_5
 set val(name), file("${newName}_transcriptome.bam")  into g_4_transcriptome_bam
 set val(name), file("${newName}Log.final.out")  into g_4_logOut_g_21

errorStrategy 'retry'

when:
(params.run_STAR && (params.run_STAR == "yes")) || !params.run_STAR

script:
    
params_STAR = params.STAR_align_Single_Best_Multimapper.params_STAR
nameAll = reads.toString()
nameArray = nameAll.split(' ')

if (nameAll.contains('.gz')) {
    newName =  nameArray[0] - ~/(\.fastq.gz)?(\.fq.gz)?$/
    file =  nameAll - '.gz' - '.gz'
    runGzip = "ls *.gz | xargs -i echo gzip -df {} | sh"
} else {
    newName =  nameArray[0] - ~/(\.fastq)?(\.fq)?$/
    file =  nameAll 
    runGzip = ''
}

"""
$runGzip
STAR ${params_STAR}  \
    --genomeDir ${params.star_index} \
    --readFilesIn $file \
    --outSAMtype BAM SortedByCoordinate \
	--outFileNamePrefix ${newName} \
	--quantMode TranscriptomeSAM \
	--outReadsUnmapped Fastx

echo "Alignment completed."

if [ -e "${newName}Aligned.sortedByCoord.out.bam" ] ; then
    mv ${newName}Aligned.sortedByCoord.out.bam ${newName}.bam
fi
samtools index ${newName}.bam

if [ -e "${newName}Aligned.toTranscriptome.out.bam" ] ; then
    mv ${newName}Aligned.toTranscriptome.out.bam ${newName}_transcriptome.bam
fi


"""


}


process STAR_Summary {

input:
 set val(name), file(alignSum) from g_4_logOut_g_21.groupTuple()

output:
 file "*.tsv"  into g_21_outputFile_g_22
 val "star_alignment_sum"  into g_21_name_g_22

shell:
'''
#!/usr/bin/env perl
use List::Util qw[min max];
use strict;
use File::Basename;
use Getopt::Long;
use Pod::Usage; 
use Data::Dumper;

my %tsv;
my @headers = ();
my $name = "!{name}";

alteredAligned();

my @keys = keys %tsv;
my $summary = "$name"."_star_sum.tsv";
my $header_string = join("\\t", @headers);
`echo "$header_string" > $summary`;
foreach my $key (@keys){
	my $values = join("\\t", @{ $tsv{$key} });
	`echo "$values" >> $summary`;
}


sub alteredAligned
{
	my @files = qw(!{alignSum});
	my $multimappedSum;
	my $alignedSum;
	my $inputCountSum;
	push(@headers, "Sample");
    push(@headers, "Total Reads");
	push(@headers, "Multimapped Reads Aligned (STAR)");
	push(@headers, "Unique Reads Aligned (STAR)");
	foreach my $file (@files){
		my $multimapped;
		my $aligned;
		my $inputCount;
		chomp($inputCount = `cat $file | grep 'Number of input reads' | awk '{sum+=\\$6} END {print sum}'`);
		chomp($aligned = `cat $file | grep 'Uniquely mapped reads number' | awk '{sum+=\\$6} END {print sum}'`);
		chomp($multimapped = `cat $file | grep 'Number of reads mapped to multiple loci' | awk '{sum+=\\$9} END {print sum}'`);
		$multimappedSum += int($multimapped);
        $alignedSum += int($aligned);
        $inputCountSum += int($inputCount);
	}
	$tsv{$name} = [$name, $inputCountSum];
	push(@{$tsv{$name}}, $multimappedSum);
	push(@{$tsv{$name}}, $alignedSum);
}
'''

}


process Merge_TSV_Files_star {

input:
 file tsv from g_21_outputFile_g_22.collect()
 val outputFileName from g_21_name_g_22.collect()

output:
 file "${name}.tsv"  into g_22_outputFileTSV_g_26

script:
name = outputFileName[0]
"""    
awk 'FNR==1 && NR!=1 {  getline; } 1 {print} ' *.tsv > ${name}.tsv
"""
}

params.single_transcript_gtf =  ""  //* @input
params.genome2bit =  ""  //* @input
params.makeavggene_main_path =  ""  //* @input
params.single_transcript_UTR_csv =  ""  //* @input
params.single_transcript_stopcodons_csv =  ""  //* @input
params.densebuilder_main_path =  ""  //* @input
params.riboseq_buildDenseTables_rpkm_path =  ""  //* @input
params.riboseq_buildDenseTables_rpkm_utr3adj_path =  ""  //* @input

thread = params.riboseq_Densebuilder.thread
control_sample_prefix = params.riboseq_Densebuilder.control_sample_prefix

//* autofill
if ($HOSTNAME == "default"){
    $CPU  = 10
    $MEMORY = 30
}
//* platform
if ($HOSTNAME == "ghpcc06.umassrc.org"){
    $TIME = 2500
    $CPU  = 10
    $MEMORY = 30
    $QUEUE = "long"
} 
//* platform
//* autofill

process riboseq_Densebuilder {

input:
 set val(name), file(bambai) from g_4_mapped_reads_g_5

output:
 set val(name), file("${dirName}")  into g_5_outputDir_g_9

when:
(params.run_riboseq_workflow && (params.run_riboseq_workflow == "yes")) || !params.run_riboseq_workflow

script:

bambai = bambai.toString()
bam = bambai.split(" ")[0]
dirName = name + "Densebuilder"
"""
#!/usr/bin/env python

import matplotlib
matplotlib.use('Agg') ### set backend
import matplotlib.pyplot as plt
plt.rcParams['pdf.fonttype'] = 42 # this keeps most text as actual text in PDFs, not outlines

import sys
import os
import argparse
import subprocess
import collections
from collections import OrderedDict
import pysam
import glob
import pandas as pd
import seaborn as sns
import numpy as np
from pathos.multiprocessing import ProcessingPool as Pool
import multiprocessing
import importlib

ftsize = [15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40]
#ftsize = [35]
### Footprint assignment
norm_type = "rpm" # either 'raw' for no normalization or 'rpm' for normalization to all mapped reads in alignment BAM file
threshold = '0'

flmin = 28
flmax = 35 # set to exact readsize that is desired
eAmin = 21
eAmax = 24
eEmin = 18
eEmax = 19

### Count tables:
cdsDenThresh = True
raw_dense_thresh = 0.01
rpm_dense_thresh = 0.001

# select insets to use to avoid start and stop peaks
inset_choice = 'default' # defaultInsets, zeroInsets, or customInsets
UTRfilestring = "${params.single_transcript_UTR_csv}"

def densebuilder_multi_raw(fiveorthreeprime, ribosome_site="0"):

    assignment = fiveorthreeprime
    totreads = 1E6
    ribosome_shift = "0"
    bampath = "${name}/${name}_star_default/"
    
    def densebuilder_one(readsize):
        readsize = str(readsize)
        riboshiftdict = "{"+readsize+':[0]'+"}"
        densitypath = "${name}/DensityUnnormalized/density"+assignment+"p_"+ribosome_shift+"shift_"+readsize+"/${name}_"+readsize+"f/"
        densityoutput = densitypath + "${name}_"+readsize+"f_"
        bamfileoutput = bampath + "${name}_"+readsize+"_match.bam"
        bamfilesorted = bampath + "${name}_"+readsize+"_match.sorted.bam"
        try:
            if not os.path.exists(densitypath): os.makedirs(densitypath)
        except OSError as err:
            print(err)
        try:
            if not os.path.exists(bampath): os.makedirs(bampath)
        except OSError as err:
            print(err)
        
        cmd = ("python ${params.densebuilder_main_path} --bamfileinput " + "${bam}" + 
            " --GTFfile " + "${params.single_transcript_gtf}" + 
            " --twobitfile " + "${params.genome2bit}" + 
            " --assignment " + assignment + 
            " --riboshiftdict " +  riboshiftdict +
            " --threshold " + str(threshold) + 
            " --totreads " + str(totreads) + 
            " --outputdata " + densityoutput +
            " --bamfileoutput " + bamfileoutput)
        print cmd
        os.system(cmd)
        ### sorting and indexing output bamfile
        sort_bamfile = 'samtools sort -o ' + bamfilesorted + " " + bamfileoutput 
        os.system(sort_bamfile)
        index_bamfile = 'samtools index %s' % (bamfilesorted)
        os.system(index_bamfile)
        clean_bamfile = 'rm %s' % (bamfileoutput)
        os.system(clean_bamfile)
    p = Pool(nodes=int(${thread}))
    p.map(densebuilder_one, ftsize)
    
    # now calculate the total number of reads from the bamfiles in ftsize
    # use pysam to calculate the number of reads present for all read lengths
    readcounter = 0
    for readlen in ftsize:
        bamfilesorted = bampath + "${name}_"+str(readlen)+"_match.sorted.bam"
        bamfile = pysam.AlignmentFile(bamfilesorted, 'rb')
        readnumber = bamfile.count()
        readcounter += readnumber
    
    readcount_outfile = "${name}/${name}_FPassigned_counts.txt"
    readcountout = open(readcount_outfile, "w")
    readcountout.write(str(readcounter))
    readcountout.close()

def densebuilder_multi_norm(self, ribosome_site="0"):
    assignment = "5" # 5' or 3' alignment
    ribosome_shift = "0"
    bampath = "${name}/${name}_star_default/"
    try:
        if not os.path.exists(bampath): os.makedirs(bampath)
    except OSError as err:
        print(err)
    
    totreads_countfile = "${name}/${name}_FPassigned_counts.txt"
    totreadcountf = open(totreads_countfile, "r")
    totreads = int(totreadcountf.read())
    totreadcountf.close()
    
    def densebuilder_one(readsize):
        readsize = str(readsize)
        riboshiftdict = "{"+readsize+':[0]'+"}"
        densitypath = "${name}/Density_rpm/density"+assignment+"p_"+ribosome_shift+"shift_"+readsize+"/${name}_"+readsize+"f/"
        densityoutput = densitypath + "${name}_"+readsize+"f_"
        bamfileoutput = bampath + "${name}_"+readsize+"_match.bam"
        bamfilesorted = bampath + "${name}_"+readsize+"_match.sorted.bam"
        try:
            if not os.path.exists(densitypath): os.makedirs(densitypath)
        except OSError as err:
            print(err)

        
        cmd = ("python ${params.densebuilder_main_path} --bamfileinput " + "${bam}" + 
            " --GTFfile " + "${params.single_transcript_gtf}" + 
            " --twobitfile " + "${params.genome2bit}" + 
            " --assignment " + assignment + 
            " --riboshiftdict " +  riboshiftdict +
            " --threshold " + str(threshold) + 
            " --totreads " + str(totreads) + 
            " --outputdata " + densityoutput +
            " --bamfileoutput " + bamfileoutput)
        print cmd
        os.system(cmd)
        ### sorting and indexing output bamfile
        sort_bamfile = 'samtools sort -o ' + bamfilesorted + " " + bamfileoutput 
        os.system(sort_bamfile)
        index_bamfile = 'samtools index %s' % (bamfilesorted)
        os.system(index_bamfile)
        clean_bamfile = 'rm %s' % (bamfileoutput)
        os.system(clean_bamfile)
    p = Pool(nodes=int(${thread}))
    p.map(densebuilder_one, ftsize)

def avggene_multi(alignposition, fiveorthreeprime):
    #Create avggenes for all readlengths defined by ftsize
	#Run all avggene scripts in parallel using pathos
	
	alignpos = alignposition # '1' for start, '2' for stop
	assignment = fiveorthreeprime # '5' for fiveprime mapping, '3' for threeprime mapping
	norm = "uneq" # 'uneq' for no normalizaiton; 'eq' to give equal weight to all genes
	ribosome_shift = "0" 
	read_length_list = ftsize
	
	if norm == 'uneq':
		equalweight = '0'
	if norm == 'eq':
		equalweight = '1'
		
	filtermodule = str(0)
	exclusionmodule = str(0)
	thresh = str(0)  # rpkm threshold
	if alignpos == '1':
		regionlength5 = '50'
		regionlength3 = '150'
		cdslength = regionlength3
	if alignpos == '2':
		regionlength5 = '150'
		regionlength3 = '100'
		cdslength = regionlength5
		
	if norm_type == "raw":
		densitystring = "DensityUnnormalized"
	elif norm_type == "rpm":
		densitystring = "Density_rpm"
	
	outfolder = "${name}/avggene"+alignpos+"_ORF"+norm_type+"_"+ribosome_shift+"shift_"+assignment+norm+cdslength
	if not os.path.exists(outfolder):   os.makedirs(outfolder)
	
	def avggene_onesize(readin):
	    # Calculate average gene for one read size
		readsize = str(readin) # set read size
		trspdictfilestring = '%s/%s/density%sp_%sshift_%s/%s_%sf/%s_%sf_' %("${name}", densitystring, assignment, ribosome_shift, readsize, "${name}", readsize, "${name}", readsize)
		outfilebase = "%s/%s_%sf" % (outfolder, "${name}", readsize)
		commandstring = 'python ${params.makeavggene_main_path} --regionlength5 %s --regionlength3 %s --trspdictfilestring %s --UTRfilestring %s --filtermodule %s --exclusionmodule %s --threshold %s --alignpos %s --equalweight %s --outfilebase %s' %(regionlength5,regionlength3,trspdictfilestring,UTRfilestring,filtermodule,exclusionmodule,thresh,alignpos,equalweight,outfilebase)
		print commandstring
		os.system(commandstring)

	p = Pool(nodes=int(${thread}))
	p.map(avggene_onesize, read_length_list)

def get_riboshift(shiftsite, alignposition):
    alignpos = alignposition
    assignment = "5"  # must use 5' start mapping for this
    norm = "uneq"
    ribosome_shift = "0"
    shiftpos = shiftsite
    if shiftpos == "P":
        shift_adjust = 0
    if shiftpos == "A":
        shift_adjust = 3 # A site is shifted 3 nt's farther than P site from 5' end
    if shiftpos == "E":
        shift_adjust = -3
    cdslength = "150"
    avggene_paths = "${name}/avggene"+alignpos+"_ORF"+norm_type+"_"+ribosome_shift+"shift_"+assignment+norm+cdslength
    frame = pd.DataFrame()
    list_ = []
    for ft in ftsize:
		infile = "%s/%s_%sf_avg_1.csv" % (avggene_paths, "${name}", ft)
		print "infile:" + infile
		df = pd.read_csv(infile, index_col="position", header=0)
		df.columns = ["%d" % (ft)]
		list_.append(df)
    frame = pd.concat(list_, axis=1)
	# Get the max (or modal read number per position) of footprints relative to start codon
	# Adding a +1 to account of +/- 1nt from digestion heterogeneity
	# Only look between 11 to 16 nts upstream of start for 5' aligned reads
	# Define full length reads from 31-34 nt's in length
	# retrieve these using 31-15 == 16, and 34-15 == 19+1, thus 16:20
    offset_fl = 50+1+shift_adjust-frame.iloc[34:39,flmin-15:flmax+1-15].idxmax()
    print "offset_fl"
    print offset_fl
	# Define emptyA reads from 21-24 nt's in length
	# retrieve these using 21-15 == 6, and 24-15 == 9+1, thus 6:10
    offset_emptyA = 50+1+shift_adjust-frame.iloc[34:39,eAmin-15:eAmax+1-15].idxmax()
    print "offset_emptyA"
    print offset_emptyA
	# Define emptyE reads from 18-19 nt's in length - need to cut close to avoid empty A
	# retrieve these using 18-15 == 3, and 19-15 == 4, thus 3:5
	# also need to only look around the start codon for these reads since the 5' end is directly over the 'A' of the AUG
    offset_emptyE = 50+1+shift_adjust-frame.iloc[49:52,eEmin-15:eEmax+1-15].idxmax()
    print "offset_emptyE"
    print offset_emptyE
    outpath = "${name}/riboshift"
    if not os.path.exists(outpath):     os.makedirs(outpath)
    outfile_fl = "%s/full_length_%soffsets.csv" % (outpath, shiftpos)
    outfile_emptyA = "%s/emptyA_%soffsets.csv" % (outpath, shiftpos)
    outfile_emptyE = "%s/emptyE_%soffsets.csv" % (outpath, shiftpos)
    offset_fl.to_csv(outfile_fl, index = True, header = False)
    offset_emptyA.to_csv(outfile_emptyA, index = True, header = False)
    offset_emptyE.to_csv(outfile_emptyE, index = True, header = False)

def densebuilder_riboshift(ribosome_site):
    #Function to map densities for ranges of read lengths to single sites of the ribosome
    #Normally map to the A site of the ribosome
    assignment = "5" # should be five prime assignment to use shift dicts
    ribosome_shift = ribosome_site # 'A' or 'P'
    bamfileinput = "${bam}"
    bamfileoutpath = "${name}/${name}_star_default/"
    if not os.path.exists(bamfileoutpath): os.makedirs(bamfileoutpath)
    densitystring = "Density_rpm"
    totreads_countfile = "${name}/${name}_FPassigned_counts.txt"
    totreadcountf = open(totreads_countfile, "r")
    totreads = int(totreadcountf.read())
    totreadcountf.close()
    ## Set paths to riboshift dicts for each species of interest
    fl_infile = "${name}/riboshift/full_length_%soffsets.csv" % (ribosome_shift)
    eA_infile = "${name}/riboshift/emptyA_%soffsets.csv" % (ribosome_shift)
    eE_infile = "${name}/riboshift/emptyE_%soffsets.csv" % (ribosome_shift)
    fl_ribo_shift = pd.read_csv(fl_infile, index_col = 0, header=None)
    fl_ribo_shift.columns = ['riboshift']
    eA_ribo_shift = pd.read_csv(eA_infile, index_col = 0, header=None)
    eA_ribo_shift.columns = ['riboshift']
    eE_ribo_shift = pd.read_csv(eE_infile, index_col = 0, header=None)
    eE_ribo_shift.columns = ['riboshift']
    fl_sizes = fl_ribo_shift.index.tolist()
    eA_sizes = eA_ribo_shift.index.tolist()
    eE_sizes = eE_ribo_shift.index.tolist()
    fl_shiftdict = collections.OrderedDict()
    eA_shiftdict = collections.OrderedDict()
    eE_shiftdict = collections.OrderedDict()
    for i in fl_sizes:
        fl_shiftdict.update( {str(i):'['+str(fl_ribo_shift.loc[i,'riboshift'])+']'})
    for i in eA_sizes:
        eA_shiftdict.update( {str(i):'['+str(eA_ribo_shift.loc[i,'riboshift'])+']'})
    for i in eE_sizes:
        eE_shiftdict.update( {str(i):'['+str(eE_ribo_shift.loc[i,'riboshift'])+']'})
    allLenghts_shiftdict = collections.OrderedDict()
    allLenghts_shiftdict.update(eE_shiftdict)
    allLenghts_shiftdict.update(eA_shiftdict)
    allLenghts_shiftdict.update(fl_shiftdict)
    fp_assign_path = "${name}"
    
    def build_fl_dense(fldict):
	    readsize_dict = fldict
	    pop = "fl"
	    minlen = str(readsize_dict.keys()[0])
	    maxlen = str(readsize_dict.keys()[-1])
	    dict_conv = dict(readsize_dict)
	    # fl_shift:
	    riboshiftdict = str(dict_conv).replace(" ","")
	    # riboshiftdict = {'33':'17','32':'17','31':'16','34':'17'} ## manually set readlengths
	    bamfileoutput = '%s/%s_%s_match.bam' % (bamfileoutpath, "${name}", pop)
	    densityfilepath = '%s/%s/density%sp_%sshift_%sto%s/%s_%sto%sf' %("${name}",densitystring, assignment, ribosome_shift, minlen, maxlen, "${name}", minlen,maxlen)
	    densityfileout = '%s/%s_%sto%sf_' % (densityfilepath, "${name}", minlen, maxlen)
	    if not os.path.exists(densityfilepath): os.makedirs(densityfilepath)
	    commandstring = 'python ${params.densebuilder_main_path} --bamfileinput %s --GTFfile %s --twobitfile %s --assignment %s --riboshiftdict %s --threshold %s --totreads %s --outputdata %s --bamfileoutput %s' % (bamfileinput, "${params.single_transcript_gtf}", "${params.genome2bit}", assignment,riboshiftdict, threshold, totreads, densityfileout, bamfileoutput)
	    print commandstring
	    os.system(commandstring)
	    ## sort and index bamfile
	    sort_bamfile = 'samtools sort -o %s/%s_%s_match.sorted.bam %s' % (bamfileoutpath, "${name}", pop, bamfileoutput)
	    os.system(sort_bamfile)
	    index_bamfile = 'samtools index %s/%s_%s_match.sorted.bam' % (bamfileoutpath,"${name}", pop)
	    os.system(index_bamfile)
	    clean_bamfile = 'rm %s' % (bamfileoutput)
	    os.system(clean_bamfile)
    def build_eA_dense(eAdict):
        readsize_dict = eAdict
        pop = "eA"
        minlen = str(readsize_dict.keys()[0])
        maxlen = str(readsize_dict.keys()[-1])
        dict_conv = dict(readsize_dict)
        # eA_shift:
        riboshiftdict = str(dict_conv).replace(" ","")
        # riboshiftdict = {'33':'17','32':'17','31':'16','34':'17'}
        bamfileoutpath = '%s/%s_star_default' % ("${name}", "${name}")
        bamfileoutput = '%s/%s_%s_match.bam' % (bamfileoutpath, "${name}", pop)
        densityfilepath = '%s/%s/density%sp_%sshift_%sto%s/%s_%sto%sf' %("${name}",densitystring, assignment, ribosome_shift, minlen, maxlen, "${name}", minlen,maxlen)
        densityfileout = '%s/%s_%sto%sf_' % (densityfilepath, "${name}", minlen, maxlen)
        if not os.path.exists(densityfilepath): os.makedirs(densityfilepath)
        commandstring = 'python ${params.densebuilder_main_path} --bamfileinput %s --GTFfile %s --twobitfile %s --assignment %s --riboshiftdict %s --threshold %s --totreads %s --outputdata %s --bamfileoutput %s' % (bamfileinput, "${params.single_transcript_gtf}", "${params.genome2bit}", assignment,riboshiftdict, threshold, totreads, densityfileout, bamfileoutput)
        print commandstring
        os.system(commandstring)
        ## sort and index bamfile
        sort_bamfile = 'samtools sort -o %s/%s_%s_match.sorted.bam %s' % (bamfileoutpath, "${name}", pop, bamfileoutput)
        os.system(sort_bamfile)
        index_bamfile = 'samtools index %s/%s_%s_match.sorted.bam' % (bamfileoutpath,"${name}", pop)
        os.system(index_bamfile)
        clean_bamfile = 'rm %s' % (bamfileoutput)
        os.system(clean_bamfile)
    def build_eE_dense(eEdict):
        readsize_dict = eEdict
        pop = "eE"
        minlen = str(readsize_dict.keys()[0])
        maxlen = str(readsize_dict.keys()[-1])
        dict_conv = dict(readsize_dict)
        # eE_shift:
        riboshiftdict = str(dict_conv).replace(" ","")
        # riboshiftdict = {'33':'17','32':'17','31':'16','34':'17'}
        bamfileoutpath = '%s/%s_star_default' % (fp_assign_path, "${name}")
        bamfileoutput = '%s/%s_%s_match.bam' % (bamfileoutpath, "${name}", pop)
        densityfilepath = '%s/%s/density%sp_%sshift_%sto%s/%s_%sto%sf' %(fp_assign_path,densitystring, assignment, ribosome_shift, minlen, maxlen, "${name}", minlen,maxlen)
        densityfileout = '%s/%s_%sto%sf_' % (densityfilepath, "${name}", minlen, maxlen)
        if not os.path.exists(densityfilepath): os.makedirs(densityfilepath)
        commandstring = 'python ${params.densebuilder_main_path} --bamfileinput %s --GTFfile %s --twobitfile %s --assignment %s --riboshiftdict %s --threshold %s --totreads %s --outputdata %s --bamfileoutput %s' % (bamfileinput, "${params.single_transcript_gtf}", "${params.genome2bit}", assignment,riboshiftdict, threshold, totreads, densityfileout, bamfileoutput)
        print commandstring
        os.system(commandstring)
        ## sort and index bamfile
        sort_bamfile = 'samtools sort -o %s/%s_%s_match.sorted.bam %s' % (bamfileoutpath, "${name}", pop, bamfileoutput)
        os.system(sort_bamfile)
        index_bamfile = 'samtools index %s/%s_%s_match.sorted.bam' % (bamfileoutpath,"${name}", pop)
        os.system(index_bamfile)
        clean_bamfile = 'rm %s' % (bamfileoutput)
        os.system(clean_bamfile)
    def build_allLengths_dense(alldict):
        readsize_dict = alldict
        pop = "aL"
        minlen = str(readsize_dict.keys()[0])
        maxlen = str(readsize_dict.keys()[-1])
        dict_conv = dict(readsize_dict)
        # aL_shift:
        riboshiftdict = str(dict_conv).replace(" ","")
        # riboshiftdict = {'33':'17','32':'17','31':'16','34':'17'}
        bamfileoutpath = '%s/%s_star_default' % (fp_assign_path, "${name}")
        bamfileoutput = '%s/%s_%s_match.bam' % (bamfileoutpath, "${name}", pop)
        densityfilepath = '%s/%s/density%sp_%sshift_%s/%s_%sf' %(fp_assign_path,densitystring, assignment, ribosome_shift, pop, "${name}", pop)
        densityfileout = '%s/%s_%sf_' % (densityfilepath, "${name}", pop)
        if not os.path.exists(densityfilepath): os.makedirs(densityfilepath)
        commandstring = 'python ${params.densebuilder_main_path} --bamfileinput %s --GTFfile %s --twobitfile %s --assignment %s --riboshiftdict %s --threshold %s --totreads %s --outputdata %s --bamfileoutput %s' % (bamfileinput, "${params.single_transcript_gtf}", "${params.genome2bit}", assignment,riboshiftdict, threshold, totreads, densityfileout, bamfileoutput)
        print commandstring
        os.system(commandstring)
        ## sort and index bamfile
        sort_bamfile = 'samtools sort -o %s/%s_%s_match.sorted.bam %s' % (bamfileoutpath, "${name}", pop, bamfileoutput)
        os.system(sort_bamfile)
        index_bamfile = 'samtools index %s/%s_%s_match.sorted.bam' % (bamfileoutpath,"${name}", pop)
        os.system(index_bamfile)
        clean_bamfile = 'rm %s' % (bamfileoutput)
        os.system(clean_bamfile)
    p1 = multiprocessing.Process(target=build_fl_dense, args=(fl_shiftdict,))
    p1.start()
    p2 = multiprocessing.Process(target=build_eA_dense, args=(eA_shiftdict,))
    p2.start()
    p3 = multiprocessing.Process(target=build_eE_dense, args=(eE_shiftdict,))
    p3.start()
    p4 = multiprocessing.Process(target=build_allLengths_dense, args=(allLenghts_shiftdict,))
    p4.start()
    p1.join()
    p2.join()
    p3.join()
    p4.join()

def avggene_riboshift(alignposition, ribosome_site, normalization, thresh):
    #Create avggenes for all readlengths defined by ftsize
    #Run all avggene scripts in parallel using pathos
    alignpos = alignposition # '1' for start, '2' for stop
    assignment = "5" # should be 5' mapped at this point
    norm = normalization # 'uneq' for no normalizaiton; 'eq' to give equal weight to allgenes
    ribosome_shift = ribosome_site # 'A' or 'P' or maybe 'E' later
    
    if norm == 'uneq':
        equalweight = '0'
    if norm == 'eq':
        equalweight = '1' 	#equalweight = '0'    # 1 gives equal weight to all genes.
    filtermodule = str(0)
    exclusionmodule = str(0)
    if alignpos == '1':
        regionlength5 = '50'
        regionlength3 = '150'
        cdslength = regionlength3
    if alignpos == '2':
        regionlength5 = '150'
        regionlength3 = '100'
        cdslength = regionlength5
    if norm_type == "raw":
        densitystring = "DensityUnnormalized"
    elif norm_type == "rpm":
        densitystring = "Density_rpm"
    else:
        print "Normalization is not set for densebuilder!"
        sys.exit()
    fp_assign_path = "${name}"
    
    outfolder =  "%s/avggene%s_ORF%s_%sshift_%s%s%s" % (fp_assign_path, alignpos, norm_type, ribosome_shift, assignment, norm, cdslength) # norm should be 'uneq' fornow
    if not os.path.exists(outfolder):   os.makedirs(outfolder)
    ### define read sizes for fl, eA, and eE - need minlen and maxlen
    fl_infile = "%s/riboshift/full_length_%soffsets.csv" % (fp_assign_path,ribosome_shift)
    eA_infile = "%s/riboshift/emptyA_%soffsets.csv" % (fp_assign_path, ribosome_shift)
    eE_infile = "%s/riboshift/emptyE_%soffsets.csv" % (fp_assign_path, ribosome_shift)
    fl_ribo_shift = pd.read_csv(fl_infile, index_col = 0, header=None)
    fl_ribo_shift.columns = ['riboshift']
    eA_ribo_shift = pd.read_csv(eA_infile, index_col = 0, header=None)
    eA_ribo_shift.columns = ['riboshift']
    eE_ribo_shift = pd.read_csv(eE_infile, index_col = 0, header=None)
    eE_ribo_shift.columns = ['riboshift']
    # print eE_ribo_shift.index.tolist()
    fl_sizes = fl_ribo_shift.index.tolist()
    eA_sizes = eA_ribo_shift.index.tolist()
    eE_sizes = eE_ribo_shift.index.tolist()
    minlen_fl = fl_sizes[0]
    maxlen_fl = fl_sizes[-1]
    minlen_eA = eA_sizes[0]
    maxlen_eA = eA_sizes[-1]
    minlen_eE = eE_sizes[0]
    maxlen_eE = eE_sizes[-1]
    def avggene_fl():
        # Calculate average gene for one read size
        # readsize = str(readin)
        # change this to read in fl reads
        pop = "fl"
        minlen = minlen_fl
        maxlen = maxlen_fl
        trspdictfilestring = '%s/%s/density%sp_%sshift_%sto%s/%s_%sto%sf/%s_%sto%sf_' % (fp_assign_path, densitystring, assignment, ribosome_shift, minlen, maxlen, "${name}", minlen, maxlen, "${name}", minlen, maxlen)
        outfilebase = "%s/%s_%s_rpkmThresh%s_%sto%sf" % (outfolder, "${name}", pop,thresh, minlen, maxlen)
        print ""
        commandstring = 'python ${params.makeavggene_main_path} --regionlength5 %s --regionlength3 %s --trspdictfilestring %s --UTRfilestring %s --filtermodule %s --exclusionmodule %s --threshold %s --alignpos %s --equalweight %s --outfilebase %s' %(regionlength5,regionlength3,trspdictfilestring,UTRfilestring,filtermodule,exclusionmodule,thresh,alignpos,equalweight,outfilebase)
        print commandstring
        os.system(commandstring)
    
    def avggene_eA():
        # Calculate average gene for one read size
        # readsize = str(readin)
        # change this to read in fl reads
        pop = "eA"
        minlen = minlen_eA
        maxlen = maxlen_eA
        trspdictfilestring = '%s/%s/density%sp_%sshift_%sto%s/%s_%sto%sf/%s_%sto%sf_' % (fp_assign_path, densitystring, assignment, ribosome_shift, minlen, maxlen, "${name}", minlen, maxlen, "${name}", minlen, maxlen)
        outfilebase = "%s/%s_%s_rpkmThresh%s_%sto%sf" % (outfolder, "${name}", pop,thresh, minlen, maxlen)
        print ""
        commandstring = 'python ${params.makeavggene_main_path} --regionlength5 %s --regionlength3 %s --trspdictfilestring %s --UTRfilestring %s --filtermodule %s --exclusionmodule %s --threshold %s --alignpos %s --equalweight %s --outfilebase %s' %(regionlength5,regionlength3,trspdictfilestring,UTRfilestring,filtermodule,exclusionmodule,thresh,alignpos,equalweight,outfilebase)
        print commandstring
        os.system(commandstring)
    
    def avggene_eE():
        # Calculate average gene for one read size
        # readsize = str(readin)
        # change this to read in fl reads
        pop = "eE"
        minlen = minlen_eE
        maxlen = maxlen_eE
        trspdictfilestring = '%s/%s/density%sp_%sshift_%sto%s/%s_%sto%sf/%s_%sto%sf_' % (fp_assign_path, densitystring, assignment, ribosome_shift, minlen, maxlen, "${name}", minlen, maxlen, "${name}", minlen, maxlen)
        outfilebase = "%s/%s_%s_rpkmThresh%s_%sto%sf" % (outfolder, "${name}", pop,thresh, minlen, maxlen)
        commandstring = 'python ${params.makeavggene_main_path} --regionlength5 %s --regionlength3 %s --trspdictfilestring %s --UTRfilestring %s --filtermodule %s --exclusionmodule %s --threshold %s --alignpos %s --equalweight %s --outfilebase %s' %(regionlength5,regionlength3,trspdictfilestring,UTRfilestring,filtermodule,exclusionmodule,thresh,alignpos,equalweight,outfilebase)
        print commandstring
        os.system(commandstring)
    def avggene_aL():
        pop = "aL"
        # minlen = minlen_eE
        # maxlen = maxlen_eE
        trspdictfilestring = '%s/%s/density%sp_%sshift_%s/%s_%sf/%s_%sf_' % (fp_assign_path, densitystring, assignment, ribosome_shift, pop, "${name}", pop,"${name}", pop)
        outfilebase = "%s/%s_rpkmThresh%s_%s_%sf" % (outfolder, "${name}", thresh, pop,pop)
        commandstring = 'python ${params.makeavggene_main_path} --regionlength5 %s --regionlength3 %s --trspdictfilestring %s --UTRfilestring %s --filtermodule %s --exclusionmodule %s --threshold %s --alignpos %s --equalweight %s --outfilebase %s' %(regionlength5,regionlength3,trspdictfilestring,UTRfilestring, filtermodule,exclusionmodule,thresh,alignpos,equalweight,outfilebase)
        print commandstring
        os.system(commandstring)
    p1 = multiprocessing.Process(target=avggene_fl)
    p1.start()
    p2 = multiprocessing.Process(target=avggene_eA)
    p2.start()
    p3 = multiprocessing.Process(target=avggene_eE)
    p3.start()
    p4 = multiprocessing.Process(target=avggene_aL)
    p4.start()
    p1.join()
    p2.join()
    p3.join()
    p4.join()

def build_count_tables(alignposition="1", fiveorthreeprime="5", normalization="uneq", ribosome_site="A"):
    # wrapper to direct building of density tables
    # count tables are based on cds, utr5, and utr3 sizes defined by UTRfilesting
	# inset_choice contains adjustments to avoid start and stop codon peaks in data tables
	# output is a csv file with densities for each region and 3'UTR occupancies

    # mappingpos = "5"
	alignpos = alignposition # '1' for start, '2' for stop
	assignment = fiveorthreeprime # '5' for fiveprime mapping, '3' for threeprime mapping
	norm = normalization # 'uneq' for no normalizaiton; 'eq' to give equal weight to allgenes
	ribosome_shift = ribosome_site # 'A' or 'P' or maybe 'E' later
	flminlen = str(flmin)
	flmaxlen = str(flmax)
	eAminlen = str(eAmin)
	eAmaxlen = str(eAmax)
	eEminlen = str(eEmin)
	eEmaxlen = str(eEmax)
	fp_assign_path = "${name}"
	outfolder =  "%s/countTables" % (fp_assign_path) # make a new folder to store counttables
	if not os.path.exists(outfolder):   os.makedirs(outfolder)
	
	### require Density_rpm for densityfilestring
	densitystring = "Density_rpm"
	totreads_countfile = "%s/%s_FPassigned_counts.txt" % (fp_assign_path, "${name}")
	totreadcountf = open(totreads_countfile, "r")
	totreads = int(totreadcountf.read())
	totreadcountf.close()
	def fl_countTable():
		trspdictfilestring = '%s/%s/density%sp_%sshift_%sto%s/%s_%sto%sf/%s_%sto%sf_' %(
				fp_assign_path, densitystring, assignment, ribosome_shift, flminlen,flmaxlen, "${name}", flminlen, flmaxlen, "${name}", flminlen, flmaxlen)
		outfilestring = "%s/%s_fl_%s_%sto%s_countTable_rpkm" % (outfolder, "${name}",norm_type, flminlen, flmaxlen)
		print ""
		# need to input trspdictfilestring, UTRfilesting, cdsDenThresh,raw_dense_thresh, rpm_dense_thresh, norm_type, inset_choice
		commandstring = 'python ${params.riboseq_buildDenseTables_rpkm_path} --trspdictfilestring %s --UTRfilestring %s --cdsDenThresh %s --norm_type %s --raw_dense_thresh %s --rpm_dense_thresh %s --inset_choice %s --outfilestring %s --totreads %s' % (
			trspdictfilestring,UTRfilestring,cdsDenThresh,norm_type,raw_dense_thresh, rpm_dense_thresh,inset_choice,outfilestring, totreads)
		print commandstring
		os.system(commandstring)
	def eA_countTable():
		trspdictfilestring = '%s/%s/density%sp_%sshift_%sto%s/%s_%sto%sf/%s_%sto%sf_' %(
				fp_assign_path, densitystring, assignment, ribosome_shift, eAminlen,eAmaxlen, "${name}", eAminlen, eAmaxlen, "${name}", eAminlen, eAmaxlen)
		outfilestring = "%s/%s_eA_%s_%sto%s_countTable_rpkm" % (outfolder, "${name}",norm_type, eAminlen, eAmaxlen)
		print ""
		# need to input trspdictfilestring, UTRfilesting, cdsDenThresh,raw_dense_thresh, rpm_dense_thresh, norm_type, inset_choice
		commandstring = 'python ${params.riboseq_buildDenseTables_rpkm_path} --trspdictfilestring %s --UTRfilestring %s --cdsDenThresh %s --norm_type %s --raw_dense_thresh %s --rpm_dense_thresh %s --inset_choice %s --outfilestring %s --totreads %s' % (
			trspdictfilestring,UTRfilestring,cdsDenThresh,norm_type,raw_dense_thresh,rpm_dense_thresh,inset_choice,outfilestring, totreads)
		print commandstring
		os.system(commandstring)
	def eE_countTable():
		trspdictfilestring = '%s/%s/density%sp_%sshift_%sto%s/%s_%sto%sf/%s_%sto%sf_' %(
				fp_assign_path, densitystring, assignment, ribosome_shift, eEminlen,eEmaxlen, "${name}", eEminlen, eEmaxlen, "${name}", eEminlen, eEmaxlen)
		outfilestring = "%s/%s_eE_%s_%sto%s_countTable_rpkm" % (outfolder, "${name}",norm_type, eEminlen, eEmaxlen)
		print ""
		# need to input trspdictfilestring, UTRfilesting, cdsDenThresh,raw_dense_thresh, rpm_dense_thresh, norm_type, inset_choice
		commandstring = 'python ${params.riboseq_buildDenseTables_rpkm_path} --trspdictfilestring %s --UTRfilestring %s --cdsDenThresh %s --norm_type %s --raw_dense_thresh %s --rpm_dense_thresh %s --inset_choice %s --outfilestring %s --totreads %s' % (
			trspdictfilestring,UTRfilestring,cdsDenThresh,norm_type,raw_dense_thresh,rpm_dense_thresh,inset_choice,outfilestring, totreads)
		print commandstring
		os.system(commandstring)
	def aL_countTable():
		trspdictfilestring = '%s/%s/density%sp_%sshift_aL/%s_aLf/%s_aLf_' %(
				fp_assign_path, densitystring, assignment, ribosome_shift, "${name}","${name}",)
		outfilestring = "%s/%s_aL_%s_countTable_rpkm" % (outfolder, "${name}", norm_type)
		print ""
		# need to input trspdictfilestring, UTRfilesting, cdsDenThresh,raw_dense_thresh, rpm_dense_thresh, norm_type, inset_choice
		commandstring = 'python ${params.riboseq_buildDenseTables_rpkm_path} --trspdictfilestring %s --UTRfilestring %s --cdsDenThresh %s --norm_type %s --raw_dense_thresh %s --rpm_dense_thresh %s --inset_choice %s --outfilestring %s --totreads %s' % (
			trspdictfilestring,UTRfilestring,cdsDenThresh,norm_type,raw_dense_thresh,rpm_dense_thresh,inset_choice,outfilestring, totreads)
		print commandstring
		os.system(commandstring)
	p1 = multiprocessing.Process(target=fl_countTable)
	p1.start()
	p2 = multiprocessing.Process(target=eA_countTable)
	p2.start()
	p3 = multiprocessing.Process(target=eE_countTable)
	p3.start()
	p4 = multiprocessing.Process(target=aL_countTable)
	p4.start()
	p1.join()
	p2.join()
	p3.join()
	p4.join()

def build_count_tables_utr3adj(alignposition="1", fiveorthreeprime="5", normalization="uneq", ribosome_site="A"):
	# wrapper to direct building of density tables
	# count tables are based on cds, utr5, and utr3 sizes defined by UTRfilesting
	# inset_choice contains adjustments to avoid start and stop codon peaks in data tables
	# output is a csv file with densities for each region and 3'UTR occupancies
	# mappingpos = "5"
	stopcodons = "${params.single_transcript_stopcodons_csv}"
	alignpos = alignposition # '1' for start, '2' for stop
	assignment = fiveorthreeprime # '5' for fiveprime mapping, '3' for threeprime mapping
	norm = normalization # 'uneq' for no normalizaiton; 'eq' to give equal weight to allgenes
	ribosome_shift = ribosome_site # 'A' or 'P' or maybe 'E' later
	flminlen = str(flmin)
	flmaxlen = str(flmax)
	eAminlen = str(eAmin)
	eAmaxlen = str(eAmax)
	eEminlen = str(eEmin)
	eEmaxlen = str(eEmax)
	fp_assign_path = '${name}'
	outfolder =  "%s/countTables" % (fp_assign_path) # make a new folder to store counttables
	if not os.path.exists(outfolder):   os.makedirs(outfolder)
	
	### require Density_rpm for densityfilestring
	densitystring = "Density_rpm"
	totreads_countfile = "%s/%s_FPassigned_counts.txt" % (fp_assign_path, '${name}')
	totreadcountf = open(totreads_countfile, "r")
	totreads = int(totreadcountf.read())
	totreadcountf.close()
	def fl_countTable():
		trspdictfilestring = '%s/%s/density%sp_%sshift_%sto%s/%s_%sto%sf/%s_%sto%sf_' %(
				fp_assign_path, densitystring, assignment, ribosome_shift, flminlen,flmaxlen, '${name}', flminlen, flmaxlen, '${name}', flminlen, flmaxlen)
		outfilestring = "%s/%s_fl_%s_%sto%s_countTable_rpkm_utr3adj" % (outfolder, "${name}", norm_type, flminlen, flmaxlen)
		print ""
		commandstring = 'python ${params.riboseq_buildDenseTables_rpkm_utr3adj_path} --trspdictfilestring %s --UTRfilestring %s --cdsDenThresh %s --norm_type %s --raw_dense_thresh %s --rpm_dense_thresh %s --inset_choice %s --outfilestring %s --totreads %s --stopcodons %s' % (
			trspdictfilestring,UTRfilestring,cdsDenThresh,norm_type,raw_dense_thresh,rpm_dense_thresh,inset_choice,outfilestring, totreads, stopcodons)
		print commandstring
		os.system(commandstring)
	def eA_countTable():
		trspdictfilestring = '%s/%s/density%sp_%sshift_%sto%s/%s_%sto%sf/%s_%sto%sf_' %(fp_assign_path, densitystring, assignment, ribosome_shift, eAminlen,eAmaxlen, '${name}', eAminlen, eAmaxlen, '${name}', eAminlen, eAmaxlen)
		outfilestring = "%s/%s_eA_%s_%sto%s_countTable_rpkm_utr3adj" % (outfolder, "${name}", norm_type, eAminlen, eAmaxlen)
		print ""
		# need to input trspdictfilestring, UTRfilesting, cdsDenThresh,raw_dense_thresh, rpm_dense_thresh, norm_type, inset_choice
		commandstring = 'python ${params.riboseq_buildDenseTables_rpkm_utr3adj_path}  --trspdictfilestring %s --UTRfilestring %s --cdsDenThresh %s --norm_type %s --raw_dense_thresh %s --rpm_dense_thresh %s --inset_choice %s --outfilestring %s --totreads %s --stopcodons %s' % (trspdictfilestring,UTRfilestring,cdsDenThresh,norm_type,raw_dense_thresh,rpm_dense_thresh,inset_choice,outfilestring, totreads, stopcodons)
		print commandstring
		os.system(commandstring)
	def eE_countTable():
		trspdictfilestring = '%s/%s/density%sp_%sshift_%sto%s/%s_%sto%sf/%s_%sto%sf_' %(fp_assign_path, densitystring, assignment, ribosome_shift, eEminlen,eEmaxlen, '${name}', eEminlen, eEmaxlen, '${name}', eEminlen, eEmaxlen)
		outfilestring = "%s/%s_eE_%s_%sto%s_countTable_rpkm_utr3adj" % (outfolder, "${name}", norm_type, eEminlen, eEmaxlen)
		print ""
		# need to input trspdictfilestring, UTRfilesting, cdsDenThresh,raw_dense_thresh, rpm_dense_thresh, norm_type, inset_choice
		commandstring = 'python ${params.riboseq_buildDenseTables_rpkm_utr3adj_path} --trspdictfilestring %s --UTRfilestring %s --cdsDenThresh %s --norm_type %s --raw_dense_thresh %s --rpm_dense_thresh %s --inset_choice %s --outfilestring %s --totreads %s --stopcodons %s' % (trspdictfilestring,UTRfilestring,cdsDenThresh,norm_type,raw_dense_thresh,rpm_dense_thresh,inset_choice,outfilestring, totreads, stopcodons)
		print commandstring
		os.system(commandstring)
	def aL_countTable():
		trspdictfilestring = '%s/%s/density%sp_%sshift_aL/%s_aLf/%s_aLf_' %(fp_assign_path, densitystring, assignment, ribosome_shift, '${name}','${name}',)
		outfilestring = "%s/%s_aL_%s_countTable_rpkm_utr3adj" % (outfolder, '${name}',norm_type)
		print ""
		# need to input trspdictfilestring, UTRfilesting, cdsDenThresh,raw_dense_thresh, rpm_dense_thresh, norm_type, inset_choice
		commandstring = 'python ${params.riboseq_buildDenseTables_rpkm_utr3adj_path} --trspdictfilestring %s --UTRfilestring %s --cdsDenThresh %s --norm_type %s --raw_dense_thresh %s --rpm_dense_thresh %s --inset_choice %s --outfilestring %s --totreads %s --stopcodons %s' % (trspdictfilestring,UTRfilestring,cdsDenThresh,norm_type,raw_dense_thresh,rpm_dense_thresh,inset_choice,outfilestring, totreads, stopcodons)
		print commandstring
		os.system(commandstring)
	p1 = multiprocessing.Process(target=fl_countTable)
	p1.start()
	p2 = multiprocessing.Process(target=eA_countTable)
	p2.start()
	p3 = multiprocessing.Process(target=eE_countTable)
	p3.start()
	p4 = multiprocessing.Process(target=aL_countTable)
	p4.start()
	p1.join()
	p2.join()
	p3.join()
	p4.join()

def main():
    fiveorthreeprime = "5"
    densebuilder_multi_raw(fiveorthreeprime)
    densebuilder_multi_norm(fiveorthreeprime)
    avggene_multi("1", fiveorthreeprime)
    avggene_multi("2", fiveorthreeprime)
    get_riboshift("P", "1") # shiftsite, alignposition
    get_riboshift("A", "1")
    densebuilder_riboshift("A")
    avggene_riboshift("1","A", 'uneq', '0') 
    avggene_riboshift("2","A", 'uneq', '0')
    avggene_riboshift("1","A", 'eq', '10')
    avggene_riboshift("2","A", 'eq', '10')
    build_count_tables()
    build_count_tables_utr3adj()
    os.system("mv ${name} ${dirName}")

if __name__ == '__main__':
	# execute only if run as a script
	main()
"""
}

params.single_transcript_mRNA_csv =  ""  //* @input
params.single_transcript_UTR_csv =  ""  //* @input
params.genome2bit =  ""  //* @input
params.rphelper_path =  ""  //* @input

//* autofill
if ($HOSTNAME == "default"){
    $CPU  = 10
    $MEMORY = 30
}
//* platform
if ($HOSTNAME == "ghpcc06.umassrc.org"){
    $TIME = 240
    $CPU  = 10
    $MEMORY = 30
    $QUEUE = "short"
} 
//* platform
//* autofill

process riboseq_avggene_cdsNorm_start {

input:
 set val(name), file(dirName) from g_5_outputDir_g_9

output:
 set val(name), file("${newdirName}")  into g_9_outputDir_g_12

script:
newdirName = name + "riboseq_avggene_cdsNorm_start_stop"
"""
#!/usr/bin/env python

import matplotlib
matplotlib.use('Agg') # set backend for matplotlib
import matplotlib.pyplot as plt 
plt.rcParams['pdf.fonttype'] = 42 # this keeps most text as actual text in PDFs, not outlines


import sys
sys.path.append("${params.rphelper_path}")
from Bio import SeqIO
import twobitreader
import GFF
import csv
import os
import struct
import pandas as pd
import rphelper as rph
import argparse
import importlib


threshold = 10
customSize = 30
pop = "fl"

inset_choice = 'default'
flmin = 28
flmax = 35 # set to exact readsize that is desired
eAmin = 21
eAmax = 24
eEmin = 18
eEmax = 19

assignment = "5"
ribosome_shift = "A"
norm_type = 'rpm'
densitystring = 'Density_rpm'

norm = 'eq'

if pop == "fl":
	minlen = str(flmin)
	maxlen = str(flmax)
elif pop == "eA":
	minlen = str(eAmin)
	maxlen = str(eAmax)
elif pop == "eE":
	minlen = str(eEmin)
	maxlen = str(eEmax)
elif pop == "custom":
	minlen = str(customSize)
	maxlen = str(customSize)
else:
	print "read lengths not set"



mRNAdf = pd.read_csv('${params.single_transcript_mRNA_csv}', sep=',', index_col=0)

def load_genomes(UTRfilestring, twobitfile):
	UTRdict= rph.readindict(open(UTRfilestring, "rU"))
	genome= twobitreader.TwoBitFile(twobitfile) 
	return UTRdict, genome

def build_avggene_first(UTRdict, region_choice, upstreamNTs, downstreamNTs):
    fp_assign_path = "${name}"
    samp = "${name}"
    trspfilestring = '%s/%s/density%sp_%sshift_%sto%s/%s_%sto%sf/%s_%sto%sf_' %(fp_assign_path, densitystring, assignment, ribosome_shift, minlen, maxlen, samp, minlen, maxlen, samp, minlen, maxlen)
    totreads_countfile = "${name}/${name}_FPassigned_counts.txt"
    totreadcountf = open(totreads_countfile, "r")
    totreads = int(totreadcountf.read())
    totreadcountf.close()
    ### This is where all the count files are loaded into a dictionary
    trspdict= rph.readcountsf(trspfilestring)
    ## create a list of 0's that is the length of the region of interest
    ## add 3 to account for the stop codon
    ## this will be added to for every transcript
    averagegene= [0 for num in range(0, 3+upstreamNTs+ downstreamNTs)] # add 3 for start or stop codon
    ## add counters and set to zero
    noUTRentry = 0 ### discard transcripts not in UTRdict
    zeroCdsdense = 0 ### discard transcripts with zero reads in CDS
    lowCdsdense = 0 ## Optional CDS density thresholding, in RPKM
    regionTooShort = 0 ### not enough 3'UTR in region of interest past first stop codon
    totalCountedTranscripts = 0 ## number included in final output
    
    ### calculate mRNA-region densities
    defaultInsets = { 'utr5Inset3' : 6, 'cdsInset5' : 18, 'cdsInset3' : 15, 'utr3Inset5' : 6 }
    zeroInsets    = { 'utr5Inset3' : 0, 'cdsInset5' : 0, 'cdsInset3' : 0, 'utr3Inset5' : 0 }
    customInsets  = { 'utr5Inset3' : 15, 'cdsInset5' : 24, 'cdsInset3' : 15, 'utr3Inset5' : 15 }
    
    if inset_choice == "default":
        insets = defaultInsets
    elif inset_choice == "zero":
        insets = zeroInsets
    elif inset_choice == "custom":
        insets = customInsets
    else:
        print "Insets were not set"
        sys.exit()
    
    ### Iterated through transcripts one at a time, retrieving counts in region of interest:
    for trsp in trspdict:
        ### Load in count file for the transcript here
        exonsplicedcounts = trspdict[trsp]
        
        if UTRdict.has_key(trsp)!=True: # check to make sure density file has an annotation in the UTR csv
            noUTRentry +=1
            continue
        
        mrnalen = int(UTRdict[trsp][3])
        cdslen = int(UTRdict[trsp][4])
        utr5len = int(UTRdict[trsp][5])
        utr3len = int(UTRdict[trsp][6])
        assert mrnalen == cdslen + utr5len + utr3len ## check that this is true
        
        ### define Coding sequence here
        cdsstart = utr5len
        cdsend = len(exonsplicedcounts) - utr3len # cdsend == first position of utr3
        if cdsstart == cdsend:
            print "Error, gene length is 0 for transcript "+ trsp
            sys.exit()
        
        ### Calculate Region Densities ###
        utr5lenMod = utr5len-insets['utr5Inset3']
        cdslenMod = cdslen-insets['cdsInset5']-insets['cdsInset3']
        utr3lenMod = utr3len-insets['utr3Inset5']
        mrnalenMod = utr5lenMod+cdslenMod+utr3lenMod
        
        utr5Counts = sum(exonsplicedcounts[:cdsstart-insets['utr5Inset3']])
        cdsCounts = sum(exonsplicedcounts[cdsstart+insets['cdsInset5']:cdsend-insets['cdsInset3']])
        utr3Counts = sum(exonsplicedcounts[cdsend+insets['utr3Inset5']:])
        mrnaCounts = utr5Counts+cdsCounts+utr3Counts
        
        ### RAW counts
        RAWutr5Counts = int(utr5Counts*(totreads/1E6))
        RAWutr3Counts = int(utr3Counts*(totreads/1E6))
        RAWcdsCounts = int(cdsCounts*(totreads/1E6))
        RAWmrnaCounts = int(mrnaCounts*(totreads/1E6))
        ### denisites
        # mrnaDensity = (mrnaCounts/mrnalenMod)
        cdsDensity = (cdsCounts/cdslenMod)
        # utr5Density = (utr5Counts/utr5lenMod)
        # utr3Density = (utr3Counts/utr3lenMod)
        #### RPKM densities
        # mrnaDensity_rpkm = (mrnaCounts/mrnalenMod) * 1000
        cdsDensity_rpkm = (cdsCounts/cdslenMod) * 1000
        # utr5Density_rpkm = (utr5Counts/utr5lenMod) * 1000
        # utr3Density_rpkm = (utr3Counts/utr3lenMod) * 1000
        ### throw out zero's
        if cdsDensity == 0:
            zeroCdsdense += 1
            continue
        if cdsDensity*float(1000)< int(threshold):	# Threshold on cds density: (thresholding on "rpkm")
            lowCdsdense += 1
            continue
        ### define vector in valid CDS region, normalize by cdsDensity
        cdsSplicedCounts = exonsplicedcounts[cdsstart+insets['cdsInset5']:cdsend-insets['cdsInset3']]
        cdsNormCounts = [rpf/cdsDensity for rpf in cdsSplicedCounts] ## just region of cds within insets for counts, sum/len == 1
        # print cdslenMod, "modified cds length"
        # print sum(cdsNormCounts), "total normalized counts, should equal length of cds"
        # print sum(cdsNormCounts)/len(cdsNormCounts) ## this should == 1
        ### define vector for whole transcript, normlaized by cdsDensity
        exonNormCounts = [rpf/cdsDensity for rpf in exonsplicedcounts] ## counts normalized by cds density, using only region within insets
        # print exonNormCounts
        # print sum(exonNormCounts)
        # print sum(exonNormCounts)/len(exonNormCounts) ## should typically be less than 1, unless greater density from utrs and start/stop codons
        ### for start codon metagenes
        if region_choice == 'start':
            ### check boundaries:
            # print mRNAdf.loc[trsp]['mRNAseqs'][cdsstart:cdsstart+3] ### start codon
            if len(exonNormCounts[cdsstart-upstreamNTs:cdsstart+3+downstreamNTs]) < upstreamNTs+3+downstreamNTs:
                regionTooShort +=1
                continue
            else:
                totalCountedTranscripts +=1
                avggene_counts = exonNormCounts[cdsstart-upstreamNTs:cdsstart+3+downstreamNTs]
                for i in range(len(avggene_counts)): ### add these counts to the running total
                    averagegene[i] += avggene_counts[i]
        
        if region_choice == 'stop':
            ### check boundaries
            # print mRNAdf.loc[trsp]['mRNAseqs'][cdsend-3:cdsend] ### stop codon
            if len(exonNormCounts[cdsend-3-upstreamNTs:cdsend+downstreamNTs]) < upstreamNTs+3+downstreamNTs:
                regionTooShort +=1
                continue
            else:
                totalCountedTranscripts +=1
                # print upstreamNTs+3+downstreamNTs
                # print len(exonNormCounts[cdsend-3-upstreamNTs:cdsend+downstreamNTs])
                avggene_counts = exonNormCounts[cdsend-3-upstreamNTs:cdsend+downstreamNTs]
                for i in range(len(avggene_counts)): ### add these counts to the running total
                    averagegene[i] += avggene_counts[i]
    
    averagegene_equal = [rpf/totalCountedTranscripts for rpf in averagegene] ### divide by total number of valid transcripts
    # print averagegene
    # print averagegene_equal
    positions = range(-upstreamNTs-1, downstreamNTs+2) # start or stop codon is [-1,0,1]
    df = pd.DataFrame({'position':positions, 'avg':averagegene_equal})
    df = df[['position', 'avg']]
    # print df['avg'].sum()/len(df)
    # print df 
    print "Avggene run compolete for sample %s" % samp
    print "Transcripts inculded %s" % totalCountedTranscripts
    print "Number of transcripts absent in UTRfile: %s" % noUTRentry
    print "Number of transcripts with zero CDS density: %s" % zeroCdsdense
    print "Number of transcripts below CDS density threshold: %s" % lowCdsdense 
    print "Number of transcripts too short avggene region: %s" % regionTooShort
    print "- - - - - - - -"
    print df
    print "- - - - - - - -"
    
    if region_choice == 'start':
        alignpos = "1"
    elif region_choice == 'stop':
        alignpos = "2"
    else:
        print "alignpos not set!!!"
    
    ### write csv file
    fp_assign_path = "${name}"
    avggene_csv_path = "%s/avggene%s_ORF%s_%sshift_%s%s150" % (fp_assign_path, alignpos, norm_type, ribosome_shift, assignment, norm)
    if not os.path.exists(avggene_csv_path):	os.makedirs(avggene_csv_path)
    # adding cdsNorm to indicate normalized densities
    csv_outfile = "%s/%s_%s_rpkmThresh%s_%sto%sf_avg_%s_cdsNorm.csv" %(avggene_csv_path, samp, pop, threshold, minlen, maxlen, alignpos)
    df.to_csv(csv_outfile, index=False)

def main():
    os.system("mv ${dirName} ${name}")
    UTRfilestring = "${params.single_transcript_UTR_csv}"
    twobitfile = "${params.genome2bit}"
    UTRdict, genomes = load_genomes(UTRfilestring, twobitfile)
    build_avggene_first(UTRdict, "start", 100, 150)
    build_avggene_first(UTRdict, "stop", 150, 100)
    os.system("mv ${name} ${newdirName}")

if __name__ == '__main__':
    main()
"""
}

params.motif_file_root_path =  ""  //* @input
params.codon_occupancy_main_path =  ""  //* @input
params.single_transcript_UTR_csv =  ""  //* @input

//* autofill
if ($HOSTNAME == "default"){
    $CPU  = 10
    $MEMORY = 30
}
//* platform
if ($HOSTNAME == "ghpcc06.umassrc.org"){
    $TIME = 240
    $CPU  = 10
    $MEMORY = 30
    $QUEUE = "short"
} 
//* platform
//* autofill

process codon_occupancies {

input:
 set val(name), file(dirName) from g_9_outputDir_g_12

output:
 set val(name), file("${newdirName}")  into g_12_outputDir_g_6

script:
newdirName = name + "codon_occupancies"
"""
#!/usr/bin/env python

import sys
import os
import numpy as np
import argparse
import importlib
from pathos.multiprocessing import ProcessingPool as Pool

motiffilerootpath = "${params.motif_file_root_path}"

motifs= " \\"['TTT','TTC','TTA','TTG','CTT','CTC','CTA','CTG','ATT','ATC','ATA','ATG','GTT','GTC','GTA','GTG','TCT','TCC','TCA','TCG','CCT','CCC','CCA','CCG','ACT','ACC','ACA','ACG','GCT','GCC','GCA','GCG','TAT','TAC','CAT','CAC','CAA','CAG','AAT','AAC','AAA','AAG','GAT','GAC','GAA','GAG','TGT','TGC','TGG','AGT','AGC','GGT','GGC','GGA','GGG','AGA','AGG','CGT','CGC','CGA','CGG','TAA','TAG','TGA']\\" "
cds5trim= '15'
cds3trim= '15'
seqwin= ' "[\\'0\\',\\'3\\']" ' # For density files shifted to the first nt of A sites, should be multiple of 3
exclusionfiles=  ' "[\\'0\\']" '
norm_type = "rpm"
flmin = 28
flmax = 35 # set to exact readsize that is desired
eAmin = 21
eAmax = 24
eEmin = 18
eEmax = 19
UTRfilestring = "${params.single_transcript_UTR_csv}"

## define handler function:

def get_codon_occ_fl(file):
	#function to get codon occupancies
	#footprint_pop: fl, eA, eE, or aL
	#ribosome_site: A, P, E, or 0 
	footprint_pop = 'fl'
	ribosome_site = 'A'
	fiveorthreeprime = '5'
	### set normalization type:
	if norm_type == "raw":
		densitystring = "DensityUnnormalized"
	elif norm_type == "rpm":
		densitystring = "Density_rpm"
	else:
		print "Normalization is not set!"
		sys.exit()

	### set footprint population:
	if footprint_pop == 'fl':
		readsizemin = str(flmin)
		readsizemax = str(flmax)
	elif footprint_pop == 'eA':
		readsizemin = str(eAmin)
		readsizemax = str(eAmax)
	elif footprint_pop == 'eE':
		readsizemin = str(eEmin)
		readsizemax = str(eEmax)
	elif footprint_pop == 'aL':
		readsizemin = 'aL'

	### density file paths:
	
	trspdictfilestring = '%s/%s/density%sp_%sshift_%sto%s/%s_%sto%sf/%s_%sto%sf_' % (
		file, densitystring, fiveorthreeprime, ribosome_site, readsizemin, 
		readsizemax, file, readsizemin, readsizemax,file,readsizemin,readsizemax)
	
	sample_name = file

	outlistfilepath = '%s/codon' % (file)
	if not os.path.exists(outlistfilepath): os.makedirs(outlistfilepath)
	outlistfile = '%s/%s_%sshift_%sto%s__%soccupancy_%scds5trim_%scds3trim_codonOccupancy.csv' % (
		outlistfilepath, file, ribosome_site, readsizemin, readsizemax, fiveorthreeprime, cds5trim, cds3trim)
	outfileparams = '%s/%s_%sshift_%sto%s__%soccupancy_%scds5trim_%scds3trim_codonOccupancy_params.txt' % (
		outlistfilepath, file, ribosome_site, readsizemin, readsizemax, fiveorthreeprime, cds5trim, cds3trim)

	commandstring= 'python ${params.codon_occupancy_main_path} \
						--motiffilerootpath %s \
						--motifs %s \
						--trspdictfilestring %s \
						--sample_name %s \
						--UTRfilestring %s \
						--cds5trim %s \
						--cds3trim %s \
						--seqwin %s \
						--exclusionfiles %s \
						--outfileparams %s \
						--outlistfile %s' % (

						motiffilerootpath, 
						motifs, 
						trspdictfilestring, 
						sample_name, 
						UTRfilestring, 
						cds5trim, 
						cds3trim, 
						seqwin, 
						exclusionfiles, 
						outfileparams, 
						outlistfile)
	print commandstring
	os.system(commandstring)


def get_codon_occ_eA(file):
	# function to get codon occupancies
	# footprint_pop: fl, eA, eE, or aL
	# ribosome_site: A, P, E, or 0 although zero is not very informative here
	print "starting sample %s" % (file)

	footprint_pop = 'eA'
	ribosome_site = 'A'
	fiveorthreeprime = '5'

	### set normalization type:
	if norm_type == "raw":
		densitystring = "DensityUnnormalized"
	elif norm_type == "rpm":
		densitystring = "Density_rpm"
	else:
		print "Normalization is not set!"
		sys.exit()

	### set footprint population:
	if footprint_pop == 'fl':
		readsizemin = str(flmin)
		readsizemax = str(flmax)
	elif footprint_pop == 'eA':
		readsizemin = str(eAmin)
		readsizemax = str(eAmax)
	elif footprint_pop == 'eE':
		readsizemin = str(eEmin)
		readsizemax = str(eEmax)
	elif footprint_pop == 'aL':
		readsizemin = 'aL'

	# what is codonocc reading these in as?
	### density file paths:
	
	trspdictfilestring = '%s/%s/density%sp_%sshift_%sto%s/%s_%sto%sf/%s_%sto%sf_' % (
		file, densitystring, fiveorthreeprime, ribosome_site, readsizemin, 
		readsizemax, file, readsizemin, readsizemax,file,readsizemin,readsizemax)
	
	sample_name = file

	outlistfilepath = '%s/codon' % (file)
	if not os.path.exists(outlistfilepath): os.makedirs(outlistfilepath)
	outlistfile = '%s/%s_%sshift_%sto%s__%soccupancy_%scds5trim_%scds3trim_codonOccupancy.csv' % (
		outlistfilepath, file, ribosome_site, readsizemin, readsizemax, fiveorthreeprime, cds5trim, cds3trim)
	outfileparams = '%s/%s_%sshift_%sto%s__%soccupancy_%scds5trim_%scds3trim_codonOccupancy_params.txt' % (
		outlistfilepath, file, ribosome_site, readsizemin, readsizemax, fiveorthreeprime, cds5trim, cds3trim)

	commandstring= 'python ${params.codon_occupancy_main_path} \
						--motiffilerootpath %s \
						--motifs %s \
						--trspdictfilestring %s \
						--sample_name %s \
						--UTRfilestring %s \
						--cds5trim %s \
						--cds3trim %s \
						--seqwin %s \
						--exclusionfiles %s \
						--outfileparams %s \
						--outlistfile %s' % (

						motiffilerootpath, 
						motifs, 
						trspdictfilestring, 
						sample_name, 
						UTRfilestring, 
						cds5trim, 
						cds3trim, 
						seqwin, 
						exclusionfiles, 
						outfileparams, 
						outlistfile)

	print commandstring
	os.system(commandstring)


def main():
    os.system("mv ${dirName} ${name}")
    samplelist = ["${name}"]
    p = Pool(nodes = ${thread})
    p.map(get_codon_occ_fl, samplelist)
    p2 = Pool(nodes = ${thread})
    p2.map(get_codon_occ_eA, samplelist)
    os.system("mv ${name} ${newdirName}")


if __name__ == '__main__':
	main()
"""
}


process aggr {

publishDir params.outdir, overwrite: true, mode: 'copy',
	saveAs: {filename ->
	if (filename =~ /${name}$/) "results/$filename"
	else if (filename =~ /${name}\/countTables\/.*.csv$/) "count_tables/$filename"
}

input:
 set val(name), file(dirName) from g_12_outputDir_g_6

output:
 file "${name}"  into g_6_outputDir_g_7, g_6_outputDir_g_10, g_6_outputDir_g_13, g_6_outputDir_g_15, g_6_outputDir_g_17
 set val(name), file("${name}/countTables/*.csv")  into g_6_outputFileTab
 file "${name}_result"  into g_6_resultsdir

"""
mv ${dirName} ${name}
mkdir ${name}_result
rsync -avzu --exclude='Density_rpm' --exclude='DensityUnnormalized' ${name}/ ${name}_result/

"""
}


process figure4B {

publishDir params.outdir, overwrite: true, mode: 'copy',
	saveAs: {filename ->
	if (filename =~ /.*.pdf$/) "reports/$filename"
}

input:
 file nameAll from g_6_outputDir_g_10.collect()

output:
 file "*.pdf"  into g_10_outputFilePdf

script:
nameAll = nameAll.collect{ '"' + it + '"'}
"""
#!/usr/bin/env python
# script is based off of riboseq_readingFrame_utr3meta
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt 
plt.rcParams['pdf.fonttype'] = 42 # this keeps most text as actual text in PDFs, not outlines

## import dependencies
import sys
import math
import matplotlib.patches as mpatches
import numpy as np
import pandas as pd
import scipy.stats as stats
pd.set_option('display.max_columns', 50)
import seaborn as sns
# from pylab import *
import argparse
import importlib


samplelist = ${nameAll}
### set colors
black = '#000000'
orange = '#ffb000'
cyan = '#63cfff'
red = '#eb4300'
green = '#00c48f'
pink = '#eb68c0'
yellow = '#fff71c'
blue = '#006eb9'

# colorList = [black, orange, red, green, blue, yellow, pink]
colorList = ['lightgray', 'gray', 'black', 'gold', 'darkorange','orangred']

def log_trans_b10(x):
	try:
		return math.log(x, 10)
	except:
		return float(-6.00)
#         return float("NaN")
def log_trans_b2(x):
	try:
		return math.log(x, 2)
	except:
		# return float("NaN")
		return float(-15.00) # set arbitrarily low value

def load_countTables():
	namelist = []
	dflist= []
	# dflu3 = []
	for samp in samplelist:
		dftemp = pd.read_csv('%s/countTables/%s_fl_rpm_28to35_countTable_rpkm_utr3adj.csv' % (samp, samp))
		dftemp['sampname'] = samp
		dftemp['cdsCountsLog2'] = dftemp['cdsCounts'].apply(log_trans_b2)
		dftemp['RAWcdsCountsLog2'] = dftemp['RAWcdsCounts'].apply(log_trans_b2)
		dftemp['utr3CountsLog2'] = dftemp['utr3Counts'].apply(log_trans_b2)
		dftemp['utr3AdjCountsLog2'] = dftemp['utr3AdjCounts'].apply(log_trans_b2)
		dftemp['utr3OccLog2'] = dftemp['utr3_occupancy'].apply(log_trans_b2)
		dftemp['RRTSlog2'] = dftemp['RRTS'].apply(log_trans_b2)

		dftemp = dftemp.loc[dftemp['RAWcdsCounts'] > 128]
		dflist.append(dftemp)
		namelist.append(samp)
		
	# ### combine into one master dataframe
	df = pd.concat(dflist, axis=0, ignore_index = True)
	return df, dflist, namelist


def plot_RRTS_boxplot(df, dflist, namelist):

	figout = "Fig4B.pdf"
	validsamps = namelist
	dfp = df.loc[df['sampname'].isin(validsamps)]
	dfp = dfp.sort_values('stopcodon', ascending=True)

	fig, ax = plt.subplots(figsize=(6, 6))
	ax = sns.boxplot(x=dfp['sampname'],  y=dfp['RRTS'], orient = 'v',
				 fliersize = 0, hue = dfp['stopcodon'], notch=True, order = validsamps, palette=colorList)

	ax = plt.gca()
	ax.set_ylim((0,1.1))
	plt.savefig(figout, format='pdf', bbox_inches = "tight")

def main():
	df, dflist, namelist = load_countTables()
	print df.head()
	plot_RRTS_boxplot(df, dflist, namelist)

if __name__ == '__main__':
	main()

"""
}


process figures {

publishDir params.outdir, overwrite: true, mode: 'copy',
	saveAs: {filename ->
	if (filename =~ /.*.pdf$/) "reports/$filename"
}

input:
 file nameAll from g_6_outputDir_g_7.collect()

output:
 file "*.pdf"  into g_7_outputFilePdf

script:
nameAll = nameAll.collect{ '"' + it + '"'}
"""
#!/usr/bin/env python
# This script is for taking input density files and plotting average genes around the stop codon
# can be used to overlay multiple samples on the same plot

## plotting
import matplotlib
matplotlib.use('Agg') ## set backend here
import matplotlib.pyplot as plt 
plt.rcParams['pdf.fonttype'] = 42 # this keeps most text as actual text in PDFs, not outlines

## import dependencies
import sys
from Bio import SeqIO
import twobitreader
import GFF
import csv
import os
import pandas as pd
import argparse
import importlib
import random

### inputs:
alignposition = "2" 	# set 1 for start codon and 2 for stop codon
pop = "fl" 	# fl, eA, or eE for density files
customSize = 30
threshold = 10
ribosome_site = "A" 	# A, P, E, or 0
flmin = 28
flmax = 35
eAmin = 21
eAmax = 24
eEmin = 18
eEmax = 19
norm_type = "rpm"

samples = ${nameAll}
print samples

sample_plot_names = samples
samples_plotted = '_vs_'.join(samples)

number_of_colors = len(samples)
colorList = ["#"+''.join([random.choice('0123456789ABCDEF') for j in range(6)])
             for i in range(number_of_colors)]
print colorList
col_list = colorList[::-1]

if pop == "fl":
	minlen = str(flmin)
	maxlen = str(flmax)
elif pop == "eA":
	minlen = str(eAmin)
	maxlen = str(eAmax)
elif pop == "eE":
	minlen = str(eEmin)
	maxlen = str(eEmax)
elif pop == "custom":
	minlen = str(customSize)
	maxlen = str(customSize)
else:
	print "whoops, something went wrong here, horribly horribly wrong!"
	# pop = ftsize


### for already shifted samples
def avggene_riboshift_plot_overlay(alignposition, ribosome_site, fiveorthreeprime='5', normalization='uneq', threshold ="0"):
	# read_length_list = ftsize
	alignpos = alignposition # '1' for start, '2' for stop
	assignment = fiveorthreeprime # should be 5' mapped at this point
	norm = normalization # 'uneq' for no normalizaiton; 'eq' to give equal weight to all genes
	ribosome_shift = ribosome_site # 'A' or 'P' or maybe 'E' later

	sample_numb = len(samples)
	sample_index = range(0,len(samples))
	sample_dict = {}
	for i in sample_index:
		sample_dict[samples[i]] = sample_index[i]

	df_fl_list = []
	df_eA_list = []
	df_eE_list = []
	df_aL_list = []
	df_custom_list = []


	for file in samples:
		fp_assign_path = '%s' % (file)
		avggene_csv_path = "%s/avggene%s_ORF%s_%sshift_%s%s150" % (fp_assign_path, alignpos, norm_type, ribosome_shift, assignment, norm) # norm should be 'uneq' for now


		fl_avggene_csv = '%s/%s_fl_rpkmThresh%s_%sto%sf_avg_%s_cdsNorm.csv' % (avggene_csv_path, file, threshold, flmin, flmax, alignpos)
		# eA_avggene_csv = '%s/%s_eA_rpkmThresh%s_%sto%sf_avg_%s_cdsNorm.csv' % (avggene_csv_path, file, threshold, eAmin, eAmax, alignpos)
		# eE_avggene_csv = '%s/%s_eE_rpkmThresh%s_%sto%sf_avg_%s_cdsNorm.csv' % (avggene_csv_path, file, threshold, eEmin, eEmax, alignpos)
		# aL_avggene_csv = '%s/%s_rpkmThresh%s_aL_aLf_avg_%s_cdsNorm.csv' % (avggene_csv_path, file, threshold, alignpos)

		if pop == "custom":
			custom_avggene_csv = '%s/%s_%s_shiftCustom_rpkmThresh0_%sto%sf_avg_%s.csv' % (avggene_csv_path, file, customSize, minlen, maxlen, alignpos)
			print custom_avggene_csv
			custom_avggene_df = pd.read_csv(custom_avggene_csv, index_col = 0, header=0)
			df_custom_list.append(custom_avggene_df)

		fl_avggene_df = pd.read_csv(fl_avggene_csv, index_col = 0, header=0)
		# eA_avggene_df = pd.read_csv(eA_avggene_csv, index_col = 0, header=0)
		# eE_avggene_df = pd.read_csv(eE_avggene_csv, index_col = 0, header=0)
		# aL_avggene_df = pd.read_csv(aL_avggene_csv, index_col = 0, header = 0)

		df_fl_list.append(fl_avggene_df)
		# df_eA_list.append(eA_avggene_df)
		# df_eE_list.append(eE_avggene_df)
		# df_aL_list.append(aL_avggene_df)

	plot_outfile = "Fig2A.pdf"


	fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(6,4))


	counter = 0 ## for making offsets
	if pop == 'fl':
		for i in range(len(df_fl_list)):
			df_fl_list[i].plot.line(x=df_fl_list[i].index+counter, y=df_fl_list[i].columns.values[0], ax=ax, color = col_list[i], lw=1, use_index=True, label=sample_plot_names[i])
			# counter +=6
			# df_fl_list[i].plot.scatter(y=df_fl_list[i].columns.values[0], ax=ax, color = col_list[i], s=2, use_index=True)
	if pop == 'eA':
		for i in range(len(df_eA_list)):
			df_eA_list[i].plot.line(y=df_eA_list[i].columns.values[0], ax=ax, color = col_list[i], lw=1, use_index=True, label=sample_plot_names[i])
	if pop == 'custom':
		for i in range(len(df_fl_list)):
			df_custom_list[i].plot.line(x=df_custom_list[i].index+counter, y=df_custom_list[i].columns.values[0], ax=ax, color = col_list[i], lw=1, use_index=True, label=sample_plot_names[i])

	plt.legend(loc=1, prop={'size': 6})
	plt.savefig(plot_outfile, format = 'pdf', bbox_inches = "tight")
	plt.close()

### for already shifted samples, for now
def avggene_riboshift_plot_overlay_zoom(alignposition, ribosome_site, fiveorthreeprime='5', normalization='uneq', threshold ="0"):
	# read_length_list = ftsize
	alignpos = alignposition # '1' for start, '2' for stop
	assignment = fiveorthreeprime # should be 5' mapped at this point
	norm = normalization # 'uneq' for no normalizaiton; 'eq' to give equal weight to all genes
	ribosome_shift = ribosome_site # 'A' or 'P' or maybe 'E' later

	sample_numb = len(samples)
	sample_index = range(0,len(samples))
	sample_dict = {}
	for i in sample_index:
		sample_dict[samples[i]] = sample_index[i]

	df_fl_list = []
	df_eA_list = []
	df_eE_list = []
	df_aL_list = []
	df_custom_list = []


	for file in samples:
		fp_assign_path = '%s' % (file)
		avggene_csv_path = "%s/avggene%s_ORF%s_%sshift_%s%s150" % (fp_assign_path, alignpos, norm_type, ribosome_shift, assignment, norm) # norm should be 'uneq' for now


		fl_avggene_csv = '%s/%s_fl_rpkmThresh%s_%sto%sf_avg_%s_cdsNorm.csv' % (avggene_csv_path, file, threshold, flmin, flmax, alignpos)
		# eA_avggene_csv = '%s/%s_eA_rpkmThresh%s_%sto%sf_avg_%s_cdsNorm.csv' % (avggene_csv_path, file, threshold, eAmin, eAmax, alignpos)
		# eE_avggene_csv = '%s/%s_eE_rpkmThresh%s_%sto%sf_avg_%s_cdsNorm.csv' % (avggene_csv_path, file, threshold, eEmin, eEmax, alignpos)
		# aL_avggene_csv = '%s/%s_rpkmThresh%s_aL_aLf_avg_%s_cdsNorm.csv' % (avggene_csv_path, file, threshold, alignpos)

		if pop == "custom":
			custom_avggene_csv = '%s/%s_%s_shiftCustom_rpkmThresh0_%sto%sf_avg_%s.csv' % (avggene_csv_path, file, customSize, minlen, maxlen, alignpos)
			print custom_avggene_csv
			custom_avggene_df = pd.read_csv(custom_avggene_csv, index_col = 0, header=0)
			df_custom_list.append(custom_avggene_df)

		fl_avggene_df = pd.read_csv(fl_avggene_csv, index_col = 0, header=0)
		# eA_avggene_df = pd.read_csv(eA_avggene_csv, index_col = 0, header=0)
		# eE_avggene_df = pd.read_csv(eE_avggene_csv, index_col = 0, header=0)
		# aL_avggene_df = pd.read_csv(aL_avggene_csv, index_col = 0, header = 0)

		df_fl_list.append(fl_avggene_df)
		# df_eA_list.append(eA_avggene_df)
		# df_eE_list.append(eE_avggene_df)
		# df_aL_list.append(aL_avggene_df)

	plot_outfile = "Fig2B.pdf"


	fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(4,4))


	ax.set_ylim(0, 0.35)
	ax.set_xlim(0,100)

	counter = 0 ## for making offsets
	if pop == 'fl':
		for i in range(len(df_fl_list)):
			df_fl_list[i].plot.line(x=df_fl_list[i].index+counter, y=df_fl_list[i].columns.values[0], ax=ax, color = col_list[i], lw=1, use_index=True, label=sample_plot_names[i])
			# counter +=6
			# df_fl_list[i].plot.scatter(y=df_fl_list[i].columns.values[0], ax=ax, color = col_list[i], s=2, use_index=True)
	if pop == 'eA':
		for i in range(len(df_eA_list)):
			df_eA_list[i].plot.line(y=df_eA_list[i].columns.values[0], ax=ax, color = col_list[i], lw=1, use_index=True, label=sample_plot_names[i])
	if pop == 'custom':
		for i in range(len(df_fl_list)):
			df_custom_list[i].plot.line(x=df_custom_list[i].index+counter, y=df_custom_list[i].columns.values[0], ax=ax, color = col_list[i], lw=1, use_index=True, label=sample_plot_names[i])

	plt.legend(loc=1, prop={'size': 6})
	plt.savefig(plot_outfile, format = 'pdf', bbox_inches = "tight")
	plt.close()

def main():
	avggene_riboshift_plot_overlay(alignposition, ribosome_site, normalization='eq', threshold=threshold)
	avggene_riboshift_plot_overlay_zoom(alignposition, ribosome_site, normalization='eq', threshold=threshold)

if __name__ == '__main__':
	main()
"""
}

//* autofill
if ($HOSTNAME == "default"){
    $CPU  = 10
    $MEMORY = 30
}
//* platform
if ($HOSTNAME == "ghpcc06.umassrc.org"){
    $TIME = 240
    $CPU  = 10
    $MEMORY = 30
    $QUEUE = "short"
} 
//* platform
//* autofill

process figure3S1C {

publishDir params.outdir, overwrite: true, mode: 'copy',
	saveAs: {filename ->
	if (filename =~ /.*.pdf$/) "reports/$filename"
}

input:
 file nameAll from g_6_outputDir_g_17.collect()

output:
 file "*.pdf"  into g_17_outputFilePdf

script:
nameAll = nameAll.collect{ '"' + it + '"'}
"""
#!/usr/bin/env python
# Script for plotting read size distributions for each region of an mRNA
# Can then look at 5'UTR vs CDS vs 3'UTR vs mRNA
# Relies on inputs from densebuilder: Unnormalized, 5'Mapped

import matplotlib
matplotlib.use('Agg') # set backend for matplotlib
import matplotlib.pyplot as plt 
plt.rcParams['pdf.fonttype'] = 42 # this keeps most text as actual text in PDFs, not outlines

import sys, os
sys.path.append("${params.rphelper_path}")
import pandas as pd
import pysam
import numpy as np
import csv
from datetime import datetime
import argparse
from pathos.multiprocessing import ProcessingPool as Pool
pd.set_option('display.max_columns', 40)
import argparse
import importlib
import rphelper as rph

samplelist = ${nameAll}
ftsize = [15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40]


### set inputs here:
defaultInsets = { 'utr5Inset3' : 6, 'cdsInset5' : 18, 'cdsInset3' : 15, 'utr3Inset5' : 6 }
zeroInsets    = { 'utr5Inset3' : 0, 'cdsInset5' : 0, 'cdsInset3' : 0, 'utr3Inset5' : 0 }
# insets= defaultInsets
insets= zeroInsets
UTRfilestring = "${params.single_transcript_UTR_csv}"
UTRdict = rph.readindict(open(UTRfilestring, "rU"))

def region_size_dist_ftsize(readsize, sample):
	### build the counts for a given region 
	fp_assign_path = '%s' % (sample)
	totreads_countfile = "%s/%s_FPassigned_counts.txt" % (fp_assign_path, sample)
	totreadcountf = open(totreads_countfile, "r")
	totreads = int(totreadcountf.read())
	totreadcountf.close()
	print "total reads for sample %s = %s" % (sample, totreads)

	# for readsize in ftsize:
	readsize = str(readsize) # convert to string
	trspdictfilestring = '%s/DensityUnnormalized/density5p_0shift_%s/%s_%sf/%s_%sf_' %(
		fp_assign_path, readsize, sample, readsize, sample, readsize)

	bamfilepath_readsize = '%s/%s_star_default/%s_%s_match.sorted.bam' % (
		fp_assign_path, sample, sample, readsize)
	bamfile = pysam.AlignmentFile(bamfilepath_readsize, 'rb')
	read_count_bam = bamfile.count()
	print "total reads in bamfile for sample: %s, read length: %s, equals == %s" % (sample, readsize, read_count_bam)

	## build the trspdict now for a given readlength:
	trspdict = rph.readcountsf(trspdictfilestring)

	## add counters 
	noUTRentry = 0
	zeroUtrlen = 0
	zeroUtrlenInsets = 0
	zeroCdsdense = 0
	lowCdsdense = 0
	lowCdsCounts = 0 # adding cds raw read counter

	### counters for output
	totUtr5Counts = 0
	totCdsCounts = 0
	totUtr3Counts = 0
	totMrnaCounts = 0

	## iterate through every transcript in the gtf file
	for trsp in trspdict:
		if UTRdict.has_key(trsp)!=True: # check to make sure density file has an annotation in the UTR csv
			noUTRentry +=1
			continue

	# define base region sizes from UTRdict
		mrnalen = int(UTRdict[trsp][3])
		cdslen = int(UTRdict[trsp][4])
		utr5len = int(UTRdict[trsp][5])
		utr3len = int(UTRdict[trsp][6])

		### Not sure if I want to keep this here... see how many have lengths of zero first
		if utr5len == 0:
			zeroUtrlen +=1
			# print("transcript has zero utr5 len %s") % trsp
			# sys.exit() 
			continue
		if utr3len == 0:
			zeroUtrlen +=1
			continue

# get counts from density file
		exonsplicedcounts = trspdict[trsp]

		# set starts and ends 
		cdsstart = utr5len
		cdsend = len(exonsplicedcounts) - utr3len
		if cdsstart == cdsend:
			print "Error, gene length is 0 for transcript %s" % trsp
			sys.exit()


		# modify region lengths using insets:
		
		utr5len = utr5len-insets['utr5Inset3']
		cdslen = cdslen-insets['cdsInset5']-insets['cdsInset3']
		utr3len = utr3len-insets['utr3Inset5']
		mrnalen = utr5len+cdslen+utr3len

		if utr5len == 0:
			zeroUtrlenInsets +=1
			# print "transcript has zero utr5 len %s" % trsp
			# sys.exit() 
			continue
		if utr3len == 0:
			zeroUtrlenInsets +=1
			continue

		utr5Counts = sum(exonsplicedcounts[:cdsstart-insets['utr5Inset3']])
		utr3Counts = sum(exonsplicedcounts[cdsend+insets['utr3Inset5']:])
		cdsCounts = sum(exonsplicedcounts[cdsstart+insets['cdsInset5']:cdsend-insets['cdsInset3']])
		mrnaCounts = utr5Counts+cdsCounts+utr3Counts

		totUtr5Counts += utr5Counts
		totCdsCounts += cdsCounts
		totUtr3Counts += utr3Counts
		totMrnaCounts += mrnaCounts
	print "UTR5total = %s, CDStotal = %s, UTR3total = %s, mRNAtotal = %s " % (totUtr5Counts, totCdsCounts, totUtr3Counts, totMrnaCounts)
	return [int(readsize), totUtr5Counts, totCdsCounts, totUtr3Counts, totMrnaCounts]



def main():
	counter = 0
	treatlist = samplelist
	for sample in samplelist:
		treat = treatlist[counter]
		sl = [sample] * len(ftsize)
		p = Pool(nodes=40) 
		results = p.map(region_size_dist_ftsize, ftsize, sl)

		dfcols = ['readLength', 'utr5', 'cds', 'utr3', 'mrna']
		df = pd.DataFrame.from_records(results, columns = dfcols)
		df = df.set_index('readLength')

		utr5reads = df['utr5'].sum()
		cdsreads = df['cds'].sum()
		utr3reads = df['utr3'].sum()
		mrnareads = df['mrna'].sum()

		df['utr5'] = df['utr5']/utr5reads
		df['cds'] = df['cds']/cdsreads
		df['utr3'] = df['utr3']/utr3reads
		df['mrna'] = df['mrna']/mrnareads

		df.drop('utr5', axis=1, inplace=True)
		df.drop('mrna', axis=1, inplace=True)
		outfig = 'Fig3S1C_%s.pdf' % (treat)

		fig, ax = plt.subplots()

		df.plot.line()
		plt.xticks([15,20,25,30,35,40])
		plt.title(sample)

		plt.savefig(outfig, format='pdf', bbox_inches="tight")
		plt.close()
		counter +=1

if __name__ == '__main__':
	main()
"""
}


process figure2S3A {

publishDir params.outdir, overwrite: true, mode: 'copy',
	saveAs: {filename ->
	if (filename =~ /.*.pdf$/) "reports/$filename"
}

input:
 file nameAll from g_6_outputDir_g_15.collect()

output:
 file "*.pdf"  into g_15_outputFilePdf

script:
nameAll = nameAll.collect{ '"' + it + '"'}
"""
#!/usr/bin/env python
# This script is for taking input density files and plotting average genes around the start codon
# can be used to overlay multiple samples on the same plot

## plotting
import matplotlib
matplotlib.use('Agg') ## set backend here
import matplotlib.pyplot as plt
plt.rcParams['pdf.fonttype'] = 42 # this keeps most text as actual text in PDFs, not outlines

## import dependencies
import sys
from Bio import SeqIO
import twobitreader
import GFF
import csv
import os
import pandas as pd
import argparse
import importlib
import random

### inputs:
alignposition = "1"     # set 1 for start codon and 2 for stop codon
pop = "fl"      # fl, eA, or eE for density files
customSize = 30
threshold = 10
ribosome_site = "A"     # A, P, E, or 0
flmin = 28
flmax = 35
eAmin = 21
eAmax = 24
eEmin = 18
eEmax = 19
norm_type = "rpm"

samples = ${nameAll}

sample_plot_names = samples
samples_plotted = '_vs_'.join(samples)
number_of_colors = len(samples)
colorList = ["#"+''.join([random.choice('0123456789ABCDEF') for j in range(6)])
    for i in range(number_of_colors)]
print colorList
col_list = colorList[::-1]

if pop == "fl":
    minlen = str(flmin)
    maxlen = str(flmax)
elif pop == "eA":
    minlen = str(eAmin)
    maxlen = str(eAmax)
elif pop == "eE":
    minlen = str(eEmin)
    maxlen = str(eEmax)
elif pop == "custom":
    minlen = str(customSize)
    maxlen = str(customSize)
else:
    print "whoops, something went wrong here, horribly horribly wrong!"
    # pop = ftsize


### for already shifted samples
def avggene_riboshift_plot_overlay(alignposition, ribosome_site, fiveorthreeprime='5', normalization='uneq', threshold ="0"):
    # read_length_list = ftsize
    alignpos = alignposition # '1' for start, '2' for stop
    assignment = fiveorthreeprime # should be 5' mapped at this point
    norm = normalization # 'uneq' for no normalizaiton; 'eq' to give equal weight to all genes
    ribosome_shift = ribosome_site # 'A' or 'P' or maybe 'E' later

    sample_numb = len(samples)
    sample_index = range(0,len(samples))
    sample_dict = {}
    for i in sample_index:
        sample_dict[samples[i]] = sample_index[i]

    df_fl_list = []
    df_eA_list = []
    df_eE_list = []
    df_aL_list = []
    df_custom_list = []


    for file in samples:
        fp_assign_path = '%s' % (file)
        avggene_csv_path = "%s/avggene%s_ORF%s_%sshift_%s%s150" % (fp_assign_path, alignpos, norm_type, ribosome_shift, assignment, norm) # norm should be 'uneq' for now


        fl_avggene_csv = '%s/%s_fl_rpkmThresh%s_%sto%sf_avg_%s_cdsNorm.csv' % (avggene_csv_path, file, threshold, flmin, flmax, alignpos)
        # eA_avggene_csv = '%s/%s_eA_rpkmThresh%s_%sto%sf_avg_%s_cdsNorm.csv' % (avggene_csv_path, file, threshold, eAmin, eAmax, alignpos)
        # eE_avggene_csv = '%s/%s_eE_rpkmThresh%s_%sto%sf_avg_%s_cdsNorm.csv' % (avggene_csv_path, file, threshold, eEmin, eEmax, alignpos)
        # aL_avggene_csv = '%s/%s_rpkmThresh%s_aL_aLf_avg_%s_cdsNorm.csv' % (avggene_csv_path, file, threshold, alignpos)

        if pop == "custom":
            custom_avggene_csv = '%s/%s_%s_shiftCustom_rpkmThresh0_%sto%sf_avg_%s.csv' % (avggene_csv_path, file, customSize, minlen, maxlen, alignpos)
            print custom_avggene_csv
            custom_avggene_df = pd.read_csv(custom_avggene_csv, index_col = 0, header=0)
            df_custom_list.append(custom_avggene_df)

        fl_avggene_df = pd.read_csv(fl_avggene_csv, index_col = 0, header=0)


        df_fl_list.append(fl_avggene_df)


    plot_outfile = "Fig2S3A.pdf"
    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(6,4))


    counter = 0 ## for making offsets
    if pop == 'fl':
        for i in range(len(df_fl_list)):
            df_fl_list[i].plot.line(x=df_fl_list[i].index+counter,
                                    y=df_fl_list[i].columns.values[0],
                                    ax=ax, color = colorList[i],
                                    lw=1,
                                    use_index=True,
                                    label=sample_plot_names[i])
            counter +=6
    ax.set_xlim(-50,100)
    plt.legend(loc=1, prop={'size': 6})
    plt.savefig(plot_outfile, format = 'pdf', bbox_inches = "tight")
    plt.close()


def main():
    avggene_riboshift_plot_overlay(alignposition, ribosome_site, normalization='eq', threshold=threshold)

if __name__ == '__main__':
    main()

"""
}


process figure2S3B {

publishDir params.outdir, overwrite: true, mode: 'copy',
	saveAs: {filename ->
	if (filename =~ /.*.pdf$/) "reports/$filename"
}

input:
 file nameAll from g_6_outputDir_g_13.collect()

output:
 file "*.pdf"  into g_13_outputFilePdf

script:
nameAll = nameAll.collect{ '"' + it + '"'}
"""
#!/usr/bin/env python
# This script plots codon occupancies relative to untreated cells
import matplotlib
matplotlib.use('Agg') ## set backend here
import matplotlib.pyplot as plt 
import matplotlib.patches as mpatches
plt.rcParams['pdf.fonttype'] = 42 # this keeps most text as actual text in PDFs, not outlines

import sys, os
import pandas as pd
import pysam
import numpy as np
import math
import argparse
import importlib
from scipy import stats

readlength = "28to35"
shift = "A"

## indexes
samplelist = ${nameAll}
ctrl_sample = '${control_sample_prefix}'
compare_list = samplelist[:]
compare_list.remove(ctrl_sample)

## colors
black = '#000000'
orange = '#ffb000'
cyan = '#63cfff'
red = '#eb4300'
green = '#00c48f'
pink = '#eb68c0'
yellow = '#fff71c'
blue = '#006eb9'

colorDict = {
	'D':cyan,
	'G':orange,
	'I':green
}

codonDict = {
	'A':['GCT', 'GCC', 'GCA', 'GCG'],
	'C':['TGT', 'TGC'],
	'D':['GAT', 'GAC'],
	'E':['GAA', 'GAG'],
	'F':['TTT', 'TTC'],
	'G':['GGT', 'GGC', 'GGA', 'GGG'],
	'H':['CAT', 'CAC'],
	'I':['ATT', 'ATC', 'ATA'],
	'K':['AAA', 'AAG'],
	'L':['CTT', 'CTC', 'CTA', 'CTG', 'TTA', 'TTG'],
	'M':['ATG'],
	'N':['AAT', 'AAC'],
	'P':['CCT', 'CCC', 'CCA', 'CCG'],
	'Q':['CAA', 'CAG'],
	'R':['AGA','AGG','CGT','CGC','CGA','CGG'],
	'S':['TCT', 'TCC', 'TCA', 'TCG'],
	'T':['ACT', 'ACC', 'ACA', 'ACG'],
	'V':['GTT', 'GTC', 'GTA', 'GTG'],
	'W':['TGG'],
	'Y':['TAT', 'TAC']
}



def log_trans_b10(x):
	try:
		return math.log(x, 10)
	except:
		return float(-6.00)
#         return float("NaN")
def log_trans_b2(x):
	try:
		return math.log(x, 2)
	except:
#         return float("NaN")
		return float(-15.00) # set arbitrarily low value

def corrfunc(x, y, **kws):
	r, _ = stats.pearsonr(x, y)
	rho, _ = stats.spearmanr(x, y)
	ax = plt.gca()
	ax.annotate("R = {:.3f}".format(r),
				xy=(.1, .9), xycoords=ax.transAxes)


def plot_codons():
	counter = 0
	for i in range(len(compare_list)):
		comp_sample = compare_list[i]
		df1 = pd.read_csv("%s/codon/%s_%sshift_%s__5occupancy_15cds5trim_15cds3trim_codonOccupancy.csv" % 
						 (ctrl_sample, ctrl_sample, shift, readlength), index_col=0)
		df2 = pd.read_csv("%s/codon/%s_%sshift_%s__5occupancy_15cds5trim_15cds3trim_codonOccupancy.csv" % 
						 (comp_sample, comp_sample, shift, readlength), index_col=0)

		dflist = [df1, df2]

		dfout = pd.concat(dflist, axis=1)
		dfout.drop(index=['TAA', 'TAG', 'TGA'],axis=0,inplace=True)
		dfout = dfout[[ctrl_sample,comp_sample]].apply(pd.to_numeric)

		for col in dfout.columns:
			dfout[col+'_log2'] = dfout[col].apply(log_trans_b2)


		dfout['ctrl_mean'] = dfout[[ctrl_sample+'_log2']].mean(axis=1)
		dfout['ctrl_std'] = dfout[[ctrl_sample+'_log2']].std(axis=1)
		dfout['ctrl_sem'] = dfout[[ctrl_sample+'_log2']].sem(axis=1)

		dfout['treat_mean'] = dfout[[comp_sample+'_log2']].mean(axis=1)
		dfout['treat_std'] = dfout[[comp_sample+'_log2']].std(axis=1)
		dfout['treat_sem'] = dfout[[comp_sample+'_log2']].sem(axis=1)



		### plotting:
		fig, ax = plt.subplots(figsize=(6,6))

		for aa in codonDict:
			for cdn in codonDict[aa]:

				plt.scatter(x=dfout.loc[cdn,'ctrl_mean'], 
							y=dfout.loc[cdn,'treat_mean'], 
							color="black", 
							edgecolor="black", 
							linewidth = 0.5, 
							s=20
						   )
				plt.errorbar(x=dfout.loc[cdn,'ctrl_mean'], 
							y=dfout.loc[cdn,'treat_mean'], 
							xerr = dfout.loc[cdn,'ctrl_std'],
							yerr = dfout.loc[cdn,'treat_std'],
							color="black", 
							ecolor="black", 
							linewidth = 0.75
							)	

		for aa in colorDict:
			for cdn in codonDict[aa]:
	#           
				plt.scatter(x=dfout.loc[cdn,'ctrl_mean'], 
							y=dfout.loc[cdn,'treat_mean'], 
							color=colorDict[aa], 
							edgecolor="black", 
							linewidth = 0.5, 
							s=20
						   )

				plt.text(x=dfout.loc[cdn,'ctrl_mean']+0.00, 
						 y=dfout.loc[cdn,'treat_mean']+0.00, 
						 color=colorDict[aa], 
						 s=aa, 
						 fontsize = 16)
				plt.text(x=dfout.loc[cdn,'ctrl_mean']+0.00, 
						 y=dfout.loc[cdn,'treat_mean']-0.00, 
						 color=colorDict[aa], #"black", 
						 s=cdn, 
						 fontsize = 10, 
						 va="top", 
						 ha="left")

		plt.title(ctrl_sample+' vs '+comp_sample, color='gray')

		lims = [
			np.min([-1.0, -1.0]),
			np.max([1.5, 1.5]),
		]

		ax.plot(lims, lims, 'k--', alpha=0.75, zorder=0)
		ax.set_aspect('equal')
		ax.set_xlim(lims)
		ax.set_ylim(lims)
		corrfunc(dfout['ctrl_mean'], dfout['treat_mean'])

		outfile = "Fig2S3B_%s_vs_%s.pdf" % (ctrl_sample, comp_sample)
		plt.savefig(outfile, format="pdf")

		counter +=1


def main():
	plot_codons()

if __name__ == '__main__':
	main()

"""
}


process Adapter_Trimmer_Quality_Module_Quality_Filtering_Summary {

input:
 file logfile from g0_20_log_file_g0_16.collect()
 val mate from g_1_mate_g0_16

output:
 file "quality_filter_summary.tsv"  into g0_16_outputFileTSV_g_26

shell:
'''
#!/usr/bin/env perl
use List::Util qw[min max];
use strict;
use File::Basename;
use Getopt::Long;
use Pod::Usage;
use Data::Dumper;

my @header;
my %all_files;
my %tsv;
my %headerHash;
my %headerText;

my $i = 0;
chomp( my $contents = `ls *.log` );
my @files = split( /[\\n]+/, $contents );
foreach my $file (@files) {
    $i++;
    my $mapper   = "";
    my $mapOrder = "1";
    if ($file =~ /(.*)\\.fastx_quality\\.log/){
        $mapper   = "fastx";
        $file =~ /(.*)\\.fastx_quality\\.log/;
        my $name = $1;    ##sample name
        push( @header, $mapper );
        my $in;
        my $out;
        chomp( $in =`cat $file | grep 'Input:' | awk '{sum+=\\$2} END {print sum}'` );
        chomp( $out =`cat $file | grep 'Output:' | awk '{sum+=\\$2} END {print sum}'` );
        $tsv{$name}{$mapper} = [ $in, $out ];
        $headerHash{$mapOrder} = $mapper;
        $headerText{$mapOrder} = [ "Total Reads", "Reads After Quality Filtering" ];
    } elsif ($file =~ /(.*)\\.trimmomatic_quality\\.log/){
        $mapper   = "trimmomatic";
        $file =~ /(.*)\\.trimmomatic_quality\\.log/;
        my $name = $1;    ##sample name
        push( @header, $mapper );
        my $in;
        my $out;
        if ( "!{mate}" eq "pair"){
            chomp( $in =`cat $file | grep 'Input Read Pairs:' | awk '{sum+=\\$4} END {print sum}'` );
            chomp( $out =`cat $file | grep 'Input Read Pairs:' | awk '{sum+=\\$7} END {print sum}'` );
        } else {
            chomp( $in =`cat $file | grep 'Input Reads:' | awk '{sum+=\\$3} END {print sum}'` );
            chomp( $out =`cat $file | grep 'Input Reads:' | awk '{sum+=\\$5} END {print sum}'` );
        }
        $tsv{$name}{$mapper} = [ $in, $out ];
        $headerHash{$mapOrder} = $mapper;
        $headerText{$mapOrder} = [ "Total Reads", "Reads After Quality Filtering" ];
    }
    
}

my @mapOrderArray = ( keys %headerHash );
my @sortedOrderArray = sort { $a <=> $b } @mapOrderArray;

my $summary          = "quality_filter_summary.tsv";
writeFile( $summary,          \\%headerText,       \\%tsv );

sub writeFile {
    my $summary    = $_[0];
    my %headerText = %{ $_[1] };
    my %tsv        = %{ $_[2] };
    open( OUT, ">$summary" );
    print OUT "Sample\\t";
    my @headArr = ();
    for my $mapOrder (@sortedOrderArray) {
        push( @headArr, @{ $headerText{$mapOrder} } );
    }
    my $headArrAll = join( "\\t", @headArr );
    print OUT "$headArrAll\\n";

    foreach my $name ( keys %tsv ) {
        my @rowArr = ();
        for my $mapOrder (@sortedOrderArray) {
            push( @rowArr, @{ $tsv{$name}{ $headerHash{$mapOrder} } } );
        }
        my $rowArrAll = join( "\\t", @rowArr );
        print OUT "$name\\t$rowArrAll\\n";
    }
    close(OUT);
}

'''
}


process Adapter_Trimmer_Quality_Module_Adapter_Removal_Summary {

input:
 file logfile from g0_18_log_file_g0_11.collect()
 val mate from g_1_mate_g0_11

output:
 file "adapter_removal_summary.tsv"  into g0_11_outputFileTSV_g_26
 file "adapter_removal_detailed_summary.tsv" optional true  into g0_11_outputFile

shell:
'''
#!/usr/bin/env perl
use List::Util qw[min max];
use strict;
use File::Basename;
use Getopt::Long;
use Pod::Usage;
use Data::Dumper;

my @header;
my %all_files;
my %tsv;
my %tsvDetail;
my %headerHash;
my %headerText;
my %headerTextDetail;

my $i = 0;
chomp( my $contents = `ls *.log` );

my @files = split( /[\\n]+/, $contents );
foreach my $file (@files) {
    $i++;
    my $mapOrder = "1";
    if ($file =~ /(.*)\\.fastx\\.log/){
        $file =~ /(.*)\\.fastx\\.log/;
        my $mapper   = "fastx";
        my $name = $1;    ##sample name
        push( @header, $mapper );

        my $in;
        my $out;
        my $tooshort;
        my $adapteronly;
        my $noncliped;
        my $Nreads;

        chomp( $in =`cat $file | grep 'Input:' | awk '{sum+=\\$2} END {print sum}'` );
        chomp( $out =`cat $file | grep 'Output:' | awk '{sum+=\\$2} END {print sum}'` );
        chomp( $tooshort =`cat $file | grep 'too-short reads' | awk '{sum+=\\$2} END {print sum}'`);
        chomp( $adapteronly =`cat $file | grep 'adapter-only reads' | awk '{sum+=\\$2} END {print sum}'`);
        chomp( $noncliped =`cat $file | grep 'non-clipped reads.' | awk '{sum+=\\$2} END {print sum}'`);
        chomp( $Nreads =`cat $file | grep 'N reads.' | awk '{sum+=\\$2} END {print sum}'` );

        $tsv{$name}{$mapper} = [ $in, $out ];
        $headerHash{$mapOrder} = $mapper;
        $headerText{$mapOrder} = [ "Total Reads", "Reads After Adapter Removal" ];
        $tsvDetail{$name}{$mapper} = [ $in, $tooshort, $adapteronly, $noncliped, $Nreads, $out ];
        $headerTextDetail{$mapOrder} = ["Total Reads","Too-short reads","Adapter-only reads","Non-clipped reads","N reads","Reads After Adapter Removal"];
    } elsif ($file =~ /(.*)\\.trimmomatic\\.log/){
        $file =~ /(.*)\\.trimmomatic\\.log/;
        my $mapper   = "trimmomatic";
        my $name = $1;    ##sample name
        push( @header, $mapper );
        
        my $in;
        my $out;

        if ( "!{mate}" eq "pair"){
            chomp( $in =`cat $file | grep 'Input Read Pairs:' | awk '{sum+=\\$4} END {print sum}'` );
            chomp( $out =`cat $file | grep 'Input Read Pairs:' | awk '{sum+=\\$7} END {print sum}'` );
        } else {
            chomp( $in =`cat $file | grep 'Input Reads:' | awk '{sum+=\\$3} END {print sum}'` );
            chomp( $out =`cat $file | grep 'Input Reads:' | awk '{sum+=\\$5} END {print sum}'` );
        }
        


        $tsv{$name}{$mapper} = [ $in, $out ];
        $headerHash{$mapOrder} = $mapper;
        $headerText{$mapOrder} = [ "Total Reads", "Reads After Adapter Removal" ];
        
    }
    
}

my @mapOrderArray = ( keys %headerHash );
my @sortedOrderArray = sort { $a <=> $b } @mapOrderArray;

my $summary          = "adapter_removal_summary.tsv";
my $detailed_summary = "adapter_removal_detailed_summary.tsv";
writeFile( $summary,          \\%headerText,       \\%tsv );
if (%headerTextDetail){
    writeFile( $detailed_summary, \\%headerTextDetail, \\%tsvDetail );  
}

sub writeFile {
    my $summary    = $_[0];
    my %headerText = %{ $_[1] };
    my %tsv        = %{ $_[2] };
    open( OUT, ">$summary" );
    print OUT "Sample\\t";
    my @headArr = ();
    for my $mapOrder (@sortedOrderArray) {
        push( @headArr, @{ $headerText{$mapOrder} } );
    }
    my $headArrAll = join( "\\t", @headArr );
    print OUT "$headArrAll\\n";

    foreach my $name ( keys %tsv ) {
        my @rowArr = ();
        for my $mapOrder (@sortedOrderArray) {
            push( @rowArr, @{ $tsv{$name}{ $headerHash{$mapOrder} } } );
        }
        my $rowArrAll = join( "\\t", @rowArr );
        print OUT "$name\\t$rowArrAll\\n";
    }
    close(OUT);
}

'''
}

g_22_outputFileTSV_g_26= g_22_outputFileTSV_g_26.ifEmpty([""]) 
g_28_outputFileTSV_g_26= g_28_outputFileTSV_g_26.ifEmpty([""]) 
g0_11_outputFileTSV_g_26= g0_11_outputFileTSV_g_26.ifEmpty([""]) 
g0_21_outputFileTSV_g_26= g0_21_outputFileTSV_g_26.ifEmpty([""]) 
g0_16_outputFileTSV_g_26= g0_16_outputFileTSV_g_26.ifEmpty([""]) 

//* autofill
//* platform
if ($HOSTNAME == "ghpcc06.umassrc.org"){
    $TIME = 30
    $CPU  = 1
    $MEMORY = 10
    $QUEUE = "short"
}
//* platform
//* autofill

process Overall_Summary {

publishDir params.outdir, overwrite: true, mode: 'copy',
	saveAs: {filename ->
	if (filename =~ /overall_summary.tsv$/) "Summary/$filename"
}

input:
 file starSum from g_22_outputFileTSV_g_26
 file rsemSum from g_28_outputFileTSV_g_26
 file adapterSum from g0_11_outputFileTSV_g_26
 file trimmerSum from g0_21_outputFileTSV_g_26
 file qualitySum from g0_16_outputFileTSV_g_26

output:
 file "overall_summary.tsv"  into g_26_outputFileTSV

shell:
'''
#!/usr/bin/env perl
use List::Util qw[min max];
use strict;
use File::Basename;
use Getopt::Long;
use Pod::Usage;
use Data::Dumper;

my @header;
my %all_rows;
my @seen_cols;
my $ID_header;

chomp(my $contents = `ls *.tsv`);
my @rawFiles = split(/[\\n]+/, $contents);
my @files = ();
# order must be in this order for chipseq pipeline: bowtie->dedup
# rsem bam pipeline: dedup->rsem, star->dedup
# riboseq ncRNA_removal->star
my @order = ("adapter_removal","trimmer","quality","extractUMI","sequential_mapping","ncRNA_removal","bowtie","star","hisat2","tophat2", "dedup","rsem");
for ( my $k = 0 ; $k <= $#order ; $k++ ) {
    for ( my $i = 0 ; $i <= $#rawFiles ; $i++ ) {
        if ( $rawFiles[$i] =~ /$order[$k]/ ) {
            push @files, $rawFiles[$i];
        }
    }
}

print Dumper \\@files;
##add rest of the files
for ( my $i = 0 ; $i <= $#rawFiles ; $i++ ) {
    push(@files, $rawFiles[$i]) unless grep{$_ == $rawFiles[$i]} @files;
}
print Dumper \\@files;

##Merge each file according to array order

foreach my $file (@files){
        open IN,"$file";
        my $line1 = <IN>;
        chomp($line1);
        ( $ID_header, my @header) = ( split("\\t", $line1) );
        push @seen_cols, @header;

        while (my $line=<IN>) {
        chomp($line);
        my ( $ID, @fields ) = ( split("\\t", $line) ); 
        my %this_row;
        @this_row{@header} = @fields;

        #print Dumper \\%this_row;

        foreach my $column (@header) {
            if (! exists $all_rows{$ID}{$column}) {
                $all_rows{$ID}{$column} = $this_row{$column}; 
            }
        }   
    }
    close IN;
}

#print for debugging
#print Dumper \\%all_rows;
#print Dumper \\%seen_cols;

#grab list of column headings we've seen, and order them. 
my @cols_to_print = uniq(@seen_cols);
my $summary = "overall_summary.tsv";
open OUT, ">$summary";
print OUT join ("\\t", $ID_header,@cols_to_print),"\\n";
foreach my $key ( keys %all_rows ) { 
    #map iterates all the columns, and gives the value or an empty string. if it's undefined. (prevents errors)
    print OUT join ("\\t", $key, (map { $all_rows{$key}{$_} // '' } @cols_to_print)),"\\n";
}
close OUT;

sub uniq {
    my %seen;
    grep ! $seen{$_}++, @_;
}

'''


}

params.run_FastQC =  "no"  //* @dropdown @options:"yes","no" @description:"FastQC provides quality control checks on raw sequence data."



process Adapter_Trimmer_Quality_Module_FastQC {

publishDir params.outdir, overwrite: true, mode: 'copy',
	saveAs: {filename ->
	if (filename =~ /.*.(html|zip)$/) "FastQC/$filename"
}

input:
 val mate from g_1_mate_g0_3
 set val(name), file(reads) from g_2_reads_g0_3

output:
 file '*.{html,zip}'  into g0_3_FastQCout

errorStrategy 'retry'
maxRetries 3

when:
(params.run_FastQC && (params.run_FastQC == "yes"))

script:
nameAll = reads.toString()
if (nameAll.contains('.gz')) {
    file =  nameAll - '.gz' - '.gz'
    runGzip = "ls *.gz | xargs -i echo gzip -df {} | sh"
} else {
    file =  nameAll 
    runGzip = ''
}
"""
${runGzip}
fastqc ${file} 
"""
}


workflow.onComplete {
println "##Pipeline execution summary##"
println "---------------------------"
println "##Completed at: $workflow.complete"
println "##Duration: ${workflow.duration}"
println "##Success: ${workflow.success ? 'OK' : 'failed' }"
println "##Exit status: ${workflow.exitStatus}"
}
