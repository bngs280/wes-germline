# #!/usr/bin/perl -w
# use strict;
# use Parse::CSV;
# use Data::Dumper;

# # convert exomedepth CNV bed file to VCF format

# my $GT_ratio = 0.25;



# my $header = <<"END_HEADER";
# ##fileformat=VCFv4.1
# ##source=ExomeDepth
# ##reference=/hpc/cog_bioinf/GENOMES/Homo_sapiens.GRCh37.GATK.illumina/genome.fa
# ##contig=<ID=1,length=249250621>
# ##contig=<ID=2,length=243199373>
# ##contig=<ID=3,length=198022430>
# ##contig=<ID=4,length=191154276>
# ##contig=<ID=5,length=180915260>
# ##contig=<ID=6,length=171115067>
# ##contig=<ID=7,length=159138663>
# ##contig=<ID=8,length=146364022>
# ##contig=<ID=9,length=141213431>
# ##contig=<ID=10,length=135534747>
# ##contig=<ID=11,length=135006516>
# ##contig=<ID=12,length=133851895>
# ##contig=<ID=13,length=115169878>
# ##contig=<ID=14,length=107349540>
# ##contig=<ID=15,length=102531392>
# ##contig=<ID=16,length=90354753>
# ##contig=<ID=17,length=81195210>
# ##contig=<ID=18,length=78077248>
# ##contig=<ID=19,length=59128983>
# ##contig=<ID=20,length=63025520>
# ##contig=<ID=21,length=48129895>
# ##contig=<ID=22,length=51304566>
# ##contig=<ID=X,length=155270560>
# ##contig=<ID=Y,length=59373566>
# ##contig=<ID=MT,length=16569>
# ##ALT=<ID=CNV,Description="Copy number variable region">
# ##FILTER=<ID=q10,Description="Quality below 10">
# ##FILTER=<ID=L10kb,Description="Length shorter than 10kb">
# ##FORMAT=<ID=GT,Number=1,Type=String, Description="inferred Genotype">
# ##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">
# ##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the variant described in this record">
# ##INFO=<ID=NEXONS,Number=1,Type=Integer,Description="Exomdepth: number of exons included in this event">
# ##INFO=<ID=BF,Number=1,Type=Integer,Description="Exomdepth: BF (Bayes Factor) value. log10 value of the likelihood ratio">
# ##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the variant described in this record">
# ##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Number of reference positions spanned by this CNV">
# ##INFO=<ID=RATIO,Number=1,Type=Float,Description="ratio of observed reads over expected. This is used to determine genotype for deletions: homozygous deletions have ratio < $GT_ratio">
# END_HEADER

# my $file = $ARGV[0];
# my $sample = $file;
# $sample =~ s{\.[^.]+$}{};
# my $outfile = "$sample\.vcf";
# open OUT, ">$outfile" or die "cannot open outfile: $!\n";

# my @files = $file;


# #add Xcall file to array if present
# my $xcallfile = $file;
# $xcallfile =~ s/_calls/_Xcalls/;

# push (@files,$xcallfile) if (-e $xcallfile);
# #    my @files = ($file,$xcallfile);

# # print out vcf header
# print OUT $header;
# my $chrom_line="#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT";
# print OUT $chrom_line,"\t$sample\n";


# # parse data and write VCF lines
# foreach my $file (@files) {
#     my $objects = Parse::CSV->new(
# 	file	=>	$file,
#         names	=>	1,
#     );

#     while ( my $object = $objects->fetch ) {
# 	my $ALT;
#         if ($$object{'type'} eq 'duplication') {
# 		$ALT="DUP";
# 	}elsif ($$object{'type'} eq 'deletion') {
# 	    $ALT="DEL";
#         }
# 	my $svlen=($$object{'end'}-$$object{'start'});
# 	my $id;
# 	if ($$object{'Conrad.hg19'} eq "NA") {
# 	    $id = '.';
# 	}else{
# 	    $id = $$object{'Conrad.hg19'};
# 	}
# 	#try to determine genotype..
# 	#only for dels
# 	my $GT="0/1";
# 	$GT = "1/1" if ($$object{'reads.ratio'}<0.25);
	
#         print OUT join("\t",$$object{'chromosome'},$$object{'start'},"$id","N","<$ALT\>","1000","PASS","SVTYPE=$ALT\;END=".$$object{'end'}.";NEXONS=".$$object{'nexons'}.";BF=".$$object{'BF'}.";SVLEN=$svlen".";RATIO=".$$object{'reads.ratio'}),"\tGT\t$GT\n";
#     }
# }

#!/usr/bin/perl -w
use strict;
use warnings;
use Parse::CSV;
use Data::Dumper;

# convert exomedepth CNV bed file to VCF format

my $GT_ratio = 0.25;

my $header = <<"END_HEADER";
##fileformat=VCFv4.1
##source=ExomeDepth
##reference=/hpc/cog_bioinf/GENOMES/Homo_sapiens.GRCh37.GATK.illumina/genome.fa
##contig=<ID=1,length=249250621>
##contig=<ID=2,length=243199373>
##contig=<ID=3,length=198022430>
##contig=<ID=4,length=191154276>
##contig=<ID=5,length=180915260>
##contig=<ID=6,length=171115067>
##contig=<ID=7,length=159138663>
##contig=<ID=8,length=146364022>
##contig=<ID=9,length=141213431>
##contig=<ID=10,length=135534747>
##contig=<ID=11,length=135006516>
##contig=<ID=12,length=133851895>
##contig=<ID=13,length=115169878>
##contig=<ID=14,length=107349540>
##contig=<ID=15,length=102531392>
##contig=<ID=16,length=90354753>
##contig=<ID=17,length=81195210>
##contig=<ID=18,length=78077248>
##contig=<ID=19,length=59128983>
##contig=<ID=20,length=63025520>
##contig=<ID=21,length=48129895>
##contig=<ID=22,length=51304566>
##contig=<ID=X,length=155270560>
##contig=<ID=Y,length=59373566>
##contig=<ID=MT,length=16569>
##ALT=<ID=CNV,Description="Copy number variable region">
##FILTER=<ID=q10,Description="Quality below 10">
##FILTER=<ID=L10kb,Description="Length shorter than 10kb">
##FORMAT=<ID=GT,Number=1,Type=String, Description="inferred Genotype">
##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">
##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the variant described in this record">
##INFO=<ID=NEXONS,Number=1,Type=Integer,Description="Exomdepth: number of exons included in this event">
##INFO=<ID=BF,Number=1,Type=Integer,Description="Exomdepth: BF (Bayes Factor) value. log10 value of the likelihood ratio">
##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the variant described in this record">
##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Number of reference positions spanned by this CNV">
##INFO=<ID=RATIO,Number=1,Type=Float,Description="ratio of observed reads over expected. This is used to determine genotype for deletions: homozygous deletions have ratio < $GT_ratio">
END_HEADER

# Check input file
unless (@ARGV) {
    die "Usage: $0 <input_csv_file>\n";
}

my $file = $ARGV[0];
unless (-e $file) {
    die "Error: Input file '$file' not found\n";
}

# Prepare output filename
my $sample = $file;
$sample =~ s{\.[^.]+$}{};
my $outfile = "$sample.vcf";

open my $out_fh, ">", $outfile or die "Cannot open output file '$outfile': $!\n";

my @files = ($file);

# Add Xcall file to array if present
my $xcallfile = $file;
$xcallfile =~ s/_calls/_Xcalls/;
push (@files, $xcallfile) if (-e $xcallfile);

# Print VCF header
print $out_fh $header;
print $out_fh "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t$sample\n";

# Process each input file
foreach my $input_file (@files) {
    unless (-e $input_file) {
        warn "Warning: Input file '$input_file' not found - skipping\n";
        next;
    }

    my $parser = Parse::CSV->new(
        file => $input_file,
        names => 1,
    );

    unless ($parser) {
        warn "Warning: Failed to create parser for '$input_file' - skipping\n";
        next;
    }

    # Process each record
    while (my $object = $parser->fetch) {
        # Skip if required fields are missing
        unless (defined $$object{'chromosome'} && defined $$object{'start'} && defined $$object{'end'} && defined $$object{'type'}) {
            warn "Warning: Skipping incomplete record in $input_file - missing required fields\n";
            next;
        }
        
        # Determine ALT type
        my $ALT;
        if ($$object{'type'} eq 'duplication') {
            $ALT = "DUP";
        } elsif ($$object{'type'} eq 'deletion') {
            $ALT = "DEL";
        } else {
            warn "Warning: Skipping record with unknown type '$$object{'type'}' in $input_file\n";
            next;
        }
        
        # Calculate SV length
        my $svlen = $$object{'end'} - $$object{'start'};
        
        # Handle ID field
        my $id = '.';
        if (exists $$object{'Conrad.hg19'} && defined $$object{'Conrad.hg19'} && $$object{'Conrad.hg19'} ne "NA") {
            $id = $$object{'Conrad.hg19'};
        }
        
        # Set default values for optional fields
        my $nexons = exists $$object{'nexons'} ? $$object{'nexons'} : 0;
        my $bf = exists $$object{'BF'} ? $$object{'BF'} : 0;
        my $ratio = exists $$object{'reads.ratio'} ? $$object{'reads.ratio'} : 1;
        
        # Determine genotype
        my $GT = "0/1";
        $GT = "1/1" if ($ALT eq "DEL" && $ratio < $GT_ratio);
        
        # Prepare INFO field
        my $info = join(";",
            "SVTYPE=$ALT",
            "END=$$object{'end'}",
            "NEXONS=$nexons",
            "BF=$bf",
            "SVLEN=$svlen",
            "RATIO=$ratio"
        );
        
        # Print VCF record
        print $out_fh join("\t",
            $$object{'chromosome'},    # CHROM
            $$object{'start'},         # POS
            $id,                      # ID
            "N",                       # REF
            "<$ALT>",                  # ALT
            "1000",                    # QUAL
            "PASS",                   # FILTER
            $info,                    # INFO
            "GT",                     # FORMAT
            $GT                       # SAMPLE
        ), "\n";
    }
}

close $out_fh or warn "Warning: Failed to properly close output file: $!\n";

print "Successfully created VCF file: $outfile\n";
exit 0;
