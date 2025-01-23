#!/usr/bin/env perl
use strict;
use warnings;
use File::Basename; 
###########################################
####Author Yangzi, wangyz.benniao@gmail.com
###########################################

###############################################################################
# This script processes a BLAST output (fmt=0) and an original VCF file.
# 1) It parses BLAST alignment results to find alignments covering position 151
#    (the "middle" base of a 301 bp query).
# 2) It outputs simplified alignment info and identifies the genotype at that
#    position for both the query and subject sequences.
# 3) It then filters and compares these genotypes with the original VCF,
#    writing a new VCF plus several summary files.
#
# USAGE:
#   perl script.pl <blast_output0> <original_vcf> <out_dir>
#
# ARGUMENTS:
#   ARGV[0] = BLAST result file (fmt=0)
#   ARGV[1] = Original VCF file
#
# OUTPUT FILES:
#   1) <blast_output0>.simplified.4ck.tmp      # intermediate
#   2) <blast_output0>.genotype.4ck.tsv        # alignment-based genotype info
#   3) <blast_output0>.4outgroup.vcf           # filtered VCF lines
#   4) <blast_output0>.final.4outgroup.tsv     # final list of query vs outgroup genotype
#   5) <blast_output0>.outgroup.fa             # outgroup sequence in FASTA
#   6) <blast_output0>.count.tsv               # summary of sites found/identical
###############################################################################

die "Usage: perl $0 <blast_output0> <original_vcf> <out_dir> \n" if (@ARGV < 3);

# Minimal aligned sequence length (for filtering) in the final step
my $min_len = 50;
my $out_dir = $ARGV[2];
chomp $out_dir;

my $input_base = basename($ARGV[0]);
###############################
# 1) Parse BLAST output (fmt=0)
#    and extract alignments that cover
#    the position 151 in a 301 bp query.
###############################

# Query= ChrS01_13311
# Length=301
# ***** No hits found *****
# Lambda      K        H
#     1.33    0.621     1.12
# Gapped
# Lambda      K        H
#     1.28    0.460    0.850
# Effective search space used: 656737064256
# Query= ChrS01_13356
# Length=301
#                                                                       Score        E
# Sequences producing significant alignments:                          (Bits)     Value
# LG03                                                                  60.2       5e-07
# >LG03
# Length=187626166
#  Score = 60.2 bits (32),  Expect = 5e-07
#  Identities = 69/86 (80%), Gaps = 6/86 (7%)
#  Strand=Plus/Minus
# Query  218       GAACAACGTAAATCCGTGTGTCTAGTTT-TAGAGCGAATTGA--CT-TTTCTTCCCTCTC  273
#                  |||| |||||||||| ||||| || | | | | | | |||||  ||  ||||||||||||
# Sbjct  13993150  GAACCACGTAAATCCTTGTGT-TATTATCTTGTGTG-ATTGATCCTCATTCTTCCCTCTC  13993093
# Query  274       TTTTCTTCTTGTCGTGCAAATTAAGA  299
#                  |||||||||||||||||  || ||||
# Sbjct  13993092  TTTTCTTCTTGTCGTGCGGATCAAGA  13993067
# Lambda      K        H
#     1.33    0.621     1.12
# Gapped
# Lambda      K        H
#     1.28    0.460    0.850
# Effective search space used: 656737064256

open my $OUT0,  "<",  $ARGV[0]
    or die "Cannot open BLAST res file $ARGV[0]: $!\n";

open my $TMP,   ">",  "$out_dir/$input_base.simplified.4ck.tmp"
    or die "Cannot write temp file: $!\n";

my ($key, $ref_chr, $bit, $E, $ID);

while (<$OUT0>) {
    chomp;
    next if (/^\s*$/);

    # Check for new query
    if (/^Query=\s*(\w+)/) {
        $key = $1;
        next;
    }

    # If line contains "*****", it often indicates "No hits found", reset $key
    if (/\*+/) {
        $key = undef;
        next;
    }

    # Capture subject chromosome or reference ID, test if it can be recognized
    if (/^>/) {
        if (/^>([\w\.\:\-]+)$/) {
            $ref_chr = $1;
            next;
        } else {
            die "reference seq id \"$_\" cannot be recognized: \\w\\.\\:\\-\n";
        }
    }

    # Capture Score and E-value lines
    if (/^\s*Score\s*=\s*([^\s]+)\s*bits[^,]+,\s*Expect\s*=\s*([^\s]+)/) {
        $bit = $1;
        $E   = $2;
        next;
    }

    # Capture the identities info (e.g., "Identities = 69/86")
    if (/^\s*Identities\s+=\s+([^\s]+)/) {
        $ID = $1;
        next;
    }

    # Lines that start with "Query <start> ... <end>"
    # We specifically check if it covers position 151 in the 301 bp query
    if (/^(Query)\s+\d+/) {
        my @fields = split /\s+/, $_;

        # If query alignment range includes 151
        if ($fields[1] <= 151 and $fields[3] >= 151) {
            # Print partial info to the temp file
            print $TMP ">$key\t$ref_chr\t$bit\t$ID\n";
            print $TMP $_,"\n";

            # Grab the next two lines (typically the alignment lines)
            my $tmp1 = <$OUT0>;
            my $tmp2 = <$OUT0>;
            chomp($tmp1);
            chomp($tmp2);
            print $TMP "$tmp1\n$tmp2\n";
        }
    }
}

close $OUT0;
close $TMP;

###############################
# 2) Analyze the extracted alignment
#    and identify the genotype (single
#    base at position 151).
###############################
open $TMP, "<", "$out_dir/$input_base.simplified.4ck.tmp"
    or die "Cannot open temp file for reading: $!\n";

open my $OUT, ">", "$out_dir/$input_base.genotype.4ck.tsv"
    or die "Cannot write genotype TSV: $!\n";

print $OUT "q_chr\tq_pos\ts_chr\ts_start\ts_end\tq_genotype\ts_genotype\tbit_score\tidentity\n";

while (<$TMP>) {
    chomp;
    # Lines that begin with '>' hold metadata: ">ChrS01_13487 LG03 60.2 69/86"
    # We then read the next three lines: the query alignment, an alignment line, and the sbjct alignment
    if (/^\>(.*)$/) {
        my $name = $1;
        $name =~ s/_/\t/;                # split query name into two parts
        my @c = split /\t+/, $name;      # [0] => q_chr, [1] => q_pos, [2] => s_chr, [3] => bit_score, [4] => identity

        my $q_line = <$TMP>;            # "Query   start  ... seq ...  end"
        my $tmp_l  = <$TMP>;            # alignment line of '|' or ' '
        my $s_line = <$TMP>;            # "Sbjct   start  ... seq ...  end"

        chomp $q_line;
        chomp $tmp_l;
        chomp $s_line;

        my @q_fields  = split /\s+/, $q_line;
        my @qseq      = split //, $q_fields[2];     # the actual query sequence
        my @s_fields  = split /\s+/, $s_line;
        my @sseq      = split //, $s_fields[2];     # the subject sequence

        # Basic format check
        die "Format error:\nq:$q_line\ns:$s_line\n" 
            if ($q_line !~ /^Query/ || $s_line !~ /^Sbjct/ || @q_fields != 4 || @s_fields != 4);

        # We want the genotype at position 151. But alignment can contain gaps.
        # The while-loop logic ensures we skip any gaps up to that "middle" base.
        my $n = 151 - ($q_fields[1] - 1);
        my ($count_seq, $count_gap, $loop) = (0, 0, $n);

      A: while (1) {
            for my $i (0..($loop-1)) {
                if ($qseq[$i] =~ /[^\-]/) {
                    $count_seq++;
                } else {
                    $count_gap++;
                }
            }
            # If not enough real bases, extend loop
            if ($count_seq < $n) {
                $loop      = $count_gap + $n;
                $count_gap = 0;
                $count_seq = 0;
            } else {
                last A;
            }
        }

        # The genotype from the query and subject at position 151
        my $q_genotype = $qseq[$loop-1];
        my $s_genotype = $sseq[$loop-1];

        # Print final line: q_chr, q_pos, s_chr, s_start, s_end, q_genotype, s_genotype, bit_score, identity
        print $OUT join("\t",
                        $c[0],          # q_chr
                        $c[1],          # q_pos
                        $c[2],          # s_start
                        $s_fields[1],   # s_chr or s_start?
                        $s_fields[3],   # s_end
                        $q_genotype,
                        $s_genotype,
                        $c[3],          # bit_score
                        $c[4]           # identity
               ), "\n";
    }
}
close $TMP;
close $OUT;

###############################
# 3) Filter based on minimal
#    alignment length, store
#    best bitscore.
###############################
my %hash;
open my $TSV, "<", "$out_dir/$input_base.genotype.4ck.tsv"
    or die "Could not open genotype.4ck.tsv: $!\n";

# Skip header line
<$TSV>;

while (<$TSV>) {
    chomp;
    my @a = split /\t+/;
    # a[0]: q_chr
    # a[1]: q_pos
    # a[2]: s_chr
    # a[3]: s_start
    # a[4]: s_end
    # a[5]: q_genotype
    # a[6]: s_genotype
    # a[7]: bit_score
    # a[8]: identity (like 69/86)

    next if ($a[6] eq "-");   # skip if outgroup genotype is gap

    # parse alignment length from identity "69/86"
    my $len;
    if ($a[8] =~ /\d+\/(\d+)/) {
        $len = $1;
    } else {
        die "Wrong format in identity field: $_\n";
    }
    next if $len < $min_len;

    # If we already have an entry for [q_chr][q_pos],
    # keep the record with the higher bitscore
    if (exists $hash{$a[0]}{$a[1]}) {
        if ($a[7] > $hash{$a[0]}{$a[1]}[2]) {
            my @b = ($a[5], $a[6], $a[7]);
            $hash{$a[0]}{$a[1]} = \@b;
        }
    } else {
        my @b = ($a[5], $a[6], $a[7]);  # store q_genotype, s_genotype, bitscore
        $hash{$a[0]}{$a[1]} = \@b;
    }
}
close $TSV;

###############################
# 4) Compare with original VCF
#    to see how many sites
#    have outgroup info,
#    write final VCF and stats.
###############################
open my $VCF, "<", $ARGV[1]
    or die "Cannot open original VCF $ARGV[1]: $!\n";

open my $OUTVCF,  ">", "$out_dir/$input_base.4outgroup.vcf"
    or die "Cannot write outgroup VCF: $!\n";

open my $OUTTSV,  ">", "$out_dir/$input_base.final.4outgroup.tsv"
    or die "Cannot write final outgroup TSV: $!\n";

open my $OUTCOUNT, ">", "$out_dir/$input_base.count.tsv"
    or die "Cannot write count file: $!\n";

open my $OUTFA,   ">", "$out_dir/$input_base.outgroup.fa"
    or die "Cannot write outgroup.fa: $!\n";
print $OUTFA ">outgroup\n";

my $vcf_sites       = 0;  # total lines of VCF with variant info
my $found_sites     = 0;  # how many found in outgroup
my $identical_sites = 0;  # outgroup genotype same as reference genotype in VCF

print $OUTTSV "q_chr\tq_pos\tq_genotype\ts_genotype\tbit_score\n";

while (<$VCF>) {
    chomp;
    if (/^#/) {
        # copy VCF header lines as is
        print $OUTVCF "$_\n";
        next;
    }
    $vcf_sites++;
    my @a = split /\t/;
    # a[0] => q_chr
    # a[1] => q_pos
    # a[3] => REF allele
    # a[4] => ALT allele
    # etc.

    if (exists $hash{$a[0]}{$a[1]}) {
        $found_sites++;
        my $char = uc($hash{$a[0]}{$a[1]}[0]);  # q_genotype from alignment
        $a[3] = uc($a[3]);                     # REF in VCF

        # Check if outgroup base agrees with VCF REF
        die "Reference mismatch for $a[0]\t$a[1]\t(Outgroup=$char,VCF REF=$a[3])\n" if ($char ne $a[3]);

        print $OUTVCF "$_\n";  # keep this VCF line

        # Append outgroup genotype sequence
        print $OUTFA uc($hash{$a[0]}{$a[1]}[1]);

        my @b = @{$hash{$a[0]}{$a[1]}};
        $b[0] = uc($b[0]);  # q_genotype
        $b[1] = uc($b[1]);  # s_genotype

        # If outgroup genotype is same as query's REF genotype => identical
        if ($b[0] eq $b[1]) {
            $identical_sites++;
        }
        # b[0], b[1], b[2] => q_genotype, s_genotype, bitscore
        my $tmp = join("\t", @b);
        print $OUTTSV "$a[0]\t$a[1]\t$tmp\n";
    }
}

print $OUTFA "\n";

# Print summary counts
my $p;
print $OUTCOUNT "Num of SNPs in VCF:\t$vcf_sites\n";
print $OUTCOUNT "Num of orthologous SNPs found in outgroup:\t$found_sites\n";
print $OUTCOUNT "Num of identical alleles between ref_sp and OG_sp:\t$identical_sites\n";

$p = $found_sites / $vcf_sites * 100;
print $OUTCOUNT "Perc. of mutations found in outgroup (\%):\t$p\n";

$p = $identical_sites / $found_sites * 100 if ($found_sites > 0);
print $OUTCOUNT "Perc. of identical alleles (\%):\t$p\n";

# $p = $identical_sites / $vcf_sites * 100 if ($vcf_sites > 0);
# print $OUTCOUNT "Perc. of identical alleles in all sites:\t$p\n";

close $OUTCOUNT;
close $VCF;
close $OUTVCF;
close $OUTFA;
close $OUTTSV;

__END__
explain of blast output6 format:
 1.      qseqid  query (e.g., unknown gene) sequence id
 2.      sseqid  subject (e.g., reference genome) sequence id
 3.      pident  percentage of identical matches
 4.      length  alignment length (sequence overlap)
 5.      mismatch        number of mismatches
 6.      gapopen         number of gap openings
 7.      qstart  start of alignment in query
 8.      qend    end of alignment in query
 9.      sstart  start of alignment in subject
 10.     send    end of alignment in subject
 11.     evalue  expect value
 12.     bitscore        bit score
