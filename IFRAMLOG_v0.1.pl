#! /usr/bin/env perl
###########################################
####Author Yangzi, wangyz.benniao@gmail.com
###########################################

# this script is used for extracting flanking 150 bp (in total 301 bp) sequences from SNPs that were stored in VCF, and output a FASTA

use strict;
use warnings;
use Getopt::Long;
use Bio::DB::Fasta;
#use Sys::CPU;

my $vcf_file;
my $ref_fasta;
my $out_base;
my $output_dir;
#check available cores 
#$cpu_count ||= 8;  # protect
#my $blast_N_CPU = int($cpu_count * 0.8) > 8 ? 8 : int($cpu_count * 0.8);
my $blast_N_CPU = 4;
my $excluded_snps = "Excluded_SNP.log";
my $outgroup_ref_fasta;
my $UPDOWN_LEN     = 150;  # default 150bp-flanking seq on each side 
my $N_CUTOFF_RATIO = 0.5;  # exclude if # of 'N' > half the % region 
my $help;

GetOptions(
    "vcf=s"       => \$vcf_file,
    "ref=s"       => \$ref_fasta,
    "OG_ref=s"    => \$outgroup_ref_fasta,
    "out_dir=s"   => \$output_dir,
    "blast_C=i"   => \$blast_N_CPU,
    "out=s"       => \$out_base,
    "FLK_len=i"   => \$UPDOWN_LEN,
    "N_perc=f"    => \$N_CUTOFF_RATIO,
    "help"        => \$help
) or die usage();

print "Using $blast_N_CPU threads for BLAST, Modeltest and RAxML.\n";

# USAGE
if ( $help or not ($vcf_file and $ref_fasta  and $out_base and $output_dir and $outgroup_ref_fasta) ) {
    die usage();
}

die "ERROR: Input VCF file $vcf_file does not exist.\n" unless -e $vcf_file;
die "ERROR: Reference FASTA file $ref_fasta does not exist.\n" unless -e $ref_fasta;
die "ERROR: Outgroup reference FASTA file $outgroup_ref_fasta does not exist.\n" unless -e $outgroup_ref_fasta;

# output dir
mkdir $output_dir unless -d $output_dir;


#############################################
######################Main logic
#############################################
{
    ################################### 1) to read in reference genome
    my $db = Bio::DB::Fasta->new($ref_fasta)
        or die "ERROR: Cannot open reference FASTA '$ref_fasta': $!\n";

    # 1a) to read in VCF and open output flanking seq fasta
    open VCF,  "<$vcf_file" or die "ERROR: Cannot open VCF '$vcf_file': $!\n";
    $excluded_snps=$out_base.$excluded_snps;
    open EXCL, ">$output_dir/$excluded_snps" or die "ERROR: Cannot write to '$output_dir/$excluded_snps': $!\n";
    open OUT_FLK,  ">$output_dir/$out_base" or die "ERROR: Cannot write to '$output_dir/$out_base': $!\n";

    print EXCL "#CHROM\tPOS\tReason\n";

    # 1b) to parse VCF header
    #parse_vcf_header(\*VCF);

    # 1c) to parse VCF
    parse_vcf_variants(
        \*VCF,       # VCF
        \*EXCL,      # excluded_SNP
        \*OUT_FLK,       # out_ref
        $db,
        $UPDOWN_LEN,
        $N_CUTOFF_RATIO
    );

    close VCF;
    close EXCL;
    close OUT_FLK;

    print <<"DONE_MSG";

------------------------------------------------------------
1st step -- generating flanking sequences from SNPs is DONE.
------------------------------------------------------------

Flanking length: $UPDOWN_LEN
N cutoff ratio:  $N_CUTOFF_RATIO
Output FASTA:    $output_dir/$out_base
Excluded SNPs:   $output_dir/$excluded_snps
DONE_MSG
   

    ################################### 2) to blast flanking sequences with OUTGROUP
    # 2a) make blast db
    my $cmd_make_db = "makeblastdb -in $outgroup_ref_fasta -dbtype nucl ".
                  "-max_file_sz 3GB -logfile $outgroup_ref_fasta.log";

    run_command($cmd_make_db, "Failed to make blast db for $outgroup_ref_fasta");

    # 2b) blast
    my $cmd_blastn = "blastn -db $outgroup_ref_fasta ".
                     "-max_target_seqs 1 ".
                     "-query $output_dir/$out_base ".
                     "-outfmt 0 ".
                     "-evalue 1e-6 ".
                     "-num_threads $blast_N_CPU ".
                     "-out $output_dir/$out_base.blast.fmt0.res";

    run_command($cmd_blastn, "Failed to run blastn on $outgroup_ref_fasta");

    print <<"DONE_MSG";

-----------------------
2nd step -- blast DONE.
-----------------------

Output:     $output_dir/$out_base.blast.fmt0.res
DONE_MSG


    ################################### 3) to parse blast result
    #to check, and get outgroup fasta file (and also corresponding SNP list)
    
    # 3a) make blast db
    my $cmd_get_OG = "perl ./get_othologous_genotype.pl ".
                    "$output_dir/$out_base.blast.fmt0.res ".
                    "$vcf_file ".
                    "$output_dir ";

    my $cmd_get_SNP_list = "grep -v \"q_chr\" $output_dir/$out_base.blast.fmt0.res.final.4outgroup.tsv ".
                    "| awk \'{print \$1\"_\"\$2}\' > ".
                    "$output_dir/$out_base.blast.fmt0.res.final.4outgroup.tsv.snp.list";

    run_command($cmd_get_OG, "Failed to run ./get_othologous_genotype.pl");

    run_command($cmd_get_SNP_list, "Failed to run cmd_get_SNP_list");

    print <<"DONE_MSG";

--------------------------------------
3rd step -- parsing blast result DONE.
--------------------------------------

Output files:
  1) $output_dir/$out_base.blast.fmt0.res.simplified.4ck.tmp            # intermediate
  2) $output_dir/$out_base.blast.fmt0.res.genotype.4ck.tsv              # alignment-based genotype info
  3) $output_dir/$out_base.blast.fmt0.res.4outgroup.vcf                 # filtered VCF lines
  4) $output_dir/$out_base.blast.fmt0.res.final.4outgroup.tsv           # final list of query vs outgroup genotype
  5) $output_dir/$out_base.blast.fmt0.res.outgroup.fa                   # outgroup sequence in FASTA
  6) $output_dir/$out_base.blast.fmt0.res.count.tsv                     # summary of sites found/identical
  7) $output_dir/$out_base.blast.fmt0.res.final.4outgroup.tsv.snp.list  # retained SNP list for downstream VCF2FASTA
DONE_MSG


    ################################### 4) to convert retained VCF to fasta, and concatenate VCF-fasta with outgroup-fasta as the final fasta for RAxML
    
    
    #4a) run VCF2FASTA

    my $run_name="$out_base.MLtree";

    my $vcf2phylip_path="./"; ####assume vcf2phylip.py, which is used for converting VCF to FASTA file, exists in the same folder with the main script
    chomp $vcf2phylip_path;
    if (! -e "$vcf2phylip_path/vcf2phylip.py" ) {
        die "$vcf2phylip_path/vcf2phylip.py doesnt exists, please check the path\n";
    } else {
        print "checking $vcf2phylip_path/vcf2phylip.py DONE.\n";
    }

    my $cmd_VCF2FASTA = "$vcf2phylip_path/vcf2phylip.py  -f -n  -i ".
                        "$output_dir/$out_base.blast.fmt0.res.4outgroup.vcf ".
                        "--output-folder $output_dir";


    print "running $vcf2phylip_path/vcf2phylip.py  -f -n  -i  $output_dir/$out_base.blast.fmt0.res.4outgroup.vcf --output-folder $output_dir\n";

    run_command($cmd_VCF2FASTA, "failed to run $vcf2phylip_path/vcf2phylip.py  -f -n  -i  $output_dir/$out_base.blast.fmt0.res.4outgroup.vcf --output-folder $output_dir\n");


    #4b) to generate fasta file that contains outgroup sequence for RAxML

    #to generate fasta file that contains only the outgroup concatenated orthologous alleles

    open SNP_LIST, "<$output_dir/$out_base.blast.fmt0.res.final.4outgroup.tsv"
        or die "$!: cannot open $output_dir/$out_base.blast.fmt0.res.final.4outgroup.tsv\n";
    
    open OUT_FA, ">$output_dir/$out_base.tmp.fa";
    print OUT_FA ">OUTGROUP\n";
    # q_chr   q_pos   q_genotype      s_genotype      bit_score
    # ChrS01  51356   C       T       252
    # ChrS01  51389   C       C       202
    while (<SNP_LIST>){
            chomp;
            next if (/^q_chr/);
            my @a=split/\t/,$_;
            
            print OUT_FA "$a[3]"; #print orthologous allele (s_genotype) from OUTGROUP species into the FASTA
            
    }
    print OUT_FA "\n";
    close OUT_FA;
    close SNP_LIST;


    #to cat two fasta files
    my $cmd_cat_two_fa = "cat $output_dir/$out_base.tmp.fa $output_dir/$out_base.blast.fmt0.res.4outgroup.min4.fasta > $output_dir/$run_name.fa";
    run_command($cmd_cat_two_fa, "Failed to cat $output_dir/$out_base.tmp.fa $output_dir/$out_base.blast.fmt0.res.4outgroup.min4.fasta > $output_dir/$run_name.fa\n");

    print <<"DONE_MSG";

--------------------------------------
4th step -- preparing fasta (contains also OUTGROUP) for RAxML DONE.
--------------------------------------

Output files:
  1) $output_dir/$out_base.tmp.fa                                   # fasta that contains only the outgroup concatenated orthologous alleles
  2) $output_dir/$out_base.blast.fmt0.res.4outgroup.min4.fasta      # fasta that contains concatenated SNPs from VCF
  3) $output_dir/$run_name.fa                                       # merged fasta from the previous two

DONE_MSG

    ################################### 5) to run MDtest and RAxML
    
    
    #5a) run MDtest
    open MDTEST, ">$output_dir/$run_name.MDtest.sh.log"
        or die "cannot open $output_dir/$run_name.MDtest.sh.log to write\n";

    print MDTEST "modeltest-ng -i $output_dir/$run_name.fa -p $blast_N_CPU -o $output_dir/$run_name.MDtest\n";
    close MDTEST;

    my $cmd_MDtest = "modeltest-ng -i $output_dir/$run_name.fa -p $blast_N_CPU  -o $output_dir/$run_name.MDtest";
    run_command($cmd_MDtest, "Failed to modeltest-ng -i $output_dir/$run_name.fa -p $blast_N_CPU -o $output_dir/$run_name.MDtest\n");

    #5b) run RAxML

    #first to extract the infered best evolution model from the Modeltest output (from BIC)
    my $cmd_extract_model = "grep \"raxml-ng\" $output_dir/$run_name.MDtest.out | head -n 1 | sed 's/>//' > $output_dir/$run_name.MDtest.BIC_modle_4raxml";
    run_command($cmd_extract_model, "Failed to grep \"raxml-ng\" $output_dir/$run_name.MDtest.out | head -n 1 | sed 's/>//' > $output_dir/$run_name.MDtest.BIC_modle_4raxml\n");

    my $raxml_cmd;
    if (! -z "$output_dir/$run_name.MDtest.BIC_modle_4raxml") {
        open MD, "<$output_dir/$run_name.MDtest.BIC_modle_4raxml"
            or die "cannot open $output_dir/$run_name.MDtest.BIC_modle_4raxml for read\n";
        $raxml_cmd = <MD>; # reading best model
        chomp $raxml_cmd;

        $raxml_cmd =~ s/^\s+|\s+$//g;  # to remove spaces
        
        $raxml_cmd = $raxml_cmd." --all --threads $blast_N_CPU --bs-metric fbp,tbe  --seed 142857 --tree pars{10},rand{10} --prefix $output_dir/$run_name --force";
    } else {
        die "$output_dir/$run_name.MDtest.BIC_modle_4raxml is empty, please check modeltest-ng step\n";
    }

    open RAXML, ">$output_dir/$run_name.RAXML.sh.log"
        or die "cannot open $output_dir/$run_name.RAXML.sh.log to write\n";

    print RAXML "$raxml_cmd\n";
    close RAXML;

    my $cmd_RAXML = $raxml_cmd;
    run_command($cmd_RAXML, "Failed to $raxml_cmd\n");

    print <<"DONE_MSG";

--------------------------------------
5th step -- running MDtest and RAxML DONE.
--------------------------------------

Output files:
  1) $output_dir/$run_name.MDtest                                   # modeltest-ng
  2) $output_dir/$run_name.raxml.supportTBE             # RAxML tree with Fast Transfer Bootstrap Expectation
  3) $output_dir/$run_name.raxml.supportFBP             # RAxML tree with Fast Bootstrap Proportions

DONE_MSG

}


#######to remove temporary files
# END {
#     unlink "$output_dir/$out_base.tmp.fa" if -e "$output_dir/$out_base.tmp.fa";
# }


########################################################################
#######################################SUBROUITNES
########################################################################

# usage
sub usage {
    return <<"USAGE";
Usage:
  The script is a pipeline used for infer ML tree based on population 
  genomics SNPs (stored in a VCF) and a outgroup species's reference 
  genome. Basically, it extracting the flanking [150] bp sequences from 
  SNPs (in total [301] bp) and blast these sequences to the outgroup 
  genome, identifying the SNPs' orthologous alleles from the outgroup 
  species. Then, modeltest and raxml were used to generate a ML tree.

  Before running the script, please makesure the pre-requists are 
  installed and can be called directly from command line: 
  
  Perl modules:
  -Bio::DB::Fasta

  Blast:
  -makeblastdb
  -blastn

  Modeltest:
  -modeltest-ng (recommended: conda install bioconda::raxml-ng)
  -raxml-ng (recommended: conda install bioconda::raxml-ng)

  perl $0 --vcf <FILE> --ref <FILE> --out_dir <PATH> --out <FILE> [options]

  --vcf      PATH   Input VCF file
  --ref      PATH   Reference FASTA file
  --OG_ref   PATH   Outgroup Reference FASTA file
  --out_dir  PATH   Output directory
  --out      PATH   Basename for output
  --blast_C  INT    Number of CPU used in blastn and RAxML [Default: 4]
  --FLK_len  INT    Flank length on each side [Default: 150]
  --N_perc   FLOAT  Fraction of region allowed to be 'N' [Default: 0.5]
  --help             Print this help message

Example:
  perl $0 \\
    --vcf input.vcf \\
    --ref genome.fa \\
    --OG_ref OG_genome.fa \\
    --out_dir /PATH/to/output/ \\
    --out snp_regions.fa \\
    --blast_C 8 \\
    --FLK_len 150 \\
    --N_perc 0.5

Contact: wangyz.benniao\@gmail.com
USAGE
}


##################################
# to read sample info ---deprecated
#############################
sub parse_vcf_header {

    my ($VCF_hf) = @_;
    my @samples;
    while (<$VCF_hf>) {
        next if /^##/;  #skip meta
        chomp;
        if (/^#CHROM/) {
            my @cols = split /\t/;
            # next if @cols < 10;

            @samples = @cols[9..$#cols]; #samp name starts from 10th
            last;
        }
    }
    die "ERROR: No #CHROM line found in VCF header\n" unless @samples;
    return @samples;
}


########################################################################
# to parse vcf
# -parse_vcf_variants($VCF_fh, $EXCL_fh, $OUT_fh, $db, $UPDOWN_LEN, $N_cut_ratio)
# -to read SNPs, extracts 301 bp, and output if passed filtration
# -record excluded SNPs to $EXCL_fh.
########################################################################
sub parse_vcf_variants {
    my ($VCF_fh, $EXCL_fh, $OUT_fh, $db, $UPDOWN_LEN, $N_cut_ratio) = @_;

    while (<$VCF_fh>) {
        chomp;
        next if (/^#/); #skip header

        # at least 10 cols exist in a normal VCF
        my @cols = split /\t/;
        next if @cols < 10;

        my ($chrom, $pos, $id, $ref, $alt) = @cols[0..4];

        # 1) remove multi-allelic SNPs
        if ($alt =~ /,/) {
            exclude_snp($EXCL_fh, $chrom, $pos, "Multi-allelic");
            next;
        }

        # 2)to get length of chrs
        my $chr_len = $db->length($chrom);
        if (!$chr_len) {
            exclude_snp($EXCL_fh, $chrom, $pos, "Chrom_not_in_reference");
            next;
        }

        # 3) to extract 301bp seq
        
        #my $start = $pos - $UPDOWN_LEN -1;
        my $start = $pos - $UPDOWN_LEN;
        my $end   = $pos + $UPDOWN_LEN;

        # Near_chrom_edge
        if ($start < 1 || $end > $chr_len) {
            exclude_snp($EXCL_fh, $chrom, $pos, "Near_chrom_edge");
            next;
        }

        # extract seq from ref fasta
        my $region_seq = $db->seq($chrom, $start => $end);
        unless (defined $region_seq && length($region_seq) == (2*$UPDOWN_LEN + 1)) {
            my $reason = !$region_seq ? "Sequence_not_found" : "Sequence_length_incorrect";
            exclude_snp($EXCL_fh, $chrom, $pos, $reason);
            next;
        }

        # 4) how many Ns in 301 seq?
        my $num_N = ($region_seq =~ tr/Nn//);
        if ($num_N > length($region_seq) * $N_cut_ratio) {
            exclude_snp($EXCL_fh, $chrom, $pos, "Too_many_N");
            next;
        }

        # output FASTA
        #my $seq_id = "$chrom-$pos";
        my $seq_id = "$chrom"."_"."$pos";
        print $OUT_fh ">$seq_id\n$region_seq\n";
    }
}

###############################################
# IUPAC rule for VCF2FASTA ----depracted
###############################################

sub iupac_heterozygous {
    my ($r, $a) = @_;
    $r = uc($r); $a = uc($a);

    foreach my $b ($r, $a) {
        return 'N' unless $b =~ /^[ACGT]$/;
    }

    my %iupac_map = (  #IUPAC rule: https://www.bioinformatics.org/sms/iupac.html
        'A_C' => 'M', 'C_A' => 'M',
        'A_G' => 'R', 'G_A' => 'R',
        'A_T' => 'W', 'T_A' => 'W',
        'C_G' => 'S', 'G_C' => 'S',
        'C_T' => 'Y', 'T_C' => 'Y',
        'G_T' => 'K', 'T_G' => 'K',
    );
    my $key = $r . '_' . $a;
    return $iupac_map{$key} // 'N';
}


###############################################
# exclude_snp($EXCL_fh, $chrom, $pos, $reason)
###############################################

sub exclude_snp {
    my ($fh, $chrom, $pos, $reason) = @_;
    print $fh "$chrom\t$pos\t$reason\n";
}

############################################
# run_command($cmd, $error_message)
############################################
sub run_command {
    my ($cmd, $err_msg) = @_;

    print "Running: $cmd\n";
    my $exit_code = system("$cmd 2>&1");

    if ($exit_code != 0) {
        # $? 
        # using $? >> 8 to get true info if command run correct
        my $status = $? >> 8;
        die "ERROR: $err_msg\nCommand was: $cmd\nExit status: $status\n";
    }
}



__END__
perl ./IFRAMLOG_v0.1.pl  --vcf ./test_data/Query_test.vcf --ref ./test_data/Query_test_ref.fa --out_dir ./out/ --out Test.flk150 --OG_ref ./test_data/Outgroup_test_ref.fa

