# Contents
[TOC]

# Hardware requirements

## RAM

When running from Bcl data, for human genome analyses, it is recommended to let Isaac use at least 50 GB of RAM on a 40-threaded 
system. See [tweaks](#tweaks) section for ways to run Isaac on limited hardware.

## IO

As a ball-park figure, if there is Y GBs of compressed BCL data, then Isaac roughly does the following:

* Reads Y GBs of BCL files
* Reads 3 GBs of reference (for human)
* Writes 3*Y GBs of Temporary data  
* Reads 3*Y GBs of Temporary Data
* Writes Y GBs of BAM Data

As a rule of thumb, given a reasonably high end modern CPU and enough memory for a lane of BCL files plus the reference, 
then the scratch storage should be able to do over 200 MB/s to avoid IO dominating the processing time and preferable over 500 MB/s.


# Sorted Reference

Although starting with Isaac-04 the path to .fa file can be supplied with -r argument to [isaac-align](#isaac-align), some processing time can be spared by pre-processing the reference and producing sorted-reference.xml metadata file using [isaac-sort-reference](#isaac-sort-reference).

The metadata file produced by the pre-processing  contains the absolute paths to the original .fa file and isaac-align uses that file. It is important to ensure the paths are valid at the time isaac-align is being run.

As the metadata uses absolute paths to reference files, manually copying or moving the sorted reference is not recommended. 
Instead, using the [isaac-pack-reference](#isaac-pack-reference)/[isaac-unpack-reference](#isaac-unpack-reference) tool pair is advised.

# Examples

**Analyze all data from a bcl run**

    $ isaac-align -r /path/to/sorted-reference.xml -b <Run_Folder>/Data/Intensities/BaseCalls -m 40

**Analyze bcl data using unprocessed reference**

    $ isaac-align -r /path/to/genome.fa -b <Run_Folder>/Data/Intensities/BaseCalls -m 40

**Analyze subset of lanes from a bcl run**

    $ isaac-align -r /path/to/sorted-reference.xml -b <Run_Folder>/Data/Intensities/BaseCalls -m 38 --tiles s_[1234]_

**Analyze paired read data from fastq.gz files**

> **NOTE:** the fastq files need to be named or symlinked in a special way so that Isaac can recognize them.

    $ ls Fastq/
    lane1_read1.fastq.gz  lane1_read2.fastq.gz  lane2_read1.fastq.gz  lane2_read2.fastq.gz

    $ isaac-align -r /path/to/sorted-reference.xml -b Fastq -m 40 --base-calls-format fastq-gz

**Analyze single-ended data from fastq.gz files**

> **NOTE:** that the fastq files need to be named or symlinked in a special way so that Isaac can recognize them.

    $ ls Fastq/
    lane1_read1.fastq.gz  lane2_read1.fastq.gz
    $ isaac-align -r /path/to/sorted-reference.xml -b Fastq -m 40 --base-calls-format fastq-gz

**Analyze data from bam file**

    $ isaac-align -r /path/to/sorted-reference.xml -b /path/to/my.bam -m 40 --base-calls-format bam

# Output folder structure

    Aligned
    |-- Projects (output data files)
    |   |-- <project name>
    |   |   |-- <sample name>
    |   |   |   |-- sorted.bam (bam file for the sample. Contains data for the project/sample from all flowcells)
    |   |   |   `-- sorted.bam.bai
    |   |   |-- ...
    |   `-- ...
    |-- Reports (navigable statistics pages)
    |   |-- gif
    |   |   |-- <flowcell id>
    |   |   |   `-- all
    |   |   |       `-- all
    |   |   |           `-- all
    |   |   |               |-- <per-tile statistic plot images>
    |   |   `-- ...
    |   `-- html
    |       `-- index.html (root html for the analysis reports)
    `-- Stats
        |-- BuildStats.xml (chromosome-level duplicate and coverage statistics)
        |-- DemultiplexingStats.xml (information about the barcode hits)
        `-- AlignmentStats.xml (tile-level yield, pair and alignment quality statistics)

# Tweaks

## Turning off data clipping

Isaac has the following data clipping mechanisms enabled by default: [--base-quality-cutoff and --clip-overlapping](#isaac-align). Both have shown to improve the consistency of variant calling, however they are not suitable for 
scenarios targeting evaluation of the quality of sequencing, library preparation and such. Please see [isaac-align](#isaac-align) command line reference manual for details.

# Algorithms

## Alignment candidate search

In order to find alignment candidates Isaac uses hash table of all K-mers found in the reference genome.
The value for K is possible to specify at startup with --seed-length command line option.
A list of ordered genome positions produced on startup for each K-mer. When sequence alignment candidates 
are needed, the reverse and forward strands of non-overlapping sequence K-mer seeds are produced and
the corresponding lists of reference positions are retrieved from the hash table. The lists are merged
so that the seed alignment positions resulting in the same sequence end alignment position are collapsed
together with the number of supporting seeds recorded for further prioritization.

The resulting alignment position lists are processed in the order of decreasing number
of supporting seeds. For each position the sequence is checked against the corresponding reference. Up to 
--repeat-threshold of alignments resulting in the least number of mismatches are retained for alignment
quality scoring. No more than two position lists are processed. For example, if the hash table search 
resulted in positions supported by 4,3 and 2 seeds, only the candidates supported by 4 and 3 seeds are 
checked. This has been demonstrated to be sufficient to accurately score alignments on simulated data.

### Seeding

1. Non-overlapping seeds are constructed from the entire sequence starting at the first cycle. Only the 
seeds having all bases with qscores above --seed-base-quality-min are used.
2. For each seed the hash table is queried for the number of corresponding forward and reverse alignment 
positions.
3. Seeds yielding lesser number of hits are evaluated first. 
4. Evaluation stops when required minimum and maximum of candidate positions reached, maximum set number 
of seeds evaluated or all seeds evaluated. Since the position list evaluation is iterative, some seeds 
might end up being unused

## Alignment quality scoring

**Probability of a Correct Read**

This probability is the product of the probability for each base to be correct, as determined by the quality score of 
the base and the alignment of the base against the reference. If Q[i] is the Phred quality score of a base at position i, 
we have:

    pBaseError[i] = 10^(-Q[i]/10)
    pBaseCorrect[i] = 1 - pBaseError[i]
    pBaseMatch[i] = pBaseCorrect[i]
    pBaseMismatch[i] = pBaseError[i]/3
    pReadCorrect = product(i=0..readLength-1, pBase[i]), where pBase[i]is
        pBaseMatch[i] if the base at position i matches the reference, or is part of an indel
        pBaseMismatch[i] otherwise

Note: using pBaseMatch[i] for indels is an arbitrary choice (doing otherwise would require a model for indels)

**Alignment Quality of a Single Read**

The alignment quality depends on the intrinsic quality of the alignment (inferred from pReadCorrect above), but also on 
the specificity of the alignment (i.e. the probability that the read aligns somewhere else). This is inferred from two 
quantities:

    pNeighbourhood = sum of pReadCorrect for all alignments in the neighborhood (all other identified alignment positions)
    rogCorrection = 2*GenomeLength/(4^ReadLength)

The rogCorrection is the "rest-of-genome" correction that gives an indication of the probability of having a random read 
aligning to the reference. This value tends to (and should) be very small and allows differentiating between the quality 
of reads with unique alignments (in which case pNeighbourhood == 0).

**Isaac Accumulates Actual Probabilities in pNeighborhood**

These estimated values are pessimistic in the sense that they assume that extending the alignment does not introduce any 
additional mismatches. The assumption is that a seed with one mismatch will lead to an alignment descriptor with exactly 
one match on the base with the worst quality (using the definition of pReadCorrect given above). Similarly, a seed with 
two mismatches will lead to an alignment descriptor with exactly two mismatches (on the bases with the two worst qualities).

    pNormalized = rogCorrection + pNeighbourhood


finally, the alignment quality is:

    alignmentQuality == -10 * log10(pNormalized/(pNormalized + pReadCorrect))

**Alignment Quality of a Pair**

This is simply the sum of the alignment quality of each fragment when there is exactly one resolved fragment. Otherwise, 
the alignment score is corrected by the total alignment score of all the resolved fragments:

alignmentQuality = -10 * log(pBestTemplateCorrect / pTotalTemplateCorrect)

where:

    pTemplateCorrect = product(pReadCorrect for all reads)
    totalRogCorrection = rogCorrection for the total length of all reads
    pTotalTemplateCorrect = totalRogCorrection + sum(pTemplateCorrect for all resolved templates)
    pBestTemplateCorrect = max(pTemplateCorrect for all resolved templates)

## Best pair selection

All pairings of best alignment candidates are analyzed. If at least one pair is possible that satisfies
the [dominant template](#dominant-template-detection) distribution, pairs that don't are ignored.

If [--rescue-shadows](#isaac-align) is enabled, a number of extra operations is executed depending on the state of
seed based alignment candidates:
* if one of the mates results in an empty seed-based alignment candidate list, the other mate is used
  as an anchor to find the best possible alignment using [shadow rescue](#shadow-rescue) mechanism
* If one of the mates results in a non-unique alignment, the unique mate is used as an anchor to search for more
  and possibly better alternatives using [shadow rescue](#shadow-rescue) mechanism.
* If best pair produced by seed alignment does not satisfy dominant template, each end is used as an anchor with 
  [shadow rescue](#shadow-rescue) mechanism to search for alternatives that do. If an alternative pair alignment is 
  found that is better than the anomalous alignment, the alternative is accepted.
* finally, if best pair produced from seed alignments has mates that have an indication of one or both of them crossing
  a breakpoint, the [shadow rescue](#shadow-rescue) mechanism is being used to detect presence of potential structural
  variant.

## Dominant template detection
For paired data Isaac uses pairs where each read aligns on the same contig within 50Kbp to collect template length
distribution statics. When statistically significant number of pairs have been identified (10,000), they are classified
into 4 different orientation models according to the orientation of each read (FR, RF, FF or RR). The histogram of sample
counts of the most common model is then scanned and the median value is determined as the length of the template that would
be in the middle of the list when sorted by the template length. The minimal and maximal insert sizes are sizes of templates
that would correspond to 3 standard deviations if the underlying distribution was normal distribution. In other words, the
confidence interval is defined as CI=erf(3/sqrt(2)) (approximately 0.9973) and the cut-off fractions of the list for the
minimal and maximal insert sizes are at (1-CI)/2 and (1+CI)/2 respectively, which is equivalent to discarding the lower and
the higher 0.135% of the fragments that match the dominant alignment model (note that the fragments that do not match that
alignment model are not used at all in the calculation of the insert size distribution). The process is repeated until
the insertion size distribution is stable. If there aren't enough pairs to reach stability, the insert size distribution is
inferred from the available pairs.

Dominant template is then used to flag 'proper pairs' and to enhance the sensitivity by attempting [shadow rescue](#shadow-rescue).

## Shadow rescue

With paired data, when one of the reads of the pair does not find an alignment candidate position or does not have enough 
evidence for being unambiguously mapped at that position, Isaac performs an exhaustive search for possible alignments of the 
read in the vicinity of its mate alignment. The search range is automatically computed based on the 
[dominant template](#dominant-template-detection) orientation and length.

Found candidates are checked for small indels using Smith-Waterman algorithm. If the best rescued shadow indicates
the possibility of crossing a breakpoint, an additional checking for [structural variant](#checking-for-structural-variants) 
is performed.

## Checking for structural variants

Search for structural variant breakpoint is done under assumption that each of PE mates will have some presence on one
or both sides of the breakpoint. The combination of mate parts that align to the same side of breakpoint therefore should 
satisfy the dominant template. In this implementation, the breakpoint search is perfromed as an extension of 
[shadow rescue](#shadow-rescue) mechanism.

The rescued shadow candidates are partitioned into those that have a low-mismatch tail and those that don't. The
ones with reasonably clean matching end of the read sequencing cycles are considered to be semi-aligned candidates for 
crossing the breakpoint. Then, the beginning of the sequencing read cycles is used to attempt a seed-based search for 
alternative locations in the genome. First, only alternative locations on the same contig are checked. If unsuccessful, the 
genome-wide search is performed.

For each of the alternatives an attempt is made to introduce a split alignment that includes the alternative's start 
cycles and and each of the semi-aligned candidate end cycles. The breakpoint cycle is chosen to minimize the number of 
mismatches in the resulting split alignment. Top-ranking split alignments are retained for best choice and scoring.

## Gap realignment

Gap realignment is done by attempting to introduce the combinations of gaps found on other alignments overlapped by a read
being realigned. The combination that yields the lowest smith-waterman score or reduces edit distance if the smith-waterman
score does not change wins. Of all the n gaps considered for realignment, n choose k combinations are attempted where k can 
be controlled by --realigned-gaps-per-fragment parameter.

In addition to gaps found automatically, known indels can be supplied as a VCF file with --known-indels command line option.
If multiple equivalent realignments are possible, the ones that contain known indels get preference.

## Duplicates marking

Isaac duplicates marking is possible only for paired data. During the alignment, individual segments of each pair are binned.
When each bin is processed, the alignment records are sorted by the criteria that would co-locate the duplicate templates.
The best quality alignments sort to the top of each group. Depending on the command line arguments, all but top alignments
are either marked duplicate or removed from the subsequent processing.

# Bam files

Isaac produces a separate bam file for each project/sample.

## Unaligned pairs

Pairs where both reads are unaligned are stored depending on the argument of [--keep-unaligned](#isaac-align) command line option.

--keep-unaligned|Behavior
:---------------|:-----------------------------------------------------------------------------------------------------
discard         |Ensures unaligned pairs are not present in the bam file
front           |Places unaligned pairs in the beginning of the bam file before the first aligned pair of the first chromosome. The Isaac-generated bam index file is specially crafted to skip those. This approach makes it easier to locate the unaligned clusters compared to the standard implementations which require reading past the last aligned pair of the last chromosome in the genome to locate the first unaligned pair. The drawback is that the standard samtools index command is unable to process such bam files. Be sure to keep the bam index files produced by Isaac.
back            |Makes unaligned pairs appear at the end of the bam file. Although this makes it somewhat difficult to extract unaligned data, this is the option that produced bam file that is compatible with samtools index command

## Singleton/Shadow pairs

Singleton/shadow pairs refer to pairs in which aligner was unable to decide on the alignment of one of the ends (shadow). 
In this case, the shadows are assigned the position of the end that does align (singleton). The shadows are stored in the 
bam file, immediately after their singleton.

## Read names

When run from Fastq or Bam, the original read names from the input data are preserved. When run from bcl, the following format is being used for read names.

    read-name     = flowcell-id "_" flowcell-idx ":" lane-number ":" tile-number ":" cluster-id ":0"
    
    flowcell-id   = ;flowcell identifier from BaseCalls/config.xml. "unknown-flowcell" if the identifier cannot be 
                    ;determined from config.xml file
    flowcell-idx  = ;unique 0-based index of the flowcell within the analysis
    lane-number   = ;Lane number 1-8
    tile-number   = ;Unpadded tile number
    cluster-id    = ;Unpadded 0-based cluster id in the order in which the clusters appear in the bcl tile. 

## Bam flags usage
Bit  |Description                                            |Isaac notes
:----|:------------------------------------------------------|:----------------------------------------------------------
0x001|template having multiple segments in sequencing        |
0x002|each segment properly aligned according to the aligner |Pair matches dominant template orientation. Single-ended templates don't have this flag set.
0x004|segment unmapped                                       |
0x008|next segment in the template unmapped                  |
0x010|SEQ being reverse complemented                         |
0x020|SEQ of the next segment in the template being reversed |
0x040|the first segment in the template                      |Read 1 (not set for single-ended data)
0x080|the last segment in the template                       |Read 2
0x100|secondary alignment                                    |Isaac does not produce secondary alignments
0x200|not passing quality controls                           |PF flag from RTA
0x400|PCR or optical duplicate                               |If [--keep-duplicates](#isaac-align) is turned off, duplicates are excluded from the bam file. If [--mark-duplicates](#isaac-align) is turned off, duplicates are not marked in the bam file.
0x800|Supplementary alignment                                |

## Extended tags

Isaac generates following tags in the output bam files. The list of tags stored can be controlled by --bam-exclude-tags command-line argument.

Tag|Isaac meaning
---|:----------------------------------------------------------------------------------------------------------
AS |Pair alignment score
BC |Barcode string.
NM |Edit distance (mismatches and gaps) including the soft-clipped parts of the read
OC |Original CIGAR before realignment or alignment splitting
OP |Original position before realignment
RG |Isaac read groups correspond to flowcell/lane/barcode. Should not be used for anything other than debugging
SM |Single read alignment score. Rescued shadows have it set to 65535 meaning that SM was not computed.
ZX |Cluster X pixel coordinate on the tile times 100 (only available when running from BCL, excluded from output by default)
ZY |Cluster Y pixel coordinate on the tile times 100 (only available when running from BCL, excluded from output by default) 

## MAPQ

Bam MAPQ for pairs that match dominant template orientation is min(max(SM, AS), 60). For reads that are not members of a 
pair matching the dominant template orientation, the MAPQ is min(SM, 60). The MAPQ might be downgraded to 0 or set to be 
unknown (255) for alignments that don't have enough evidence to be correctly scored. This behavior depends on the 
[--dodgy-alignment-score](#isaac-align) argument.

# Multiplexed data

Isaac supports processing of multiplexed data. A sample sheet is required to describe the required
output. Below is an example of a valid sample sheet:

    FCID,Lane,SampleID,SampleRef,Barcode,Description,Control,Recipe,Operator,Project
    A805CKABXX,1,AR005,human,ACAGTG,Library testing,N,101+7,RP,Demo
    A805CKABXX,1,AR008,human,ACTTGA,Library testing,N,101+7,RP,Demo
    A805CKABXX,1,PhiX,phix,TTAGGC,Library testing,Y,101+7,RP,Demo
    A805CKABXX,1,,unknown,Undetermined,Ignored clusters with unmatched barcodes for lane 1,N,101+7,RP,Demo

A sample sheet can be used to split the data even if it does not have barcodes. Such data can be
configured to produce a separate output per lane or a combination of lanes. Leave the Barcode column
empty if the barcode cycles are not present in the data. For multi-component barcodes use '-' to
separate the components.

# Toolkit Reference

## isaac-align

**Usage**

isaac-align -r <reference> -b <base calls> -m <memory limit> [optional arguments]

**Options**

    --allow-empty-flowcells arg (=0)                Avoid failure when some of the --base-calls contain no data
    --anchor-mate arg (=1)                          Allow entire pair to be anchored by only one read if it has not 
                                                    been realigned. If not set, each read is anchored individually and 
                                                    does not affect anchoring of its mate.
    --anomalous-pair-handicap arg (=240)            When deciding between an anomalous pair and a rescued pair, this is
                                                    proportional to the number of mismatches anomalous pair needs to 
                                                    have less in order to be accepted instead of a rescued pair.
    --bam-exclude-tags arg (=ZX,ZY)                 Comma-separated list of regular tags to exclude from the output BAM
                                                    files. Allowed values are: all,none,AS,BC,NM,OC,RG,SM,ZX,ZY
    --bam-gzip-level arg (=1)                       Gzip level to use for BAM
    --bam-header-tag arg                            Additional bam entries that are copied into the header of each 
                                                    produced bam file. Use '\t' to represent tab separators.
    --bam-pessimistic-mapq arg (=0)                 When set, the MAPQ is computed as MAPQ:=min(60, min(SM, AS)), 
                                                    otherwise MAPQ:=min(60, max(SM, AS))
    --bam-produce-md5 arg (=1)                      Controls whether a separate file containing md5 checksum is 
                                                    produced for each output bam.
    --bam-pu-format arg (=%F:%L:%B)                 Template string for bam header RG tag PU field. Ordinary characters
                                                    are directly copied. The following placeholders are supported:
                                                      - %F             : Flowcell ID
                                                      - %L             : Lane number
                                                      - %B             : Barcode
    --barcode-mismatches arg (=1)                   Multiple entries allowed. Each entry is applied to the 
                                                    corresponding base-calls. Last entry applies to all the 
                                                    bases-calls-directory that do not have barcode-mismatches 
                                                    specified. Last component mismatch value applies to all subsequent 
                                                    barcode components should there be more than one. Examples: 
                                                      - 1:0             : allow one mismatch for the first barcode 
                                                    component and no mismatches for the subsequent components.
                                                      - 1               : allow one mismatch for every barcode 
                                                    component.
                                                      - 0               : no mismatches allowed in any barcode 
                                                    component. This is the default.
    -b [ --base-calls ] arg                         full path to the base calls. Multiple entries allowed. Path should 
                                                    point either to a directory or a file depending on 
                                                    --base-calls-format
    -f [ --base-calls-format ] arg                  Multiple entries allowed. Each entry is applied to the 
                                                    corresponding base-calls. Last entry is applied to all --base-calls
                                                    that don't have --base-calls-format specified.
                                                      - bam             : --base-calls points to a Bam file. All data 
                                                    found in bam file is assumed to come from lane 1 of a single 
                                                    flowcell.
                                                      - bcl             : --base-calls points to RunInfo.xml file. Data
                                                    is made of uncompressed bcl files.
                                                      - bcl-gz          : --base-calls points to RunInfo.xml file. Bcl 
                                                    cycle tile files are individually compressed and named 
                                                    s_X_YYYY.bcl.gz
                                                      - bcl-bgzf        : --base-calls points to RunInfo.xml file. Bcl 
                                                    data is stored in cycle files that are named CCCC.bcl.bgzf
                                                      - fastq           : --base-calls points to a directory containing
                                                    one fastq per lane/read named lane<X>_read<Y>.fastq. Use 
                                                    lane<X>_read1.fastq for single-ended data.
                                                      - fastq-gz        : --base-calls points to a directory containing
                                                    one compressed fastq per lane/read named lane<X>_read<Y>.fastq.gz. 
                                                    Use lane<X>_read1.fastq.gz for single-ended data.
    --base-quality-cutoff arg (=15)                 3' end quality trimming cutoff. Value above 0 causes low quality 
                                                    bases to be soft-clipped. 0 turns the trimming off.
    --bcl-tiles-per-chunk arg (=1)                  Increase this number when the tiles are too small for the 
                                                    processing to be efficient. In particular, collecting the template 
                                                    length statistics requires several tens of thousands clusters to 
                                                    work. If tiles are small and data is heavily multiplexed, there 
                                                    might be not enough clusters in a single tile to collect the tls 
                                                    for a sample
    --bin-regex arg (=all)                          Define which bins appear in the output bam files
                                                    all                   : Include all bins in the bam and all contig 
                                                    entries in the bam header.
                                                    skip-empty             : Include only the contigs that have aligned
                                                    data.
                                                    REGEX                 : Is treated as comma-separated list of 
                                                    regular expressions. Bam files will be filtered to contain only the
                                                    bins that match by the name.
    --candidate-matches-max arg (=800)              Maximum number of candidate matches to be considered for finding 
                                                    the best alignment. If seeds yield a greater number, the alignment 
                                                    generally is not performed. Other mechanisms such as shadow rescue 
                                                    may still place the fragment.
    --cleanup-intermediary arg (=0)                 When set, Isaac will erase intermediate input files for the stages 
                                                    that have been completed. Notice that this will prevent resumption 
                                                    from the stages that have their input files removed. --start-from 
                                                    Last will still work.
    --clip-overlapping arg (=1)                     When set, the pairs that have read ends overlapping each other will
                                                    have the lower-quality end soft-clipped.
    --clip-semialigned arg (=0)                     When set, reads have their bases soft-clipped on either sides until
                                                    a stretch of 5 matches is found
    -c [ --cluster ] arg                            Restrict the alignment to the specified cluster Id (multiple 
                                                    entries allowed)
    --clusters-at-a-time arg (=8000000)             Bam and fastq only. When not set, number of clusters to process 
                                                    together when input is bam or fastq is computed automatically based
                                                    on the amount of available RAM. Set to non-zero value to force 
                                                    deterministic behavior.
    --decoy-regex arg (=decoy)                      Contigs that have matching names are marked as decoys and enjoy 
                                                    reduced effort. In particular: 
                                                      - Smith waterman is not used for alignments
                                                      - Suspicious alignments are marked dodgyFor example, to mark 
                                                    everything that does not begin with chr as decoy use the following 
                                                    regex: ^(?!chr.*)
    --default-adapters arg                          Multiple entries allowed. Each entry is associated with the 
                                                    corresponding base-calls. Flowcells that don't have 
                                                    default-adapters provided, don't get adapters clipped in the data. 
                                                    Each entry is a comma-separated list of adapter sequences written 
                                                    in the direction of the reference. Wildcard (* character) is 
                                                    allowed only on one side of the sequence. Entries with * apply only
                                                    to the alignments on the matching strand. Entries without * apply 
                                                    to all strand alignments and are matched in the order of appearance
                                                    in the list.
                                                    Examples:
                                                      ACGT*,*TGCA       : Will clip ACGT and all subsequent bases in 
                                                    the forward-strand alignments and mirror the behavior for the 
                                                    reverse-strand alignments.
                                                      ACGT,TGCA         : Will find the following sequences in the 
                                                    reads: ACGT, TGCA, ACGTTGCA  (but not TGCAACGT!) regardless of the 
                                                    alignment strand. Then will attempt to clip off the side of the 
                                                    read that is shorter. If both sides are roughly equal length, will 
                                                    clip off the side that has less matches.
                                                      Standard          : Standard protocol adapters. Same as 
                                                    AGATCGGAAGAGC*,*GCTCTTCCGATCT
                                                      Nextera           : Nextera standard. Same as 
                                                    CTGTCTCTTATACACATCT*,*AGATGTGTATAAGAGACAG
                                                      NexteraMp         : Nextera mate-pair. Same as 
                                                    CTGTCTCTTATACACATCT,AGATGTGTATAAGAGACAG
    --description arg                               Free form text to be stored in the Isaac @PG DS bam header tag
    --detect-template-block-size arg (=10000)       Number of pairs to use as a single block for template length 
                                                    statistics detection
    --disable-resume arg (=0)                       If eanbled, Isaac does not persist the state of the analysis on 
                                                    disk. This might save noticeable amount of runtime at the expense 
                                                    of not being able to use --start-from option.
    --dodgy-alignment-score arg (=0)                Controls the behavior for templates where alignment score is 
                                                    impossible to assign:
                                                     - Unaligned        : marks template fragments as unaligned
                                                     - 0-254            : exact MAPQ value to be set in bam
                                                     - Unknown          : assigns value 255 for bam MAPQ. Ensures SM 
                                                    and AS are not specified in the bam
    --enable-numa [=arg(=1)] (=0)                   Replicate static data across NUMA nodes, lock threads to their NUMA
                                                    nodes, allocate thread private data on the corresponding NUMA node
    --expected-bgzf-ratio arg (=1)                  compressed = ratio * uncompressed. To avoid memory overallocation 
                                                    during the bam generation, Isaac has to assume certain compression 
                                                    ratio. If Isaac estimates less memory than is actually required, it
                                                    will fail at runtime. You can check how far you are from the 
                                                    dangerous zone by looking at the resident/swap memory numbers for 
                                                    your process during the bam generation. If you see too much showing
                                                    as 'swap', it is safe to reduce the --expected-bgzf-ratio.
    --expected-coverage arg (=60)                   Expected coverage is required for Isaac to estimate the efficient 
                                                    binning of the aligned data.
    --fastq-q0 arg (=!)                             Character to serve as base quality 0 in fastq input.
    --gap-scoring arg (=bwa)                        Gapped alignment algorithm parameters:
                                                     - eland            : equivalent of 2:-1:-15:-3:-25
                                                     - bwa              : equivalent of 0:-3:-11:-4:-20
                                                     - bwa-mem          : equivalent of 1:-4:-6:-1:-20
                                                     - m:mm:go:ge:me:gl : colon-delimited string of values where:
                                                         m              : match score
                                                         mm             : mismatch score
                                                         go             : gap open score
                                                         ge             : gap extend score
                                                         me             : min extend score (all gaps reaching this 
                                                    score will be treated as equal)
    --hash-table-buckets arg (=4294967296)          Number of buckets to use for reference hash table. Larger number of
                                                    buckets requires more RAM but it tendsto speed up the execution.
    -h [ --help ]                                   produce help message and exit
    --help-defaults                                 produce tab-delimited list of command line options and their 
                                                    default values
    --help-md                                       produce help message pre-formatted as a markdown file section and 
                                                    exit
    --ignore-missing-bcls arg (=0)                  When set, missing bcl files are treated as all clusters having N 
                                                    bases for the corresponding tile cycle. Otherwise, encountering a 
                                                    missing bcl file causes the analysis to fail.
    --ignore-missing-filters arg (=0)               When set, missing filter files are treated as if all clusters pass 
                                                    filter for the corresponding tile. Otherwise, encountering a 
                                                    missing filter file causes the analysis to fail.
    --input-concurrent-load arg (=64)               Maximum number of concurrent file read operations for --base-calls
    -j [ --jobs ] arg (=40)                         Maximum number of compute threads to run in parallel
    --keep-duplicates arg (=1)                      Keep duplicate pairs in the bam file (with 0x400 flag set in all 
                                                    but the best one)
    --keep-unaligned arg (=back)                    Available options:
                                                     - discard          : discard clusters where both reads are not 
                                                    aligned
                                                     - front            : keep unaligned clusters in the front of the 
                                                    BAM file
                                                     - back             : keep unaligned clusters in the back of the 
                                                    BAM file
    --known-indels arg                              path to a VCF file containing known indels fore realignment.
    --lane-number-max arg (=8)                      Maximum lane number to look for in --base-calls-directory (fastq 
                                                    only).
    --mapq-threshold arg (=-1)                      If any fragment alignment in template is below the threshold, 
                                                    template is not stored in the BAM.
    --mark-duplicates arg (=1)                      If not set and --keep-duplicates is set, the duplicates are not 
                                                    discarded and not flagged.
    --match-finder-shadow-split-repeats arg (=100000)
                                                    Maximum number of seed candidate matches to be considered for 
                                                    finding a possible alignment split.
    --match-finder-too-many-repeats arg (=4000)     Maximum number of seed matches to be looked at for each attempted 
                                                    seed
    --match-finder-way-too-many-repeats arg (=100000)
                                                    Maximum number of seed matches to be looked at if in a pair one 
                                                    read has candidate alignments and the otherhas gone over 
                                                    match-finder-too-many-repeats on all seeds or over 
                                                    candidate-matches-max when seed position merge was attempted 
    --memory-control arg (=off)                     Define the behavior in case unexpected memory allocations are 
                                                    detected: 
                                                      - warning         : Log WARNING about the allocation.
                                                      - off             : Don't monitor dynamic memory usage.
                                                      - strict          : Fail memory allocation. Intended for 
                                                    development use.
    -m [ --memory-limit ] arg (=0)                  Limits major memory consumption operations to a set number of 
                                                    gigabytes. 0 means no limit, however 0 is not allowed as in such 
                                                    case Isaac will most likely consume all the memory on the system 
                                                    and cause it to crash. Default value is taken from ulimit -v.
    --neighborhood-size-threshold arg (=0)          Threshold used to decide if the number of reference 32-mers sharing
                                                    the same prefix (16 bases) is small enough to justify the 
                                                    neighborhood search. Use large enough value e.g. 10000 to enable 
                                                    alignment to positions where seeds don't match exactly.
    --output-concurrent-save arg (=120)             Maximum number of concurrent file write operations for 
                                                    --output-directory
    -o [ --output-directory ] arg (=./Aligned)      Directory where the final alignment data be stored
    --per-tile-tls arg (=0)                         Forces template length statistics(TLS) to be recomputed for each 
                                                    tile. When not set, the first tile that produces stable TLS will 
                                                    determine TLS for the rest of the tiles of the lane. Notice that as
                                                    the tiles are not guaranteed to be processed in the same order 
                                                    between different runs, some pair alignments might vary between two
                                                    runs on the same data unless --per-tile-tls is set. It is not 
                                                    recommended to set --per-tile-tls when input data is not randomly 
                                                    distributed (such as bam) as in such cases, the shadow rescue range
                                                    will be biased by the input data ordering.
    --pf-only arg (=1)                              When set, only the fragments passing filter (PF) are generated in 
                                                    the BAM file
    --pre-allocate-bins arg (=0)                    Use fallocate to reduce the bin file fragmentation. Since bin files
                                                    are pre-allocated based on the estimation of their size, it is 
                                                    recommended to turn bin pre-allocation off when using RAM disk as 
                                                    temporary storage.
    --pre-sort-bins arg (=1)                        Unset this value if you are working with references that have many 
                                                    contigs (1000+)
    --read-name-length arg (=0)                     Maximum read name length (fastq and bam only). Value of 0 causes 
                                                    the read name length to be determined by reading the first records 
                                                    of the input data. Shorter than needed read names can cause 
                                                    duplicate names in the output bam files.
    --realign-dodgy arg (=0)                        If not set, the reads without alignment score are not realigned 
                                                    against gaps found in other reads.
    --realign-gaps arg (=sample)                    For reads overlapping the gaps occurring on other reads, check if 
                                                    applying those gaps reduces mismatch count. Significantly reduces 
                                                    number of false SNPs reported around short indels.
                                                      - no              : no gap realignment
                                                      - sample          : realign against gaps found in the same sample
                                                      - project         : realign against gaps found in all samples of 
                                                    the same project
                                                      - all             : realign against gaps found in all samples
    --realign-mapq-min arg (=60)                    Gaps from alignments with lower MAPQ will not be used as candidates
                                                    for gap realignment
    --realign-vigorously arg (=0)                   If set, the realignment result will be used to search for more gaps
                                                    and attempt another realignment, effectively extending the 
                                                    realignment over multiple deletions not covered by the original 
                                                    alignment.
    --realigned-gaps-per-fragment arg (=4)          Maximum number of gaps the realigner can introduce into a fragment.
                                                    For 100 bases long DNA it is reasonable to keep it no bigger than 
                                                    2. RNA reads can overlap multiple introns. Therefore a larger 
                                                    number is probably required for RNA. Notice that bigger values can 
                                                    significantly slow down the bam generation as there is a n choose k
                                                    combination try with n being the number of gaps detected by all 
                                                    other fragment alignments that overlap the fragment being 
                                                    realigned.
    -r [ --reference-genome ] arg                   Full path to the reference genome XML descriptor or .fa file.
    -n [ --reference-name ] arg (=default)          Unique symbolic name of the reference. Multiple entries allowed. 
                                                    Each entry is associated with the corresponding --reference-genome 
                                                    and will be matched against the 'reference' column in the sample 
                                                    sheet. 
                                                    Special names:
                                                      - unknown         : default reference to use with data that did 
                                                    not match any barcode.
                                                      - default         : reference to use for the data with no 
                                                    matching value in sample sheet 'reference' column.
    --remap-qscores arg                             Replace the base calls qscores according to the rules provided.
                                                     - identity   : No remapping. Original qscores are preserved
                                                     - bin:8      : Equivalent of 0-1:0,2-9:7,10-19:11,20-24:22,25-29:2
                                                    7,30-34:32,35-39:37,40-63:40
    --repeat-threshold arg (=100)                   Threshold used to decide if matches must be discarded as too 
                                                    abundant (when the number of repeats is greater or equal to the 
                                                    threshold)
    --rescue-shadows arg (=1)                       Scan within dominant template range off an orphan, for a possible 
                                                    shadow alignment
    --response-file arg                             file with more command line arguments
    -s [ --sample-sheet ] arg                       Multiple entries allowed. Each entry is applied to the 
                                                    corresponding base-calls.
                                                      - none            : process flowcell as if there is no sample 
                                                    sheet
                                                      - default         : use <base-calls>/SampleSheet.csv if it 
                                                    exists. This is the default behavior.
                                                      - <file path>     : use <file path> as sample sheet for the 
                                                    flowcell.
    --scatter-repeats arg (=1)                      When set, extra care will be taken to scatter pairs aligning to 
                                                    repeats across the repeat locations 
    --seed-base-quality-min arg (=3)                Minimum base quality for the seed to be used in alignment candidate
                                                    search.
    --seed-length arg (=16)                         Length of the seed in bases. Only 10 11 12 13 14 15 16 17 18 19 20 
                                                    are allowed. Longer seeds reduce sensitivity on noisy data but 
                                                    improve repeat resolution and run time.
    --shadow-scan-range arg (=-1)                   -1     - scan for possible mate alignments between template min and
                                                    max
                                                    >=0    - scan for possible mate alignments in range of template 
                                                    median += shadow-scan-range
    --single-library-samples arg (=1)               If set, the duplicate detection will occur across all read pairs in
                                                    the sample. If not set, different lanes are assumed to originate 
                                                    from different libraries and duplicate detection is not performed 
                                                    across lanes.
    --smith-waterman-gap-size-max arg (=16)         Maximum length of gap detectable by smith waterman algorithm.
    --smith-waterman-gaps-max arg (=4)              Maximum number of gaps that can be introduced into an alignment by 
                                                    Smith-Waterman algorithm. If the optimum alignment has more gaps, 
                                                    it is simply ignored as an alignment candidate.
    --split-alignments arg (=1)                     When set, alignments crossing a structural variant are allowed to 
                                                    be split with SA tag.
    --split-gap-length arg (=10000)                 Maximum length of insertion or deletion allowed to exist in a read.
                                                    If a gap exceeds this limit, the read gets broken up around the gap
                                                    with SA tag introduced
    --start-from arg (=Start)                       Start processing at the specified stage:
                                                      - Start            : don't resume, start from beginning
                                                      - Align            : same as Start
                                                      - AlignmentReports : regenerate alignment reports and bam
                                                      - Bam              : resume at bam generation
                                                      - Finish           : Same as Bam.
                                                      - Last             : resume from the last successful step
                                                    Note that although Isaac attempts to perform some basic validation,
                                                    the only safe option is 'Start' The primary purpose of the feature 
                                                    is to reduce the time required to diagnose the issues rather than 
                                                    be used on a regular basis.
    --stats-image-format arg (=none)                Format to use for images during stats generation
                                                     - gif        : produce .gif type plots
                                                     - none       : no stat generation
    --stop-at arg (=Finish)                         Stop processing after the specified stage is complete:
                                                      - Start            : perform the first stage only
                                                      - Align            : same as Start
                                                      - AlignmentReports : don't perform bam generation
                                                      - Bam              : finish when bam is done
                                                      - Finish           : stop at the end.
                                                      - Last             : perform up to the last successful step only
                                                    Note that although Isaac attempts to perform some basic validation,
                                                    the only safe option is 'Finish' The primary purpose of the feature
                                                    is to reduce the time required to diagnose the issues rather than 
                                                    be used on a regular basis.
    --target-bin-size arg (=0)                      Isaac will attempt to bin temporary data so that each bin is close 
                                                    to targetBinSize in megabytes (1024 * 1024 bytes). Value of 0 will 
                                                    cause Isaac to compute the target bin size automatically based on 
                                                    the available memory.
    --temp-concurrent-load arg (=4)                 Maximum number of concurrent file read operations for 
                                                    --temp-directory
    --temp-concurrent-save arg (=680)               Maximum number of concurrent file write operations for 
                                                    --temp-directory
    -t [ --temp-directory ] arg (=./Temp)           Directory where the temporary files will be stored (matches, 
                                                    unsorted alignments, etc.)
    --tiles arg                                     Comma-separated list of regular expressions to select only a subset
                                                    of the tiles available in the flow-cell.
                                                    - to select all the tiles ending with '5' in all lanes: --tiles 
                                                    [0-9][0-9][0-9]5
                                                    - to select tile 2 in lane 1 and all the tiles in the other lanes: 
                                                    --tiles s_1_0002,s_[2-8]
                                                    Multiple entries allowed, each applies to the corresponding 
                                                    base-calls.
    --tls arg                                       Template-length statistics in the format 
                                                    'min:median:max:lowStdDev:highStdDev:M0:M1', where M0 and M1 are 
                                                    the numeric value of the models (0=FFp, 1=FRp, 2=RFp, 3=RRp, 4=FFm,
                                                    5=FRm, 6=RFm, 7=RRm)
    --trim-pe arg (=1)                              Trim overhanging ends of PE alignments
    --use-bases-mask arg                            Conversion mask characters:
                                                      - Y or y          : use
                                                      - N or n          : discard
                                                      - I or i          : use for indexing
                                                    
                                                    If not given, the mask will be guessed from the config.xml file in 
                                                    the base-calls directory.
                                                    
                                                    For instance, in a 2x76 indexed paired end run, the mask 
                                                    I<Y76,I6n,y75n> means:
                                                      use all 76 bases from the first end, discard the last base of the
                                                    indexing read, and use only the first 75 bases of the second end.
    --use-smith-waterman arg (=smart)               One of the following:
                                                     - always           : Use smith-waterman to reduce the amount of 
                                                    mismatches in aligned reads
                                                     - smart            : apply heuristics to avoid executing costly 
                                                    smith-waterman on sequences that are unlikely to produce gaps
                                                     - never            : Don't use smith-waterman
    --variable-read-length arg                      Unless set, Isaac will fail if the length of the sequence changes 
                                                    between the records of a fastq or a bam file.
    --verbosity arg (=2)                            Verbosity: FATAL(0), ERRORS(1), WARNINGS(2), INFO(3), DEBUG(4) (not
                                                    supported yet)
    -v [ --version ]                                print program version information


## isaac-merge-references

**Usage**

isaac-merge-references [options]

**Options**

    -h [ --help ]                                         Print this message
    -v [ --version ]                                      Only print version information

    -i [ --input-file ] arg                               Path to sorted-reference.xml to be merged. 
                                                          Multiple entries allowed.
    -o [ --output-directory ] arg (./IsaacIndex.20161121) Location where the results are stored


## isaac-pack-reference

**Usage**

isaac-pack-reference [options]

**Options**

    -h [ --help ]                                         Print this message
    -n [ --dry-run ]                                      Don't actually run any commands; just print them
    -v [ --version ]                                      Only print version information
    -j [ --jobs ] arg (=40)                               Maximum number of parallel operations

    -r [ --reference-genome ] arg                         Path to sorted-reference.xml 
    -o [ --output-file ] arg (./packed-reference.tar.gz)  Archive path


## isaac-reorder-reference

**Usage**

isaac-reorder-reference

**Options**

    -h [ --help ]                 produce help message and exit
    --help-defaults               produce tab-delimited list of command line options and their default values
    --help-md                     produce help message pre-formatted as a markdown file section and exit
    --order arg                   Comma-separated list of contig names in the order in which they will appear in the 
                                  new .fa file.
    -d [ --output-directory ] arg Path for the reordered fasta and annotation files.
    -x [ --output-xml ] arg       Path for the new xml file.
    -r [ --reference-genome ] arg Full path to the reference genome XML descriptor.
    -v [ --version ]              print program version information


## isaac-sort-reference

**Usage**

isaac-sort-reference [options]

**Options**

    -g [ --genome-file ] arg                              Path to fasta file containing the reference contigs 
    -h [ --help ]                                         Print this message
    -n [ --dry-run ]                                      Don't actually run any commands; just print them
    -o [ --output-directory ] arg (./IsaacIndex.20161121) Location where the results are stored
    -q [ --quiet ]                                        Avoid excessive logging
    -v [ --version ]                                      Only print version information
    --target arg (all)                                    Individual target to make

## isaac-unpack-reference

**Usage**

isaac-unpack-reference [options]

**Options**

    -h [ --help ]                                        Print this message
    -i [ --input-file ] arg                              Archive path
    --make-movable                                       Store relative paths in sorted-reference.xml so that the entire
                                                         folder can be copied elsewhere
    -n [ --dry-run ]                                     Don't actually run any commands; just print them
    -v [ --version ]                                     Only print version information
