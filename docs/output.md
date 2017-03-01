# Output #

Neptune's output directory contains the following items:

* **candidates**: directory containing signature candidates
* **filtered**: directory containing filtered candidates in extracted order
* **sorted**: directory containing filtered signatures in sorted order
* **consolidated**: directory containing the consolidate signatures
* **database**: directory containing Neptune's constructed databases
* **aggregate.kmers**: file containing all observed \textit{k}-mers
* **receipt.txt**: file containing Neptune's run receipt

A file with the same name as each reference will be placed in each output directory, corresponding to the reference file from which it was derived.

## Candidate Signatures ##

The candidate signatures are the sequences produced from the signature extraction step. These signatures will relatively sensitive, but not necessarily specific. This is because signature extraction is done using exact *k*-mer matches. The candidate signatures are guaranteed to contain no more exact matches with any exclusion *k*-mer than specified by the --exhits parameter. However, there may be inexact matches with exclusion targets.

## Filtered Signatures ##

The filtering step is designed to remove signatures which are not interesting enough to warrant further investigation, because the negative component of their score is prohibitively large. The filtering step removes signatures that align sufficiently with any exclusion target. The filtered signatures are a subset of the candidate signatures.

## Sorted Signatures ##

The sorted signatures files are organized as FASTA records containing the same signatures as their filtered signatures counterparts. However, the signatures are listed in descending order by their signature score. Signatures are assigned a score corresponding to their highest-scoring BLAST alignments with all inclusion and exclusion targets, which is a sum of a positive inclusion component and a negative exclusion component. This score is maximized when all inclusion targets contain a region exactly matching the entire signature and there exists no exclusion targets that match the signature. The signatures have the following format:

    >[ID] [SCORE] [IN SCORE] [EX SCORE] [LENGTH] [REF] [POS]
    [SEQUENCE]

The following is an example:

    >425 score=0.86 in=0.98 ex=-0.13 len=31 ref=ecoli pos=160
    TGTCATTCTCCTGTTCTGCCTGTATCACTGC

Where:

* **[ID]**: an __arbitrary__, run-unique ID assigned to the signature
* **[SCORE]**: the total signature score
* **[IN SCORE]**: positive inclusion component of signature score
* **[EX SCORE]**: negative exclusion component of signature score
* **[LENGTH]**: signature length in bases
* **[REF]**: name of the contig from which the signature was extracted
* **[POS]**: starting position of the signature in the reference
* **[SEQUENCE]**: sequence content of the signature

## Consolidated Signatures ##

The sorted signatures from all references are combined into a single "consolidated.fasta" file, located within the "consolidated" directory. Signatures are added to the consolidated signatures file in a greedy manner by selecting the next highest scoring signature available from all references. While effort is taken to prevent signatures from overlapping entirely, it is possible for consolidate signatures to have a small amount of overlap. In many circumstances, this output might be considered the final output of Neptune.

## Databases ##

The databases directory contains BLAST databases constructed from the inclusion and exclusion files.

## Aggregate k-mers ##

The aggregated *k*-mers file, aggregated.kmers, contains a list of all *k*-mers observed in the inclusion and exclusion groups. These *k*-mers are sorted and followed by two integers: the number of inclusion and exclusion targets the *k*-mer appears in, respectively.

## Run Receipt ##

The run receipt contains information about the Neptune execution. It contains a list of all the files in the inclusion and exclusion group, and the command line parameters used for the execution.
