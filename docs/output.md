# Output

Neptune's output directory contains the following items:

| Item | Type | Description |
|---|---|---|
| candidates | directory | The directory containing signature candidates in extracted order. |
| filtered | directory | The directory containing filtered signature candidates in extracted order. |
| sorted | directory | The directory containing filtered signatures in signature-score sorted order. |
| consolidated | directory | The directory containing the consolidate signatures from multiple sorted-signature reference files. |
| database | directory | The directory containing Neptune's BLAST constructed databases. |
| aggregate.kmers | file | The *k*-mer file containing all observed k-mers. |
| receipt.txt | file | The file containing Neptune's run receipt. |

A file with the same name as each reference will be placed in each output directory (candidates, filtered, sorted), corresponding to the reference file from which it was derived.

## Candidates

The candidate signatures are the sequences produced from the signature extraction step. These signatures will relatively sensitive, but not necessarily specific. This is because signature extraction is done using exact *k*-mer matches. The candidate signatures are guaranteed to contain no more exact matches with any exclusion *k*-mer than specified by the `--exhits` parameter. However, there may be inexact matches with exclusion targets.

## Filtered

The filtering step is designed to remove signatures which are not interesting enough to warrant further investigation, because the negative component of their score is prohibitively large. The filtering step removes signatures that align sufficiently with any exclusion target. The filtered signatures are a subset of the candidate signatures.

## Sorted

The sorted signatures files are organized as FASTA records containing the same signatures as their filtered signatures counterparts. However, the signatures are listed in descending order by their signature score. Signatures are assigned a score corresponding to their highest-scoring BLAST alignments with all inclusion and exclusion targets, which is a sum of a positive inclusion component and a negative exclusion component. This score is maximized when all inclusion targets contain a region exactly matching the entire signature and there exists no exclusion targets that match the signature.

## Consolidated

The sorted signatures from all references are combined into a single "consolidated.fasta" file, located within the "consolidated" directory. Signatures are added to the consolidated signatures file in a greedy manner by selecting the next highest scoring signature available from all references. While effort is taken to prevent signatures from overlapping entirely, it is possible for consolidate signatures to have a small amount of overlap. In many circumstances, this output might be considered the final output of Neptune.

## Databases

The databases directory contains BLAST databases constructed from the inclusion and exclusion files.

## Aggregate k-mers

The aggregated *k*-mers file, aggregated.kmers, contains a list of all *k*-mers observed in the inclusion and exclusion groups. These *k*-mers are sorted and followed by two integers: the number of inclusion and exclusion targets the *k*-mer appears in, respectively.

## Run Receipt

The run receipt contains information about the Neptune execution. It contains a list of all the files in the inclusion and exclusion group, and the command line parameters used for the execution.
