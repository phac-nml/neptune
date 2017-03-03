# Signature Format #

The signatures produced by Neptune are output in FASTA format with additional information in the description line. Signatures are output in the following format:

    >[ID] [SCORE] [IN SCORE] [EX SCORE] [LENGTH] [REF] [POS]
    [SEQUENCE]

The following is an example:

    >425 score=0.86 in=0.98 ex=-0.13 len=31 ref=ecoli pos=160
    TGTCATTCTCCTGTTCTGCCTGTATCACTGC

Where:

| Item | Full Name | Description |
|---|---|---|
| [ID] | ID | An __arbitrary__, run-unique ID assigned to the signature. |
| [SCORE] | Score | The total signature score. This is the sum of the inclusion (sensitivity) and exclusion (specificity) scores. |
| [IN SCORE] | Inclusion Score | The positive inclusion component of signature score (sensitivity). |
| [EX SCORE] | Exclusion Score | The negative exclusion component of signature score (specificity). |
| [LENGTH] | Length | The signature length in bases. |
| [REF] | Reference | The unique identifier of the contig from which the signature was extracted. |
| [POS] | Position | The starting position of the signature in the reference. |
| [SEQUENCE] | Sequence | The sequence content of the signature. |

## ID ##

The signature ID is an __arbitrary__, run-unique ID assigned to the signature. The signatures within the same FASTA file will have unique IDs, relative to each other. However, signatures within multiple output files will have overlapping signature IDs. This will be the case when using multiple references or not specifying any reference files. The signatures within the `consolidated.fasta` output will have unique signature IDs.

## Total Score ##

Signatures are assigned a score corresponding to their highest-scoring BLAST alignments with all inclusion and exclusion targets, which is a sum of the positive inclusion score (sensitivity) and the negative exclusion component (specificity). This score is maximized when all inclusion targets contain a region exactly matching the entire signature and there exists no exclusion targets that match the signature.

## Inclusion Score ##

The inclusion score is a non-negative number between 0.00 and 1.00 and relates to the signature's sensitivity. This score is determined by the signature's highest-scoring BLAST alignments with all inclusion targets. The inclusion score is maximized (good) when the signature is found exactly and completely in all inclusion targets and minimized (bad) when the signature is not found whatsoever in any inclusion targets.

## Exclusion Score ##

The exclusion score is a non-positive number between -1.00 and 0.00 and relates to the signature's specificity. This score is determined by the signature's highest-scoring BLAST alignments with all exclusion targets. The exclusion score is maximized (bad) when the signature is found exactly and completely in all exclusion targets and minimized (good) when the signature is not found whatsoever in any exclusion targets.

## Length ##

The length describes the length of the signature in bases. Although this can be calculated from the sequence, it is included in the FASTA description to accommodate other tools.

## Reference ##

The reference describes the sequence identifier of the contig the signature was extracted from. This is useful for determining where the signature lies and what sequence surrounds it.

## Position ##

The position describes the base position of the signature within the contig reference it was extracted from. This is useful for determining where the signature lies and what sequence surrounds it.

## Sequence ##

The sequence describes the sequence content of the signature and follows the specifications of FASTA format. However, the sequence will not contain line breaks, regardless of the sequence length.

