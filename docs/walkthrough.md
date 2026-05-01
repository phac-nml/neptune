# Walkthrough

## Overview

The purpose of this walkthrough will be to illustrate a simple, but complete example of using Neptune to locate discriminatory sequences. We will identity signature sequences within an artificial data set containing three inclusion sequences and three exclusion sequences. The output will be a list of signatures, sorted by score, for each inclusion input, and one consolidated signatures file, sorted by signature score, containing signatures from all inclusion inputs.

## Input Data

We will be using very small, artificial genomes for this walkthrough. However, these small genomes will be sufficient to illustrate the operation of Neptune. The artificial genome sequence content is derived from *Escherichia coli* and has been modified to introduce simple variation between genomes.

The example inclusion genomes are located in the following location:

```bash
neptune/tests/data/example/inclusion/
```

The example exclusion genomes are located in the following location:

```bash
neptune/tests/data/example/exclusion/
```

The inclusion and exclusion directories each contain three FASTA-format genomes. The genomes all have some insertions and deletions that differentiate them from each other. However, the three inclusion genomes primarily differ from the three exclusion genomes in that they share large sequences that are absent from all exclusion genomes.

## Running Neptune

Neptune will automatically calculate many of the parameters that might otherwise be specified by the user, such as the minimum number of inclusion inputs the signature sequence must be present within for it to be considered shared sequence. At minimum, Neptune requires the user specify the inclusion sequences, exclusion sequences, and an output directory. We will provide Neptune inclusion and exclusion sequences in the form of FASTA file genomes located within directories. The following command will run Neptune on the example data and output to the specified directory:

```bash
neptune
    --inclusion tests/data/example/inclusion/
    --exclusion tests/data/example/exclusion/
    --output output/
```

## Output

### Standard Output

After running Neptune, very similar output will be printed to standard output, indicating that Neptune is starting and completing different stages of operation:

```text
Neptune v2.0.0

Estimating k-mer size ...
k = 15

k-mer Counting...
Submitted 6 jobs.
0.04893639404326677 seconds

k-mer Aggregation...
Submitted 65 jobs.
0.0374402878805995 seconds

Signature Extraction...
Submitted 3 jobs.
0.024856548057869077 seconds

Signature Filtering...
Submitted 2 jobs.
Submitted 3 jobs.
1.965036110021174 seconds

Consolidate Signatures...
Submitted 1 jobs.
0.21438432089053094 seconds

Complete!
```

### Consolidated Signatures

As we did not specify references from which to extract signatures, Neptune will automatically investigate all inclusion genomes for signatures and consolidate those signatures into a single consolidated signature file. The `output/consolidated/consolidated.fasta` file contains these consolidated signatures. This file may be understood as the final output of the application. The following FASTA-format output is from the consolidated signatures file produced from this example:

```text
>0.0 score=1.0000 in=1.0000 ex=0.0000 len=103 ref=inclusion3 pos=99
TAGTCTCCAGGATTCCCGGGGCGGTTCAGATAATCTTAGCATTGACCGCCTTTATATAGAAGCTGTTATTCAAGAAGCATTTTCAAGCAGTGATGTAAGAAAA
>1.1 score=0.9979 in=0.9979 ex=0.0000 len=640 ref=inclusion2 pos=3494
CGCGGGCGATATTTTCACAGCCATTTTCAGGAGTTCAGCCATGAACGCTTATTACATTCAGGATCGTCTTGAGGCTCAGAGCTGGGAGCGTCACTACCAGCAGATCGCCCGTGAAGAGAAAGAGGCAGAACTGGCAGACACATGGAAAAAGGCCTGCCCCAGCACCTGTTTTGAATCGCTATGCATCGATCATTTGCAACGCCACGGGGCCAGCAAAAAAGCCATTACCCGTGCGTTTGATGACGATGTTGAGTTTCAGGAGCGCATGGCAGAACACATCCGGTACTGGTTAAACCATTGCTCACCACCAGGTTGATATTGATTCAGAGGTATAAAACGAATGAGTACAGCACTCGCAACGCTGGCAGGGAAGCTGGCTGAACGTGTCGGCATGGATTCTGTCGACCCACAGGAACTGATCACCACTCTTCGCCAGACGGCATTTAAAGGTGATGCCAGCGATGCGCAGTTCATCGCATTGCTGATCGTCGCCAACCAGTACGGTCTTAATCCGTGGACGAAAGAAATTTACGCCTTTCCTGATAAGCAGAACGGCATCGTTCCGGTGGTGGGCGTTGATGGCTGTCCCGTATCATCAATGAAAACCAGCAGTTTGAGGCATGGTACTTTGAGCAGGACA
>2.2 score=0.9966 in=0.9966 ex=0.0000 len=98 ref=inclusion1 pos=5209
GCGAGTTTTGCGAGATGGTGCCGGAGTTCATCGAAAAAATGGACGAGGCACTGCTGAAATTGGTTTTGTATTTGGGGAGCAATGGCGATGAAGCATCC
```

The FASTA header contains information relavent to the identified signature. A detailed explanation of this information is located within the [output section](../format) of this documentation. The `score` is the sum of the `in` (inclusion/sensitivity) and `ex` (exclusion/specificity) scores, and represents a combined measure of sensitivity and specificity. The `length` describes the length of the signature in bases. The `ref` (reference) and `pos` (position) describe the location of the signature within the reference FASTA record it was extracted from.

In this example, Neptune identified three signatures: `0.0`, `1.1`, and `2.2` of lengths 103, 640, and 98, respectively. These IDs (`0.0`, `1.1`, `2.2`) are arbitrarily assigned by Neptune and are only used to uniquely identify each signature. They have no other meaning.

We can see that the `0.0` signature is derived from the `inclusion3` reference, the `1.1` from the `inclusion2` reference, and the `2.2` from the `inclusion1` reference. These signatures were located at positions 99, 3497, and 5209 within the inclusion1 reference. These signatures are of very high quality, within the context of our data set, with scores of 1.0000, 0.9979, and 0.9969, within the possible range of score values from -1.00 to +1.00.

### Sorted Signatures

If we're interested in looking at the signatures produced from each individual inclusion input, we need to investigate the output in the `output/sorted` directory. The following are the signatures extracted exclusively from the `inclusion1.fasta` input:

```text
>0 score=1.0000 in=1.0000 ex=0.0000 len=103 ref=inclusion1 pos=99
TAGTCTCCAGGATTCCCGGGGCGGTTCAGATAATCTTAGCATTGACCGCCTTTATATAGAAGCTGTTATTCAAGAAGCATTTTCAAGCAGTGATGTAAGAAAA
>1 score=0.9979 in=0.9979 ex=0.0000 len=640 ref=inclusion1 pos=3497
CGCGGGCGATATTTTCACAGCCATTTCAGGAGTTCAGCCATGAACGCTTATTACATTCAGGATCGTCTTGAGGCTCAGAGCTGGGAGCGTCACTACCAGCAGATCGCCCGTGAAGAGAAAGAGGCAGAACTGGCAGACACATGGAAAAAGGCCTGCCCCAGCACCTGTTTGAATCGCTATGCATCGATCATTTGCAACGCCACGGGGCCAGCAAAAAAGCCATTACCCGTGCGTTTGATGACGATGTTGAGTTTCAGGAGCGCATGGCAGAACACATCCGGTACTGGTTGAAACCATTGCTCACCACCAGGTTGATATTGATTCAGAGGTATAAAACGAATGAGTACAGCACTCGCAACGCTGGCAGGGAAGCTGGCTGAACGTGTCGGCATGGATTCTGTCGACCCACAGGAACTGATCACCACTCTTCGCCAGACGGCATTTAAAGGTGATGCCAGCGATGCGCAGTTCATCGCATTGCTGATCGTCGCCAACCAGTACGGTCTTAATCCGTGGACGAAAGAAATTTACGCCTTTCCTGATAAGCAGAACGGCATCGTTCCGGTGGTGGGCGTTGATAGGCTGTCCCGTATCATCAATGAAAACCAGCAGTTTGAGGCATGGTACTTTGAGCAGGACA
>2 score=0.9966 in=0.9966 ex=0.0000 len=98 ref=inclusion1 pos=5209
GCGAGTTTTGCGAGATGGTGCCGGAGTTCATCGAAAAAATGGACGAGGCACTGCTGAAATTGGTTTTGTATTTGGGGAGCAATGGCGATGAAGCATCC
```

The following are the signatures extracted exclusively from the `inclusion2.fasta` input:

```text
>0 score=1.0000 in=1.0000 ex=0.0000 len=103 ref=inclusion2 pos=99
TAGTCTCCAGGATTCCCGGGGCGGTTCAGATAATCTTAGCATTGACCGCCTTTATATAGAAGCTGTTATTCAAGAAGCATTTTCAAGCAGTGATGTAAGAAAA
>1 score=0.9979 in=0.9979 ex=0.0000 len=640 ref=inclusion2 pos=3494
CGCGGGCGATATTTTCACAGCCATTTTCAGGAGTTCAGCCATGAACGCTTATTACATTCAGGATCGTCTTGAGGCTCAGAGCTGGGAGCGTCACTACCAGCAGATCGCCCGTGAAGAGAAAGAGGCAGAACTGGCAGACACATGGAAAAAGGCCTGCCCCAGCACCTGTTTTGAATCGCTATGCATCGATCATTTGCAACGCCACGGGGCCAGCAAAAAAGCCATTACCCGTGCGTTTGATGACGATGTTGAGTTTCAGGAGCGCATGGCAGAACACATCCGGTACTGGTTAAACCATTGCTCACCACCAGGTTGATATTGATTCAGAGGTATAAAACGAATGAGTACAGCACTCGCAACGCTGGCAGGGAAGCTGGCTGAACGTGTCGGCATGGATTCTGTCGACCCACAGGAACTGATCACCACTCTTCGCCAGACGGCATTTAAAGGTGATGCCAGCGATGCGCAGTTCATCGCATTGCTGATCGTCGCCAACCAGTACGGTCTTAATCCGTGGACGAAAGAAATTTACGCCTTTCCTGATAAGCAGAACGGCATCGTTCCGGTGGTGGGCGTTGATGGCTGTCCCGTATCATCAATGAAAACCAGCAGTTTGAGGCATGGTACTTTGAGCAGGACA
>2 score=0.9933 in=0.9933 ex=0.0000 len=99 ref=inclusion2 pos=5206
GCGAGTTTTGACGAGATGGTGCCGGAGTTCATCGAAAAAATGGACGAGGCACTGCTGAAATTGGTTTTGTATTTGGGGAGCAATGGCGATGAAGCATCC
```

The following are the signatures extracted exclusively from the `inclusion3.fasta` input:

```text
>0 score=1.0000 in=1.0000 ex=0.0000 len=103 ref=inclusion3 pos=99
TAGTCTCCAGGATTCCCGGGGCGGTTCAGATAATCTTAGCATTGACCGCCTTTATATAGAAGCTGTTATTCAAGAAGCATTTTCAAGCAGTGATGTAAGAAAA
>2 score=0.9833 in=0.9833 ex=0.0000 len=100 ref=inclusion3 pos=5203
GCGAGTTTTAACGAGATGGTGCCGGAGTTCATCGAAAAAATGGACCGAGGCACTGCTGAAATTGGTTTTGTATTTGGGGAGCAATGGCGATGAAGCATCC
>1 score=0.9792 in=0.9979 ex=0.0187 len=640 ref=inclusion3 pos=3492
CGCGGGCGATATTTTCACAGCCATTTTCAGGAGTTCAGCCATGAACGCTTATTACATTCAGGATCGTCTTGAGGCTCAGAGCTGGGAGCGTCACTACCAGCAGATCGCCCGTGAAGAAAAGAGGCAGAACTGGCAGACACATGGAAAAAGGCCTGCCCCAGCACCTGTTTGAATCGCTATGCATCGATCATTTGCAACGCCACGGGGCCAGCAAAAAATGCCATTACCCGTGCGTTTGATGACGATGTTGAGTTTCAGGAGCGCATGGCAGAACACATCCGGTACTGGTTGAAACCATTGCTCACCACCAGGTTGATATTGATTCAGAGGTATAAAACGAATGAGTACAGCACTCGCAACGCTGGCAGGGAAGCTGGCTGAACGTGTCGGCATGGATTCTGTCGACCCACAGGAACTGATCACCACTCTTCGCCAGACGGCATTTAAAGGTGATGCCAGCGATGCGCAGTTCATCGCATTGCTGATCGTCGCCAACCAGTACGGTCTTAATCCGTGGACGAAAGAAATTTACGCCTTTCCTGATAAGCAGAACGGCATCGTTCCGGTGGTGGGCGTTGATGGCTGTCCCGTATCATCAATGAAAACCAGCAGTTTGAGGCATGGTACTTTGAGCAGGACA
```

The output from these files appears very similar, as is expected when Neptune identifies highly discriminatory signatures from a homogeneous data set. However, there are some slight differences between some of these signatures. For example, signature ID `2` in `inclusion1.fasta` corresponds to `2` in `inclusion2.fasta`, but corresponds to ID `1` in `inclusion3.fasta`. The lengths (98, 99, 100) and scores (0.9966, 0.9933, 0.9833) are slightly different. Another slight difference is the sequence similarity of signature ID `1` in `inclusion3.fasta` with exclusion sequence:

```text
>1 score=0.9792 in=0.9979 ex=0.0187 len=640 ref=inclusion3 pos=3492
CGCGGGCGATATTTTCACAGCCATTTTCAGGAGTTCAGCCATGAACGCTTATTACATTCAGGATCGTCTTGAGGCTCAGAGCTGGGAGCGTCACTACCAGCAGATCGCCCGTGAAGAAAAGAGGCAGAACTGGCAGACACATGGAAAAAGGCCTGCCCCAGCACCTGTTTGAATCGCTATGCATCGATCATTTGCAACGCCACGGGGCCAGCAAAAAATGCCATTACCCGTGCGTTTGATGACGATGTTGAGTTTCAGGAGCGCATGGCAGAACACATCCGGTACTGGTTGAAACCATTGCTCACCACCAGGTTGATATTGATTCAGAGGTATAAAACGAATGAGTACAGCACTCGCAACGCTGGCAGGGAAGCTGGCTGAACGTGTCGGCATGGATTCTGTCGACCCACAGGAACTGATCACCACTCTTCGCCAGACGGCATTTAAAGGTGATGCCAGCGATGCGCAGTTCATCGCATTGCTGATCGTCGCCAACCAGTACGGTCTTAATCCGTGGACGAAAGAAATTTACGCCTTTCCTGATAAGCAGAACGGCATCGTTCCGGTGGTGGGCGTTGATGGCTGTCCCGTATCATCAATGAAAACCAGCAGTTTGAGGCATGGTACTTTGAGCAGGACA
```

This signature had some similarity with exclusion sequence, represented by the `ex=0.0187`, and indicates a very small amount of imprecise discrimination in this signature. This example further illustrates that the `score` (0.9792) is the sum of the `in` (0.9979) and `ex` (0.0187) values.

These differences in signatures from each inclusion target are a consequence of input sequence variation. The user's discretion will be required in determining which of these are most appropriate. Nonetheless, as described above, Neptune will attempt to consolidate these signatures into a single output file, if a single answer is desirable.
