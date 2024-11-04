# Walkthrough

## Overview

The purpose of this walkthrough will be to illustrate a simple, but complete example of using Neptune to locate discriminatory sequences. We will identity signature sequences within an artificial data set containing three inclusion sequences and three exclusion sequences. The output will be a list of signatures, sorted by score, for each inclusion target, and one consolidated signatures file, sorted by signature score, containing signatures from all inclusion targets.

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

The inclusion and exclusion directories each contain three FASTA format genomes. The genomes all have some insertions and deletions that differentiate them from each other. However, the three inclusion genomes primarily differ from the three exclusion genomes in that they share large sequences that are absent from all exclusion genomes.

## Running Neptune

Neptune will automatically calculate many of the parameters that might otherwise be specified by the user, such as the minimum number of targets signature sequence must be present within for it to be considered shared sequence. At minimum, Neptune requires the user specify the inclusion sequences, exclusion sequences, and an output directory. We will provide Neptune inclusion and exclusion sequences in the form of FASTA file genomes located within directories. The following command will run Neptune on the example data and output to the specified directory:

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
Neptune v1.2.4

Estimating k-mer size ...
k = 15

k-mer Counting...
Submitted 12 jobs.
0.002063 seconds

k-mer Aggregation...
Submitted 65 jobs.
0.010473 seconds

Signature Extraction...
Submitted 6 jobs.
0.000771 seconds

Signature Filtering...
Submitted 2 jobs.
Submitted 6 jobs.
0.002498 seconds

Consolidate Signatures...
Submitted 1 jobs.
0.000411 seconds

Complete!
```

### Consolidated Signatures

As we did not specify references from which to extract signatures, Neptune will automatically investigate all inclusion genomes for signatures and consolidate those signatures into a single consolidated signature file. The `output/consolidated/consolidated.fasta` file contains these consolidated signatures. This file may be understood as the final output of the application. The following FASTA output is from the consolidated signatures file produced from this example:

```text
>1.0 score=1.0000 in=1.0000 ex=0.0000 len=103 ref=inclusion1 pos=99
TAGTCTCCAGGATTCCCGGGGCGGTTCAGATAATCTTAGCATTGACCGCCTTTATATAGAAGCTGTTATTCAAGAAGC...
>1.1 score=0.9979 in=0.9979 ex=0.0000 len=640 ref=inclusion1 pos=3497
CGCGGGCGATATTTTCACAGCCATTTCAGGAGTTCAGCCATGAACGCTTATTACATTCAGGATCGTCTTGAGGCTCAG...
>1.2 score=0.9966 in=0.9966 ex=0.0000 len=98 ref=inclusion1 pos=5209
GCGAGTTTTGCGAGATGGTGCCGGAGTTCATCGAAAAAATGGACGAGGCACTGCTGAAATTGGTTTTGTATTTGGGGA...
```

The FASTA header contains information relavent to the identified signature. A detailed explanation of this information is located within the [output section](output/#sorted-signatures) of this documentation. The `score` is the sum of the `in` (inclusion/sensitivity) and `ex` (exclusion/specificity) scores, and represents a combined measure of sensitivity and specificity. The `length` describes the length of the signature in bases. The `ref` (reference) and `pos` (position) describe the location of the signature within the reference FASTA record it was extracted from.

In this example, Neptune identified three signatures: 1.0, 1.1, and 1.2 of lengths 103, 640, and 98, respectively. We see that all of these signatures originated from the inclusion1 reference. These signatures were located at positions 99, 3497, and 5209 within the inclusion1 reference. These signatures are of very high quality, within the context of our data set, with scores of 1.0000, 0.9979, and 0.9969, within the possible range of score values from -1.00 to +1.00.

### Sorted Signatures

If we're interested in looking at the signatures produced from each individual inclusion target, we need to investigate the output in the `output/sorted` directory. The following are the signatures extracted exclusively from the inclusion1.fasta target:

```text
>0 score=1.0000 in=1.0000 ex=0.0000 len=103 ref=inclusion1 pos=99
TAGTCTCCAGGATTCCCGGGGCGGTTCAGATAATCTTAGCATTGACCGCCTTTATATAGAAGCTGTTATTCAAGAAGC...
>1 score=0.9979 in=0.9979 ex=0.0000 len=640 ref=inclusion1 pos=3497
CGCGGGCGATATTTTCACAGCCATTTCAGGAGTTCAGCCATGAACGCTTATTACATTCAGGATCGTCTTGAGGCTCAG...
>2 score=0.9966 in=0.9966 ex=0.0000 len=98 ref=inclusion1 pos=5209
GCGAGTTTTGCGAGATGGTGCCGGAGTTCATCGAAAAAATGGACGAGGCACTGCTGAAATTGGTTTTGTATTTGGGGA...
```

The following are the signatures extracted exclusively from the \textit{inclusion2.fasta} target:

```text
>0 score=1.0000 in=1.0000 ex=0.0000 len=103 ref=inclusion2 pos=99
TAGTCTCCAGGATTCCCGGGGCGGTTCAGATAATCTTAGCATTGACCGCCTTTATATAGAAGCTGTTATTCAAGAAGC...
>1 score=0.9979 in=0.9979 ex=0.0000 len=640 ref=inclusion2 pos=3494
CGCGGGCGATATTTTCACAGCCATTTTCAGGAGTTCAGCCATGAACGCTTATTACATTCAGGATCGTCTTGAGGCTCA...
>2 score=0.9933 in=0.9933 ex=0.0000 len=99 ref=inclusion2 pos=5206
GCGAGTTTTGACGAGATGGTGCCGGAGTTCATCGAAAAAATGGACGAGGCACTGCTGAAATTGGTTTTGTATTTGGGG...
```

The following are the signatures extracted exclusively from the \textit{inclusion3.fasta} target:

```text
>0 score=1.0000 in=1.0000 ex=0.0000 len=103 ref=inclusion3 pos=99
TAGTCTCCAGGATTCCCGGGGCGGTTCAGATAATCTTAGCATTGACCGCCTTTATATAGAAGCTGTTATTCAAGAAGC...
>2 score=0.9833 in=0.9833 ex=0.0000 len=100 ref=inclusion3 pos=5203
GCGAGTTTTAACGAGATGGTGCCGGAGTTCATCGAAAAAATGGACCGAGGCACTGCTGAAATTGGTTTTGTATTTGGG...
>1 score=0.9792 in=0.9979 ex=0.0187 len=640 ref=inclusion3 pos=3492
CGCGGGCGATATTTTCACAGCCATTTTCAGGAGTTCAGCCATGAACGCTTATTACATTCAGGATCGTCTTGAGGCTCA...
```

The output from these files appears very similar, as is expected when Neptune identifies highly discriminatory signatures from a homogeneous data set. However, there are some slight differences between some of these signatures. For example, the signatures in each of these output files have corresponding ID numbers and some of these signatures have slight differences. However, because Neptune assigns signature IDs arbitrarily, this correspondence will usually never happen when using real data. Nonetheless, we see that signature ID 2 is slightly different sizes in all three inclusion targets (5209, 5206, and 5203) with slightly different scores (0.9966, 0.9933, 0.9833). Another slight difference between the signatures is the sequence similarity of signature ID 1 in inclusion3.fasta with exclusion sequence:

```text
>1 score=0.9792 in=0.9979 ex=0.0187 len=640 ref=inclusion3 pos=3492
CGCGGGCGATATTTTCACAGCCATTTTCAGGAGTTCAGCCATGAACGCTTATTACATTCAGGATCGTCTTGAGGCTCA...
```

This signature had some similarity with exclusion sequence, represented by the `ex=0.0187`, and indicates a small amount of imprecision in this signature. This example also illustrates that the `score` (0.9792) is the sum of the `in` (0.9979) and `ex` (0.0187) values.

These differences in signatures from each inclusion target are a consequence of sequence differences. The user's discretion will be required in determining which of these are most appropriate. Nonetheless, as described above, Neptune will attempt to consolidate these signatures into a single output file, if a single answer is desirable.
