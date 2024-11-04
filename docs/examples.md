# Examples

## Basic Execution

The following basic example will report all of the signatures that are sufficiently shared by the (FASTA) sequences in the inclusion directory and sufficiently absent from the (FASTA) sequences in the exclusion directory. Neptune will automatically calculate many of the parameters used in this execution.

```bash
neptune
    --inclusion inclusion_directory/
    --exclusion exclusion_directory/
    --output output_directory/
```

The output of immediate interest will be located in the follow file:

    output_directory/consolidated/consolidated.fasta

This file will contain a consolidated list of signatures, sorted by their Neptune score, which is a combined estimate of sensitivity and specificity. The signatures with higher scores, near the top of the file, are considered the most discriminatory signatures.

## Faster Execution

The following example highlights options that allow Neptune to run faster when running in parallel mode (default). It will attempt to run Neptune on 16 parallel processes (`--parallelization`) and parallelize *k*-mer counting and aggregation into 64 tasks (`--organization`) distributed over the 16 parallel processes available.

```bash
neptune
    --inclusion inclusion_directory/
    --exclusion exclusion_directory/
    --output output_directory/
    --parallelization 16
    --organization 3
```

## Specifying File Locations

You may wish to specify particular files used in signature discovery. This may be important when specifying references for signature extraction:

```bash
neptune
    --inclusion inclusion_dir/ in1.fasta in2.fasta
    --exclusion exclusion_dir/ ex1.fasta ex2.fasta
    --reference in1.fasta in2.fasta
    --output output/
```