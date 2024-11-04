# Neptune

A genomic signature is a genomic locus that is sufficiently represented in an inclusion group, and sufficiently absent from a background, or exclusion group. A signature might correlate genomic features with phenotypic traits, such as the presence of a gene with increased organism pathogenicity.

Neptune locates genomic signatures using an exact *k*-mer matching strategy while accommodating *k*-mer mismatches. The software identifies sequences that are sufficiently represented within inclusion targets and sufficiently absent from exclusion targets. The signature discovery process is accomplished using probabilistic models instead of heuristic strategies. Neptune may be leveraged to reveal discriminatory signature sequences to uniquely delineate one group of organisms, such as isolates associated with a disease cluster or event, from unrelated sporadic or environmental microbes.

## Neptune v2.0.0

2024-10-21

This release updates Neptune to Python3, removes DRMAA support, fixes a crash when no signatures are produced, and updates the installation process.

## Resources

* **Source**: [https://github.com/phac-nml/neptune](https://github.com/phac-nml/neptune)
* **Installation**: [https://phac-nml.github.io/neptune/install/](https://phac-nml.github.io/neptune/install/)
* **Walkthrough**: [https://phac-nml.github.io/neptune/walkthrough/](https://phac-nml.github.io/neptune/walkthrough/)

## Contact

* **Eric Marinier**: eric.marinier@phac-aspc.gc.ca
* **Gary van Domselaar**: gary.vandomselaar@phac-aspc.gc.ca
