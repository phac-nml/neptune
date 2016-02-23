# Change Log

All notable changes to Neptune will be documented in this file.

## [1.1.1] - 2016-02-23

This release of Neptune updates the installation instructions to be more
informative.

### Changed
- Updated README and manual installation instructions.

## [1.1.0] - 2016-01-19

This release of Neptune introduces a simple signature consolidation step, which
consolidates signatures produced from multiple files into a single file.
Furthermore, the software has been updated to be compatible with the Slurm
scheduler.

### Added
- Neptune now automatically consolidates signatures into a single file.
- DRMAA job names.
- Neptune now maintains DRMAA log files.
- Added the ability to specify the BLAST seed size.

### Changed
- The run receipt has been reorganized.
- Removed some unneeded output files.
- Removed some unneeded print statements.
- Removed the --verbose parameter. There was no functionality.

### Fixed
- Neptune is now compatible with the Slurm scheduler.
- Updated PEP8/Flake8 code compliance (W503).

## [1.0.0] - 2015-11-18

This is the initial release of Neptune.
