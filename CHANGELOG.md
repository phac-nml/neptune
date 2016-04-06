# Change Log

All notable changes to Neptune will be documented in this file.

## [1.2.2] - 2016-04-06

This release includes some Galaxy improvements and fixes a signature scoring
problem.

### Changed
- Galaxy XML files have been updated to use different packages of Python.

### Fixed
- A bug confusing inclusion and exclusion has been fixed.

## [1.2.1] - 2016-03-23

This release of Neptune adds support for Galaxy.

### Added
- Galaxy-related files: capsules, XML files.

### Changed
- Neptune.py and Execution.py are now compatible with Galaxy.

## [1.2.0] - 2016-03-18

This release of Neptune allows for execution on a single machine without
requiring DRMAA. Furthermore, several command line parameters have been
modified.

### Added
- Neptune may be run in parallel on a single machine without DRMAA.
- "--version" command line option.

### Changed
- Several command-line parameters have been changed.
- The "--parallelization" / "-p" parameter effects all parallelization.
- The exclusion score is now displayed as a positive number.

## [1.1.1] - 2016-02-24

This release of Neptune updates the installation instructions to be more
informative.

### Changed
- Updated README and manual installation instructions.
- Modified the style of code examples in the manual.

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
