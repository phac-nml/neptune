# Change Log

All notable changes to Neptune will be documented in this file.

## 2.0.0 ##

2024-10-21

This release updates Neptune to Python3, removes DRMAA support, fixes a crash when no signatures are produced, and updates the installation process.

### Changed ###

- Python3 has replaced Python2.
- Improved and updated the installation process.

### Fixed ###

- Fixed a crash that occurred when candidate signatures were of such low quality (as a consequence of ambiguous sequence characters) that these regions could not be aligned with themselves using BLAST.

### Removed ###

- DRMAA support.

## 1.2.5 ##

2017-05-03

This release provides fixes for ambiguous crashes and improvements to the code quality.

### Changed ###

- We have made an effort to improve the readability of function comments in the source code.

### Fixed ###

- When running Neptune in parallel (non-DRMAA mode), runtime errors in forked jobs now correctly inform the calling process instead of hanging forever with no meaningful error message. Additionally, the runtime error message is reported to the user. This relates the a known error in Python 2.7 (https://bugs.python.org/issue9400).
- Inputs containing no A, C, G, or T characters will now cause an appropriate runtime error with an informative message about this problem.
- Lowercase characters are no longer ingnored when calculating the GC content of inputs.

## 1.2.4 ##

2017-02-27

This release makes several small improvements, including: reducing the standard output clutter, adding timings to stages, and updating the documentation.

### Added ###

- Links in the README to the manual.
- Walkthrough to the manual.
- Example data to test the software.
- Timings for stages.

### Changed ###

- Improved clarity in manual.
- Codeblocks in the manual.

### Removed ###

- Considerable clutter has been removed from standard output.

## 1.2.3 ##

2016-07-11

This release simplifies the installation process.

### Added ###

- A script for automatically installing Debian dependencies.

### Changed ###

- The dependencies have changed. Several are now installed as part of Neptune.
- The Neptune installation no longer requires security privilages.
- Neptune may be installed multiple times in multiple locations.
- NumPy and SciPy are now installed using pip.

### 1.2.2 ###

2016-04-06

This release includes some Galaxy improvements and fixes a signature scoring problem.

### Changed ###

- Galaxy XML files have been updated to use different packages of Python.

### Fixed ###

- A bug confusing inclusion and exclusion has been fixed.

## 1.2.1 ##

2016-03-23

This release of Neptune adds support for Galaxy.

### Added ###

- Galaxy-related files: capsules, XML files.

### Changed ###

- Neptune.py and Execution.py are now compatible with Galaxy.

## 1.2.0 ##

2016-03-18

This release of Neptune allows for execution on a single machine without requiring DRMAA. Furthermore, several command line parameters have been modified.

### Added ###

- Neptune may be run in parallel on a single machine without DRMAA.
- "--version" command line option.

### Changed ###

- Several command-line parameters have been changed.
- The "--parallelization" / "-p" parameter effects all parallelization.
- The exclusion score is now displayed as a positive number.

## 1.1.1 ##

2016-02-24

This release of Neptune updates the installation instructions to be more informative.

### Changed ###

- Updated README and manual installation instructions.
- Modified the style of code examples in the manual.

## 1.1.0 ##

2016-01-19

This release of Neptune introduces a simple signature consolidation step, which consolidates signatures produced from multiple files into a single file. Furthermore, the software has been updated to be compatible with the Slurm scheduler.

### Added ###

- Neptune now automatically consolidates signatures into a single file.
- DRMAA job names.
- Neptune now maintains DRMAA log files.
- Added the ability to specify the BLAST seed size.

### Changed ###

- The run receipt has been reorganized.
- Removed some unneeded output files.
- Removed some unneeded print statements.
- Removed the --verbose parameter. There was no functionality.

### Fixed ###

- Neptune is now compatible with the Slurm scheduler.
- Updated PEP8/Flake8 code compliance (W503).

## 1.0.0 ##

2015-11-18

This is the initial release of Neptune.
