# Changelog

All notable changes to this project will be documented in this file. The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/), and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [2.9.2] - 2024-06-03
### Updated
- `baseline_correction()` function:
  - Fixed an issue where the `Random=False` part was yielding an error.

## [2.9.1] - 2024-03-08
### Updated
- `filter_PSMs()` and `filter_peptides()` functions:
  - Changed the order of filtering; removal of any TMT channels with at least one NA value is now performed at the end.

## [2.9.0] - 2024-03-06
### Updated
- `filter_PSMs()` and `filter_peptides()` functions with several enhancements:
  - Updates in `filter_PSMs`:
    - Removed filtering for TMT channels with at least one NA value.
    - Other filtering criteria remain unchanged.
  - Updates in `filter_peptides`:
    - Removed filtering for TMT channels with at least one NA value.

## [2.8.5] - 2024-02-06
### Added
- `PSMs_to_Peptide()` function:
  - Separates the functionality from `baseline_correction` for merging PSMs into peptides using specified columns.
  - Auto-detects 'Theo. MH+ [Da]' column using regular expressions even if the name does not fully match.

## [2.8.4] - 2024-01-18
### Updated
- The script now supports PSMs and Peptide file processing for both MS2 and MS3 measurements across both `PD_input` and `plain_text_input` classes.
- `filter_peptides()` function now specifically works for peptide files.
- Added `filter_PSMs()` function for PSM files, removing isolation interference > 50%.
- Updated `extract_heavy()` and `extract_light()` functions to support a broader range of labels.
- Major updates to `baseline_correction()` function:
  - New usage parameters and functionality adjustments.
  - Improved baseline correction process and integration with PSM and Peptides file identification.

## [2.7.0] - 2023-12-05
### Updated
- `filter_peptides()` function:
  - Enhanced to remove sum of intensities equal to 0.
  - Now converts NaN values to 0 for cleaner data handling.