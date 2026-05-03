The main script is `tz_enterovirus_tree_manipulation.R` which takes as inputs:

* vp1_tree_with_dates.nex
* full_length_tree_with_dates.nex
* all_microreact_format_metadata.csv

## The metadata file

There are two metadata files, `all_minimal_metadata.csv` (307 records) and `coxsackie_a24_full_length_minimal_metadata.csv` (67 records). Each of these contains metadata fields
Accession, Country and Collection_Date and can be reformated into Microreact format with `format_for_microreact.py`.

