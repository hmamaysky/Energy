# Energy paper (Calomiris, &#199;akir Melek, Mamaysky 2021) code

Step 1: Construct variables
* `TextProcessing/` -- Reuters archive analytics
* `VariableConstruction/` -- code to construct dependent and explanatory variables
* `Analysis/` -- SDF construction
  * other parts of this (OOS runs tests and and outliers plots) largely deprecated (as of 7/2/2024)
* `Reconciliation/` -- code to compare new topic model (circa 2024) to the prior topic model (from the KC Fed working paper version)
  * the code in the `Reconciliation/` directory compares some files from `2023/` to files in `2020-11-16/`
  * see description in the data section below

Step 2: Run analysis, first in-sample, then OOS
* `InSample/` -- all the in-sample forward selection and bootstrapping code
* `OutOfSample/` -- various out-of-sample tests

Partial output of analysis:
* `data/` -- data generated after all the processing steps, which is then used in the IS and OOS analysis (after Step 1 is done)
* Intermediate data generated during Step 1 are stored in `/shared/share_mamaysky-glasserman/energy_drivers/`
  *  `README` file on the grid contains additional information about directory structure
  * `2023/` -- contains the intermediate data for round 2 of MS submission
    * contains the rolling topic model intermediate data
    * contains the final OOS analysis results
    * also contains a new version (circa 2023-2024) of the full-sample topic model from `2020-11-16/` (KC Fed WP version)
  * `2020-11-16/DataProcessing` -- contains intermediate data for the initial MS submission
    * `info/` contains the monthly lists of articles with text stats, like sentiment, entropy, topic allocations
    * `Louvain/` contains the topic model, i.e., words in topics
    * `topic_allocation/` shows headlines of articles and their topic allocations in monthly files
  * Any other directory in the shared grid folder is deprecated and no longer used (as of 7/2/2024)
