# Energy paper (Calomiris, &#199;akir Melek, Mamaysky 2021) code

* `Analysis/` -- directory contains the SDF construction, OOS runs tests, and outliers plots.
* `InSample/` -- all the in-sample forward selection and bootstrapping code
* `OutOfSample/` -- various out-of-sample tests
* `Reconciliation` -- code to comapre new topic model (circa 2024) to the prior topic model (from the KC Fed working paper version)
* `TextProcessing/` -- Reuters archive analytics
* `VariableConstruction/` -- code to construct dependent and explanatory variables
* `data/` -- data output of the project
  * the text pull data for the original MS submission are in `/shared/share_mamaysky-glasserman/energy_drivers/2020-11-16/DataProcessing`
    * `info/` contains the monthly lists of articles with text stats, like sentiment, entropy, topic allocations
    * `Louvain/` contains the topic model, i.e., words in topics
    * `topic_allocation/` shows headlines of articles and their topic allocations in monthly files

Each directory (for the most part) contains its own set of documentation.
