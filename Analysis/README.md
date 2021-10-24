# This code generates some analysis for the paper

* *0714 Subperiod details.xls* is a file from Hongyu about the persistence of OOS performance in our
  7 (?) subperiods.  This details the number of runs, defined as a model that works in period t and
  that continues to work in period t+1.
  
 ## OOS Analysis
 
 * Use the *driver_energy.py* file
 * Run the top line to load in the *energy.py* code as *en*
 * Then to generate the results, run:
    ```
    oos = en.OOSResults()
    oos.calc()
    ```
 * The code snippet below this in the *driver_energy.py* file shows an example of the *p-value* calculation.