# Data notes

As of 5/31/2021, these are the latest data files with the new timing conventions. Each series has a suffix to indicate when it's calculated. The general philosophy is to take a variable from the prior day if it's not clear that it won't overlap with a day t variable, and to take it from day t if we are sure there is no look-ahead bias. For example, we calculate text-based series to end at 2:30pm on each day, so there is no overlap with these at either the 2:30pm energy market closes or the 4pm stock market closes.

Some of the independent series are from the prior month, and these are indicated with a _monthly suffix in the variable names. For example, we calculate Momentum (Mom_monthly) for month m as the cumulative return from month m-11 to m-1. After construction, we use its month m-1 value for all weeks in month m in order to avoid potential look-ahead bias.

The physical regressions use levels of industrial production and inventories that are reported on Wednesday at 10am. The RHS variables for these are all from Tues or earlier.

The physical variables are on occasion reported on either Thurs or even Friday. For this reason, when we are forecasting price-based series that use the physical oil market variables as forecasting variables, we use Friday and forward market-based observations to make sure there is no overlap with the production or inventory series in that week.

For rdsa (Royal Dutch Shell) we use its Friday or Tuesday value as an independent variable, and use its Monday week t+1 value when it is a dependent variable, because it trades in Europe and has a 10 or 11am NY close.

For hedging pressure (HedgPres_Tue), we use Tuesday value for both price and physical regressions because CFTC provides long/short hedge positions data for “CRUDE OIL, LIGHT SWEET - NEW YORK MERCANTILE EXCHANGE” only on a weekly basis. Although it usually reports Tuesday values, there are some weeks where the values are from Wednesday, Monday or Prior Friday. Unlike Monday and prior Friday values, there is potential look-ahead bias for Wednesday value. So, when week t value is from Wednesday, we replace it with week t-1's Tuesday value.
