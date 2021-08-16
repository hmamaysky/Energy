# Data Key

Artcount – weekly average of daily article counts

Entropy – weekly average of daily entropy; daily counts use word-weighted averages

ftopic[N] – for N=1,…,7 these are the weekly average of daily topic frequencies; daily measures use word-weighted averages.  Some articles don’t have any energy words, and these articles have all 7 frequencies set to zero.

stopic[N] – for N=1,…,7, the weekly average topical sentiment; daily topical sentiment is calculated using word-weighted averaging

N – topic number

1. Company (Co)
1. Global Oil Markets (Gom)
1. Environment (Env)
1. Energy/Power Generation (Epg)
1. Crude Oil Physical (Bbl)
1. Refining & Petrochemicals (Rpc)
1. Exploration & Production (Ep)

ftopic[N]_4wk – 4-week averages of ftopic[N]

stopic[N]_4wk – 4-week averages of stopic[N]

f[topic name]_Tue – 4-week average of frequency for [topic name] (shown above), calculated as of 2:30pm on the given Tuesday; e.g., gRpc_Tue

s[topic name]_Tue – 4-week average of topical sentiment for [topic name] (shown above), calculated as of 2:30pm on the given Tuesday; e.g., sGom_Tue

# File names
The _prices_ file contains variables for the Friday regressions.  The _physical_ file contains variables for the Tuesday regression.
