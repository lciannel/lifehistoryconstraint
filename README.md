# lifehistoryconstraint
This respository contains code and data used to estimate habitat constraint metrics for walleye pollock in the Bering Sea and for simulated data.  Code for extimating habitat constraints for the remaining 9 species examined in Ciannelli et al. 2021 follow the pollock template, and is avilable upon request.
This repository also contains code to generate spatially correlated random fields and to estimate habitat constraints on these fields; and data to estimate relationshps between habitat constraint and ontogenetic index acorss all 10 case studies examined in Ciannelli et al. 2021.

**Code**//
BS_Pollock_Constraint4_online.R: estimate habitat constraints on Bering sea pollock data
distance.function.R: fuction to calculate great circle distance, used in pollock script
Habitat_Constraint_4_online.Rmd: estimate relationships between habitat constraint and species ontogenetic stage across all 20 case studies examined
HC_Simulations_5_online.R: Simulate spatial random fields and estimate habitat constraint metrics on each field

\\
**Data**\\
BeringDepth.txt: Bering Sea bathymetry
Constraint_Analysis_4.csv: habitat constraint metrics across all 10 species examined
pollock_eggs.csv: pollock egg data in Bering Sea
pollock_juvenile_adult_catch.csv: pollock juvenile and adult catch data in the Bering Sea
pollock_juvenile_adult_length.csv: pollock juvenile and adult length data in the Bering Sea
pollock_larvae.csv: pollock larval data in the Berig Sea




