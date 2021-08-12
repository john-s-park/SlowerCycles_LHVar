## ABSTRACT

Populations in nature are comprised of individual life histories, whose variation underpins ecological and evolutionary processes. Yet the forces of environmental selection that shape intrapopulation life-history variation are still not well understood, and efforts have largely focused on random (stochastic) fluctuations of the environment. However, a ubiquitous mode of environmental fluctuation in nature is cyclical, whose periodicities can change independently of stochasticity. Here we test theoretically-based hypotheses for whether shortened (‘Fast’) or lengthened (‘Slow’) environmental cycles should generate higher intrapopulation variation of life history phenotypes. We show, through a combination of agent-based modelling and a multi-generational laboratory selection experiment using the tidepool copepod Tigriopus californicus, that slower environmental cycles maintain higher levels of intrapopulation variation. Surprisingly, the effect of environmental periodicity on variation was much stronger than that of stochasticity. Thus, our results show that periodicity is an important facet of fluctuating environments for life-history variation.


## File Descriptions

Two data files “EvolExp_xxxx_data.csv” contain life history trait measurements from Tigriopus californicus laboratory experiment. The two share the first four common columns: “Source_pop”, “Treatment”, “Treat_rep”, and “Ind_rep”.

“Source_pop” denotes arbitrary incubator containers that were used for pre-experiment rearing of stock populations, but was not used for analyses.

“Treatment” is categorial disturbance regime administered to replicate populations. During the data collection phase, “H” stood for High frequency but is subsequently indexed as “Fast” for a more intuitive term; “L” stood for Low frequency but is subsequently indexed as “Slow” for a more intuitive term; “S” is “Stochastic”. Treatments indexed with “_C” were arbitrary ‘controls’ but were irrelevant for analyses in this paper; they are included for completeness.

“Treat_rep” is numerical replicate population index under Treatments.

“Ind_rep” is numerical replicate individual index measured in a population.

Subsequent columns, for both “EvolExp_xxxx_data.csv” files, contain event timing records. “EvolExp_fecundity_data.csv” contains Date+AM/PM records of productions of clutches of juveniles. “EvolExp_maturity_data.csv” contains Date records of appearances of gravid females in mating wells. See main text for detailed descriptions of measurement methods, and associated .R scripts for how analyses using event records were coded.

 
 

R Code:

“EvolExp_results_dataclean.R” reproduces how life history event record data were cleaned, processed, and summarized.

“EvolExp_results_means.R” reproduces summary of intrapopulation means of life history traits. Dependent on “~/dataclean.R”.

“EvolExp_results_variances.R” reproduces summary of intrapopulation variances of life history traits. Dependent on “~/dataclean.R”.

“permutation_intrapopvar_pvals.R” reproduces Monte Carlo permutation method used to test for differences in intrapopulation variances of life history traits. Dependent on “~/dataclean.R”.

“ABM_full_final.R” reproduces full Agent-Based Model simulations, including plots of mean and variance trajectories of realizations, statistical analyses on endpoints, and plots of demographic dynamics under the four environmental regime scenarios.