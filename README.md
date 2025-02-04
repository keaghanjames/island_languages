NOTE: Since publishing this code, the authors have been made awayre of an error in Dow (2008) (see preprint by Alexander Koplenig for details https://osf.io/preprints/osf/xbu2k_v1). Reanalyysus using an alternative method (semiparametric eigenvector spatial filterting) suggests that this issue does not have any substantive impact on the conclusions of the study. 


# island_languages
Code, data and output results pertaining to the manuscript 'Islands are engines of language diversity'.

Please note that four distance matrices are required to replicate these analyses, the geographic and phylogenetic distances among languages and the geographic and phyloSor distances between islands and island language communities. These matrices exceed the maximum allowable filesize of Github. They can be accessed at dropbox via this link: https://www.dropbox.com/scl/fi/1toxqgzjfswapbqgkhdxy/distance_matrices.zip?rlkey=eurdc6n5oecn5at1tdrcr0ltg&dl=0 

Additionally, the endangerment analysis requires an R workspace used in Bromham et al.'s 2022 study of language endangerment which contains their best-fitting model of language endangerment. This workspace also exceeds the allowable limit of Github but can be accessed via https://www.dropbox.com/scl/fi/vbyyvqtwk8gj4q63i424f/bestmodel?rlkey=s3yqzhlz8shnvks5lg2d8ea6g&dl=0  

Code/

Cookie_Cutter.R - code used to generate a null model where language diversity of islands is compared to a sample of a mainland/continent of the same size and shape and at a similar latitude.

endangerment.R - code used to assess the impact of island endemism on language endangerment. The code is based on Bromham et al.'s 2022 study of global predictors of language endangerment and uses the best-fitting model of language endangerment from that study. 

gram_complexity.R - code used to implement a phylospatial analysis of the influence of island endemism on Shcherbakova et al. 2023's indices of grammatical complexity, fusion/boundedness and informativity.

island_diversity_driversLR/IE/WE.R - code used to implement a phylospatial analysis of the predictors of language richness, number of island endemic languages and language-weighted endemism across all inhabited islands in our dataset. Note that there is a separate analysis script for each diversity index.

diversity_model_analysis.R - code for jointly summarising the output of the three analyses of the island language diversity indices.

phoneme_count.R - code used to implement a phylospatial analysis of the influence of island endemism on phoneme counts using both the Phoible and Creanza et al. (2019) databases.

Please not that directories will need be updated by the user in all instances. 

Data/

Continents/mainlands_25km2/island_polygons_united.rds - spatial polygon objects of the two non-island treatments used in the study as well as the islands analysed. These maps are called by various scripts included in Code/.

island_summary.csv - dataset used for island language diversity analyses. Includes the language richness, number of island endemic languages and language-weighted endemism for each island in the dataset, as well as all observations for predictor variables. Called by island_diversity_driversLR/IE/WE.R to implement phylospatial analysis of island language diversity. 

grambank/phoneme_data_by_lang.csv - Data used to implement phylospatial analyses of grammatical complexity and phonement counts. Dataset is called by gram_complexity.R and phoneme_count.R.

languoid_data_for_analysis.csv - Supplementary data from Bromham et al. 2022 including L1 population sizes for languages. Called by various scripts in Code/.

Output/

cookie_cutter_summary.csv - Summary output of 10,000 iterations of the null model using both a non-island treatment of mainland (25k km2) and continent. The dataset includes the language richness, weighted endemism and proportion of range scores for all pairs of island and mainland samples, along with the polygon area and centroid coordinate. 

decent_modelsComplexity.csv -  Summaries of reasonable models (delta BIC < 6) for the grammatical complexity and phoneme count analyses. The table includes each model's BIC, delta, weight, pseudo r2 (pR2), number of parameters (Nparam), sample size (N) and coefficient estimates for each predictor including 95%CI. For each coefficient, significance is indicated by the number of astrices (* = p < 0.05, ** = p < 0.01, *** = p < 0.001).

decent_models_diversity.csv - Summaries of reasonable models (delta BIC < 6) for the island diversity analyses. The table includes each model's BIC, delta, weight, pseudo r2 (pR2), number of parameters (Nparam) and coefficient estimates for each predictor including 95%CI. For each coefficient, significance is indicated by the number of astrices (* = p < 0.05, ** = p < 0.01, *** = p < 0.001).

model_summaries_Gram_boundness/informativity.csv - these tables show the BIC, delta, weights, number of parameters (Nparam) and parameter combinations but not coefficient values, for all 32 possible models fit for the grammatical complexity analyses. 

model_summaries_Gram_PC_Creanza/Phoible.csv - these tables show the BIC, delta, weights, number of parameters (Nparam) and parameter combinations but not coefficient values, for all 32 possible models fit for the phoneme count analyses. 


model_summaries_Gram_LR/IE/WE.csv - these tables show the BIC, delta, weights, number of parameters (Nparam) and parameter combinations but not coefficient values, for all 128 possible models fit for island diversity analyses. 


