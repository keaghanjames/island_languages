# island_languages
Code, data and output results pertaining to the manuscript 'Islands are engines of language diversity'.

Please note that four distance matrices are required to replicate this analysis, the geographic and phylogenetic distances among languages and the geographic and phyloSor distances between islands. These matrices exceed the maximum allowable filesize of Github. They can be accessed at dropbox via this link: https://www.dropbox.com/scl/fi/1toxqgzjfswapbqgkhdxy/distance_matrices.zip?rlkey=eurdc6n5oecn5at1tdrcr0ltg&dl=0 

Additionally, the endangerment analysis requires a R workspace used in Bromham et al.'s 2022 study of language endangerment which contains their best fitting model of language endangerment. This workspace also exceeds the allowable limit of Github but can be accessed via https://www.dropbox.com/scl/fi/vbyyvqtwk8gj4q63i424f/bestmodel?rlkey=s3yqzhlz8shnvks5lg2d8ea6g&dl=0  

Code/

Cookie_Cutter.R - code used to generate a null model where language diversity of islands is compared to a sample of a mainland/continent of the same size and shape and at a similar latitude.

endangerment.R - code used to assess the impact of island endemism on language endangerment. The code is based on Bromham et al.'s 2022 study of global predictors of language endangerment and uses the best-fitting model of language endangerment from that study. 

gram_complexity.R - code used to implement a phylospatial analysis of the influence of island endemism on Shcherbakova et al. 2023's indices of grammatical complexity, fusion/boundedness and informativity.

island_diversity_driversLR/IE/WE.R - code used to implement a phylospatial analysis of the predictors of language richness, number of island endemic languages and language-weighted endemism across all inhabited islands in our dataset. Note that there is a separate analysis script for each diversity index.

diversity_model_analysis.R - code for jointly summarising the output of the three analyses of the island language diversity indices.

phoneme_count.R - code used to implement a phylospatial analysis of the influence of island endemism on phoneme counts using both the Phoible and Creanza et al. (2019) databases.

Data/


