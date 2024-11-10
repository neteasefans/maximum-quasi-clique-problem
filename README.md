# An opposition-based memetic algorithm for the maximum quasi-clique problem (MQCP)
This repository includes the source code of the proposed OBMA algorithm published in an EJOR paper titled with "An opposition-based memetic algorithm for the maximum quasi-clique problem".

The four sets of 122 benchmark instances are from the classic DIMACS, BHOSLIB, University of Florida Sparse Matrix Collection and Stanford Large Network Dataset Collection libraries (refer to [3]). To facilitate the further research, we upload the instances here. Due to the too large size of the instances from Stanford Large Network Dataset Collection, we do not upload them.

We made comparisons between OBMA and some state-of-the-art methods from the following MQCP works:

[1] Oliveira, A. B. , Plastino, A. , & Ribeiro, C. C. (2013). Construction heuristics for the maximum cardinality quasi-clique problem. In Abstracts of the 10th metaheuris- tics international conference (MIC 2013) (p. 84). Singapore.

[2]. Pinto, B. Q. , Plastino, A. , Ribeiro, C. C. , & Rosseti, I. (2015). A biased random key genetic algorithm to the maximum cardinality quasi-clique problem. In Abstracts of the 11th metaheuristics international conference (MIC 2015). agadir.

[3]. Pinto, B. Q. , Ribeiro, C. C. , Rosseti, I. , & Plastino, A. (2018). A biased random-key genetic algorithm for the maximum quasi-clique problem. European Journal of Operational Research, 271 (3), 849–865 .

Please cite our work as:

Zhou, Q., Benlic, U., & Wu, Q. (2020). An opposition-based memetic algorithm for the maximum quasi-clique problem. European Journal of Operational Research, 286(1), 63-83.

** Instructions to use the source code of OBMA

*** To compile:

q.zhou$ make

q.zhou$

*** To run:

q.zhou$ ./OBMA_MQCP ./input_file gamma time seed ./output_stat_file ./output_sol_file

(where input_file is the instance name, gamma is the required density, time is the cutoff time, seed is the random seed, such as 1, 2, ...,\
			output_stat_file is a file used to store the running information, output_sol_file stores the solution information)

q.zhou$

*** To clean

q.zhou$ make clean

q.zhou$
