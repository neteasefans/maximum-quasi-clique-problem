#ifndef _SOLUTION_H
#define _SOLUTION_H
#include"instance.h"
class solution
{
public:
	int *sol_arr;				//sol_arr[i] =0 indicates that vertex i is in the current solution
	int *flag_arr;
	int penalty_f;				//the number of edges induced by current solution
	int sol_len;
	bool greedy_create(solution &, instance ins, int fixed_k);
	bool random_create(solution &sol, instance ins, int fixed_k);
	bool frequency_based_initialize(solution &sol, instance ins, int fixed_k, int *freq);
};
#endif
