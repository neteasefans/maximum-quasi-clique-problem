#ifndef _MEMETIC_H
#define _MEMETIC_H
#include"solution.h"
#include"instance.h"
bool oppo_initial_population(solution &sol_cur, solution &sol_best, solution &sol_ms_best, instance ins, int &iter, int iter_max, int tabu_length, int tabu_tenure,
	int fixed_k, int &best_k, double gamma, double edge_density, double &end_time, solution *pop, int pop_size);
bool initial_population(solution &sol_cur, solution &sol_best, solution &sol_ms_best, instance ins, int &iter, int iter_max, int tabu_length, int tabu_tenure, int fixed_k,
	int &best_k, double gamma, double edge_density, double &end_time, solution *pop, int pop_size);
void prob_same_cross_over(solution &sol_cur, solution *pop, int pop_size, instance ins, int fixed_k);
void common_shared_cross_over(solution &sol_cur, solution*pop, int pop_size, instance ins, int fixed_k);
void mutation(solution &sol_cur, instance ins);
void update_population(solution *pop, solution offspring, instance ins, int pop_size);
void generate_opposition_solution(solution sol_cur, solution &sol_cur_oppo, instance ins, int fixed_k);
#endif
