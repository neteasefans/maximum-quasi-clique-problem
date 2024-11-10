#ifndef _COMMON_FUNC_H
#define _COMMON_FUNC_H


#include"solution.h"
#include"instance.h"

int compute_penalty_f(solution sol, instance ins);
bool is_feasible_solution(instance ins, int fixed_k, int necessary_edge, solution sol);
void copy_solution(solution source, solution &des, instance ins);
void free_memory(instance);
bool greedy_create_noK(solution &sol, instance ins, double gamma, int &fixed_k);
double observe_sim(solution *pop, int pop_size);
int compInc(const void *a, const void *b);
void perturbation(solution &sol_cur, const solution sol_best, instance ins, int fixed_k);
void perturbation_greedy(solution &sol_cur, const solution sol_best, instance ins, int fixed_k);
#endif