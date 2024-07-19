#ifndef _LOCAL_SEARCH_H
#define _LOCAL_SEARCH_H
#include"solution.h"
#include"instance.h"
void update_deg_matrix(int deg_matrix[], instance ins, int v_out, int v_in);
bool tabu_search_based_swap(solution &sol_cur, solution &sol_best, instance ins, int &iter, int iter_max, int tabu_length, int tabu_tenure, int necessary_edge, int fixed_k, double edge_density, int *freq);
bool tabu_search_constrained_swap(solution &sol_cur, solution &sol_best, instance ins, int &iter_ms, int iter_max, int tabu_length, int tabu_tenure, int necessary_edge, int fixed_k, double edge_density, int *freq);
bool multi_start(solution &sol_cur, solution &sol_best, solution &sol_ms_best, instance ins, int iter_max, int tabu_length, int tabu_tenure, int fixed_k, int &best_k, double gamma, double edge_density, double &end_time);
bool local_search(solution &sol_cur, solution &sol_best, solution &sol_ms_best, instance ins, int &iter, int iter_max, int tabu_length, int tabu_tenure, int fixed_k, int &best_k, double gamma, double edge_density, double &end_time);
#endif