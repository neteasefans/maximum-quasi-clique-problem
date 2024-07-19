#include<iostream>
#include<time.h>
#include<stdlib.h>
#include<string.h>
#include"solution.h"
#include"common_func.h"
#include"local_search.h"
#include"memetic.h"
using namespace std;

//初始化种群
bool oppo_initial_population(solution &sol_cur, solution &sol_best, solution &sol_ms_best, instance ins, int &iter, int iter_max, int tabu_length, int tabu_tenure, int fixed_k,
	int &best_k, double gamma, double edge_density, double &end_time, solution *pop, int pop_size)
{	
	solution sol_cur_oppo;
	sol_cur_oppo.flag_arr = new int[ins.num_node];
	sol_cur_oppo.sol_arr = new int[ins.num_node];
	bool succ=false;
	for (int i = 0; i < pop_size; i++)
	{		
		sol_cur.random_create(sol_cur, ins, fixed_k);
		generate_opposition_solution(sol_cur, sol_cur_oppo, ins, fixed_k);
		succ = local_search(sol_cur, sol_best, sol_ms_best, ins, iter, iter_max, tabu_length, tabu_tenure, fixed_k, best_k, gamma, edge_density, end_time);
		copy_solution(sol_best, sol_cur, ins);
		if(succ)
			break;
		succ = local_search(sol_cur_oppo, sol_best, sol_ms_best, ins, iter, iter_max, tabu_length, tabu_tenure, fixed_k, best_k, gamma, edge_density, end_time);
		copy_solution(sol_best, sol_cur_oppo, ins);
		if(succ)
			break;
		if (sol_cur.penalty_f>sol_cur_oppo.penalty_f)
			copy_solution(sol_cur, pop[i], ins);
		else
			copy_solution(sol_cur_oppo, pop[i], ins);
	}
	delete [] sol_cur_oppo.flag_arr; 
	delete [] sol_cur_oppo.sol_arr; 
	return succ;
}

//初始化种群
bool initial_population(solution &sol_cur, solution &sol_best, solution &sol_ms_best, instance ins, int &iter, int iter_max, int tabu_length, int tabu_tenure, int fixed_k,
	int &best_k, double gamma, double edge_density, double &end_time, solution *pop, int pop_size)
{	
	bool succ=false;
	for (int i = 0; i < pop_size; i++)
	{
		pop[i].random_create(pop[i], ins, fixed_k);
		succ = local_search(pop[i], sol_best, sol_ms_best, ins, iter, iter_max, tabu_length, tabu_tenure, fixed_k, best_k, gamma, edge_density, end_time);
		copy_solution(sol_best, pop[i], ins);
		if(succ)
		   return succ;
		pop[i].random_create(pop[i], ins, fixed_k);
		bool succ = local_search(pop[i], sol_best, sol_ms_best, ins, iter, iter_max, tabu_length, tabu_tenure, fixed_k, best_k, gamma, edge_density, end_time);
		if(sol_best.penalty_f>sol_cur.penalty_f)
	          copy_solution(sol_best, pop[i], ins);
	        else 
		  copy_solution(sol_cur, pop[i], ins);
		if(succ)
		   return succ;
	}
	return succ;
}

void uniform_cross_over(solution &sol_cur, solution*pop, int pop_size, instance ins, int fixed_k)
{	 
	int rand_x;
	for (int i = 0; i < ins.num_node; i++)
		sol_cur.flag_arr[i] = 0;
	int p1 = rand() % pop_size;
	int p2 = rand() % pop_size;
	while (p1 == p2)
		p2 = rand() % pop_size;

	int min_len = pop[p1].sol_len > pop[p2].sol_len ? pop[p2].sol_len : pop[p1].sol_len;
	for (int i = 0; i < min_len; i++)
	{
		rand_x = rand() % 2;	
		if (rand_x == 0)
			sol_cur.flag_arr[i] = pop[p1].flag_arr[i];
		else
			sol_cur.flag_arr[i] = pop[p2].flag_arr[i];
	}	
	sol_cur.sol_len = 0;
	for (int i = 0; i < ins.num_node; i++)
	{
		if (sol_cur.flag_arr[i] == 1)
			sol_cur.sol_arr[sol_cur.sol_len++] = i;
	}	
	while (sol_cur.sol_len < fixed_k)
	{
		rand_x = rand() % ins.num_node;
		while (sol_cur.flag_arr[rand_x] == 1)	
			rand_x = rand() % ins.num_node;
		sol_cur.sol_arr[sol_cur.sol_len++] = rand_x;
		sol_cur.flag_arr[rand_x] = 1;
	}
	sol_cur.penalty_f = compute_penalty_f(sol_cur, ins);
}

void common_shared_cross_over(solution &sol_cur, solution*pop, int pop_size, instance ins, int fixed_k)
{
	int common_node_len = 0;
	int index1 = 0, index2 = 0;
	int p1 = rand() % pop_size;
	int p2 = rand() % pop_size;
	while (p1 == p2)
		p2 = rand() % pop_size;
	sol_cur.sol_len = 0;
	memset(sol_cur.flag_arr, 0, sizeof(int)*ins.num_node);
	qsort(pop[p1].sol_arr, pop[p1].sol_len, sizeof(pop[p1].sol_arr[0]), compInc);
	qsort(pop[p2].sol_arr, pop[p2].sol_len, sizeof(pop[p2].sol_arr[0]), compInc);
	while (index1 < pop[p1].sol_len && index2 < pop[p2].sol_len)
	{
		if (pop[p1].sol_arr[index1] > pop[p2].sol_arr[index2])
			index2++;
		else if (pop[p1].sol_arr[index1] == pop[p2].sol_arr[index2])
		{
			common_node_len++;		
			int ele = pop[p1].sol_arr[index1];
			sol_cur.sol_arr[sol_cur.sol_len++] = ele;
			sol_cur.flag_arr[ele] = 1;
			index1++;
			index2++;
		}
		else
			index1++;
	}
	//cout << "sol_cur.sol_len=" << sol_cur.sol_len <<",fixed_k="<<fixed_k<< endl;
	while (sol_cur.sol_len < fixed_k)
	{
		int rand_x = rand() % ins.num_node;
		while (sol_cur.flag_arr[rand_x] == 1)
			rand_x = rand() % ins.num_node;
		sol_cur.sol_arr[sol_cur.sol_len++] = rand_x;
		sol_cur.flag_arr[rand_x] = 1;
	}
	//mutation(sol_cur, ins);
	sol_cur.penalty_f = compute_penalty_f(sol_cur, ins);
}

void mutation(solution &sol_cur,instance ins)
{
	double mu = 0.15;
	if ((ins.num_node - sol_cur.sol_len) < mu*ins.num_node)
		return;
	int *address = new int[ins.num_node];
	for (int i = 0; i < sol_cur.sol_len; i++)
	{
		int ele = sol_cur.sol_arr[i];
		address[ele] = i;
	}
	int len = 0;
	while (len < mu*ins.num_node)
	{
		int u = sol_cur.sol_arr[rand() % sol_cur.sol_len];
		int u_address = address[u];
		int v = rand() % ins.num_node;
		while (sol_cur.flag_arr[v])
			v = rand() % ins.num_node;
		sol_cur.sol_arr[u_address] = v;
		sol_cur.flag_arr[u] = 0;
		sol_cur.flag_arr[v] = 1;
		address[v] = u_address;
		len++;
	}
	delete [] address;
}


void generate_opposition_solution(solution sol_cur, solution &sol_cur_oppo, instance ins, int fixed_k)
{
	int len1 = 0;
	int *flag_node_n = new int[ins.num_node];
	int *rand_node_n = new int[ins.num_node];
	memset(flag_node_n, 0, sizeof(int)*ins.num_node);
	sol_cur_oppo.sol_len = 0;
	memset(sol_cur_oppo.flag_arr, 0, sizeof(int)*ins.num_node);
	while (len1 < ins.num_node)
	{
		int rn = rand() % ins.num_node;
		if (flag_node_n[rn] == 0)
		{
			rand_node_n[len1++] = rn;
			flag_node_n[rn] = 1;
		}
	}
	for (int i = 0; i < ins.num_node; i++)
	{
		int node = rand_node_n[i];
		if (sol_cur.flag_arr[node] == 0 && sol_cur_oppo.sol_len < fixed_k)
		{
			sol_cur_oppo.sol_arr[sol_cur_oppo.sol_len++] = node;
			sol_cur_oppo.flag_arr[node] = 1;
		}
	}
	while (sol_cur_oppo.sol_len < fixed_k)
	{
		int rand_x = rand() % ins.num_node;
		while (sol_cur_oppo.flag_arr[rand_x] == 1)
			rand_x = rand() % ins.num_node;
		sol_cur_oppo.sol_arr[sol_cur_oppo.sol_len++] = rand_x;
		sol_cur_oppo.flag_arr[rand_x] = 1;
	}
	sol_cur_oppo.penalty_f = compute_penalty_f(sol_cur_oppo, ins);
	delete [] flag_node_n; 
	delete [] rand_node_n; 
}

void update_population(solution *pop, solution offspring, instance ins, int pop_size)
{	
	int worst_cost = ins.num_node*(ins.num_node + 1);
	int sel_idx;
	for (int i = 0; i < pop_size;i++)
	{
		if (pop[i].penalty_f < worst_cost)
		{
			worst_cost = pop[i].penalty_f;
			sel_idx = i;
		}
	}
	if (offspring.penalty_f > worst_cost)	
		copy_solution(offspring, pop[sel_idx], ins);	
}
