#include"common_func.h"
#include<memory.h>
#include<iostream>
#include<string.h>
#include<stdlib.h>
#include"solution.h"
#include"instance.h"
#include"local_search.h"
using namespace std;
#define MINVALUE -99999999
int compute_penalty_f(solution sol, instance ins)
{
	int edge_sol = 0;
	for (int i = 0; i < sol.sol_len; i++)
	{
		int ele1 = sol.sol_arr[i];
		for (int j = i+1; j < sol.sol_len; j++)
		{
			int ele2 = sol.sol_arr[j];
			if (ins.adjacent_matrix[ele1][ele2])
				edge_sol++;
		}
	}
	return edge_sol;
}

bool is_feasible_solution(instance ins, int fixed_k, int necessary_edge, solution sol)
{
	int sol_len = 0;
	bool flag = true;
	for (int i = 0; i < ins.num_node; i++)
	{
		if (sol.flag_arr[i])
			sol_len++;
	}
	if (sol_len != fixed_k || sol_len != sol.sol_len)
	{
		cout << "an error is check 1, in is_feasible_solution method" <<",sol_len="<<sol_len<<",fixed_k="<<fixed_k<<",sol.sol_len="<<sol.sol_len<< endl;
		flag = false;
		//getchar();
	}
	/*for (int i = 0; i < sol.sol_len; i++)
	{
		int ele1 = sol.sol_arr[i];
		for (int j = i + 1; j < sol.sol_len; j++)
		{
			int ele2 = sol.sol_arr[j];
			if (ele1 == ele2)
			{
				cout << "an error is check 2 ,in is_feasible_solution method" << ",i=" << i << ",j=" << j << ",ele=" << ele1 << endl;
				flag = false;
				//getchar();
			}
		}
	}*/
	int f = compute_penalty_f(sol, ins);
	if (f < necessary_edge)
	{
		//cout << "an error is check 3, in is_feasible_solution method" << ",f=" << f << ",necessary_edge=" << necessary_edge << endl;
		//getchar();
		flag = false;
	}
	return flag;
}

void copy_solution(solution source, solution &des, instance ins)
{
	for (int i = 0; i < ins.num_node; i++)
	{
		des.sol_arr[i] = source.sol_arr[i];
		des.flag_arr[i] = source.flag_arr[i];
	}
	des.sol_len = source.sol_len;
	des.penalty_f = source.penalty_f;
	
}

void free_memory(instance ins)
{
	delete ins.v_edge_cnt; ins.v_edge_cnt = NULL;
	for (int i = 0; i < ins.num_node; i++)
	{
		delete ins.adjacent_matrix[i]; ins.adjacent_matrix[i] = NULL;
		delete ins.v_adj_matrix[i]; ins.v_adj_matrix[i] = NULL;
	}
}

int compInc(const void *a, const void *b)
{
	return *(int *)a - *(int *)b;
}

//与greedy_create 思路一样，只是没有了fixed_k这个先验知识，并且生成的是feasible solution, 对后续tabu search计算时间可能有较大增加。
bool greedy_create_noK(solution &sol, instance ins, double gamma, int &fixed_k)
{
	double sg_den = 1;
	int *flag_in_C = new int[ins.num_node];
	memset(flag_in_C, 0, sizeof(int)*ins.num_node);
	int rand_v = rand() % ins.num_node;
	int max_deg, deg;
	int sel_v, candidate_len;
	int *candidate_v = new int[ins.num_node];
	sol.sol_len = 0;
	sol.sol_arr[sol.sol_len++] = rand_v;
	memset(sol.flag_arr, 0, sizeof(int)*ins.num_node);
	sol.flag_arr[rand_v] = 1;
	flag_in_C[rand_v] = 1;
	int edge_temp_clique = 0;
	while (sg_den >= gamma)
	{
		max_deg = -1;
		sel_v = -1;
		candidate_len = 0;
		for (int idx = 0; idx < ins.num_node; idx++)
		{
			int ele1 = idx;
			if (!flag_in_C[ele1])
			{
				deg = 0;
				for (int jdx = 0; jdx < sol.sol_len; jdx++)
				{
					int ele2 = sol.sol_arr[jdx];
					if (ins.adjacent_matrix[ele1][ele2])
						deg++;
				}
				if (deg > max_deg)
				{
					max_deg = deg;
					candidate_len = 0;
					candidate_v[candidate_len++] = ele1;
				}
				else if (deg == max_deg)
					candidate_v[candidate_len++] = ele1;
			}
		}
		sg_den = 1.0*(edge_temp_clique + max_deg) / (sol.sol_len*(sol.sol_len + 1) / 2);
		if (sg_den >= gamma)
		{
			sel_v = candidate_v[rand() % candidate_len];
			sol.sol_arr[sol.sol_len++] = sel_v;
			flag_in_C[sel_v] = 1;
			sol.flag_arr[sel_v] = 1;
			edge_temp_clique += max_deg;
		}
		else
			break;
	}
	sol.penalty_f = compute_penalty_f(sol, ins);
	fixed_k = sol.sol_len;
	//cout << "in greedy_create_noK method" << ", sol.sol_len=" << sol.sol_len << ",gamma=" << gamma << endl;
	//getchar();
	delete [] flag_in_C; 
	delete [] candidate_v;
	return true;
}

//观察种群相似性
double observe_sim(solution *pop,int pop_size)
{
	double sim = 0;
	for (int i = 0; i < pop_size; i++)
	{
		for (int j = 0; j < pop[i].sol_len; j++)
		{
			for (int k = j + 1; k < pop[i].sol_len; k++)
			{
				if (pop[i].sol_arr[j] > pop[i].sol_arr[k])
				{
					int temp = pop[i].sol_arr[j];
					pop[i].sol_arr[j] = pop[i].sol_arr[k];
					pop[i].sol_arr[k] = temp;
				}
			}
		}
	}
	//cout<<"pop_size="<<pop_size;
	//for(int i=0;i<pop_size;i++)
	//	cout<<"i="i<<",pop_len="<<pop[i].sol_len<<",penalty_f="<<pop[i].penalty_f<<endl;
	for (int i = 0; i < pop_size; i++)
	{
		for (int j = i + 1; j < pop_size; j++)
		{
			int index1 = 0;
			int index2 = 0;
			int len_same = 0;
			while (index1 < pop[i].sol_len && index2 < pop[j].sol_len)
			{
				if (pop[i].sol_arr[index1] > pop[j].sol_arr[index2])
					index2++;
				else if (pop[i].sol_arr[index1] == pop[j].sol_arr[index2])
				{
					len_same++;					
					index1++;
					index2++;
				}
				else
					index1++;	
			}
			//cout << "i=" << i << ",j=" << j << ",pop[i].len=" << pop[i].sol_len << ",pop[j].len=" << pop[j].sol_len << ",len_same=" << len_same << endl;
			sim += 1.0*(len_same * 2) / (pop[i].sol_len + pop[j].sol_len);			
		}
	}
	sim /= (pop_size*(pop_size - 1) / 2);
	return sim;
}

void perturbation(solution &sol_cur, const solution sol_best, instance ins, int fixed_k)
{
	copy_solution(sol_best, sol_cur, ins);
	int *address = new int[ins.num_node];
	while (sol_cur.sol_len < fixed_k)
	{
		int v = rand() % ins.num_node;
		while (sol_cur.flag_arr[v])
			v = rand() % ins.num_node;
		sol_cur.flag_arr[v] = 1;
		sol_cur.sol_arr[sol_cur.sol_len++] = v;
	}
	for (int i = 0; i < sol_cur.sol_len; i++)
	{
		int ele = sol_cur.sol_arr[i];
		address[ele] = i;
	}
	double pert_str = 0.10;
	int len = 0;
	while (len < pert_str*ins.num_node)
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
	sol_cur.penalty_f = compute_penalty_f(sol_cur, ins);
	delete [] address; 
}

void perturbation_greedy(solution &sol_cur, const solution sol_best, instance ins, int fixed_k)
{
	copy_solution(sol_best, sol_cur, ins);
	int *address = new int[ins.num_node];
	int *deg_matrix = new int[ins.num_node];
	memset(deg_matrix, 0, sizeof(int)*ins.num_node);
	for (int i = 0; i < ins.num_node; i++)
	{
		for (int j = 0; j < sol_cur.sol_len; j++)
		{
			int ele = sol_cur.sol_arr[j];
			if (ins.adjacent_matrix[i][ele])
				deg_matrix[i]++;
		}
	}
	while (sol_cur.sol_len < fixed_k)
	{
		int v = rand() % ins.num_node;
		while (sol_cur.flag_arr[v])
			v = rand() % ins.num_node;
		sol_cur.flag_arr[v] = 1;
		sol_cur.sol_arr[sol_cur.sol_len++] = v;
	}
	for (int i = 0; i < sol_cur.sol_len; i++)
	{
		int ele = sol_cur.sol_arr[i];
		address[ele] = i;
	}
	int *can_v = new int[ins.num_node];
	int v_len = 0;
	int len = 0;
	double pert_str = 0.10;
	while (len < pert_str*ins.num_node)
	{
		int u = sol_cur.sol_arr[rand() % sol_cur.sol_len];
		int u_address = address[u];
		int deg_max = MINVALUE;

		for (int i = 0; i < ins.num_node; i++)
		{
			if (sol_cur.flag_arr[i] == 0)
			{
				if (deg_matrix[i]>deg_max)
				{
					deg_max = deg_matrix[i];
					v_len = 0;
					can_v[v_len++] = i;
				}
				else if (deg_matrix[i] == deg_max)
				{
					can_v[v_len++] = i;
				}
			}
		}
		int v = can_v[rand() % v_len];
		sol_cur.sol_arr[u_address] = v;
		sol_cur.flag_arr[u] = 0;
		sol_cur.flag_arr[v] = 1;
		address[v] = u_address;
		update_deg_matrix(deg_matrix, ins, u, v);
		len++;
	}
	sol_cur.penalty_f = compute_penalty_f(sol_cur, ins);
	delete [] address; 
	delete [] deg_matrix; 
	delete [] can_v; 
}
