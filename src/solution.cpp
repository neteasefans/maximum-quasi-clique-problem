#include"solution.h"
#include"instance.h"
#include"common_func.h"
#include<iostream>
#include<stdlib.h>
#include<memory.h>
#include<stdio.h>
using namespace std;
#define MAXVALUE 99999999

bool solution::greedy_create(solution &sol, instance ins, int fixed_k)
{
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
	while (sol.sol_len < fixed_k)
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
		if (candidate_len == 0)
		{
			cout << "error is checked 1, in greedy_create method " << endl;
			getchar();
		}
		else
		{
			sel_v = candidate_v[rand() % candidate_len];
			sol.sol_arr[sol.sol_len++] = sel_v;
			flag_in_C[sel_v] = 1;
			sol.flag_arr[sel_v] = 1;
		}
	}
	sol.penalty_f = compute_penalty_f(sol, ins);
	delete [] flag_in_C;
	delete [] candidate_v; 
	return true;
}

bool solution::random_create(solution &sol, instance ins, int fixed_k)
{
	sol.sol_len = 0;
	memset(sol.flag_arr, 0, sizeof(int)*ins.num_node);
	while (sol.sol_len < fixed_k)
	{
		int rand_v = rand() % ins.num_node;
		while (sol.flag_arr[rand_v])
			rand_v = rand() % ins.num_node;
		sol.sol_arr[sol.sol_len++] = rand_v;
		sol.flag_arr[rand_v] = 1;
	}
	sol.penalty_f = compute_penalty_f(sol, ins);
	return true;
}

bool solution::frequency_based_initialize(solution &sol, instance ins, int fixed_k, int *freq)
{
	int min_freq = MAXVALUE;
	int seed_v;
	for (int i = 0; i < ins.num_node; i++)
	{
		if (freq[i] < min_freq)
		{
			min_freq = freq[i];
			seed_v = i;
		}
	}

	int *flag_in_C = new int[ins.num_node];
	memset(flag_in_C, 0, sizeof(int)*ins.num_node);
	int max_deg, deg;
	int sel_v, candidate_len, candidate_len2;
	int *candidate_v = new int[ins.num_node];
	int *candidate_v2 = new int[ins.num_node];
	sol.sol_len = 0;
	sol.sol_arr[sol.sol_len++] = seed_v;
	memset(sol.flag_arr, 0, sizeof(int)*ins.num_node);
	sol.flag_arr[seed_v] = 1;
	flag_in_C[seed_v] = 1;

	while (sol.sol_len < fixed_k)
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
		if (candidate_len == 0)
		{
			cout << "error is checked 1, in frequency_based_initialize method " << endl;
			getchar();
		}
		else
		{
			min_freq = MAXVALUE;
			for (int i = 0; i < candidate_len; i++)
			{
				int v = candidate_v[i];
				if (freq[v] < min_freq)
				{
					min_freq = freq[v];
					candidate_len2 = 0;
					candidate_v2[candidate_len2++] = v;
				}
				else if (freq[v] == min_freq)
					candidate_v2[candidate_len2++] = v;
			}
			sel_v = candidate_v2[rand() % candidate_len2];
			sol.sol_arr[sol.sol_len++] = sel_v;
			flag_in_C[sel_v] = 1;
			sol.flag_arr[sel_v] = 1;
		}
	}
	sol.penalty_f = compute_penalty_f(sol, ins);
	delete [] flag_in_C;
	delete [] candidate_v; 
	delete [] candidate_v2; 
	return true;
}

