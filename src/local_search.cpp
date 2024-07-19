#include<iostream>
#include<stdlib.h>
#include<stdio.h>
#include<math.h>
#include<string.h>
#include"local_search.h"
#include"common_func.h"
#include<time.h>
using namespace std;

#define MINVALUE -99999999
#define MAXVALUE 99999999
#define MAXNUMCAN 100			//最多候选顶点对

struct t_uv{
	int u;
	int v;
};

void update_deg_matrix(int deg_matrix[], instance ins, int v_out, int v_in)
{
	for (int i = 0; i < ins.num_node; i++)
		if (ins.adjacent_matrix[i][v_out])
			deg_matrix[i]--;		
	for (int i = 0; i < ins.num_node; i++)
		if (ins.adjacent_matrix[i][v_in])
			deg_matrix[i]++;

	
}

void check_move(instance ins, solution sol)
{
	int f = compute_penalty_f(sol, ins);
	if (f != sol.penalty_f)
	{
		cout << "an error is checked 1, in check_move method" <<",f="<<f<<", sol.penalty_f="<<sol.penalty_f<< endl;
		getchar();
	}
	else
	{
		cout << "this iteration is ok, in check_move method" << ",f=" << f << ", sol.penalty_f=" << sol.penalty_f << endl;
	}
}

void perturbation(int *deg_matrix, int *address, solution &sol_cur, double per_str, instance ins)
{
	int per_len = per_str*sol_cur.sol_len;
	int len = 0;
	int *can_v = new int[ins.num_node];
	int v_len = 0;
	while (len < per_len)
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
	delete [] can_v;
}

bool tabu_search_based_swap(solution &sol_cur, solution &sol_best, instance ins, int &iter_ms, int iter_max, int tabu_length, int tabu_tenure, int necessary_edge, int fixed_k, double edge_density, int *freq)
{
	int tabu_iter = 0;
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
	int *tabu_list = new int[ins.num_node];
	memset(tabu_list, 0, sizeof(int)*ins.num_node);	
	int *address_ins = new int[ins.num_node];
	for (int i = 0; i < sol_cur.sol_len; i++)
	{
		int ele = sol_cur.sol_arr[i];
		address_ins[ele] = i;
	}
	int diff, delta, delta_tabu, uv_len, uv_len_tabu, u, v;
	int target_u[MAXNUMCAN], target_v[MAXNUMCAN];
	int target_u_tabu[MAXNUMCAN], target_v_tabu[MAXNUMCAN];

	while (tabu_iter < tabu_length && iter_ms < iter_max)
	{
		delta = MINVALUE;
		delta_tabu = MINVALUE;
		uv_len = 0;
		uv_len_tabu = 0;
		for (int i = 0; i < ins.num_node; i++)
		{
			if (sol_cur.flag_arr[i])			//S
			{
				for (int j = 0; j < ins.num_node; j++)
				{
					if (sol_cur.flag_arr[j] == 0)	//V\S
					{
						if (ins.adjacent_matrix[i][j])
							diff = deg_matrix[j] - deg_matrix[i] - 1;
						else
							diff = deg_matrix[j] - deg_matrix[i];						
						if (tabu_list[i] <= tabu_iter && tabu_list[j] <= tabu_iter)
						{
							if (diff > delta)
							{
								uv_len = 0;
								delta = diff;
								target_u[uv_len] = i;
								target_v[uv_len] = j;
								uv_len++;
							}
							else if (diff == delta && uv_len < MAXNUMCAN)
							{
								target_u[uv_len] = i;
								target_v[uv_len] = j;
								uv_len++;
							}
						}
						else
						{
							if (diff > delta_tabu)
							{
								uv_len_tabu = 0;
								delta_tabu = diff;
								target_u_tabu[uv_len_tabu] = i;
								target_v_tabu[uv_len_tabu] = j;
								uv_len_tabu++;
							}
							else if (diff == delta_tabu &&uv_len_tabu<MAXNUMCAN)
							{
								target_u_tabu[uv_len_tabu] = i;
								target_v_tabu[uv_len_tabu] = j;
								uv_len_tabu++;
							}
						}
					}
				}
			}
		}
		int le = fixed_k*(fixed_k - 1) / 2 - sol_cur.penalty_f;
		double pro_p = 1.0*(le + 2) / ins.num_node;
		if (pro_p > 0.05)
			pro_p = 0.05;
		double fa = rand() / (double)(RAND_MAX);
		//cout << "pro_p=" << pro_p << ",fa=" << fa << endl;
		u = -1, v = -1;
		if (fa >= pro_p)
		{
			if (((uv_len_tabu > 0) && (delta_tabu > delta) && ((sol_cur.penalty_f + delta_tabu) > sol_best.penalty_f))
				|| ((uv_len == 0) && (uv_len_tabu > 0)))
			{
				int x = rand() % uv_len_tabu;
				u = target_u_tabu[x];
				v = target_v_tabu[x];
				sol_cur.penalty_f += delta_tabu;
			}
			else if (uv_len>0)
			{
				int x = rand() % uv_len;
				u = target_u[x];
				v = target_v[x];
				sol_cur.penalty_f += delta;
			}
		}
		if (fa<pro_p)
		{			
			u = sol_cur.sol_arr[rand() % sol_cur.sol_len];			
			int flag_v = 1;
			while (flag_v)
			{
				v = rand() % ins.num_node;				
				if (sol_cur.flag_arr[v] == 0 )
					flag_v = 0;
				//cout << "v=" << v << endl;
			}
			sol_cur.penalty_f += deg_matrix[v] - deg_matrix[u];
			if (ins.adjacent_matrix[u][v])
				sol_cur.penalty_f += -1;
		}
		//cout << "u=" << u << ",v=" << v << endl;
		//交换u,v
		int u_address = address_ins[u];
		sol_cur.sol_arr[u_address] = v;
		sol_cur.flag_arr[u] = 0;
		sol_cur.flag_arr[v] = 1;
		address_ins[v] = u_address;
		freq[u] += 1;
		freq[v] += 1;
		int freq_len = 0;
		for (int i = 0; i < ins.num_node; i++)
		{
			if (freq[i] >= fixed_k)
				freq_len++;
			else
				break;
		}
		if (freq_len == ins.num_node)
			memset(freq, 0, sizeof(int)*ins.num_node);
		//更新禁忌表
		if (le > 10)
			le = 10;
		int lc = fixed_k / 40;
		if (lc<6)
			lc = 6;
		
		//tabu_list[u] = tabu_iter + le + rand() % lc;
		//tabu_list[v] = tabu_iter + 0.6*le + rand() % (int(0.6*lc));
		tabu_list[u] = tabu_iter + tabu_tenure;
		tabu_list[v] = tabu_iter + tabu_tenure;
		update_deg_matrix(deg_matrix, ins, u, v);
		if (sol_cur.penalty_f > sol_best.penalty_f)
			copy_solution(sol_cur, sol_best, ins);		
		if (sol_cur.penalty_f >= necessary_edge)					
			return true;		
		//cout << "tabu_iter=" << tabu_iter << ",iter_ms=" << iter_ms <<",sol_cur.penalty_f="<<sol_cur.penalty_f<<",sol_best.penalty_f="<<sol_best.penalty_f<<",fixed_k="<<fixed_k<< endl;
		//check_move(ins, sol_cur);
		tabu_iter++;
		iter_ms++;
	}

	delete [] deg_matrix; 
	delete [] tabu_list; 
	delete [] address_ins; 
	return false;
}


void free_memory(int *deg_matrix, int *tabu_list, int *address_ins, int *subset_a, int *subset_b, t_uv *t_uv_ins,
	int *target_u_tabu, int *target_v_tabu, int *target_u, int *target_v)
{
	delete [] deg_matrix; 
	delete [] tabu_list; 
	delete [] address_ins; 
	delete [] subset_a; 
	delete [] subset_b; 
	delete [] t_uv_ins; 
	delete [] target_u_tabu;
	delete [] target_v_tabu; 
	delete [] target_v;
	delete [] target_u;
}

bool tabu_search_constrained_swap(solution &sol_cur, solution &sol_best, instance ins, int &iter_ms, int iter_max, int tabu_length, int tabu_tenure, int necessary_edge, int fixed_k, double edge_density, int *freq)
{
	
	int tabu_iter = 0;
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
	int *tabu_list = new int[ins.num_node];
	memset(tabu_list, 0, sizeof(int)*ins.num_node);
	int *address_ins = new int[ins.num_node];
	for (int i = 0; i < sol_cur.sol_len; i++)
	{
		int ele = sol_cur.sol_arr[i];
		address_ins[ele] = i;
	}
	int diff, u, v;
	//int target_u_tabu, target_v_tabu;
	int min_in_s, max_out_s;
	int *subset_a = new int[ins.num_node];
	int *subset_b = new int[ins.num_node];
	t_uv *t_uv_ins = new t_uv[ins.num_node];
	int t_uv_len;
	int subset_a_len, subset_b_len;

	int *target_u_tabu = new int[ins.num_node];
	int *target_v_tabu = new int[ins.num_node];
	int target_u_tabu_len = 0, target_v_tabu_len = 0;
	int *target_u = new int[ins.num_node];
	int *target_v = new int[ins.num_node];
	int target_u_len = 0, target_v_len = 0;
	int iter_no_imp = 0;
	int pert_length = 1000;
	double pert_str = 0.2;
	while (tabu_iter < tabu_length)
	{
		int deg_min = MAXVALUE, deg_max = MINVALUE;
		for (int i = 0; i < ins.num_node; i++)
		{
			if (sol_cur.flag_arr[i] && deg_matrix[i] < deg_min && tabu_list[i] <= tabu_iter)
				deg_min = deg_matrix[i];
			if (sol_cur.flag_arr[i] == 0 && deg_matrix[i] > deg_max && tabu_list[i] <= tabu_iter)
				deg_max = deg_matrix[i];
		}
		min_in_s = deg_min;
		max_out_s = deg_max;
		subset_a_len = 0;
		subset_b_len = 0;
		diff = MINVALUE;
		for (int i = 0; i < ins.num_node; i++)
		{
			if (sol_cur.flag_arr[i] && tabu_list[i] <= tabu_iter && deg_matrix[i] == min_in_s)
				subset_a[subset_a_len++] = i;
			if (sol_cur.flag_arr[i] == 0 && tabu_list[i] <= tabu_iter && deg_matrix[i] == max_out_s)
				subset_b[subset_b_len++] = i;
		}
		t_uv_len = 0;
		for (int i = 0; i < subset_a_len; i++)
		{
			int u = subset_a[i];
			for (int j = 0; j < subset_b_len; j++)
			{
				int v = subset_b[j];
				if (ins.adjacent_matrix[u][v] == 0 && t_uv_len < ins.num_node)
				{
					t_uv_ins[t_uv_len].u = u;
					t_uv_ins[t_uv_len].v = v;
					t_uv_len++;
				}
			}
		}
		int k1 = -1, k2 = -1;
		int le = fixed_k*(fixed_k - 1) / 2 - sol_cur.penalty_f;
		double pro_p = 1.0*(le + 2) / ins.num_node;
		if (pro_p > 0.05)
			pro_p = 0.05;
		double fa = rand() / (double)(RAND_MAX);
		if (fa >= pro_p)
		{
			if (t_uv_len == 0 && subset_a_len > 0 && subset_b_len > 0)		//(u,v), u \in A, v \in B,全部有边相连，即T为空集
			{
				int rand_a = rand() % subset_a_len;
				k1 = subset_a[rand_a];
				int rand_b = rand() % subset_b_len;
				k2 = subset_b[rand_b];
				diff = deg_matrix[k2] - deg_matrix[k1] - 1;

			}
			else if (t_uv_len>0)      //从t中挑选
			{
				int rand_id = rand() % t_uv_len;
				k1 = t_uv_ins[rand_id].u;
				k2 = t_uv_ins[rand_id].v;
				diff = deg_matrix[k2] - deg_matrix[k1];
			}
		}
		if (fa<pro_p)
		{
			k1 = sol_cur.sol_arr[rand() % sol_cur.sol_len];
			int flag_v = 1;
			while (flag_v)
			{
				k2 = rand() % ins.num_node;
				if (sol_cur.flag_arr[k2] == 0)
					flag_v = 0;
			}	
			diff = deg_matrix[k2] - deg_matrix[k1];
			if (ins.adjacent_matrix[k1][k2])
				diff += -1;
		}
	
		
		//aspiration criterion， 这里有较大问题， 要注意
		deg_min = MAXVALUE;
		deg_max = MINVALUE;
		for (int i = 0; i < ins.num_node; i++)
		{
			if (sol_cur.flag_arr[i] && deg_matrix[i] < deg_min && tabu_list[i] > tabu_iter)
				deg_min = deg_matrix[i];
			if (sol_cur.flag_arr[i] == 0 && deg_matrix[i] > deg_max && tabu_list[i] > tabu_iter)
				deg_max = deg_matrix[i];
		}
		min_in_s = deg_min;
		max_out_s = deg_max;
		target_u_tabu_len = 0, target_v_tabu_len = 0;
		target_u_len = 0, target_v_len = 0;
		for (int i = 0; i < ins.num_node; i++)
		{
			if (sol_cur.flag_arr[i] && tabu_list[i] > tabu_iter && deg_matrix[i] == min_in_s)
				target_u_tabu[target_u_tabu_len++] = i;
			if (sol_cur.flag_arr[i] == 0 && tabu_list[i] > tabu_iter && deg_matrix[i] == max_out_s)
				target_v_tabu[target_v_tabu_len++] = i;
			if (sol_cur.flag_arr[i] && deg_matrix[i] == min_in_s)
				target_u[target_u_len++] = i;
			if (sol_cur.flag_arr[i] == 0 && deg_matrix[i] == max_out_s)
				target_v[target_v_len++] = i;
		}
		int merge_con1 = 1;
		if (subset_a_len == 0 || subset_b_len == 0)
			merge_con1 = 0;
		int merge_con2 = 1;
		if (target_u_tabu_len == 0 && target_v_tabu_len == 0)
			merge_con2 = 0;
		int diff_asp = MINVALUE;
		if (merge_con2 == 1)
		{
			int rand_u, rand_v;
			if (target_u_tabu_len > 0)
			{
				rand_u = target_u_tabu[rand() % target_u_tabu_len];
				rand_v = target_v[rand() % target_v_len];
			}
			else if (target_v_tabu_len > 0)
			{
				rand_v = target_v_tabu[rand() % target_v_tabu_len];
				rand_u = target_u[rand() % target_u_len];
			}
			if (ins.adjacent_matrix[rand_u][rand_v])
				diff_asp = deg_matrix[rand_v] - deg_matrix[rand_u] - 1;
			else
				diff_asp = deg_matrix[rand_v] - deg_matrix[rand_u];
			if ((merge_con1 == 0) || ((diff_asp > diff) && (sol_cur.penalty_f + diff_asp > sol_best.penalty_f)))
			{
				k1 = rand_u;
				k2 = rand_v;
				diff = diff_asp;
				//cout << "the aspiration is encourted" <<",tabu_iter="<<tabu_iter<< endl;
			}
		}

		u = k1;
		v = k2;
		sol_cur.penalty_f += diff;	

		//cout << "u=" << u << ",v=" << v << endl;
		//交换u,v
		int u_address = address_ins[u];
		sol_cur.sol_arr[u_address] = v;
		sol_cur.flag_arr[u] = 0;
		sol_cur.flag_arr[v] = 1;
		address_ins[v] = u_address;
		freq[u] += 1;
		freq[v] += 1;
		int freq_len = 0;
		for (int i = 0; i < ins.num_node; i++)
		{
			if (freq[i] >= fixed_k)
				freq_len++;
			else
				break;
		}
		if (freq_len == ins.num_node)
			memset(freq, 0, sizeof(int)*ins.num_node);
		//更新禁忌表
		if (le > 10)
			le = 10;
		int lc = fixed_k / 40;
		if (lc<6)
			lc = 6;

		//tabu_list[u] = tabu_iter + le + rand() % lc;
		//tabu_list[v] = tabu_iter + 0.6*le + rand() % (int(0.6*lc));
		tabu_list[u] = tabu_iter + tabu_tenure;
		tabu_list[v] = tabu_iter + tabu_tenure;
		update_deg_matrix(deg_matrix, ins, u, v);
		if (sol_cur.penalty_f > sol_best.penalty_f)
		{
			copy_solution(sol_cur, sol_best, ins);
			iter_no_imp = 0;
		}
		else
			iter_no_imp++;
		if (sol_cur.penalty_f >= necessary_edge)
		{			
			free_memory(deg_matrix, tabu_list, address_ins, subset_a, subset_b, t_uv_ins, target_u_tabu, target_v_tabu, target_u, target_v);
			return true;
		}
		//if (iter_no_imp == pert_length)
		//	perturbation(deg_matrix, address_ins, sol_cur, pert_str, ins);
		//cout << "tabu_iter=" << tabu_iter << ",iter_ms=" << iter_ms << ",sol_cur.penalty_f=" << sol_cur.penalty_f << ",sol_best.penalty_f=" << sol_best.penalty_f <<",fixed_k="<<fixed_k<< endl;
		//check_move(ins, sol_cur);
		tabu_iter++;
		iter_ms++;
	}

	free_memory(deg_matrix, tabu_list, address_ins, subset_a, subset_b, t_uv_ins, target_u_tabu, target_v_tabu, target_u, target_v);
	return false;
}

bool multi_start(solution &sol_cur, solution &sol_best, solution &sol_ms_best, instance ins, int iter_max, int tabu_length, int tabu_tenure, int fixed_k, int &best_k, double gamma, double edge_density, double &end_time)
{
	int iter_ms = 0;
	int edge_k_clique = fixed_k*(fixed_k - 1) / 2;
	int necessary_edge = ceil(gamma*edge_k_clique);
	int *freq = new int[ins.num_node];
	memset(freq, 0, sizeof(int)*ins.num_node);

	sol_cur.greedy_create(sol_cur, ins, fixed_k);
	copy_solution(sol_cur, sol_best, ins);
	while (iter_ms < iter_max)
	{
		bool feasible = is_feasible_solution(ins, fixed_k, necessary_edge, sol_best);
		if (feasible)
		{
			if (fixed_k > best_k)
				best_k = fixed_k;
			//cout << "succeed: fixed_k is solved" << ",fixed_k=" << fixed_k << ",best_k=" << best_k  << "#################"<<endl;
			copy_solution(sol_best, sol_ms_best, ins);
			end_time = clock();
			return true;
		}
		//bool found = tabu_search_based_swap(sol_cur, sol_best, ins, iter_ms, iter_max, tabu_length, tabu_tenure, necessary_edge, fixed_k, edge_density, freq);
		bool found = tabu_search_constrained_swap(sol_cur, sol_best, ins, iter_ms, iter_max, tabu_length, tabu_tenure, necessary_edge, fixed_k, edge_density, freq);
		if (found == true)
		{

		}
		else
		{
			sol_cur.greedy_create(sol_cur, ins, fixed_k);
			//sol_cur.frequency_based_initialize(sol_cur, ins, fixed_k, freq);
		}
	}
	delete [] freq; 
	return false;
}

bool local_search(solution &sol_cur, solution &sol_best, solution &sol_ms_best, instance ins, int &iter, int iter_max, int tabu_length, int tabu_tenure, int fixed_k, int &best_k, double gamma, double edge_density, double &end_time)
{	
	//int iter_ms = 0;
	int edge_k_clique = fixed_k*(fixed_k - 1) / 2;
	int necessary_edge = ceil(gamma*edge_k_clique);
	int *freq = new int[ins.num_node];
	bool found = false;
	memset(freq, 0, sizeof(int)*ins.num_node);		
	copy_solution(sol_cur, sol_best, ins);			
	bool feasible = is_feasible_solution(ins, fixed_k, necessary_edge, sol_best);	//sol_best不一定合法，不过没有关系，tabu search 考虑的就包括不合法解
	if (!feasible)
	{
		found = tabu_search_constrained_swap(sol_cur, sol_best, ins, iter, iter_max, tabu_length, tabu_tenure, necessary_edge, fixed_k, edge_density, freq);		
		//found = tabu_search_based_swap(sol_cur, sol_best, ins, iter, iter_max, tabu_length, tabu_tenure, necessary_edge, fixed_k, edge_density, freq);		
	}	
	if (found || feasible)
	{
		if (fixed_k > best_k)
			best_k = fixed_k;
		//cout << "succeed: fixed_k is solved" << ",fixed_k=" << fixed_k << ",best_k=" << best_k  << "#################"<<endl;
		copy_solution(sol_best, sol_ms_best, ins);
		end_time = clock();
		delete freq; freq = NULL;
		return true;
	}
	delete [] freq; 
	return false;
}
