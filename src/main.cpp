#include<iostream>
#include<fstream>
#include<time.h>
#include<stdlib.h>
#include"instance.h"
#include"solution.h"
#include"local_search.h"
#include"common_func.h"
#include"memetic.h"
using namespace std;

#define IS_LINUX
#define MEMETIC
#define MEMETIC_OPPO11
#define MULTI_START11
#define ITERATED_LOCAL11
#define POP_SIZE 10



#ifdef IS_LINUX
int main(int argc, char **argv)
{
	srand(unsigned(time(NULL)));
	instance ins;
	if (argc <4)
	{
		cout << "usage: TS.exe, input_file, gamma, time_limit ";
		exit(-1);
	}
	char *file_in = argv[1];	
	double gamma = atof(argv[2]);
	double time_limit=atof(argv[3]);
	time_limit/=0.98;
	int fixed_k=-1;
	bool open_success = ins.read_instance_file(file_in, ins);
	if (!open_success)
	{
		cout << "open file_in error" << endl;
		exit(-1);
	}
	double edge_density = -1;				//未使用
	int best_k = -1;						//所求目标值
	int global_best_k = -1;
	int iter_max = 100000;
	int tabu_length = 10000;				//论文中是5000， 实际上这个参数大一点结果应该会更好
	int tabu_tenure = 10;
	double sim;
	solution *pop = new solution[POP_SIZE];
	solution sol_cur;
	solution sol_best;
	solution global_best;
	solution sol_ms_best;
	solution sol_cur_oppo;
	solution sol_cur_temp;
	sol_cur.sol_arr = new int[ins.num_node], sol_cur.flag_arr = new int[ins.num_node];
	sol_cur_oppo.sol_arr = new int[ins.num_node], sol_cur_oppo.flag_arr = new int[ins.num_node];
	sol_cur_temp.sol_arr = new int[ins.num_node], sol_cur_temp.flag_arr = new int[ins.num_node];
	sol_best.sol_arr = new int[ins.num_node], sol_best.flag_arr = new int[ins.num_node];
	global_best.sol_arr = new int[ins.num_node], global_best.flag_arr = new int[ins.num_node];
	sol_ms_best.sol_arr = new int[ins.num_node], sol_ms_best.flag_arr = new int[ins.num_node];
	for (int i = 0; i < POP_SIZE; i++)
	{
		pop[i].flag_arr = new int[ins.num_node];
		pop[i].sol_arr = new int[ins.num_node];
	}
	int runs = 1;
	//int runs = 100;
	double *run_time_size = new double[runs];
	double average_time, average_size;
	double start_time, run_time, end_time;
	ofstream resultsFile;
	ofstream valuesFile;
	resultsFile.open("/home/fudama/wu/zhouqing/MQCP_1130/results/oppo_bb_crossover_stat.txt", ofstream::app);
	valuesFile.open("/home/fudama/wu/zhouqing/MQCP_1130/results/oppo_bb_crossover_bestsol.txt", ofstream::app);	
	//resultsFile.open("E:\\MQCP_output_sum\\OPPO_penalty_cn_2_stat_infile.txt", ofstream::app);
	//valuesFile.open("E:\\MQCP_output_sum\\OPPO_penalty_cn_2_best_sol_infile.txt", ofstream::app);
	if (!resultsFile)
	{
		cout << "resultsFile open error" << endl;
		exit(-1);
	}
	if (!valuesFile)
	{
		cout << "valuesFile open error" << endl;
		exit(-1);
	}
	average_size = 0;
	average_time = 0;
	global_best.sol_len = -1;
	int step_forward = 1;
	int step_backward = 1;
	for (int i = 0; i < runs; i++)
	{
		start_time = clock();		
		bool res_val;		
		fixed_k = -1;		//不使用best_known 这个先验知识，保证公平。
		best_k = -1, sol_ms_best.sol_len = -1;		
		int iter = 0;
#ifdef MEMETIC
		greedy_create_noK(sol_cur, ins, gamma, fixed_k);
		if (fixed_k > best_k)
		{
			best_k = fixed_k;
			copy_solution(sol_cur, sol_ms_best, ins);
			end_time = clock();
		}
		res_val = true;
		//cout << "fixed_k=" << fixed_k << endl;
		while (res_val && fixed_k <= ins.num_node)
		{
			iter = 0;
			res_val=false;	
			bool succ = initial_population(sol_cur, sol_best, sol_ms_best, ins, iter, iter_max, tabu_length, tabu_tenure, fixed_k, best_k, gamma, edge_density, end_time, pop, POP_SIZE);		
			if(succ)
			{	
			    res_val = true;
			    fixed_k += step_forward;
			    continue;
			}
			//while (iter < iter_max)
			while(1.0*(clock()-start_time)/CLOCKS_PER_SEC < time_limit)
			{
				common_shared_cross_over(sol_cur, pop, POP_SIZE, ins, fixed_k);
				res_val = local_search(sol_cur, sol_best, sol_ms_best, ins, iter, iter_max, tabu_length, tabu_tenure, fixed_k, best_k, gamma, edge_density, end_time);
				update_population(pop, sol_best, ins, POP_SIZE);
				if (res_val)
				{
					fixed_k += step_forward;
					break;
				}					
			}
		}
#endif
#ifdef MEMETIC_OPPO
		greedy_create_noK(sol_cur, ins, gamma, fixed_k);
		if (fixed_k > best_k)
		{
			best_k = fixed_k;
			copy_solution(sol_cur, sol_ms_best, ins);
			end_time = clock();
		}
		res_val = true;
		while (res_val && fixed_k <= ins.num_node)
		{
			iter = 0;
			res_val=false;
			bool succ = oppo_initial_population(sol_cur, sol_best, sol_ms_best, ins, iter, iter_max, tabu_length, tabu_tenure, fixed_k, best_k, gamma, edge_density, end_time, pop, POP_SIZE);		
			if(succ)
			{	
			    res_val = true;
			    fixed_k += step_forward;
			    continue;
			}
			//while (iter < iter_max)
			while(1.0*(clock()-start_time)/CLOCKS_PER_SEC < time_limit)
			{
				common_shared_cross_over(sol_cur, pop, POP_SIZE, ins, fixed_k);
				//prob_same_cross_over(sol_cur, pop, POP_SIZE, ins, fixed_k);
				copy_solution(sol_cur, sol_cur_temp, ins);
				res_val = local_search(sol_cur, sol_best, sol_ms_best, ins, iter, iter_max, tabu_length, tabu_tenure, fixed_k, best_k, gamma, edge_density, end_time);
				update_population(pop, sol_best, ins, POP_SIZE);
				if (res_val)
				{
					fixed_k += step_forward;
					break;
				}
				generate_opposition_solution(sol_cur_temp, sol_cur_oppo, ins, fixed_k);
				res_val = local_search(sol_cur_oppo, sol_best, sol_ms_best, ins, iter, iter_max, tabu_length, tabu_tenure, fixed_k, best_k, gamma, edge_density, end_time);
				update_population(pop, sol_best, ins, POP_SIZE);
				if (res_val)
				{
					fixed_k += step_forward;
					break;
				}
			}
		}
#endif
#ifdef MULTI_START	
		
		greedy_create_noK(sol_cur, ins, gamma, fixed_k);	//得到fixed_k值				
		if (fixed_k > best_k)
		{
			best_k = fixed_k;
			copy_solution(sol_cur, sol_ms_best, ins);
			end_time = clock();
		}
		res_val = true;
		while(res_val && fixed_k <= ins.num_node)
		{
			sol_cur.random_create(sol_cur, ins, fixed_k);
			iter = 0;
			res_val=false;
			//while (iter < iter_max)
			while(1.0*(clock()-start_time)/CLOCKS_PER_SEC < time_limit)
			{					
				res_val = local_search(sol_cur, sol_best, sol_ms_best, ins, iter, iter_max, tabu_length, tabu_tenure, fixed_k, best_k, gamma, edge_density, end_time);
				if (res_val)
				{
					fixed_k += step_forward;
					break;
				}					
				else
					sol_cur.random_create(sol_cur, ins, fixed_k);
				//cout<<"iter="<<iter<<",fixed_k="<<fixed_k<<endl;
			}
		}
#endif
#ifdef ITERATED_LOCAL
		//迭代式local search效果最差
		greedy_create_noK(sol_cur, ins, gamma, fixed_k);		//得到fixed_k值				
		res_val = true;
		while (res_val && fixed_k <= ins.num_node)
		{
			sol_cur.random_create(sol_cur, ins, fixed_k);
			copy_solution(sol_cur, sol_cur_temp, ins);
			iter = 0;
			//while (iter < iter_max)
			while(1.0*(clock()-start_time)/CLOCKS_PER_SEC < time_limit)
			{					
				perturbation(sol_cur, sol_cur_temp, ins, fixed_k);
				//perturbation_greedy(sol_cur, sol_cur_temp, ins, fixed_k);
				res_val = local_search(sol_cur, sol_best, sol_ms_best, ins, iter, iter_max, tabu_length, tabu_tenure, fixed_k, best_k, gamma, edge_density, end_time);
				if (sol_best.penalty_f > sol_cur_temp.penalty_f)
					copy_solution(sol_best, sol_cur_temp, ins);
				if (res_val)
				{
					fixed_k += step_forward;
					break;
				}
				cout<<"iter="<<iter<<",fixed_k="<<fixed_k<<endl;
			}
		}
#endif
		sim = observe_sim(pop, POP_SIZE);
		cout << "instance_name=" << file_in << ", i="<< i <<", fixed_k=" << fixed_k << ", best_k=" << best_k <<", sim="<<sim<< endl;
		if (sol_ms_best.sol_len>global_best.sol_len)
			copy_solution(sol_ms_best, global_best, ins);
		run_time = (end_time - start_time) / CLOCKS_PER_SEC;
		run_time_size[i] = run_time;
		average_size += best_k;
		average_time += run_time;
	}
	average_size /= runs;
	average_time /= runs;		
	cout << "instance_name=" << file_in << ",gamma=" << gamma << ",best.size=" <<
		global_best.sol_len << ",avg_size=" << average_size << ",avg_time=" << average_time << endl;
	resultsFile << "instance_name=" << file_in << ",gamma=" << gamma << ",best_size=" <<
		global_best.sol_len << ",avg_size=" << average_size << ",avg_time=" << average_time << endl;
	valuesFile << "instance_name=" << file_in << "---------------------------------cut line" << endl;
	for (int i = 0; i < global_best.sol_len; i++)
		valuesFile << global_best.sol_arr[i] << " ";
	valuesFile << endl;
	delete [] run_time_size; 
	delete [] sol_cur.sol_arr; 
	delete [] sol_cur.flag_arr;
	delete [] sol_best.sol_arr;
	delete [] sol_best.flag_arr; 
	delete [] global_best.sol_arr; 
	delete [] global_best.flag_arr;
	delete [] sol_ms_best.flag_arr; 
	delete [] sol_ms_best.sol_arr; 
	delete [] sol_cur_oppo.sol_arr;
	delete [] sol_cur_oppo.flag_arr;
	delete [] sol_cur_temp.sol_arr; 
	delete [] sol_cur_temp.flag_arr;
	for (int i = 0; i < POP_SIZE; i++)
	{
		delete [] pop[i].flag_arr; 
		delete [] pop[i].sol_arr; 
	}
	delete[]pop;
	delete ins.v_edge_cnt; 
	for(int i=0;i < ins.num_node; i++)
	{
		delete [] ins.v_adj_matrix[i];
		delete [] ins.adjacent_matrix[i];
	}
	delete []ins.v_adj_matrix; 
	delete []ins.adjacent_matrix;	
	resultsFile.close();
	valuesFile.close();	
	cout << "finished" << endl;
	getchar();
}
#endif
