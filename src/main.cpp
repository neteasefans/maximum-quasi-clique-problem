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

int main(int argc, char *argv[])
{
	if (argc < 7)
	{
		cout << "TS.exe usage: ./input_file gamma time seed ./output_stat_file ./output_sol_file";
		cout << "(where input_file is the instance name, gamma is the required density, time is the cutoff time, seed is the random seed, such as 1, 2, ...,\
			output_stat_file is a file used to store the running information, output_sol_file stores the solution information)" << endl;
		exit(-1);
	}
	instance ins;
	char *file_in = argv[1];	
	double gamma = atof(argv[2]);
	double time_limit = atof(argv[3]);
	time_limit /= 0.98;
	int seed = atoi(argv[4]);
	srand(seed);
	char *output_stat_file = argv[5];
	char *output_sol_file = argv[6];
	
	//the following are used parameters
	int pop_size = 10;							//the population size
	int tabu_length = 5000; 					//the search depth of tabu search
	int tabu_tenure = 10;						//the tabu tenure

	
	int fixed_k=-1;
	bool open_success = ins.read_instance_file(file_in, ins);
	if (!open_success)
	{
		cout << "open file_in error" << endl;
		exit(-1);
	}
	int best_k = -1;							//所求目标值
	int runs = 1;

	double sim;
	solution *pop = new solution[pop_size];
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
	for (int i = 0; i < pop_size; i++)
	{
		pop[i].flag_arr = new int[ins.num_node];
		pop[i].sol_arr = new int[ins.num_node];
	}
	//int runs = 100;
	double *run_time_size = new double[runs];
	double average_time, average_size;
	double start_time, run_time, end_time;
	ofstream resultsFile;
	ofstream valuesFile;
	resultsFile.open(output_stat_file, ofstream::app);
	valuesFile.open(output_sol_file, ofstream::app);

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
	
	for (int i = 0; i < runs; i++)
	{
		start_time = clock();		
		bool res_val;		
		fixed_k = -1;											//不使用best_known 这个先验知识，保证公平。
		best_k = -1, sol_ms_best.sol_len = -1;		
		int iter = 0;

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
			bool succ = oppo_initial_population(sol_cur, sol_best, sol_ms_best, ins, iter, tabu_length, tabu_tenure, fixed_k, best_k, gamma, end_time, pop, pop_size);		
			if(succ)
			{	
			    res_val = true;
			    fixed_k += step_forward;
			    continue;
			}
			while(1.0*(clock()-start_time)/CLOCKS_PER_SEC < time_limit)
			{
				common_shared_cross_over(sol_cur, pop, pop_size, ins, fixed_k);
				copy_solution(sol_cur, sol_cur_temp, ins);
				res_val = local_search(sol_cur, sol_best, sol_ms_best, ins, iter, tabu_length, tabu_tenure, fixed_k, best_k, gamma, end_time);
				update_population(pop, sol_best, ins, pop_size);
				if (res_val)
				{
					fixed_k += step_forward;
					break;
				}
				generate_opposition_solution(sol_cur_temp, sol_cur_oppo, ins, fixed_k);
				res_val = local_search(sol_cur_oppo, sol_best, sol_ms_best, ins, iter, tabu_length, tabu_tenure, fixed_k, best_k, gamma, end_time);
				update_population(pop, sol_best, ins, pop_size);
				if (res_val)
				{
					fixed_k += step_forward;
					break;
				}
			}
		}

		sim = observe_sim(pop, pop_size);
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
	valuesFile << "instance_name=" << file_in << ",best_k = " << global_best.sol_len << "---------------------cut line" << endl;
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
	for (int i = 0; i < pop_size; i++)
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
}
