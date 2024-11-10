#include<iostream>
#include<fstream>
#include<string>
#include<string.h>
#include<assert.h>
#include"instance.h"
using namespace std;
#define DEBUG_READ_FILE

bool instance::read_instance_file(char *in_file, instance &ins)
{
	ifstream infile(in_file);
	char line[1024];
	char tmps1[1024];
	char tmps2[1024];
	if (!infile)
	{
		fprintf(stderr, "Can not find file %s\n", in_file);
		return 0;
	}
	infile.getline(line, 1024);
	while (line[0] != 'p')	
		infile.getline(line, 1024);
	sscanf(line, "%s %s %d %d", tmps1, tmps2, &ins.num_node, &ins.num_edge);
	infile.getline(line, 1024);		//读取p edge下面的一行。
	ins.adjacent_matrix = new int*[ins.num_node];
	for (int i = 0; i < ins.num_node; i++)
	{
		ins.adjacent_matrix[i] = new int[ins.num_node];
		memset(ins.adjacent_matrix[i], 0, sizeof(int) * ins.num_node);
	}
	int ecnt = 0;
	ins.v_edge_cnt = new int[ins.num_node];
	memset(ins.v_edge_cnt, 0, sizeof(int) * ins.num_node);
	while (infile.getline(line, 1024))
	{
		int v1, v2;
		if (strlen(line) == 0)
			continue;
		if (line[0] != 'e')
			fprintf(stderr, "ERROR in line %d\n", ecnt + 1);
		sscanf(line, "%s %d %d", tmps1, &v1, &v2);
		v1--, v2--;
		ins.adjacent_matrix[v1][v2] = 1;
		ins.adjacent_matrix[v2][v1] = 1; 
		if (v1 != v2)
		{
			ins.v_edge_cnt[v1]++;
			ins.v_edge_cnt[v2]++;
		}
		else
			ins.v_edge_cnt[v1]++;
		ecnt++;
	}
	//cout << "num_edge=" << ins.num_edge << ",ecnt=" << ecnt << endl;
	assert(ins.num_edge == ecnt);	
	ins.v_adj_matrix = new int*[ins.num_node];

	for (int i = 0; i < ins.num_node; i++){
		int adj_cnt = 0;
		ins.v_adj_matrix[i] = new int[ins.v_edge_cnt[i]];
		
		for (int j = 0; j < ins.num_node; j++){
			if (ins.adjacent_matrix[i][j] == 1)
				ins.v_adj_matrix[i][adj_cnt++] = j;
		}
		//cout << "i=" << i << ",v_edge_cnt=" << ins.v_edge_cnt[i] << ",adj_cnt=" << adj_cnt << endl;
		if (ins.v_edge_cnt[i] != adj_cnt)
		{
			cout << "i=" << i << ",v_edge_cnt="<<ins.v_edge_cnt[i]<<",adj_cnt="<<adj_cnt<< endl;
		}
		assert(ins.v_edge_cnt[i] == adj_cnt);		
	}
	cout << "instance_name=" << in_file << ",node=" << ins.num_node << ",edge=" << ins.num_edge << endl;
	return true;
}

void instance::print_org_graph(instance ins){
	cout << "node="<<ins.num_node << ", " << "edge="<<ins.num_edge << " " << endl;
	for (int v1 = 0; v1 < ins.num_node; v1++)
	{
		for (int v2 = 0; v2 < ins.num_node; v2++)
		{
			if (ins.adjacent_matrix[v1][v2] == 1 && v1>v2)
				cout << "v1=" << v1+1 << ",v2=" << v2+1 << endl;
		}
		//cout << endl;
	}
}


