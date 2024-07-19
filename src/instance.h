#ifndef _INSTANCE_H
#define _INSTANCE_H

class instance
{		
public:
	int num_node;
	int num_edge;
	int **adjacent_matrix;	//adjacent_matrix[i][j]: the vertices i and j are adjacent
	int **v_adj_matrix;		//v_adj_matrix[v]: the set of vertices that are adjacent to vertex v
	int *v_edge_cnt;		//v_edge_cnt[v]: the number of vertices that are adjacent to vertex v
	double edge_density;
	bool read_instance_file(char *, instance &);
	void print_org_graph(instance ins);
};
#endif