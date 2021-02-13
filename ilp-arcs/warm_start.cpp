#include <bits/stdc++.h>
#include <sys/resource.h>

using namespace std;

void get_warm_start(string instance_name, vector<int> *warm_solution)
{
	FILE *solution_file;
	vector<string> solutions_paths(4);
	string solution_path, warm_solution_chosen;
	int best_sol_heur, sol_value, best_sol_value, vertex_id;
	best_sol_heur = best_sol_value = INT_MAX;

	solutions_paths[0] = "solutions_gr/" + instance_name + ".sol";
	solutions_paths[1] = "solutions_wgr/" + instance_name + ".sol";
	solutions_paths[2] = "solutions_rg/" + instance_name + ".sol";
	solutions_paths[3] = "solutions_sg/" + instance_name + ".sol";

	for(int i = 0; i < 4; i++)
	{
		char buffer[1000];
		solution_path = solutions_paths[i];
		cout << solution_path << endl;
		solution_file = fopen(solution_path.c_str(), "r");
		if (!fscanf(solution_file, "%1000[^\n]\n", buffer)) printf("Erro de leitura\n"); // Skip first line
		if (!fscanf(solution_file, "%1000[^\n]\n", buffer)) printf("Erro de leitura\n"); // Skip second line
		if (!fscanf(solution_file, "%d", &sol_value)) printf("Erro de leitura\n");
		if(sol_value < best_sol_value)
		{
			best_sol_value = sol_value;
			best_sol_heur = i;
			(*warm_solution).resize(sol_value);

			for(int j = 0; j < sol_value; j++)
			{
				if(!fscanf(solution_file, "%d", &vertex_id)) printf("Erro de leitura\n");
				(*warm_solution)[j] = vertex_id;
			}
		}

		fclose(solution_file);
	}

	printf("Warm start solution: %d", best_sol_value);
	if(best_sol_heur == 0)
		printf(" (GR)\n\n");
	else if(best_sol_heur == 1)
		printf(" (WGR)\n\n");
	else if(best_sol_heur == 2)
		printf(" (RG)\n\n");
	else if(best_sol_heur == 3)
		printf(" (SG)\n\n");
}