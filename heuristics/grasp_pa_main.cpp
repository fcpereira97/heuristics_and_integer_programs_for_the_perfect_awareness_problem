#include <bits/stdc++.h>
#include <sys/resource.h>
#include "my_struct.h"
#include "load_input.h"
#include "preprocessing_by_collapse.h"
#include "preprocessing_by_components.h"
#include "grasp.h"
#include "propagation.h"

using namespace std;

// Macro for getting runtime execution
// this macro produces a time equal to the one produced by clock(),
// but does not suffer from the wraparound problem of clock()
extern int getrusage();
#define CPUTIME(ruse) (getrusage(RUSAGE_SELF,&ruse),ruse.ru_utime.tv_sec + ruse.ru_stime.tv_sec + 1e-6 * (ruse.ru_utime.tv_usec + ruse.ru_stime.tv_usec))
struct rusage grb_ruse_main; // Global variable used on time counter

// Time variables
double time_begin, total_time;
int iterations_limit, min_time, time_limit;

// Input/output variables
string input_path, output_path, solutions_directory_path, instance_name; 
FILE *input_file;
FILE *output_file;

// Graph variables
int n_vertices, n_edges, n_components, original_n_vertices, original_n_edges;
double threshold_rate;
vector<Vertex*> vertices;
vector<vector<Vertex*>> components;

//Grasp variables
long int seed_rng;
double construction_phase_parameter, local_search_phase_parameter;
string construction_phase_strategie, local_search_phase_strategie;
int n_seeds_per_insertion, target, heap_of_seed_sets_size;

void parse_args(int argc, char **argv) {
  for (int i = 1; i < argc; i++)
  {
    string str = argv[i];
    if(str[0] == '-')
    {
        if(str.substr(1) ==  "tr") 
        	threshold_rate = atof(argv[i+1]);
        else if(str.substr(1) ==  "cps")
            construction_phase_strategie = argv[i+1];
        else if(str.substr(1) ==  "cpp")
            construction_phase_parameter = atof(argv[i+1]);
        else if(str.substr(1) ==  "nspi")
            n_seeds_per_insertion = atoi(argv[i+1]);
        else if(str.substr(1) ==  "lsps")
            local_search_phase_strategie = argv[i+1];
        else if(str.substr(1) ==  "lspp")
            local_search_phase_parameter = atof(argv[i+1]);
        else if(str.substr(1) ==  "tmin")
            min_time = atoi(argv[i+1]);
        else if(str.substr(1) ==  "tmax")
            time_limit = atoi(argv[i+1]);	
        else if(str.substr(1) ==  "itmax")
            iterations_limit = atoi(argv[i+1]);
        else if(str.substr(1) ==  "target")
            target = atoi(argv[i+1]);
        else if(str.substr(1) ==  "ip")
            input_path = argv[i+1];
        else if(str.substr(1) ==  "op")
            output_path = argv[i+1];
        else if(str.substr(1) ==  "sdp")
            solutions_directory_path = argv[i+1];
        else if(str.substr(1) ==  "srng")
            seed_rng = strtol(argv[i+1], NULL, 10);
        else if(str.substr(1) ==  "hs")
            heap_of_seed_sets_size = atoi(argv[i+1]);
        else
            cout << "Invalid option: -" << str.substr(1) << "\n";

        i++;
    }
  }
}

void print_log_header(int argc, char *argv[])
{
    cout << "Command:";
    for (int i = 0; i < argc; i++)
        cout << " " << argv[i];
    cout << endl << "Seed for random number generator: " << seed_rng << endl << endl;

    cout << "Instance name: " << instance_name << endl;
    cout << "#Vertices: " << n_vertices << endl;
    cout << "#Edges: " << n_edges << endl;
    cout << "#Components: " << n_components << endl << endl;
}

void print_sol_file_header(int argc, char *argv[], FILE* sol_file)
{
    fprintf(sol_file, "Command:");
    for (int i = 0; i < argc; i++)
        fprintf(sol_file, " %s", argv[i]);
    fprintf(sol_file, "\nSeed for random number generator: %ld\n", seed_rng);
}

int main (int argc, char *argv[])
{	
	time_begin = CPUTIME(grb_ruse_main); //Start time counter 

	threshold_rate = 0.5;
	n_seeds_per_insertion = 1;
	local_search_phase_strategie = "OBO";
	min_time = 0;
	time_limit = INT_MAX;
	iterations_limit = INT_MAX;
	target = -1;
	output_path = "results.csv";
	solutions_directory_path = "solutions/";
	seed_rng = -1;
    heap_of_seed_sets_size = 0;

	parse_args(argc, argv);

    if(seed_rng == -1)
        seed_rng = chrono::steady_clock::now().time_since_epoch().count();
    mt19937 rng(seed_rng);

	// Get instance name
	instance_name = input_path.substr(input_path.find_last_of("/") + 1, input_path.find_last_of(".") - (input_path.find_last_of("/") + 1));

	// Solution file
	FILE *sol_file;
	string sol_file_path = solutions_directory_path + instance_name + ".sol";
	sol_file = fopen(sol_file_path.c_str(), "w");
	print_sol_file_header(argc, argv, sol_file);

	// Set input and output files
	input_file = fopen(input_path.c_str(), "r");
	output_file = fopen(output_path.c_str(), "a");

	// Load graph from input file
	load_input(input_file, n_vertices, n_edges, vertices, threshold_rate);

	original_n_vertices = n_vertices;
	original_n_edges = n_edges;

	preprocessing_by_collapse(n_vertices, n_edges, vertices, threshold_rate, instance_name, solutions_directory_path);
	preprocessing_by_components(n_vertices, n_components, vertices, components);

    print_log_header(argc, argv);

	time_limit += time_begin;
	min_time += time_begin;

	//print_graph(n_vertices, &vertices);

	vector<vector<Vertex*>> best_seed_set_comp(n_components);
	grasp(original_n_vertices, original_n_edges, n_components, &components,
		construction_phase_strategie, construction_phase_parameter, n_seeds_per_insertion,
		local_search_phase_strategie, local_search_phase_parameter, time_begin, min_time, time_limit,
		iterations_limit, target, instance_name, output_file, sol_file, rng,
        &best_seed_set_comp, heap_of_seed_sets_size);


    vector<Vertex*> seed_set;
    for(int i = 0; i < n_components; i++)
        for(int j = 0; j < (int)best_seed_set_comp[i].size(); j++)
            seed_set.push_back(best_seed_set_comp[i][j]);

    /*
    printf("\nChecking if the seed set is perfect...\n");
    printf("Seed set size = %d\n", (int)seed_set.size());
    int n_aware, n_rounds;
    propagate_from_initial_state(n_vertices, vertices, seed_set, n_aware, n_rounds);

    if(n_aware == n_vertices)
        printf("The seed set is perfect.\n");
    else if(n_aware < n_vertices)
        printf("The seed set is NOT perfect.\n");

    printf("\n --------------------------------------- \n");
    */

    fclose(sol_file);
	fclose(input_file);
	fclose(output_file);

    total_time = CPUTIME(grb_ruse_main) - time_begin;

    //printf("\n%d %.2f", (int)seed_set.size(), total_time);

	return 0;
}