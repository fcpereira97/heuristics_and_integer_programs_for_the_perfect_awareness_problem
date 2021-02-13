#include <bits/stdc++.h>
#include <sys/resource.h>
#include <ilcplex/ilocplex.h>
#include "my_struct.h"
#include "load_input.h"
#include "preprocessing_by_collapse.h"
#include "nodal_heuristic.h"
#include "warm_start.h"

using namespace std;

#define EPSILON 0.000000001

// Macro for getting runtime execution
// this macro produces a time equal to the one produced by clock(),
// but does not suffer from the wraparound problem of clock()
extern int getrusage();
#define CPUTIME(ruse) (getrusage(RUSAGE_SELF,&ruse),ruse.ru_utime.tv_sec + ruse.ru_stime.tv_sec + 1e-6 * (ruse.ru_utime.tv_usec + ruse.ru_stime.tv_usec))
struct rusage grb_ruse; // Global variable used on time counter

string itos(int i) {stringstream s; s << i; return s.str(); } //Macro for converting int into str

vector<Vertex*> vertices;
IloInt simplex_root_iter;
double **vars_s_aux;
int n_vertices, n_edges, original_n_vertices, original_n_edges;
int max_tau;

double obj_first_rl = -1;
double dual_root = -1;
int total_node = -1;

int node_best_primal = -1;
int best_primal = -1;
double best_dual = -1;
int cut_iter = 0;
double time_root_relax = 0;
double time_bc = 0;
double total_time = 0;

// Initialize status and data of vertices
void initialize_vertices(int n_vertices, vector<Vertex*> *vertices, double threshold_ratio)
{
	for(int i = 0; i < n_vertices; i++)
	{
		(*vertices)[i]-> degree = (*vertices)[i]-> neighbors.size();
		if((*vertices)[i]-> degree > 0)
			(*vertices)[i]-> threshold = ceil((double)(*vertices)[i]-> degree * threshold_ratio);
		else
			(*vertices)[i]-> threshold = 1;
	}
}

void print_graph(int n_vertices, vector<Vertex*> *vertices){
	for(int i = 0; i < n_vertices; i++)
	{
		cout << "Vertex " << (*vertices)[i]-> id << endl;
		cout << "Degree: " << (*vertices)[i]-> degree << endl;
		cout << "Threshold: " << (*vertices)[i]-> threshold << endl;
		cout << "Neighs: ";
		for (list<Vertex*>::iterator it = (*vertices)[i]-> neighbors.begin(); it != (*vertices)[i]-> neighbors.end(); ++it)
    		cout << (*it)-> id << " ";
    	cout << endl << endl;
	}
}

void print_relax_info()
{
	for(int j = 0; j < max_tau; j++)
	{
		for(int i = 0; i < n_vertices; i++)
		{
			printf("vars_s_aux_%d_%d = %.30lf\n", i, j, vars_s_aux[i][j]);
		}
		printf("\n");
	}
	printf("\n");
}

ILOSIMPLEXCALLBACK0(Simplex)
{
	if(cut_iter == 0)
		simplex_root_iter = getNiterations();
}

ILOHEURISTICCALLBACK1(PrimalHeur, IloArray<IloBoolVarArray>, vars_s)
{
	int sol_value = 0;
	IloNumVarArray vars(getEnv(), n_vertices*max_tau);
	IloNumArray vals(getEnv(), n_vertices*max_tau);

	for(int i = 0; i < n_vertices; i++)
		vertices[i]-> var_value = getValue(vars_s[i][0]);

	nodal_heuristic(n_vertices, &vertices, &sol_value);

	int k = 0;
	for(int i = 0; i < n_vertices; i++)
	{
		for(int j = 0; j < max_tau; j++)
		{
			vars[k] = vars_s[i][j];

			if(vertices[i]-> round_became_spreader == -1)
				vals[k] = 0;
			else if(j < vertices[i]-> round_became_spreader)
				vals[k] = 0;
			else
				vals[k] = 1;

			k++;
		}
	}

	//cout << "SOL_VALUE = " << sol_value << endl;
	setSolution(vars, vals, sol_value);
}

ILOUSERCUTCALLBACK1(Cuts, IloArray<IloBoolVarArray>, vars_s)
{
	double current_dual = getObjValue();

 	static bool first_node = true; /* Indica se é o primeiro nó */

	/* guarda o valor da função objetivo no primeiro nó */
	if (first_node && isAfterCutLoop()) {

		first_node = false;
		dual_root = current_dual;
		time_root_relax = CPUTIME(grb_ruse) - time_root_relax;
		time_bc = CPUTIME(grb_ruse);
	}

	/* guarda o valor da função objetivo da primeira relaxação */
	if (cut_iter == 0)
	{
		dual_root = obj_first_rl = current_dual;

		/* Pega a solução do LP. */
		for(int i = 0; i < n_vertices; i++)
			for(int j = 0; j < max_tau; j++)
				vars_s_aux[i][j] = getValue(vars_s[i][j]);
	}

	cut_iter++;
}

void run_cplex(string instance_name, int time_limit, vector<Vertex*> *vertices, string conf,
	bool flag_cplex_heur, bool flag_warm_start, bool flag_nodal_heur, int n_vertices,
	int n_edges, FILE* output_file)
{
	int max_tau_new_const = 0;

	if(conf == "RSP")
		max_tau_new_const = 0;
	else if(conf == "RSP-IMP2")
		max_tau_new_const = 1;
	else if(conf == "RSP-IMP1")
		max_tau_new_const = max_tau;

	IloEnv env; 	/* ambiente do cplex */
	IloModel model(env); 	  /* objeto que representa o modelo */
	IloExpr constr_expr_1(env);
	IloExpr constr_expr_2(env);
	// VARIABLES

	IloArray<IloBoolVarArray> vars_s(env, n_vertices);

 	for (int i = 0; i < n_vertices; i++)
 		vars_s[i] = IloBoolVarArray (env, max_tau + 1);

 	for (int i = 0; i < n_vertices; i++)
 	{
 	 	for(int j = 0; j <= max_tau; j++)
 	 	{
 	 		string var_name = "x_" + itos(i) + "_" + itos(j);
 			vars_s[i][j].setName(var_name.c_str());
 	 	}
 	}

	// OBJECTIVE

	// Objective function (1)
	IloExpr obj(env);
	for(int i = 0; i < n_vertices; i++)
		obj += vars_s[i][0];
	model.add(IloMinimize(env, obj));

	// CONSTRAINTS

	// Constraints (2)
	for(int i = 0; i < n_vertices; i++)
		for(int j = 1; j <= max_tau; j++)
			model.add(vars_s[i][j-1] <= vars_s[i][j]);

	// Constraints (3)
	for(int i = 0; i < n_vertices; i++)
	{
		for (int j = 1; j <= max_tau; j++)
		{
			for (list<Vertex*>::iterator neigh = (*vertices)[i]-> neighbors.begin(); neigh != (*vertices)[i]-> neighbors.end(); ++neigh)
			{
				constr_expr_1 += vars_s[(*neigh)->id][j-1];
			}
			
			model.add((vars_s[i][j] - vars_s[i][0]) * (*vertices)[i]->threshold <= constr_expr_1);
			constr_expr_1.clear();
		}
	}
	
	// Constraints (4)
	for(int i = 0; i < n_vertices; i++)
	{
		for (int j = 1; j <= max_tau; j++)
		{
			for (list<Vertex*>::iterator neigh = (*vertices)[i]-> neighbors.begin(); neigh != (*vertices)[i]-> neighbors.end(); ++neigh)
			{
				constr_expr_1+= vars_s[(*neigh)->id][j-1];
			}
			model.add(constr_expr_1 <= ((*vertices)[i]->threshold - 1) + ((*vertices)[i]->degree - (*vertices)[i]->threshold + 1) * vars_s[i][j]);
			constr_expr_1.clear();
		}
	}

	// Constraints (5)
	
	for(int i = 0; i < n_vertices; i++)
	{
		int j = max_tau;
		for (list<Vertex*>::iterator neigh = (*vertices)[i]-> neighbors.begin(); neigh != (*vertices)[i]-> neighbors.end(); ++neigh)
		{
			constr_expr_1 += vars_s[(*neigh)->id][j-1];
		}
		model.add(1 <= vars_s[i][0] + constr_expr_1);
		constr_expr_1.clear();
	}
	
	// Constraints (6)
	for(int i = 1; i <= max_tau_new_const; i++)
	{
		for(int j = 0; j < n_vertices; j++)
		{
			constr_expr_1 += vars_s[j][i] - vars_s[j][i-1];
		}

		for(int j = 0; j < n_vertices; j++)
		{
			for (list<Vertex*>::iterator neigh = (*vertices)[j]-> neighbors.begin(); neigh != (*vertices)[j]-> neighbors.end(); ++neigh)
			{
				constr_expr_2 += vars_s[(*neigh)->id][i];
			}
			constr_expr_2 += vars_s[j][i];
			model.add(1 <= constr_expr_1 + constr_expr_2);
			constr_expr_2.clear();
		}

		constr_expr_1.clear();
	}

	// CPLEX PARAMETERS
	IloCplex cplex(env);	/* cria objeto do cplex */
	cplex.extract(model); 	/* carrega o modelo */
	cplex.setParam(IloCplex::Param::TimeLimit, time_limit); /* limita o tempo de execucao */
	cplex.setParam(IloCplex::Param::MIP::Strategy::Search, CPX_MIPSEARCH_TRADITIONAL);	/* Turn on traditional search for use with control callbacks */
	cplex.setParam(IloCplex::Param::Threads, 1);	/* Desabilita paralelismo  */
	cplex.setParam(IloCplex::Param::Simplex::Tolerances::Feasibility, 1e-9); /* Specifies the feasibility tolerance, that is, the degree to which the basic variables of a model may violate their bounds.*/
	cplex.setParam(IloCplex::Param::MIP::Display, 4);
	//cplex.setParam(IloCplex::Param::MIP::Limits::CutsFactor, 1.0); /* Desabilita cortes do CPLEX*/
	//cplex.setParam(IloCplex::Param::MIP::Limits::CutPasses, 1000); /* Sets the upper limit on the number of cutting plane passes CPLEX performs when solving the root node of a MIP model*/	
	//cplex.setParam(IloCplex::Param::MIP::Limits::Nodes, 1); /* Limita o numero de nós visitados pelo CPLEX*/
	//cplex.setParam(IloCplex::Param::MIP::Tolerances::Integrality, 1e-16); /* Diminui a tolerância de integralidade*/
	//cplex.setParam(IloCplex::Param::Feasopt::Tolerance, 1e-16); /* Controls the amount of relaxation for the routine CPXfeasopt in the C API or for the method feasOpt in the object-oriented APIs*/
	//cplex.setParam(IloCplex::Param::MIP::Tolerances::Linearization, 1e-16); /* Diminui a tolerância na linearização*/
	if(!flag_cplex_heur)
	{
		cplex.setParam(IloCplex::Param::Preprocessing::Presolve, CPX_OFF); /* Desabilita o PRESOLVE*/
		cplex.setParam(IloCplex::Param::MIP::Strategy::HeuristicFreq, -1); /* Desabilita heurísticas */
		cplex.setParam(IloCplex::Param::MIP::Strategy::RINSHeur, -1); /* Desabilita heurísticas */
		cplex.setParam(IloCplex::Param::MIP::Strategy::FPHeur, -1); /* Desabilita heurísticas */
	}

	cplex.use(Cuts(env, vars_s));
	cplex.use(Simplex(env));

	if(flag_nodal_heur)
		cplex.use(PrimalHeur(env, vars_s));
    //cplex.exportModel("model.lp");

	// Warm start
	int sol_warm_value = -1;
	if(flag_warm_start)
	{
		cout << "Getting Warm Solution" << endl;
		
		IloNumArray vals_warm(env, n_vertices);
		IloNumVarArray vars_warm(env, n_vertices);

		vector<int> warm_solution;
		get_warm_start(instance_name, &warm_solution);
		sol_warm_value = (int)warm_solution.size();

		for(int i = 0; i < n_vertices; i++)
		{
			vars_warm[i] = vars_s[i][0];
			vals_warm[i] = 0;
		}
		for(int i = 0; i < (int)warm_solution.size(); i++)
		{
			vals_warm[warm_solution[i]] = 1;
		}

		cplex.addMIPStart(vars_warm, vals_warm);
	}

    time_root_relax = CPUTIME(grb_ruse);

    // Solve model
    bool found_incumbent = cplex.solve();

    // Export solution file
    string solution_path = "solutions_pli/" + instance_name + ".mst";
    cplex.writeSolution(solution_path.c_str());

    if(time_bc != 0)
    	time_bc = CPUTIME(grb_ruse) - time_bc;
    else
    	time_root_relax = CPUTIME(grb_ruse) - time_root_relax;

    total_time = time_bc + time_root_relax;
	total_node = cplex.getNnodes();
	best_dual = cplex.getBestObjValue();

	if (found_incumbent) {

	    best_primal = cplex.getObjValue();
	    node_best_primal = cplex.getIncumbentNode();
	    if(best_primal == best_dual)
   			printf("\n\nMain: programa achou uma solução ótima!\n\n");
   		else
   				printf("\n\nMain: programa terminou e achou solucao inteira!\n\n");

	    printf("best_primal: %d\n", best_primal);
	    printf("node_best_primal: %d\n", node_best_primal);
	}
	else
	{
    	printf("\n\nMain: programa terminou sem achar solucao inteira!\n\n");
   	}

   	printf("obj_first_rl: %.8lf\n", obj_first_rl);
	printf("dual_root: %.8lf\n", dual_root);
    printf("best_dual: %.8lf\n", best_dual);
	printf("total_node: %d\n", total_node);
	printf("simplex_root_iter: %d\n", int(simplex_root_iter));
	printf("time_root_relax: %.2lf\n", time_root_relax);
	printf("time_bc: %.2lf\n", time_bc);
	printf("total_time: %.2lf\n", total_time);
   	printf("\n\n");

   	fprintf(output_file, "%s, %d, %d, %s, %d, %d, %d, %.8lf, %.8lf, %.8lf, %d, %d, %.8lf, %.8lf, %.8lf\n",
   		instance_name.c_str(), original_n_vertices, original_n_edges, conf.c_str(), sol_warm_value, best_primal, node_best_primal, obj_first_rl, dual_root,
   		best_dual, total_node, int(simplex_root_iter), time_root_relax, time_bc, total_time);

	//print_relax_info();
}

int main (int argc, char *argv[])
{
	FILE *input_file, *output_file;
	string input_path, output_path, conf, instance_name, solutions_directory_path;
	bool flag_nodal_heur, flag_warm_start, flag_cplex_heur;
	int time_limit;
	double threshold_rate;

	input_path = argv[1];
	threshold_rate = atof(argv[2]);
	conf = argv[3];
	flag_cplex_heur = atoi(argv[4]);
	flag_warm_start = atoi(argv[5]);
	flag_nodal_heur = atoi(argv[6]);
	time_limit = atoi(argv[7]);
	output_path = argv[8];
	input_file = fopen(input_path.c_str(), "r");
	output_file = fopen(output_path.c_str(), "a");

		// Get instance name
	instance_name = input_path.substr(input_path.find_last_of("/") + 1, input_path.find_last_of(".") - (input_path.find_last_of("/") + 1));
	solutions_directory_path = "solutions_pli/";

	// Load graph from input file
	load_input(input_file, n_vertices, n_edges, vertices, threshold_rate);
	original_n_vertices = n_vertices;
	original_n_edges = n_edges;

	preprocessing_by_collapse(n_vertices, n_edges, vertices, threshold_rate, instance_name, solutions_directory_path);

	cut_iter = 0;
	simplex_root_iter = 0;
	max_tau = n_vertices;

	vars_s_aux = new double*[n_vertices];
	for (int i = 0; i < n_vertices; i++)
		vars_s_aux[i] = new double[n_vertices];

	for(int i = 1; i < argc; i++)
    	printf("%s ", argv[i]);
    printf("\n");

    print_graph(n_vertices, &vertices);

	run_cplex(instance_name, time_limit, &vertices, conf, flag_cplex_heur, flag_warm_start, flag_nodal_heur,
		n_vertices, n_edges, output_file);
	
	fclose(input_file);
	fclose(output_file);
	return 0;
}