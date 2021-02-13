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
double *vars_s_aux;
double **vars_a_aux;
double *vars_s_rlx;
double **vars_a_rlx;
int n_vertices, n_edges, original_n_vertices, original_n_edges;

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
int itersep = 0;


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
	for(int i = 0; i < n_vertices; i++)
	{
		printf("vars_s_rlx[%d] = %.30lf\n", i, vars_s_rlx[i]);
	}
	printf("\n");

	for(int i = 0; i < n_vertices; i++)
	{
		for(int j = 0; j < n_vertices; j++)
		{
			printf("vars_a_rlx[%d][%d] = %.30lf\n", i, j, vars_a_rlx[i][j]);
		}
	}
	printf("\n");
}

void dfs(int cur_index, vector<Vertex*>* vertices, IloEnv *env,
IloArray<IloBoolVarArray> *vars_a, queue<IloExpr> *cuts_lhs, queue<IloExpr> *cuts_rhs)
{
	//Exploration started
	(*vertices)[cur_index]-> color = GRAY;

	//For each neigh of curr index
	for (list<Vertex*>::iterator neigh = (*vertices)[cur_index]-> neighbors.begin(); neigh != (*vertices)[cur_index]-> neighbors.end(); ++neigh)
	{
		//If not visited
		if((*neigh)-> color == WHITE && vars_a_aux[cur_index][(*neigh)-> index] > 0)
		{
			//Set parent and visit
			(*neigh)-> parent = cur_index;
			dfs((*neigh)-> index, vertices, env, vars_a, cuts_lhs, cuts_rhs);
		}
		else if ((*neigh)-> color == GRAY && vars_a_aux[cur_index][(*neigh)-> index] > 0) // If neigh is not finished yet so a cycled was detected
		{
			// Backtracking cycle to build inequality
			IloExpr constr_lhs;
			IloExpr constr_rhs;
			constr_lhs = IloExpr(*env);
			constr_rhs = IloExpr(*env);

			int size = 0;
			double sum = 0;
			int u = cur_index;
			int v = (*neigh)-> index;

			sum += vars_a_aux[u][v];
			constr_lhs += (*vars_a)[u][v];
			size++;

			while(u != (*neigh)-> index)
			{
				v = u;
				u = (*vertices)[u]-> parent;

				sum += vars_a_aux[u][v];
				constr_lhs += (*vars_a)[u][v];
				size++;
			} 

			constr_rhs += size;
			if(sum - size + 1 > EPSILON) // If violates constraint (6)
			{
				// Add inequality
				//printf("sum = %.16lf size = %d\n", sum, size - 1);
				(*cuts_lhs).push(constr_lhs);
				(*cuts_rhs).push(constr_rhs - 1);
			}
		}
	}

	//Exploration finished
	(*vertices)[cur_index]-> color = BLACK;
}

//Check violation of constr6
void sep_cycle(vector<Vertex*>* vertices, IloEnv *env, IloArray<IloBoolVarArray> *vars_a, queue<IloExpr> *cuts_lhs, queue<IloExpr> *cuts_rhs)
{

	for(int i = 0; i < n_vertices; i++)
	{
		(*vertices)[i]-> color = WHITE;
		(*vertices)[i]-> parent = -1;
	}
	for(int i = 0; i < n_vertices; i++)
		if((*vertices)[i]-> color == WHITE)
			dfs(i, vertices, env, vars_a, cuts_lhs, cuts_rhs);
}


// Callback de corte
ILOUSERCUTCALLBACK3(Cuts_Cycle, IloBoolVarArray, vars_s, IloArray<IloBoolVarArray>, vars_a, vector<Vertex*>*, vertices)
{
	//cout << "Itersep = " << itersep << endl;

	//Recupera Ambiente
	IloEnv env = getEnv();

	for(int i = 0; i < n_vertices; i++)
		vars_s_aux[i] = getValue(vars_s[i]);

	for(int i = 0; i < n_vertices; i++)
		for (list<Vertex*>::iterator neigh = (*vertices)[i]-> neighbors.begin(); neigh != (*vertices)[i]-> neighbors.end(); ++neigh)
			vars_a_aux[i][(*neigh)-> index] = getValue(vars_a[i][(*neigh)-> index]);

	double lpobjval = getObjValue();
 	static bool first_node = true; // Indica se é o primeiro nó 

	// guarda o valor da função objetivo no primeiro nó
	if (first_node && isAfterCutLoop()) {
		first_node = false;
		dual_root = lpobjval;
	}

	queue<IloExpr> cuts_lhs;
	queue<IloExpr> cuts_rhs;

	sep_cycle(vertices, &env, &vars_a, &cuts_lhs, &cuts_rhs);

	while(!cuts_lhs.empty())
	{
		//cout << cuts_lhs.front() << " <= " << cuts_rhs.front() << endl;
		add(cuts_lhs.front() <= cuts_rhs.front());
		cuts_lhs.pop();
		cuts_rhs.pop();
	}

	itersep++;
}


// Calback de corte com lazy constraints
ILOLAZYCONSTRAINTCALLBACK3(Lazy_Cuts_Cycle, IloBoolVarArray, vars_s, IloArray<IloBoolVarArray>, vars_a, vector<Vertex*>*, vertices)
{
	//cout << "Itersep = " << itersep << endl;

	//Recupera Ambiente
	IloEnv env = getEnv();

	for(int i = 0; i < n_vertices; i++)
		vars_s_aux[i] = getValue(vars_s[i]);

	for(int i = 0; i < n_vertices; i++)
		for (list<Vertex*>::iterator neigh = (*vertices)[i]-> neighbors.begin(); neigh != (*vertices)[i]-> neighbors.end(); ++neigh)
			vars_a_aux[i][(*neigh)-> index] = getValue(vars_a[i][(*neigh)-> index]);

	queue<IloExpr> cuts_lhs;
	queue<IloExpr> cuts_rhs;

	sep_cycle(vertices, &env, &vars_a, &cuts_lhs, &cuts_rhs);

	while(!cuts_lhs.empty())
	{
		//cout << cuts_lhs.front() << " <= " << cuts_rhs.front() << endl;
		add(cuts_lhs.front() <= cuts_rhs.front());
		cuts_lhs.pop();
		cuts_rhs.pop();
	}

	itersep++;
}

ILOSIMPLEXCALLBACK0(Simplex)
{
	if(cut_iter == 0)
		simplex_root_iter = getNiterations();
}

ILOUSERCUTCALLBACK3(Root_Relaxation, IloBoolVarArray, vars_s, IloArray<IloBoolVarArray>, vars_a, vector<Vertex*>*, vertices)
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
			vars_s_rlx[i] = getValue(vars_s[i]);

	for(int i = 0; i < n_vertices; i++)
		for (list<Vertex*>::iterator neigh = (*vertices)[i]-> neighbors.begin(); neigh != (*vertices)[i]-> neighbors.end(); ++neigh)
			vars_a_rlx[i][(*neigh)-> index] = getValue(vars_a[i][(*neigh)-> index]);
	}

	cut_iter++;
}

void run_cplex(string instance_name, int time_limit, vector<Vertex*> *vertices,
	bool flag_cplex_heur, bool flag_warm_start, bool flag_nodal_heur, int n_vertices,
	int n_edges, FILE* output_file)
{

	IloEnv env; 	/* ambiente do cplex */
	IloModel model(env); 	  /* objeto que representa o modelo */
	IloExpr constr_expr_1(env);
	IloExpr constr_expr_2(env);

	// VARIABLES

	IloBoolVarArray vars_s(env, n_vertices);
 	IloArray<IloBoolVarArray> vars_a(env, n_vertices);

 	for (int i = 0; i < n_vertices; i++)
 	{
 		vars_a[i] = IloBoolVarArray (env, n_vertices);

 		string var_name = "s_" + itos(i);
 		vars_s[i].setName(var_name.c_str());
 	}

 	for (int i = 0; i < n_vertices; i++)
 	{
 		for (int j = 0; j < n_vertices; j++)
 		{
 			string var_name = "a_" + itos(i) + "_" + itos(j);
 			vars_a[i][j].setName(var_name.c_str());
 		}
 	}

	// OBJECTIVE

	// Objective function (1)
	IloExpr obj(env);
	for(int i = 0; i < n_vertices; i++)
		obj += vars_s[i];
	model.add(IloMinimize(env, obj));

	// CONSTRAINTS

	// Constraints (2)
	for(int i = 0; i < n_vertices; i++)
	{
		for (list<Vertex*>::iterator neigh = (*vertices)[i]-> neighbors.begin(); neigh != (*vertices)[i]-> neighbors.end(); ++neigh)
		{
			constr_expr_1 += vars_a[(*neigh)-> index][i];
		}
		model.add(constr_expr_1 - (*vertices)[i]->degree * (1 - vars_s[i]) <= 0);
		constr_expr_1.clear();
	}

	// Constraints (3)
	for(int i = 0; i < n_vertices; i++)
	{
		for (list<Vertex*>::iterator neigh = (*vertices)[i]-> neighbors.begin(); neigh != (*vertices)[i]-> neighbors.end(); ++neigh)
		{
			constr_expr_1 += vars_a[(*neigh)-> index][i];
		}

		for (list<Vertex*>::iterator neigh = (*vertices)[i]-> neighbors.begin(); neigh != (*vertices)[i]-> neighbors.end(); ++neigh)
		{
			model.add(constr_expr_1 - (*vertices)[i]->threshold * (vars_a[i][(*neigh)-> index] - vars_s[i]) >= 0);
		}
		constr_expr_1.clear();
	}

	// Constraints (4)
	for(int i = 0; i < n_vertices; i++)
	{
		for (list<Vertex*>::iterator neigh = (*vertices)[i]-> neighbors.begin(); neigh != (*vertices)[i]-> neighbors.end(); ++neigh)
		{
			constr_expr_1 += vars_a[(*neigh)-> index][i];
		}
		model.add(constr_expr_1 + vars_s[i] >= 1);
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

	cplex.use(Lazy_Cuts_Cycle(env, vars_s, vars_a, vertices));
	//cplex.use(Cuts_Cycle(env, vars_s, vars_a, vertices));
	cplex.use(Root_Relaxation(env, vars_s, vars_a, vertices));
	cplex.use(Simplex(env));

	//if(flag_nodal_heur)
		//cplex.use(PrimalHeur(env, vars_s));
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
			vars_warm[i] = vars_s[i];
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

   	fprintf(output_file, "%s, %d, %d, %d, %d, %d, %.8lf, %.8lf, %.8lf, %d, %d, %.8lf, %.8lf, %.8lf\n",
   		instance_name.c_str(), original_n_vertices, original_n_edges, sol_warm_value, best_primal, node_best_primal, obj_first_rl, dual_root,
   		best_dual, total_node, int(simplex_root_iter), time_root_relax, time_bc, total_time);

	print_relax_info();
}

int main (int argc, char *argv[])
{
	FILE *input_file, *output_file;
	string input_path, output_path, instance_name, solutions_directory_path;
	bool flag_nodal_heur, flag_warm_start, flag_cplex_heur;
	int time_limit;
	double threshold_rate;

	input_path = argv[1];
	threshold_rate = atof(argv[2]);
	flag_cplex_heur = atoi(argv[3]);
	flag_warm_start = atoi(argv[4]);
	flag_nodal_heur = atoi(argv[5]);
	time_limit = atoi(argv[6]);
	output_path = argv[7];
	input_file = fopen(input_path.c_str(), "r");
	output_file = fopen(output_path.c_str(), "a");

		// Get instance name
	instance_name = input_path.substr(input_path.find_last_of("/") + 1, input_path.find_last_of(".") - (input_path.find_last_of("/") + 1));
	solutions_directory_path = "solutions_pli/";

	// Load graph from input file
	load_input(input_file, n_vertices, n_edges, vertices, threshold_rate);
	original_n_vertices = n_vertices;
	original_n_edges = n_edges;

	//preprocessing_by_collapse(n_vertices, n_edges, vertices, threshold_rate, instance_name, solutions_directory_path);

	cut_iter = 0;
	simplex_root_iter = 0;

	vars_s_aux = new double[n_vertices];
	vars_a_aux = new double*[n_vertices];
	for (int i = 0; i < n_vertices; i++)
		vars_a_aux[i] = new double[n_vertices];


	vars_s_rlx = new double[n_vertices];
	vars_a_rlx = new double*[n_vertices];
	for (int i = 0; i < n_vertices; i++)
		vars_a_rlx[i] = new double[n_vertices];

	for(int i = 1; i < argc; i++)
    	printf("%s ", argv[i]);
    printf("\n");

    for(int i = 0; i < n_vertices; i++)
	{
		vertices[i]-> index = i;
	}

    print_graph(n_vertices, &vertices);
    
    itersep = 0;
	run_cplex(instance_name, time_limit, &vertices, flag_cplex_heur, flag_warm_start, flag_nodal_heur,
		n_vertices, n_edges, output_file);
	
	fclose(input_file);
	fclose(output_file);
	return 0;
}