/**
 * An example of benchmarking random search on a COCO suite. A grid search optimizer is also
 * implemented and can be used instead of random search.
 *
 * Set the global parameter BUDGET_MULTIPLIER to suit your needs.
 */
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <assert.h>

#include "coco.h"
#define max(a,b) ((a) > (b) ? (a) : (b))
#define N 100
#define REPAIR 0//0:lamarckian,1:baldwinian
#define TITLE "20D-lamarckian-.txt"
#define D 20
int instance_cnt = 0;
/**
 * The maximal budget for evaluations done by an optimization algorithm equals dimension * BUDGET_MULTIPLIER.
 * Increase the budget multiplier value gradually to see how it affects the runtime.
 */
static const unsigned int BUDGET_MULTIPLIER = 10000;

/**
 * The maximal number of independent restarts allowed for an algorithm that restarts itself.
 */
static const long INDEPENDENT_RESTARTS = 1e5;

/**
 * The random seed. Change if needed.
 */
static const uint32_t RANDOM_SEED = 0xdeadbeef;

/**
 * A function type for evaluation functions, where the first argument is the vector to be evaluated and the
 * second argument the vector to which the evaluation result is stored.
 */
typedef void (*evaluate_function_t)(const double *x, double *y);

/**
 * A pointer to the problem to be optimized (needed in order to simplify the interface between the optimization
 * algorithm and the COCO platform).
 */
static coco_problem_t *PROBLEM; 

/**
 * Calls coco_evaluate_function() to evaluate the objective function
 * of the problem at the point x and stores the result in the vector y
 */
static void evaluate_function(const double *x, double *y) {
  coco_evaluate_function(PROBLEM, x, y);
}

/**
 * Calls coco_evaluate_constraint() to evaluate the constraints
 * of the problem at the point x and stores the result in the vector y
 */
static void evaluate_constraint(const double *x, double *y) {
  coco_evaluate_constraint(PROBLEM, x, y);
}

/* Declarations of all functions implemented in this file (so that their order is not important): */
void example_experiment(const char *suite_name,
                        const char *suite_options,
                        const char *observer_name,
                        const char *observer_options,
                        coco_random_state_t *random_generator);

void my_random_search(evaluate_function_t evaluate_func,
                      evaluate_function_t evaluate_cons,
                      const size_t dimension,
                      const size_t number_of_objectives,
                      const size_t number_of_constraints,
                      const double *lower_bounds,
                      const double *upper_bounds,
                      const size_t number_of_integer_variables,
                      const size_t max_budget,
                      coco_random_state_t *random_generator);

void my_grid_search(evaluate_function_t evaluate_func,
                    evaluate_function_t evaluate_cons,
                    const size_t dimension,
                    const size_t number_of_objectives,
                    const size_t number_of_constraints,
                    const double *lower_bounds,
                    const double *upper_bounds,
                    const size_t number_of_integer_variables,
                    const size_t max_budget);

void my_de_nopcm(evaluate_function_t evaluate_func,
                    evaluate_function_t evaluate_cons,
                    const size_t dimension,
                    const size_t number_of_objectives,
                    const size_t number_of_constraints,
                    const double *lower_bounds,
                    const double *upper_bounds,
                    const size_t number_of_integer_variables,
                    const size_t max_budget,
                    coco_random_state_t *random_generator);

/* Structure and functions needed for timing the experiment */
typedef struct {
	size_t number_of_dimensions;
	size_t current_idx;
	char **output;
	size_t previous_dimension;
	size_t cumulative_evaluations;
	time_t start_time;
	time_t overall_start_time;
} timing_data_t;
static timing_data_t *timing_data_initialize(coco_suite_t *suite);
static void timing_data_time_problem(timing_data_t *timing_data, coco_problem_t *problem);
static void timing_data_finalize(timing_data_t *timing_data);

/**
 * The main method initializes the random number generator and calls the example experiment on the
 * bbob suite.
 */
int main(void) {

  coco_random_state_t *random_generator = coco_random_new(RANDOM_SEED);

  /* Change the log level to "warning" to get less output */
  coco_set_log_level("info");

  printf("Running the example experiment... (might take time, be patient)\n");
  fflush(stdout);

  /**
   * Start the actual experiments on a test suite and use a matching logger, for
   * example one of the following:
   *   bbob                 24 unconstrained noiseless single-objective functions
   *   bbob-biobj           55 unconstrained noiseless bi-objective functions
   *   [bbob-biobj-ext       92 unconstrained noiseless bi-objective functions]
   *   [bbob-constrained*   48 constrained noiseless single-objective functions]
   *   bbob-largescale      24 unconstrained noiseless single-objective functions in large dimension
   *   bbob-mixint          24 unconstrained noiseless single-objective functions with mixed-integer variables
   *   bbob-biobj-mixint    92 unconstrained noiseless bi-objective functions with mixed-integer variables
   *
   * Suites with a star are partly implemented but not yet fully supported.
   *
   * Adapt to your need. Note that the experiment is run according
   * to the settings, defined in example_experiment(...) below.
   */
  coco_set_log_level("info");

  /**
   * For more details on how to change the default suite and observer options, see
   * http://numbbo.github.io/coco-doc/C/#suite-parameters and
   * http://numbbo.github.io/coco-doc/C/#observer-parameters. */

  if(REPAIR == 0){
    example_experiment("bbob-mixint", "", "bbob-mixint", "result_folder:lamarckian", random_generator);
  }
  else{
    example_experiment("bbob-mixint", "", "bbob-mixint", "result_folder:baldwinian", random_generator);
  }

  printf("Done!\n");
  fflush(stdout);

  coco_random_free(random_generator);

  return 0;
}

/**
 * A simple example of benchmarking random search on a given suite with default instances
 * that can serve also as a timing experiment.
 *
 * @param suite_name Name of the suite (e.g. "bbob" or "bbob-biobj").
 * @param suite_options Options of the suite (e.g. "dimensions: 2,3,5,10,20 instance_indices: 1-5").
 * @param observer_name Name of the observer matching with the chosen suite (e.g. "bbob-biobj"
 * when using the "bbob-biobj-ext" suite).
 * @param observer_options Options of the observer (e.g. "result_folder: folder_name")
 * @param random_generator The random number generator.
 */
void example_experiment(const char *suite_name,
                        const char *suite_options,
                        const char *observer_name,
                        const char *observer_options,
                        coco_random_state_t *random_generator) {

  size_t run;
  coco_suite_t *suite;
  coco_observer_t *observer;
  timing_data_t *timing_data;

  /* Initialize the suite and observer. */
  suite = coco_suite(suite_name, "", suite_options);
  observer = coco_observer(observer_name, observer_options);

  /* Initialize timing */
  timing_data = timing_data_initialize(suite);

  /* Iterate over all problems in the suite */
  while ((PROBLEM = coco_suite_get_next_problem(suite, observer)) != NULL) {
    //size_t function_id = coco_problem_get_suite_dep_index(PROBLEM);
    // if(function_id == 0){
      const char *function_name = coco_problem_get_name(PROBLEM);
      size_t dimension = coco_problem_get_dimension(PROBLEM);
      
      // for(int m = 0; m < dimention; m++){
      //   printf("%lf ",bsf[i]);
      // }
      //printf("%s\n",function_name);
      if(strstr(function_name, "f001") != NULL && dimension == D) {
        /* Run the algorithm at least once */
        for (run = 1; run <= 1 + INDEPENDENT_RESTARTS; run++) {
          long evaluations_done = (long) (coco_problem_get_evaluations(PROBLEM) +
                coco_problem_get_evaluations_constraints(PROBLEM));
          long evaluations_remaining = (long) (dimension * BUDGET_MULTIPLIER) - evaluations_done;

          /* Break the loop if the target was hit or there are no more remaining evaluations */
          if ((coco_problem_final_target_hit(PROBLEM) &&
              coco_problem_get_number_of_constraints(PROBLEM) == 0)
              || (evaluations_remaining <= 0))
            break;

          /* Call the optimization algorithm for the remaining number of evaluations */
          my_de_nopcm(evaluate_function,
                          evaluate_constraint,
                          dimension,
                          coco_problem_get_number_of_objectives(PROBLEM),
                          coco_problem_get_number_of_constraints(PROBLEM),
                          coco_problem_get_smallest_values_of_interest(PROBLEM),
                          coco_problem_get_largest_values_of_interest(PROBLEM),
                          coco_problem_get_number_of_integer_variables(PROBLEM),
                          (size_t) evaluations_remaining,
                          random_generator);

          /* Break the loop if the algorithm performed no evaluations or an unexpected thing happened */
          if (coco_problem_get_evaluations(PROBLEM) == evaluations_done) {
            printf("WARNING: Budget has not been exhausted (%lu/%lu evaluations done)!\n",
                (unsigned long) evaluations_done, (unsigned long) dimension * BUDGET_MULTIPLIER);
            break;
          }
          else if (coco_problem_get_evaluations(PROBLEM) < evaluations_done)
            coco_error("Something unexpected happened - function evaluations were decreased!");
        }
      }
      /* Keep track of time */
      timing_data_time_problem(timing_data, PROBLEM);
    //}
    // else{
    //   break;
    // }
  }
  // }

  /* Output and finalize the timing data */
  timing_data_finalize(timing_data);

  coco_observer_free(observer);
  coco_suite_free(suite);

}

/**
 * A random search algorithm that can be used for single- as well as multi-objective optimization. The
 * problem's initial solution is evaluated first.
 *
 * @param evaluate_func The function used to evaluate the objective function.
 * @param evaluate_cons The function used to evaluate the constraints.
 * @param dimension The number of variables.
 * @param number_of_objectives The number of objectives.
 * @param number_of_constraints The number of constraints.
 * @param lower_bounds The lower bounds of the region of interested (a vector containing dimension values).
 * @param upper_bounds The upper bounds of the region of interested (a vector containing dimension values).
 * @param number_of_integer_variables The number of integer variables (if > 0, all integer variables come
 * before any continuous ones).
 * @param max_budget The maximal number of evaluations.
 * @param random_generator Pointer to a random number generator able to produce uniformly and normally
 * distributed random numbers.
 */
void my_random_search(evaluate_function_t evaluate_func,
                      evaluate_function_t evaluate_cons,
                      const size_t dimension,
                      const size_t number_of_objectives,
                      const size_t number_of_constraints,
                      const double *lower_bounds,
                      const double *upper_bounds,
                      const size_t number_of_integer_variables,
                      const size_t max_budget,
                      coco_random_state_t *random_generator) {
  double *x = coco_allocate_vector(dimension);
  double *functions_values = coco_allocate_vector(number_of_objectives);
  double *constraints_values = NULL;
  double range;
  size_t i, j;

  if (number_of_constraints > 0 )
    constraints_values = coco_allocate_vector(number_of_constraints);

  coco_problem_get_initial_solution(PROBLEM, x);
  // for(int tmp = 0; tmp < dimension; tmp++){
  //   printf("%lf ",upper_bounds[tmp]);
  // }
  // printf("\n");
  evaluate_func(x, functions_values);

  for (i = 1; i < max_budget; ++i) {

    /* Construct x as a random point between the lower and upper bounds */
    for (j = 0; j < dimension; ++j) {
      range = upper_bounds[j] - lower_bounds[j];
      x[j] = lower_bounds[j] + coco_random_uniform(random_generator) * range;
      /* Round the variable if integer */
      if (j < number_of_integer_variables)
        x[j] = floor(x[j] + 0.5);
    }

    /* Evaluate COCO's constraints function if problem is constrained */
    if (number_of_constraints > 0 )
      evaluate_cons(x, constraints_values);

    /* Call COCO's evaluate function where all the logging is performed */
    evaluate_func(x, functions_values);

  }

  coco_free_memory(x);
  coco_free_memory(functions_values);
  if (number_of_constraints > 0 )
    coco_free_memory(constraints_values);
}

/**
 * A grid search optimizer that can be used for single- as well as multi-objective optimization.
 *
 * @param evaluate_func The evaluation function used to evaluate the solutions.
 * @param evaluate_cons The function used to evaluate the constraints.
 * @param dimension The number of variables.
 * @param number_of_objectives The number of objectives.
 * @param number_of_constraints The number of constraints.
 * @param lower_bounds The lower bounds of the region of interested (a vector containing dimension values).
 * @param upper_bounds The upper bounds of the region of interested (a vector containing dimension values).
 * @param number_of_integer_variables The number of integer variables (if > 0, all integer variables come
 * before any continuous ones).
 * @param max_budget The maximal number of evaluations.
 *
 * If max_budget is not enough to cover even the smallest possible grid, only the first max_budget
 * nodes of the grid are evaluated.
 */
void my_grid_search(evaluate_function_t evaluate_func,
                    evaluate_function_t evaluate_cons,
                    const size_t dimension,
                    const size_t number_of_objectives,
                    const size_t number_of_constraints,
                    const double *lower_bounds,
                    const double *upper_bounds,
                    const size_t number_of_integer_variables,
                    const size_t max_budget) {


  double *x = coco_allocate_vector(dimension);
  double *func_values = coco_allocate_vector(number_of_objectives);
  double *cons_values = NULL;
  long *nodes = (long *) coco_allocate_memory(sizeof(long) * dimension);
  double *grid_step = coco_allocate_vector(dimension);
  size_t i, j;
  size_t evaluations = 0;
  long *max_nodes = (long *) coco_allocate_memory(sizeof(long) * dimension);
  long integer_nodes = 1;

  /* Initialization */
  for (j = 0; j < dimension; j++) {
    nodes[j] = 0;
    if (j < number_of_integer_variables) {
      grid_step[j] = 1;
      max_nodes[j] = (long) floor(upper_bounds[j] + 0.5);
      assert(fabs(lower_bounds[j]) < 1e-6);
      assert(max_nodes[j] > 0);
      integer_nodes *= max_nodes[j];
    }
    else {
      max_nodes[j] = (long) floor(pow((double) max_budget / (double) integer_nodes,
          1 / (double) (dimension - number_of_integer_variables))) - 1;
      /* Take care of the borderline case */
      if (max_nodes[j] < 1)
        max_nodes[j] = 1;
      grid_step[j] = (upper_bounds[j] - lower_bounds[j]) / (double) max_nodes[j];
    }
  }

  if (number_of_constraints > 0 )
    cons_values = coco_allocate_vector(number_of_constraints);

  while (evaluations < max_budget) {

    /* Stop if there are no more nodes */
    if ((number_of_integer_variables == dimension) && (evaluations >= integer_nodes))
      break;

    /* Construct x and evaluate it */
    for (j = 0; j < dimension; j++) {
      x[j] = lower_bounds[j] + grid_step[j] * (double) nodes[j];
    }

    /* Evaluate COCO's constraints function if problem is constrained */
    if (number_of_constraints > 0 )
      evaluate_cons(x, cons_values);

    /* Call COCO's evaluate function where all the logging is performed */
    evaluate_func(x, func_values);
    evaluations++;

    /* Inside the grid, move to the next node */
    if (nodes[0] < max_nodes[0]) {
      nodes[0]++;
    }

    /* At an outside node of the grid, move to the next level */
    else if (max_nodes[0] > 0) {
      for (j = 1; j < dimension; j++) {
        if (nodes[j] < max_nodes[j]) {
          nodes[j]++;
          for (i = 0; i < j; i++)
            nodes[i] = 0;
          break;
        }
      }

      /* At the end of the grid, exit */
      if ((j == dimension) && (nodes[j - 1] == max_nodes[j - 1]))
        break;
    }
  }

  coco_free_memory(x);
  coco_free_memory(func_values);
  if (number_of_constraints > 0 )
    coco_free_memory(cons_values);
  coco_free_memory(nodes);
  coco_free_memory(grid_step);
  coco_free_memory(max_nodes);
}

void round_vec(double population[], size_t dimention_size, const double lower_bounds[], const double upper_bounds[]){
  double y[16] = {0};
  double y_star;

  for(int i = 0; i < dimention_size; i++){
    if(upper_bounds[i] != 5){
      // 補助値を計算
      for(int  j = 0; j <= (int)upper_bounds[i]; j++){
          y[j] = j;
      }

      //連続値に最も近い補助値 y^* を求める
      double min_dist = fabs(y[0] - population[i]);
      y_star = y[0];
      for (int j = 1; j <= (int)upper_bounds[i]; j++) {
          double dist = fabs(y[j] - population[i]);
          if (dist <= min_dist) {
              min_dist = dist;
              y_star = y[j];
          }
      }
      //printf("%lf->%lf\n",population[i],y_star);
      population[i] = y_star;
    }
  }
}

/**
 * A random search algorithm that can be used for single- as well as multi-objective optimization. The
 * problem's initial solution is evaluated first.
 *
 * @param evaluate_func The function used to evaluate the objective function.
 * @param evaluate_cons The function used to evaluate the constraints.
 * @param dimension The number of variables.
 * @param number_of_objectives The number of objectives.
 * @param number_of_constraints The number of constraints.
 * @param lower_bounds The lower bounds of the region of interested (a vector containing dimension values).
 * @param upper_bounds The upper bounds of the region of interested (a vector containing dimension values).
 * @param number_of_integer_variables The number of integer variables (if > 0, all integer variables come
 * before any continuous ones).
 * @param max_budget The maximal number of evaluations.
 * @param random_generator Pointer to a random number generator able to produce uniformly and normally
 * distributed random numbers.
 */
void my_de_nopcm(evaluate_function_t evaluate_func,
                      evaluate_function_t evaluate_cons,
                      const size_t dimension,
                      const size_t number_of_objectives,
                      const size_t number_of_constraints,
                      const double *lower_bounds,
                      const double *upper_bounds,
                      const size_t number_of_integer_variables,
                      const size_t max_budget,
                      coco_random_state_t *random_generator) {
  size_t population_size = N; // 通常、10倍のサイズを使用
  double F = 0.5; // スケーリングファクター
  double CR = 0.9; // クロスオーバー確率
  double **population = (double**)malloc(population_size * sizeof(double*));
  double *functions_values = coco_allocate_vector(number_of_objectives);
  double *constraints_values = NULL;
  double **trial = (double**)malloc(population_size * sizeof(double*));
  double *mutate = coco_allocate_vector(dimension);
  double *rnd_vals = coco_allocate_vector(dimension);
  double *cons_values = NULL;
  size_t evaluation = 0;
  size_t i, j;
  int vector[3];
  double value_population[population_size];
  double value_trial = 0;
  FILE *fp;
  char titlestr[30] = TITLE;
  char num[30];
  sprintf(num, "%d", instance_cnt);
  strcat(titlestr,num);
  fp = fopen(titlestr, "w");
  // double bsf = 1000000000;

  for (i = 0; i < population_size; i++) {
        population[i] = coco_allocate_vector(dimension);
        trial[i] = coco_allocate_vector(dimension);
        if (population[i] == NULL) {
            printf("メモリの確保に失敗しました。\n");
        }
  }

  if (number_of_constraints > 0 )
    constraints_values = coco_allocate_vector(number_of_constraints);

  //initialization
  for (i = 0; i < population_size; i++) {
    // printf("bound:");
    // for (j = 0; j < dimension; j++) {
    //   printf("%lf ",upper_bounds[j]);
    // }
    // printf("\n");
    // printf("beforepop:");
    for (j = 0; j < dimension; j++) {
      double range = upper_bounds[j] - lower_bounds[j];
      //printf("%f:%f\n",upper_bounds[j],lower_bounds[j]);
      population[i][j] = lower_bounds[j] + coco_random_uniform(random_generator) * range;
      if (j < number_of_integer_variables){
        // if(REPAIR == 0){
        //   population[i][j] = floor(population[i][j] + 0.5);
        // }
      }
      // printf("%lf ",population[i][j]);
    }
    // printf("\n");
    // printf("afterpop:");
    if(REPAIR == 0){
      round_vec(population[i],dimension,lower_bounds,upper_bounds);
    }
    // for (j = 0; j < dimension; j++) {
    //   printf("%lf ",population[i][j]);
    // }
    // printf("\n");
    evaluate_func(population[i], functions_values);
    evaluation++;
    value_population[i] = functions_values[0];
    // if(fabs(value_population) < bsf){
    //   bsf = fabs(value_population);
    // }
  }
  // printf("\n");
  // for (i = 0; i < population_size; i++) {
  //   for (j = 0; j < dimension; j++) {
  //     printf("%lf ",population[i][j]);
  //   }
  //   printf("%s\n",coco_problem_get_name(PROBLEM));
  //   evaluate_func(population[i], functions_values);
  //   evaluation++;
  //   value_population = functions_values[0];
  //   printf(" value:%lf\n",functions_values[0]);
  //   printf("\n");
  // }

  while(evaluation  < max_budget){
    for(i = 0; i < population_size; i++){
      for(j = 0; j < dimension; j++){
        fprintf(fp,"%lf ",population[i][j]);
      }
      fprintf(fp,"\n");
    }
    fprintf(fp,"\n");
    for(i = 0; i < population_size; i++){
      //selection
      vector[0] = (int)(coco_random_uniform(random_generator)*N);
      do {
          vector[1] = (int)(coco_random_uniform(random_generator)*N);
      } while (vector[1] == vector[0]);

      do {
          vector[2] = (int)(coco_random_uniform(random_generator)*N);
      } while (vector[2] == vector[0] || vector[2] == vector[1]);
      //printf("vector:%d,%d,%d\n",vector[0],vector[1],vector[2]);
      //mutation
      for (j = 0; j < dimension; j++) {
        mutate[j] = population[vector[0]][j] + F * (population[vector[1]][j] - population[vector[2]][j]);
        if (mutate[j] < lower_bounds[j]){
            //v[i] = x[vector[0]][i] + ((double)rand() / RAND_MAX)*(0 - x[vector[0]][i]);  
            mutate[j] = (lower_bounds[j] + population[i][j]) / 2.0;
        }
        else if(mutate[j] > upper_bounds[j]){
            //v[i] = x[vector[0]][i] + ((double)rand() / RAND_MAX)*(l[cnt]-1 - x[vector[0]][i]);
            mutate[j] = (upper_bounds[j] + population[i][j]) / 2.0;
        }
      }

      //crossover
      int j_rand = (int)(coco_random_uniform(random_generator)*(int)dimension);// j_rand is a random index between 0 and dim-1

      // Generate random values between 0 and 1
      for (j = 0; j < dimension; j++) {
          rnd_vals[j] = coco_random_uniform(random_generator);
      }

      // Set rnd_vals[j_rand] to 0.0
      rnd_vals[j_rand] = 0.0;

      // Perform binomial crossover
      for (j = 0; j < dimension; j++) {
          if (rnd_vals[j] <= CR) {
              trial[i][j] = mutate[j];
          } else {
              trial[i][j] = population[i][j];
          }
      }
      if(REPAIR == 0){
        round_vec(trial[i],dimension,lower_bounds,upper_bounds);
      }
    }
    if(number_of_constraints > 0 ){
      evaluate_cons(population[i], cons_values);
    }
    //evaluation population
    for(i = 0; i < population_size; i++){
      evaluate_func(trial[i], functions_values);
      evaluation++;
      value_trial = functions_values[0];
      // if(value_trial < bsf){
      //   bsf = value_trial;
      // }
      // if(value_trial < 0){
      //   printf("bug\n");
      //   break;
      // }
      if(value_trial <= value_population[i]){
        for (j = 0; j < dimension; j++) {
          population[i][j] = trial[i][j];
          value_population[i] = value_trial;
          //printf("%lf ",population[i][j]);
        }
        //printf(" value:%lf\n",value_trial);
      }
      //printf("\n");
    }
    // if(bsf < 1e-8){
    //   printf("bsf:%lf\n",bsf);
    //   break;
    // }
  }
  instance_cnt++;
  fclose(fp);
  //printf("bsf:%lf\n",bsf);
  // printf("evaluation%ld\n",evaluation);
  for (i = 0; i < population_size; ++i) {
    coco_free_memory(population[i]);
  }
  free(population);
  coco_free_memory(functions_values);
  coco_free_memory(trial);
  coco_free_memory(mutate);
  coco_free_memory(rnd_vals);
  if (number_of_constraints > 0 )
    coco_free_memory(constraints_values);
}


/**
 * Allocates memory for the timing_data_t object and initializes it.
 */
static timing_data_t *timing_data_initialize(coco_suite_t *suite) {

	timing_data_t *timing_data = (timing_data_t *) coco_allocate_memory(sizeof(*timing_data));
	size_t function_idx, dimension_idx, instance_idx, i;

	/* Find out the number of all dimensions */
	coco_suite_decode_problem_index(suite, coco_suite_get_number_of_problems(suite) - 1, &function_idx,
			&dimension_idx, &instance_idx);
	timing_data->number_of_dimensions = dimension_idx + 1;
	timing_data->current_idx = 0;
	timing_data->output = (char **) coco_allocate_memory(timing_data->number_of_dimensions * sizeof(char *));
	for (i = 0; i < timing_data->number_of_dimensions; i++) {
		timing_data->output[i] = NULL;
	}
	timing_data->previous_dimension = 0;
	timing_data->cumulative_evaluations = 0;
	time(&timing_data->start_time);
	time(&timing_data->overall_start_time);

	return timing_data;
}

/**
 * Keeps track of the total number of evaluations and elapsed time. Produces an output string when the
 * current problem is of a different dimension than the previous one or when NULL.
 */
static void timing_data_time_problem(timing_data_t *timing_data, coco_problem_t *problem) {

	double elapsed_seconds = 0;

	if ((problem == NULL) || (timing_data->previous_dimension != coco_problem_get_dimension(problem))) {

		/* Output existing timing information */
		if (timing_data->cumulative_evaluations > 0) {
			time_t now;
			time(&now);
			elapsed_seconds = difftime(now, timing_data->start_time) / (double) timing_data->cumulative_evaluations;
			timing_data->output[timing_data->current_idx++] = coco_strdupf("d=%lu done in %.2e seconds/evaluation\n",
					timing_data->previous_dimension, elapsed_seconds);
		}

		if (problem != NULL) {
			/* Re-initialize the timing_data */
			timing_data->previous_dimension = coco_problem_get_dimension(problem);
			timing_data->cumulative_evaluations = coco_problem_get_evaluations(problem);
			time(&timing_data->start_time);
		}

	} else {
		timing_data->cumulative_evaluations += coco_problem_get_evaluations(problem);
	}
}

/**
 * Outputs and finalizes the given timing data.
 */
static void timing_data_finalize(timing_data_t *timing_data) {

	/* Record the last problem */
	timing_data_time_problem(timing_data, NULL);

  if (timing_data) {
  	size_t i;
  	double elapsed_seconds;
		time_t now;
		int hours, minutes, seconds;

		time(&now);
		elapsed_seconds = difftime(now, timing_data->overall_start_time);

  	printf("\n");
  	for (i = 0; i < timing_data->number_of_dimensions; i++) {
    	if (timing_data->output[i]) {
				printf("%s", timing_data->output[i]);
				coco_free_memory(timing_data->output[i]);
    	}
    }
  	hours = (int) elapsed_seconds / 3600;
  	minutes = ((int) elapsed_seconds % 3600) / 60;
  	seconds = (int)elapsed_seconds - (hours * 3600) - (minutes * 60);
  	printf("Total elapsed time: %dh%02dm%02ds\n", hours, minutes, seconds);

    coco_free_memory(timing_data->output);
    coco_free_memory(timing_data);
  }
}
