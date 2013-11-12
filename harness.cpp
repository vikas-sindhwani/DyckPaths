/**
 * @author pkambadu
 */

#include <iostream>
#include <fstream>
#include <functional>
#include <cstring>
#include <vector>
#include <utility>
#include <cstdio>
#include <omp.h>

#include <boost/math/distributions/chi_squared.hpp>

#include "dyck_options.hpp"
#include "generate_Galton_Watson.hpp"
#include "dyck_path.hpp"

typedef generate_Galton_Watson_t::graph_type graph_type;
typedef generate_Galton_Watson_t::vertex_name_map_t vertex_name_map_t;
typedef generate_Galton_Watson_t::vertex_distance_map_t vertex_weight_map_t;
typedef generate_Galton_Watson_t::vertex_distance_map_t edge_weight_map_t;
typedef boost::graph_traits<graph_type>::vertex_iterator vertex_iter_t;

static const double PI = 
    3.1415926535897932384626433832795028841971693993751058209749445923078164;

/**
 * The one annoying problem is that OutputIterator's do not define
 * a value type, so it may not always work! Hrmp!
 */ 
int generate_seed () {
  int r_int;
#pragma omp critical
  {
    std::ifstream rf("/dev/random", std::ios::binary);
    if(rf.is_open())rf.read(reinterpret_cast<char*>(&r_int),sizeof(r_int));
    else r_int = std::time(NULL);
  }

  return r_int;
}

graph_type generate_one_graph (const dyck_options_t& options,
                               int& gen_seed) {
  graph_type G;

  gen_seed = (0>options.seed)?generate_seed():options.seed;
  if (0 == strcmp("Poisson", options.gen_method.c_str())) {
    G = generate_Galton_Watson_t::generate_Poisson (options.n_max,
                                                    gen_seed,
                                                    options.lambda,
                                                    options.use_random_weights,
                                                    options.verbosity);
  } else if (0 == strcmp("Binomial", options.gen_method.c_str())) {
    G = generate_Galton_Watson_t::generate_Binomial (options.n_max,
                                                     gen_seed,
                                                     options.k,
                                                     options.p,
                                                     options.use_random_weights,
                                                     options.verbosity);
  } else if (0 == strcmp("Geometric", options.gen_method.c_str())) {
    G = generate_Galton_Watson_t::generate_Geometric(options.n_max,
                                                     gen_seed,
                                                     options.p,
                                                     options.use_random_weights,
                                                     options.verbosity);
  } else if (0 == strcmp("Binary-0-2", options.gen_method.c_str())) {
    G = generate_Galton_Watson_t::generate_Binary_0_2
                                             (options.n_max,
                                              gen_seed,
                                              options.use_random_weights,
                                              options.verbosity);
  } else if (0 == strcmp("Binary-0-1-2", options.gen_method.c_str())) {
    G = generate_Galton_Watson_t::generate_Binary_0_1_2
                                              (options.n_max,
                                               gen_seed,
                                               options.use_random_weights,
                                               options.verbosity);
  } else {
    std::cout << "Chosen method of graph generation not supported" 
              << std::endl;
  }
  return G;
}

void run_one_experiment (const dyck_options_t& options,
                         int& gen_seed,
                         int& mle_seed,
                         int& num_nodes_in_G,
                         int& vertex_to_find,
                         double& vertex_height,
                         double& mean_height) {
  graph_type G;
  vertex_name_map_t name_map;
  vertex_weight_map_t weight_map;

  do {
    /** Generate a graph */
    G = generate_one_graph (options, gen_seed);

    /** Get the number of nodes in the graph */
    vertex_iter_t v_begin, v_end;
    boost::tie(v_begin, v_end) = boost::vertices(G);
    num_nodes_in_G = std::distance(v_begin, v_end);
  } while (num_nodes_in_G < (options.n_min));

  name_map = boost::get(boost::vertex_name, G);
  weight_map = boost::get(boost::vertex_distance, G);

  /** 
   * If we want to identify the height of a random (non-root) vertex in the
   * tree, then we want generate a random number between [1,n), where n is 
   * the number of nodes in the generated tree.
   */
  vertex_to_find = -1;
  mle_seed = (0>options.mle_seed)?generate_seed():options.mle_seed;
  boost::mt19937 engine (mle_seed);
  boost::uniform_int<int> dist (1, num_nodes_in_G-1);
  vertex_to_find = dist(engine);

  /** Now, compute the Dyck path using depth-first search */
  std::vector<double> y_axis;
  vertex_height = dyck_path (G, name_map, weight_map, 
                             std::back_inserter(y_axis), 
                             vertex_to_find, 
                             options.verbosity);

  /** Compute the mean height of the tree */
  assert (y_axis.size() == 2*num_nodes_in_G);
  mean_height = 0.0;
  for (size_t i=0; i<y_axis.size(); ++i) mean_height += y_axis[i];
  mean_height /= (2*pow(num_nodes_in_G,1.5));

#pragma omp critical
  if (options.print_tree) {
    std::ostream* output; 
    if (0==strcmp("stdout", options.dot_out.c_str())) output = &(std::cout);
    else output = new std::ofstream(options.dot_out.c_str());

    boost::write_graphviz(*output, G, 
           boost::make_label_writer(boost::get(boost::vertex_name, G)),
           boost::make_label_writer(boost::get(boost::edge_weight, G)));

    if (0!=strcmp("stdout", options.dot_out.c_str())) delete output;
  }

#pragma omp critical
  if (options.print_path) {
    /** Write stuff out where needed */
    /** We need to ensure that the x-axis is such that the slope is always 45*/
    std::vector<double> x_axis(y_axis.size());
    x_axis[0] = 0;
    for (size_t i=1; i<x_axis.size(); ++i)
      x_axis[i] = x_axis[i-1] + std::abs(y_axis[i] - y_axis[i-1]);

    /** Remember, DO NOT PRINT THE LAST INDEX -- this is an artifact of DFS */
    std::ostream* output; 
    if (0==strcmp("stdout", options.dyck_out.c_str())) output = &(std::cout);
    else output = new std::ofstream(options.dyck_out.c_str());
    for (size_t i=0; i<(x_axis.size()-1); ++i) {
      *output << x_axis[i] << "  " << y_axis[i] << std::endl;
    }
    if (0!=strcmp("stdout", options.dyck_out.c_str())) delete output;
  }

}

int main (int argc, char** argv) {
  dyck_options_t options (argc, argv);
  if (options.exit_on_return) { return -1; }

  if (0<options.verbosity) options.pretty_print();

  std::vector<int> gen_seed_vec(options.num_trials);
  std::vector<int> mle_seed_vec(options.num_trials);
  std::vector<int> num_nodes_vec(options.num_trials);
  std::vector<int> vertex_to_find_vec(options.num_trials);
  std::vector<double> vertex_height_vec(options.num_trials);
  std::vector<double> mean_height_vec(options.num_trials);

#pragma omp parallel for num_threads(options.num_threads)
  for (int i=0; i<options.num_trials; ++i) {
    run_one_experiment (options, gen_seed_vec[i], mle_seed_vec[i], 
                        num_nodes_vec[i], vertex_to_find_vec[i], 
                        vertex_height_vec[i], mean_height_vec[i]);
  }

  double true_var;
  if (0 == strcmp("Poisson", options.gen_method.c_str())) 
    true_var=options.lambda;
  else if (0 == strcmp("Binomial", options.gen_method.c_str())) 
    true_var=options.k*options.p*(1-options.p);
  else if (0 == strcmp("Geometric", options.gen_method.c_str()))
    true_var=(1-options.p)/(options.p*options.p);
  else if (0 == strcmp("Binary-0-2", options.gen_method.c_str()))
    true_var=1;
  else if (0 == strcmp("Binary-0-1-2", options.gen_method.c_str()))
    true_var=2./3.;

  if (options.dump_numbers) {
    for (int i=0; i<options.num_trials; ++i) {
      double x_i = (1.0/(sqrt((double)num_nodes_vec[i])))*vertex_height_vec[i];
      printf (" %-8.3f", x_i);
    }
    printf ("\n");
    for (int i=0;i<options.num_trials;++i) {
      printf(" %-8.3f", mean_height_vec[i]);
    }
    printf ("\n");
  }

  if (options.measure_mle) {
    if (options.verbosity) {
      char header_string[1024];
      int h_count = sprintf (header_string, 
            "%9s %2s %10s %2s %7s %7s %7s %7s %12s %12s %12s", 
                "dist-type", "k", "p", "lm", "n-max", "n-min", 
                "n-tru", "node-#", "node-H-raw", "node-H-nrm" , "mean-H");
      if (1<options.verbosity)
   sprintf (&(header_string[h_count]), " %10s %10s", "gen-seed", "MLE-seed");
      std::cout << header_string << std::endl;
    }

    double sum_of_x_i_sqr = 0.0;
    double mean_of_x_i_sqr = 0.0;
    for (int i=0; i<options.num_trials; ++i) {
      double x_i = (1.0/(sqrt((double)num_nodes_vec[i])))*vertex_height_vec[i];
      sum_of_x_i_sqr += x_i*x_i;
      mean_of_x_i_sqr += x_i;

      if (options.verbosity) {
        char value_string[1024];
        int v_count = sprintf (value_string, 
               "%9s %2d %.4e %2d %7d %7d %7d %7d %.6e %.6e %.6e",
                               options.gen_method.c_str(),
                               options.k,
                               options.p,
                               options.lambda,
                               options.n_max,
                               options.n_min,
                               num_nodes_vec[i],
                               vertex_to_find_vec[i],
                               vertex_height_vec[i],
                               x_i,
                               mean_height_vec[i]);
        if (1<options.verbosity)
          sprintf (&(value_string[v_count]), " %10u %10u", 
                          gen_seed_vec[i], mle_seed_vec[i]);
        std::cout << value_string << std::endl;
      }
    }
    mean_of_x_i_sqr /= options.num_trials;
    mean_of_x_i_sqr *= mean_of_x_i_sqr;

    double mle_estimate_var = (2*options.num_trials)/sum_of_x_i_sqr;
    double ano_estimate_var = PI/(2*mean_of_x_i_sqr);

    std::cout << "True variance = " << true_var
              << " MLE = " << mle_estimate_var 
              << " Another = " << ano_estimate_var << std::endl;
    std::cout << " Mean of x_i^2 =  " << mean_of_x_i_sqr << std::endl;
  }

  if (options.test_confidence) {
    int num_batches = options.num_trials / options.batch_size;
    int num_contained = 0;

    int num_nodes = num_nodes_vec[0];
    bool experimental_error = false;
    for (int i=1; i<options.num_trials; ++i) {
      if (num_nodes != num_nodes_vec[i]) {
        experimental_error = true;
        break;
      }
    }

    if (experimental_error) {
      std::cout << "Run this experiment again, please" << std::endl;
      goto END;
    }

    std::vector<double> other_variances;
    if (0 == strcmp("Poisson", options.gen_method.c_str())) {
      other_variances.push_back(1 - 1./2.);
      other_variances.push_back(1 - 1./3.);
      other_variances.push_back(1 - 1./4.);
      other_variances.push_back(1 - 1./5.);
      other_variances.push_back(2);
    } else if (0 == strcmp("Binomial", options.gen_method.c_str())) {
      other_variances.push_back(1);
      other_variances.push_back(2);
      other_variances.push_back(2.0/3.0);
    } else if (0 == strcmp("Geometric", options.gen_method.c_str())) {
      other_variances.push_back(1);
      other_variances.push_back(1 - 1./2.);
      other_variances.push_back(1 - 1./3.);
      other_variances.push_back(1 - 1./4.);
      other_variances.push_back(1 - 1./5.);
    } else if (0 == strcmp("Binary-0-2", options.gen_method.c_str())) {
      other_variances.push_back(1 - 1./2.);
      other_variances.push_back(1 - 1./4.);
      other_variances.push_back(1 - 1./5.);
      other_variances.push_back(2);
      other_variances.push_back(2.0/3.0);
    } else if (0 == strcmp("Binary-0-1-2", options.gen_method.c_str())) {
      other_variances.push_back(1 - 1./2.);
      other_variances.push_back(1 - 1./3.);
      other_variances.push_back(1 - 1./4.);
      other_variances.push_back(1 - 1./5.);
      other_variances.push_back(1);
      other_variances.push_back(2);
    }
    std::vector<double> power_against_others (other_variances.size(), 0.0);

    boost::math::chi_squared chi_dist (2*options.batch_size);
    const double a = boost::math::quantile(chi_dist, 0.025);
    const double b = boost::math::quantile(chi_dist, 0.975);
    
    for (int i=0; i<options.num_trials; i+=options.batch_size) {
      double sum_of_x_j_sqr = 0.0;
      for (int j=0; j<options.batch_size; ++j) {
        double x_j = (1.0/(sqrt((double)num_nodes_vec[i+j])))*
                                        vertex_height_vec[i+j];
        sum_of_x_j_sqr += x_j*x_j;
      }
      
      const double A = a/sum_of_x_j_sqr;
      const double B = b/sum_of_x_j_sqr;
      if (A <= true_var && B >= true_var) ++num_contained;
      for (int l=0; l<other_variances.size(); ++l)
        if (A <= other_variances[l] && B >= other_variances[l]) 
          power_against_others[l]+=1;

      if (1<options.verbosity) {
        std::cout << "Lower = " << A << " Upper = " << B  
                  << " sum_of_x_j_sqr = " << sum_of_x_j_sqr << std::endl;
      }
    }
    const double contained = 
       100*static_cast<double>(num_contained)/static_cast<double>(num_batches);

    printf ("%12s %5d %3d %7.3f%%(tru-var=%4.3f) ", options.gen_method.c_str(),
                                                    num_batches,
                                                    options.batch_size,
                                                    contained,
                                                    true_var);
    for (int l=0; l<other_variances.size(); ++l) 
      printf (" %7.3f%%(var=%4.3f)", 100*(power_against_others[l]/num_batches),
                                     other_variances[l]);
    printf ("\n");
  }

END:

  return 0;
}
