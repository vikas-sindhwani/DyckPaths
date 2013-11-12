#ifndef GENERATE_GALTON_WATSON_HPP
#define GENERATE_GALTON_WATSON_HPP

#include <cstdlib>
#include <string>
#include <queue>
#include <algorithm>
#include <numeric>
#include <cmath>

#ifdef _OPENMP
#include <omp.h>
#endif

#include <boost/random.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/graphviz.hpp>
#include <boost/math/special_functions/binomial.hpp>
#include <boost/math/special_functions/factorials.hpp>

/**
 * Given:
 * (1) Distribution: Poisson, Negative Geometric, Binomial
 * (2) Number of nodes for the conditioned Galton-Watson process.
 * (3) Parameters needed for the distributions:
 *    (a) Probability of success (p)
 *    (b) Lambda if using Poisson distribution
 *    (c) 'k', the total number of trials if we are using Binomial
 *    (d) Seed value for the random number
 */
struct generate_Galton_Watson_t {
  typedef boost::mt19937 engine_type;
  typedef boost::binomial_distribution<int, double> binomial_type;
  typedef boost::uniform_real<double> uni_real_type;
  typedef boost::uniform_int<int> uni_int_type;

  typedef boost::adjacency_list
                   <boost::listS,
                    boost::vecS,
                    boost::directedS,
                    boost::property<boost::vertex_name_t, std::string,
                    boost::property<boost::vertex_distance_t, double> >,
                    boost::property<boost::edge_weight_t, double> > graph_type;
    typedef boost::property_map<graph_type, boost::vertex_name_t>::type 
                                                vertex_name_map_t;
    typedef boost::property_map<graph_type, boost::vertex_distance_t>::type 
                                                vertex_distance_map_t;
    typedef boost::property_map<graph_type, boost::edge_weight_t>::type
                                                edge_weight_map_t;
    typedef boost::graph_traits<graph_type>::vertex_descriptor vertex_t;

  private:

  struct always_1_t {
    typedef double result_type;
    result_type operator()() { return 1.0; }
  };

  struct Poisson_computer_t {
    int lambda;
    Poisson_computer_t (int lambda) : lambda(lambda) {}

    double operator()(int i) const { 
      return (pow(lambda,i)*exp(-lambda)) / boost::math::factorial<double>(i);
    }
  };

  struct geometric_computer_t {
    double p;
    geometric_computer_t (double p) : p(p) {}

    double operator()(int i) const { return pow((1-p), i) * p; }
  };

  struct binomial_computer_t {
    int k;
    double p;
    binomial_computer_t (int k, double p) : k(k), p(p) {}

    double operator()(int i) const { 
      return
      boost::math::binomial_coefficient<double>(k,i)*pow(p,i)*pow((1-p),(k-i));
    }
  };

  struct binary_0_2_computer_t {
    double operator()(int i) const { 
      if (0==i || 2==i) return 1./2.;
      else return 0;
    }
  };

  struct binary_0_1_2_computer_t {
    double operator()(int i) const { 
      if (0==i || 1==i || 2==i) return 1./3.;
      else return 0;
    }
  };

  /**
   * Thin wrapper around boost::uniform<> so that it conforms to random_shuffle.
   */
  struct random_functor_t : std::unary_function<int, int> {
    engine_type engine;
    
    random_functor_t (engine_type engine) : engine(engine) {}

    int operator()(const int& ULIMIT) {
      uni_int_type dist (0, ULIMIT-1);
      return dist(engine);
    }
  };

  static void rotate_to_get_tree (std::vector<int>& psi_vec) {
    /** Get the vector S for the queues in the tree */
    size_t n = psi_vec.size();
    std::vector<int> S(n);

    for (size_t i=0; i<n; ++i) S[i] = ((0==i)?1:(S[i-1]))+(psi_vec[i]-1);

    int min_ele=std::numeric_limits<int>::max();
    int min_index=-1;
    for (size_t i=0; i<n; ++i)
      if (min_ele>S[i]) { min_ele=S[i]; min_index=i; }

    std::rotate(psi_vec.begin(), psi_vec.begin()+min_index+1, psi_vec.end());
  }

  template <typename WDistribution>
  static graph_type generate_impl (int n, 
                                   const std::vector<int>& psi_vec,
                                   WDistribution w_prng,
                                   int verbosity) {
    graph_type G;
    vertex_name_map_t name_map = boost::get (boost::vertex_name, G);
    vertex_distance_map_t dist_map = boost::get (boost::vertex_distance, G);
    edge_weight_map_t weight_map = boost::get (boost::edge_weight, G);

    /** Set up a queue for the children and add the root vertex */
    vertex_t root = boost::add_vertex(G);
    boost::put (name_map, root, boost::lexical_cast<std::string>(0));
    boost::put (dist_map, root, 0);

    std::queue<vertex_t> vertex_queue;
    vertex_queue.push(root);

    /** Count the number of children left */
    int n_left = (n-1);

    /** Iterate till we are left with nothing or we have maxed out children */
    while ((0<=n_left) && (false==vertex_queue.empty())) {
      const vertex_t parent = vertex_queue.front();
      vertex_queue.pop();

#pragma omp critical
      if (1<verbosity) std::cout << "Popped " << parent << std::endl;

      int num_children = psi_vec[parent];
#pragma omp critical
      if (1<verbosity) std::cout << "num children of " << parent << " = " 
                               << num_children << std::endl;

      for (int child_number=0; child_number<num_children; ++child_number) {
        const vertex_t child = boost::add_vertex(G);
        vertex_queue.push(child);
        boost::put (name_map, 
                    child, 
                    boost::get (name_map, parent) +
                    boost::lexical_cast<std::string>(child_number));

        /** 
         * Duplicating a map here, but I think it's fine. We need the vertex
         * weights to easily be able to get the Dyck/Harris paths.
         */ 
        double edge_weight = w_prng();
        boost::put (dist_map, child, edge_weight);
        boost::put (weight_map, 
                    (boost::add_edge (parent, child, G)).first,
                    edge_weight);

        --n_left;
      }

#pragma omp critical
      if (1<verbosity) std::cout << "N Left = " << n_left << std::endl;
    }

    return G;
  }

  template <typename ProbabilityComputer>
  static graph_type generate_tree (int n,
                                   int seed,
                                   ProbabilityComputer compute,
                                   bool use_random_weights,
                                   int verbosity) {

    std::vector<int> N_j_vec;
    engine_type engine(seed);
    int j_times_N_j, j;

    do {

      j_times_N_j = j = 0;
      int sum_of_N_i = 0;
      int i = 0;
      double sum_of_p_i = 0;
      N_j_vec.resize(0);
      while (0 != (n-sum_of_N_i)) {
        const double p_i = compute(i);
        engine_type new_engine(engine());
        binomial_type dist((n-sum_of_N_i), p_i/(1-sum_of_p_i));
        const int N_j = dist(new_engine);
        N_j_vec.push_back (N_j);
        sum_of_N_i += N_j;
        sum_of_p_i += p_i;
        j_times_N_j += j*N_j;
        ++i;
        ++j;
      }

      if (2 < verbosity) std::cout << "j*N_j = " << j_times_N_j << std::endl;

    } while (j_times_N_j != (n-1));

    std::vector<int> psi_vec (n, 0);
    int offset = 0;
    for (int j=0; j<N_j_vec.size(); ++j) {
      for (int i=0; i<N_j_vec[j]; ++i) psi_vec[offset+i] = j;
      offset += N_j_vec[j];
    }

    random_functor_t shuffle_prng (engine);
    std::random_shuffle (psi_vec.begin(), psi_vec.end(), shuffle_prng);

    if (2<verbosity) {
      for (int i=0; i<n; ++i) std::cout<<psi_vec[i]<<" "; 
      std::cout<<std::endl;
    }

    rotate_to_get_tree (psi_vec);

    if (2<verbosity) {
      for (int i=0; i<n; ++i) std::cout<<psi_vec[i]<<" "; 
      std::cout<<std::endl;
    }

    if (use_random_weights) {
      uni_real_type weight_dist;
      boost::variate_generator<engine_type, uni_real_type> 
                                 w_prng (engine, weight_dist);
      return generate_impl (n, psi_vec, w_prng, verbosity);
    } else {
      return generate_impl (n, psi_vec, always_1_t(), verbosity);
    }
  }

  public:

  static graph_type generate_Poisson (int n,
                                      int seed,
                                      int lambda,
                                      bool use_random_weights,
                                      int verbosity) {
    return generate_tree (n, seed, Poisson_computer_t(lambda), 
                          use_random_weights, verbosity);
  }

  static graph_type generate_Geometric (int n,
                                        int seed,
                                        double p,
                                        bool use_random_weights,
                                        int verbosity) {
    return generate_tree (n, seed, geometric_computer_t(p), 
                          use_random_weights, verbosity);
  }

  static graph_type generate_Binomial (int n,
                                       int seed,
                                       int k,
                                       double p,
                                       bool use_random_weights,
                                       int verbosity) {
    return generate_tree (n, seed, binomial_computer_t(k, p), 
                          use_random_weights, verbosity);
  }

  static graph_type generate_Binary_0_2 (int n,
                                         int seed,
                                         bool use_random_weights,
                                         int verbosity) {
    return generate_tree (n, seed, binary_0_2_computer_t(), 
                          use_random_weights, verbosity);
  }

  static graph_type generate_Binary_0_1_2 (int n,
                                           int seed,
                                           bool use_random_weights,
                                           int verbosity) {
    return generate_tree (n, seed, binary_0_1_2_computer_t(), 
                          use_random_weights, verbosity);
  }
};

#endif // GENERATE_GALTON_WATSON_HPP
