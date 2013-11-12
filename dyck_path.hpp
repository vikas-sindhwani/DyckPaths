#ifndef DYCK_PATH_HPP
#define DYCK_PATH_HPP

#ifdef _OPENMP
#include <omp.h>
#endif

#include <iterator>
#include <boost/graph/visitors.hpp>
#include <boost/graph/depth_first_search.hpp>
#include <boost/graph/graph_utility.hpp>

template<typename Tag, 
         typename NameMap,
         typename WeightMap,
         typename OutputIterator>
struct vertex_visitor_t : public boost::base_visitor<
             vertex_visitor_t<Tag, NameMap, WeightMap, OutputIterator> > {
  typedef Tag event_filter;

  vertex_visitor_t(std::string label,
                   NameMap& name_map,
                   WeightMap& weight_map,
                   OutputIterator y_axis,
                   double& y_val,
                   int vertex_to_find,
                   double& vertex_height,
                   int verbosity):label(label),
                                  name_map(name_map),
                                  weight_map(weight_map),
                                  y_axis(y_axis),
                                  y_val(y_val),
                                  vertex_to_find(vertex_to_find),
                                  vertex_height(vertex_height),
                                  verbosity(verbosity) {}

  double bin_op (double a,double b,boost::on_discover_vertex)const{return a+b;}
  double bin_op (double a,double b,boost::on_finish_vertex)const{return a-b;}

  void set_vertex_height (int vertex, 
                          double height, 
                          boost::on_discover_vertex){
    if (vertex_to_find == vertex) vertex_height = height;
  }
  void set_vertex_height (int vertex, 
                          double height, 
                          boost::on_finish_vertex){}

  /**
   * FIXME:
   * The Graph is being copied every single time. This is not correct.
   * Change the argument type to be "Graph&" or "const Graph&".
   */ 
  template <class Vertex, class Graph>
  void operator()(const Vertex& v, const Graph& G) {
#pragma omp critical
    if (1<verbosity)
      std::cout << label << " vertex : " << v 
                << " (label: " << boost::get(name_map, v)  << ")" << std::endl;

    y_val = bin_op (y_val, boost::get(weight_map, v), event_filter());
    *y_axis = y_val;
    set_vertex_height (v, y_val, event_filter());
    ++y_axis;
  }

  private:
  std::string label;
  NameMap& name_map;
  WeightMap& weight_map;
  OutputIterator y_axis;
  double& y_val;
  int vertex_to_find;
  double& vertex_height;
  int verbosity;
};

template <typename Graph, 
          typename NameMap,
          typename WeightMap,
          typename OutputIterator>
double dyck_path (Graph& G, 
                  NameMap& name_map,
                  WeightMap& weight_map,
                  OutputIterator y_axis,
                  int vertex_to_find,
                  int verbosity) {
  typedef vertex_visitor_t<boost::on_discover_vertex,
                           NameMap,
                           WeightMap,
                           OutputIterator> discover_visitor_t;
  typedef vertex_visitor_t<boost::on_finish_vertex,
                           NameMap,
                           WeightMap,
                           OutputIterator> finish_visitor_t;

  double y_init = 0;
  double vertex_height = -1;

  discover_visitor_t discover_visitor ("Discovered", 
                                       name_map,
                                       weight_map,
                                       y_axis, 
                                       y_init, 
                                       vertex_to_find, 
                                       vertex_height, 
                                       verbosity);
  finish_visitor_t finish_visitor ("Finished", 
                                    name_map,
                                    weight_map,
                                    y_axis, 
                                    y_init, 
                                    vertex_to_find, 
                                    vertex_height, 
                                    verbosity);
  boost::depth_first_search(G, 
         boost::visitor(
           boost::make_dfs_visitor(
             boost::make_list(discover_visitor, finish_visitor))));

  return vertex_height;
}

#endif // DYCK_PATH_HPP
