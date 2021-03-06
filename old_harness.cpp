#include <iostream>
#include <cstdlib>
#include <vector>
#include <cmath>
#include <algorithm>
#include <fstream>
#include <string>

const int INVALID_NODE = -1;

/**
 * Errors:
 * (1) The root always has a left child.
 * (2) The contour plot repeats a zero at the end.
 * (3) Change it to handle n-ary trees.
 * (4) Include edge weights from an exponential random variable.
 * (5) Add a parameter that takes in the number of trees to generate.
 *
 * A function to generate a random binary tree using rand() function. This 
 * function is NOT reentrant (i.e., NOT thread safe).
 * @param[inout] preorder The vector containing the preorder encoding.
 * @param[inout] inorder The vector containing the inorder encoding.
 * @param[in] n Number of nodes in the tree
 *
 * This function is based on the method described by Devroye and Kruszewski 
 * (1996). Generate random binary trees with n nodes by generating a
 * real-valued random variable x in the unit interval [0,1), assigning the
 * first xn nodes (rounded down to an integer number of nodes) to the left
 * subtree, the next node to the root, and the remaining nodes to the right
 * subtree, and continuing recursively in each subtree. If x is chosen
 * uniformly at random in the interval, the result is the same as the random
 * binary search tree generated by a random permutation of the nodes, as any
 * node is equally likely to be chosen as root; however, this formulation
 * allows other distributions to be used instead. For instance, in the
 * uniformly random binary tree model, once a root is fixed each of its two
 * subtrees must also be uniformly random, so the uniformly random model may
 * also be generated by a different choice of distribution for x. As Devroye
 * and Kruszewski show, by choosing a beta distribution on x and by using an
 * appropriate choice of shape to draw each of the branches, the mathematical
 * trees generated by this process can be used to create realistic-looking
 * botanical trees.
 * Source: http://en.wikipedia.org/wiki/Random_binary_tree
 */
void generate_random_tree (std::vector<int>& preorder,
                           std::vector<int>& inorder,
                           std::vector<int>& contour,
                           const int n, 
                           const int start = 0) { 

  switch (n) {
    case 0: 
      {
        // Do nothing to the contour, inorder, and preorder. 
        // There is nothing to do.
        break; 
      }
    case 1: 
      { 
        // There was a leaf node here, so we increase the contour function.
        preorder.push_back (start); 
        inorder.push_back (start); 
        break;
      }
    default:
    {
  // Generate a random number 'x' between [0,1).
  const double x = static_cast<double>(rand())/static_cast<double>(RAND_MAX);

  // Assign floor (x*n) to the left subtree.
  const int n_left = floor (x*n);
  const int left_start = start;
  const int root = start + n_left;
  const int n_right = n - n_left - 1;
  const int right_start = start + n_left + 1;

  preorder.push_back (root);

  if (0<n_left) {
    contour.push_back (contour[contour.size()-1]+1);
    generate_random_tree (preorder, inorder, contour, n_left, left_start);
    contour.push_back (contour[contour.size()-1]-1);
  }

  inorder.push_back (root);

  if (0<n_right) {
    contour.push_back (contour[contour.size()-1]+1);
    generate_random_tree (preorder, inorder, contour, n_right, right_start);
    contour.push_back (contour[contour.size()-1]-1);
  }

  break;
    }
  }
}

/**
 * Reconstruct the tree from the preorder and inorder traversals.
 * @param[inout] current_preorder_index Position of the current root.
 * @param[in] inorder_first Iterator pointing to the first element.
 * @param[in] inorder_last Iterator pointing to the last element.
 * @param[inout] dot_file_stream Output stream.
 */
template <typename InputIterator>
int print_tree (InputIterator& current_preorder_index,
                InputIterator inorder_first,
                InputIterator inorder_last,
                std::ostream& dot_file_stream) {
  if (inorder_first != inorder_last) {

    // Get the current root
    const int current_root = *current_preorder_index++;

    // Find the position of the root to help process left and right subtrees.
    InputIterator root_position = 
      std::find (inorder_first, inorder_last, current_root);

    // Iterate the left subtree
    const int left_child = print_tree (current_preorder_index,
                                       inorder_first,
                                       root_position,
                                       dot_file_stream);

    // Iterate the right subtree
    const int right_child = print_tree (current_preorder_index,
                                        (root_position+1),
                                        inorder_last,
                                        dot_file_stream);

    // Print the edge
    if (INVALID_NODE != left_child) 
      dot_file_stream << current_root 
                      << "->" << left_child << "[color=red];" << std::endl;
    if (INVALID_NODE != right_child) 
      dot_file_stream << current_root 
                      << "->" << right_child << "[color=blue];" << std::endl;

    // return this node
    return current_root;
  } else {
    return INVALID_NODE;
  }
}


int main (int argc, char** argv) {
  if (4 > argc) {
    std::cout << "./harness <tree-size> <prng-seed> <pretty-print> [<filename>]" 
              << std::endl;
    return -3;
  }

  const int n = atoi (argv [1]);
  const int s = atoi (argv [2]);
  const int pretty_print = (0==strcmp(argv[3], "true")) ? true: false;

  std::vector<int> preorder;
  std::vector<int> inorder;
  std::vector<int> contour;

  srand (s);
  contour.push_back (0);
  generate_random_tree (preorder, inorder, contour, n);
  contour.push_back (0);

  if (pretty_print) { 
    if (5 != argc) {
      std::cout << "Print option used without output file name" << std::endl;
      return -4;
    }
    std::string dot_file_name (argv[4]);
    std::string plot_file_name (argv[4]);
    dot_file_name = dot_file_name + ".dot";
    plot_file_name = plot_file_name + ".plt";
    std::ofstream dot_file_stream (dot_file_name.c_str());
    std::ofstream plot_file_stream (plot_file_name.c_str());

    std::vector<int>::iterator preorder_iter = preorder.begin ();
    dot_file_stream << "digraph " << argv[4] << " { " << std::endl;
    print_tree 
      (preorder_iter, inorder.begin (), inorder.end (), dot_file_stream);
    dot_file_stream << "}" << std::endl;
    dot_file_stream.close();

    for (int i=0; i<contour.size(); ++i) 
      plot_file_stream << i << " " << contour[i] << std::endl;
    plot_file_stream.close();
  }

  return 0;
}
