#ifndef LCS_SOLVER_STRUCTURES_SIMPLE_DAG_HPP
#define LCS_SOLVER_STRUCTURES_SIMPLE_DAG_HPP

#include <cstddef>
#include <stack>
#include <vector>

namespace lcs_solver::structures {

/*******************************************************************************
 * @brief structure for working with directed acyclic graphs
 ******************************************************************************/
class SimpleDAG {
 public:
  using uint = std::size_t;
  using AdjacencyList = std::vector<std::vector<uint>>;

  explicit SimpleDAG(uint numVertices); ///< Constructor

  void addEdge(uint u, uint v); ///< Method to add edge to the graph
  void reverseAllEdges();       ///< Method to reverse all edges
  [[nodiscard]] std::vector<uint> topologicalSort() const; ///< Function to perform a topological sort on the DAG

  [[nodiscard]] uint getNumOfVertices() const;  // Getter for #Vertices
  [[nodiscard]] uint getNumOfEdges() const;     // Getter for #Edgess
  [[nodiscard]] const AdjacencyList& getAdjacencyList() const;

 private:
  uint numOfVertices;             // Number of vertices in the graph
  uint numOfEdges;                // Number of edges in the graph
  AdjacencyList adjList;          // Main adjacency list
  AdjacencyList reverseAdjList;   // Reversed adjacency list
};

/*******************************************************************************
 * @brief Class for traversing a `SimpleDAG` in a DFS
 ******************************************************************************/
class SimpleDAGIterator_DFS {
  public:
  using uint = std::size_t;

  SimpleDAGIterator_DFS(const SimpleDAG &g, uint startVertex);
  uint *next();

 private:
  const SimpleDAG &graph;
  std::vector<bool> visited;
  std::stack<uint> dfsStack;
};

} //end of namespace
#endif /* LCS_SOLVER_STRUCTURES_SIMPLE_DAG_HPP */
