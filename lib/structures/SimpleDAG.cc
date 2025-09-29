/******************************************************************************
 * @file SimpleDAG.cc
 * @author Steinkopp:Felix
 * @version 0
 * @brief Simple Directed Acyclic Graph Class. Supports topological sorting
 ******************************************************************************/

#include "structures/SimpleDAG.h"

#include <algorithm>
#include <stdexcept>

namespace {
using uint = ::lcs_solver::structures::SimpleDAG::uint;
}

namespace lcs_solver::structures {
// constructor creates DAG with numVertices vertices 0, 1, ..., numVertices and with zero edges
SimpleDAG::SimpleDAG(uint numVertices)
    : numOfVertices(numVertices),
      numOfEdges(0),
      adjList(AdjacencyList(numOfVertices)),
      reverseAdjList(AdjacencyList(numOfVertices)) {

}

// method to add edge to the graph
void SimpleDAG::addEdge(uint u, uint v) {
  adjList[u].push_back(v);
  reverseAdjList[v].push_back(u);
  numOfEdges++;
}

// method to reverse all edges
void SimpleDAG::reverseAllEdges() {
  std::swap(adjList, reverseAdjList);
}

// getter for number of vertices
SimpleDAG::uint SimpleDAG::getNumOfVertices() const {
  return numOfVertices;
}

// getter for number of edges
SimpleDAG::uint SimpleDAG::getNumOfEdges() const {
  return numOfEdges;
}

// Function to perform a topological sort on the DAG
std::vector<SimpleDAG::uint> SimpleDAG::topologicalSort() const {
  std::vector<uint> sortedVertices;
  std::vector<bool> visited(numOfVertices, false);

  // Helper DFS function using SimpleDAGIterator_DFS
  auto dfsVisit = [&](uint startVertex) {
    SimpleDAGIterator_DFS dfsIterator(*this, startVertex);

    uint *vertexPtr;
    while ((vertexPtr = dfsIterator.next()) != nullptr) {
      uint vertex = *vertexPtr;
      if (!visited[vertex]) {
        visited[vertex] = true;
        sortedVertices.push_back(vertex);
      }
    }
  };

  // Perform DFS on unvisited vertices
  for (uint vertex = 0; vertex < numOfVertices; ++vertex) {
    if (!visited[vertex]) {
      dfsVisit(vertex);
    }
  }

  // Reverse the order to get the topological sorting
  std::reverse(sortedVertices.begin(), sortedVertices.end());

  // Ensure the graph is a DAG (no cycles)
  if (sortedVertices.size() != numOfVertices) {
    throw std::runtime_error(
        "The graph contains cycles and cannot be topologically sorted.");
  }

  return sortedVertices;
}
const SimpleDAG::AdjacencyList &SimpleDAG::getAdjacencyList() const {
  return adjList;
}

//=== SimpleDAGIterator_DFS ===================================================
SimpleDAGIterator_DFS::SimpleDAGIterator_DFS(const SimpleDAG &g,
                                             uint startVertex)
    : graph(g), visited(g.getNumOfVertices(), false) {
  dfsStack.push(startVertex);
  visited[startVertex] = true;
}

SimpleDAG::uint *SimpleDAGIterator_DFS::next() {
  if (dfsStack.empty()) {
    return nullptr;
  }

  uint &currentVertex = dfsStack.top();
  dfsStack.pop();

  // Add unvisited neighbors to the stack
  auto adjList = graph.getAdjacencyList();
  for (uint neighbor : adjList[currentVertex]) {
    if (!visited[neighbor]) {
      dfsStack.push(neighbor);
      visited[neighbor] = true;
    }
  }

  return &currentVertex;
}
}