///******************************************************************************
// * @file LCS_STD_DAG.cpp
// * @author Steinkopp:Felix
// * @version 0
// * @brief Implementation of an algorithm for LCS2_Classic
// * @details Time Complexity:
// *          Space Complexity:
// *****************************************************************************/
//
//
//#include <iostream>
//#include <algorithm>        // Use of std::reverse_copy (a deque is slower than reversing a vector)
//#include <unordered_set>    // Set S represent visited nodes: S = {(i,j) | DFS was there}
//
//#include "../../../include/algorithms/LCS/LCS_STD_2_DAG.hpp"
//
//namespace lcs_solver::algorithms::lcs {
//LCS_FL_DAG::LCS_FL_DAG(std::string &seq1, std::string &seq2, int &len1, int &len2)
//        :v(seq1),w(seq2),m(len1),n(len2)
//{
//    m = v.size();
//    n = w.size();
//}
//
///**
// * @brief Perform preprocessing steps for finding the Longest Common Subsequence (LCS). Time Complexity: O(n*m)
// *
// * This function performs the necessary preprocessing steps to find the LCS between two sequences, v and w.
// * The steps include filling the dynamic programming table (2D matrix), building a Directed Acyclic Graph (DAG)
// * representing the LCS paths, trimming the graph to remove unimportant nodes, and initializing the necessary
// * variables for generating LCS strings.
// */
//void LCS_FL_DAG::doPreprocessing()
//{
//    dp = std::vector<std::vector<int>>(m + 1, std::vector<int>(n + 1, 0));
//
//    // Step 1: Fill the dynamic programming table (2D matrix)
//    for (int i = 1; i <= m; ++i) {
//        for (int j = 1; j <= n; ++j) {
//            if (v[i - 1] == w[j - 1]) {
//                dp[i][j] = dp[i - 1][j - 1] + 1;
//            }
//            else {
//                dp[i][j] = std::max(dp[i - 1][j], dp[i][j - 1]);
//            }
//        }
//    }
//
//    // Step 2: Build the Directed Acyclic Graph (DAG) representing the LCS paths
//    // Add edges to the graph based on the values in the dp table.
//    // Each position (i, j) in the dp table corresponds to a node in the graph.
//    // Open Task: Optimizations for the average case
//    // for (int i = 1; i <= m; ++i) {
//    //     for (int j = 1; j <= n; ++j) {
//    //         std::pair<int, int> start = std::make_pair(i, j);
//    //         if (dp[i][j] == dp[i][j-1]) {
//    //             std::pair<int, int> end   = std::make_pair(i, j-1);
//    //             graph[start].push_back(end);
//    //         }
//    //         if (dp[i][j] == dp[i-1][j]) {
//    //             std::pair<int, int> end   = std::make_pair(i-1, j);
//    //             graph[start].push_back(end);
//    //         }
//    //         if (v[i-1] == w[j-1]) {
//    //             std::pair<int, int> end   = std::make_pair(i-1, j-1);
//    //             graph[start].push_back(end);
//    //         }
//    //     }
//    // }
//    for (int i = m; 1 < i; --i) {
//        for (int j = n; 1 < j; --j) {
//            std::pair<int, int> start = std::make_pair(i, j);
//            graph[start].push_back(std::make_pair(i-1, j  ));
//            graph[start].push_back(std::make_pair(i  , j-1));
//            graph[start].push_back(std::make_pair(i-1, j-1));
//        }
//    }
//    //graph[std::make_pair(-1, -1)].push_back(std::make_pair(n,m));
//
//    // Step 2.1 Save Source Nodes
//    // Open Task: Set sources to (n,m) if v[m-1]==w[n-1] or search for important nodes of (m,n)
//    printDataStructure("Before Step 3");
//
//    // Step 3: Trim the graph by removing unimportant nodes
//    // Use a DFS search and keep nodes in the set { (i,j) ⊂ [m]x[n] : v(i)==w(j) }
//    stack = std::stack<std::pair<int, int>>();                              // New (empty) stack for DFS
//    auto visited = std::unordered_set<std::pair<int, int>, pair_hash>();    // New (empty) visited set for DFS
//
//    // Start the DFS from the bottom-right node (position (m, n))
//    std::pair<int, int> start = std::make_pair(m, n);
//    //std::pair<int, int> start = std::make_pair(-1, -1);
//    stack.push(start);
//    while (!stack.empty())
//    {
//        std::pair<int,int> node = stack.top();
//        stack.pop();
//
//        // Ensure that each node is trimmed only once
//        if (visited.find(node) == visited.end()) {
//            trim(node, graph); // trims unimportant nodes and modifies graph[node]
//            visited.insert(node);
//        }
//
//        // Get all important neighbors of the current node and push them to the stack for further processing
//        for(auto i = graph[node].begin(); i != graph[node].end();i++)
//            if(visited.find(*i)==visited.end())
//                stack.push(*i);
//    }
//
//    // Step 4: Initialize the stacks and solution variables for generating LCS strings (in linear time)
//    resetQuerying();
//}
//
///**
// * @brief Resets `stack`, `nRightRightPos`, `currentLCS` and `solution` for new queries after preprocessing
// */
//void LCS_FL_DAG::resetQuerying()
//{
//    // Clean up variables used in pre-processing and genNextSol()
//    stack = std::stack<std::pair<int, int>>();
//    nRightRightPos = std::stack<int>();
//    currentLCS = "";
//    // Setup DFS for the DAG (without a visited array)
//    stack.push(std::make_pair(m,n));
//    //stack.push(std::make_pair(-1,-1));
//    nRightRightPos.push(0);
//}
//
//int LCS_FL_DAG::getSolLength() const
//{
//    return (dp.empty()) ? -1 : dp[m][n];
//}
//
///**
// * @brief Generate the next LCS solution.
// *
// * This function generates the next solution in the Longest Common Subsequence (LCS) problem. It pops a node from
// * the stack representing the current position in the LCS graph and appends the corresponding character to the current
// * LCS string. It continues this process until a complete LCS solution is generated, or until all possible solutions
// * have been generated.
// *
// * @return The next LCS solution as a string, or an empty string if no more solutions are available.
// */
//std::string LCS_FL_DAG::genNextSol()
//{
//    // Generate the next LCS solution (LCS::initStacksAndSol was called in the past)
//    while(!stack.empty()){
//        std::pair<int, int> & node = stack.top();
//        int k = nRightRightPos.top(); // Number of right chars form right to left in the final solution for the next lcs string
//        stack.pop();
//        nRightRightPos.pop();
//
//        if(v[node.first-1] == w[node.second-1]){ // Because (m,n) could be an unimportant node, but still be in the DGA.
//            currentLCS.resize(k);           // CurrentLCS = w[m_lcs] .. w[m_k] with m_lcs >= m_k
//            currentLCS += v[node.first-1];  // Add char corresponding to a node with v[i]==w[j]
//            k++;
//
//            // Found a termination node and a new lcs
//            if(graph[node].empty() && (int) currentLCS.size()==dp[m][n]){
//                std::string lcs;
//                std::reverse_copy(currentLCS.begin(), currentLCS.end(), std::back_inserter(lcs));
//                return lcs;
//            }
//        }
//
//        // Add children of node for getting the next symbols in the future
//        for(auto i : graph[node]){
//            stack.push(i);
//            nRightRightPos.push(k);
//        }
//    }
//
//    return std::string();
//}
//
///**
// * @brief Print the data structure used in the LCS algorithm.
// *
// * This function prints the data structure used in the Longest Common Subsequence (LCS) algorithm. It outputs the
// * dynamic programming table and the graph adjacency list, providing a visual representation of the underlying
// * data used in the LCS computation.
// */
//void LCS_FL_DAG::printDataStructure(std::string msg)
//{
//    if(!msg.empty()){
//        std::cout << "[Info]: " << msg << std::endl;
//    }
//
//    if(!dp.empty()){
//        // Print dp row by row
//        for (const auto& row : dp) {
//            for (const auto& element : row) {
//                std::cout << element << " ";
//            }
//            std::cout << std::endl;
//        }
//    }
//    else{
//        std::cout << "[Info]: dp is empty" << std::endl;
//    }
//
//    if(!graph.empty()){
//        for (const auto& pair : graph) {
//            std::cout << "(" << pair.first.first << ", " << pair.first.second << ") : ";
//            for (const auto& edge : pair.second) {
//                std::cout << "(" << edge.first << ", " << edge.second << "), ";
//            }
//            std::cout << std::endl;
//        }
//    }
//    else{
//        std::cout << "[Info]: adjacency list is empty" << std::endl;
//    }
//    std::cout << std::endl;
//}
//
///**
// * @brief Simplifies the directed acyclic graph by trimming nodes in such a way that querying in O(n) time is is possible
// *
// * This function trims the Directed Acyclic Graph (DAG) representing the LCS paths by removing nodes that do not
// * contribute to the LCS solution (meaning nodes (i,j) with v[i]!=w[j]). It performs a Depth-First Search (DFS)
// * traversal starting from the specified start node and keeps track of the visited nodes. In this DFS are
// * contributing nodes are leavees.
// *
// * @param start The start node for the DFS traversal.
// * @param graph The graph representing the LCS paths.
// */
//void LCS_FL_DAG::trim(std::pair<int, int>& start, std::unordered_map<std::pair<int, int>, std::list<std::pair<int, int>>, pair_hash>& graph)
//{
//    std::stack<std::pair<int, int>> st;  // Stack for DFS
//    st.push(start);                      // Push the start node
//    std::unordered_set<std::pair<int, int>, pair_hash> visited_trim;
//
//    // Decl. of resulting neighbors for the `start` node: { (i,j) ⊂ [m]x[n] : v(i)==w(j)}
//    std::list<std::pair<int, int>> candidates = std::list<std::pair<int, int>> ();
//
//    // DFS search
//    while (!st.empty()) {
//        std::pair<int, int> & node = st.top();
//        st.pop();
//
//        // Skip nodes representing an empty string
//        if(node.first==0 || node.second==0)
//            continue;
//
//        // Check, if `node` is a candidate to simplify the adjacency list
//        bool nodeWasVisited = visited_trim.find(node) == visited_trim.end();
//        bool isNotStartNode = node!=start;
//        bool isOnLCSPath1 = isMatch(start) && isMatch(node)
//                            && dp[node.first][node.second] + 1 == dp[start.first][start.second];
//        bool isOnLCSPath2 = !(isMatch(start)) && isMatch(node)
//                            && dp[node.first-1][node.second-1] + 1 == dp[node.first][node.second];
//        if (nodeWasVisited && isNotStartNode) {    // first time dfs visit?
//            if(isOnLCSPath1 || isOnLCSPath2){
//                //??? below this line code is very bad
//                for(int i=node.first-1;dp[start.first][start.second]-1==dp[i][node.second];i--)
//                    if(v[i-1]==w[node.second-1])
//                        candidates.push_back(std::make_pair(i,node.second));
//                for(int j=node.second-1;dp[start.first][start.second]-1==dp[node.second][j];j--)
//                    if(v[node.first-1]==w[j-1])
//                        candidates.push_back(std::make_pair(node.second,j));
//                //???? above this line code is very bad
//                candidates.push_back(node);
//                visited_trim.insert(node);
//                continue;
//            }
//            visited_trim.insert(node);
//        }
//
//        // Get all neighbors of node and push them onto the stack `st` for further processing
//        for (auto i : graph[node])
//            if (visited_trim.find(i) == visited_trim.end()) // if the neighbor hasn't been visited yet
//                st.push(i);
//    }
//
//    // Overwrite the list with the found neighbors of `start`
//    graph[start] = candidates;
//}
//
//bool LCS_FL_DAG::isMatch(int i, int j){
//    return v[i-1]==w[j-1];
//}
//
//bool LCS_FL_DAG::isMatch(std::pair<int,int> &p){
//    return v[p.first-1]==w[p.second-1];
//}
//
//// void LCS::solve()
//// {
////     // Reconstruct the LCS from the table
////     std::vector<std::string> currentLCS;
////     std::function<void(int, int)> backtrack = [&](int i, int j) {
////         if (i == 0 || j == 0) {
////             std::string lcs;
////             for (const std::string& str : currentLCS) {
////                 lcs += str;
////             }
////             std::reverse(lcs.begin(), lcs.end());
////             solution.push_back(lcs);
////             return;
////         }
////         if (v[i - 1] == w[j - 1]) {
////             currentLCS.push_back(std::string(1, v[i - 1]));
////             backtrack(i - 1, j - 1);
////             currentLCS.pop_back();
////         } else {
////             if (dp[i - 1][j] >= dp[i][j - 1]) {
////                 backtrack(i - 1, j);
////             }
////             if (dp[i][j - 1] >= dp[i - 1][j]) {
////                 backtrack(i, j - 1);
////             }
////         }
////     };
////     backtrack(m, n);
//// }
//}