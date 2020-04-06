#include <Rcpp.h>
using namespace Rcpp;

#include <iostream>
#include <queue>

using namespace std;

// [[Rcpp::export]]
bool isBipartite(IntegerMatrix G,int V, int src)
{

  IntegerVector colorArr(V);
  for (int i = 0; i < V; ++i)
    colorArr[i] = -1;

  // Assign first color to source
  colorArr[src] = 1;

  // Create a queue (FIFO) of vertex
  // numbers and enqueue source vertex
  // for BFS traversal
  queue <int> q;
  q.push(src);

  // Run while there are vertices
  // in queue (Similar to BFS)
  while (!q.empty())
  {
    int u = q.front();
    q.pop();

    // Return false if there is a self-loop
    if (G(u,u) == 1)
      return false;

    // Find all non-colored adjacent vertices
    for (int v = 0; v < V; ++v)
    {
      // An edge from u to v exists and
      // destination v is not colored
      if ((G(u,v)==1) & (colorArr[v] == -1))
      {
        // Assign alternate color to this adjacent v of u
        colorArr[v] = 1 - colorArr[u];
        q.push(v);
      }

      // An edge from u to v exists and destination
      // v is colored with same color as u
      else if ((G(u,v)==1) & (colorArr[v] == colorArr[u]))
        return false;
    }
  }

  // If we reach here, then all adjacent vertices can be colored with alternate color
  return true;
}

