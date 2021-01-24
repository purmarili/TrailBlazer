/******************************************************************************
 * File: Trailblazer.cpp
 *
 * Implementation of the graph algorithms that comprise the Trailblazer
 * assignment.
 */

#include "Trailblazer.h"
#include "TrailblazerGraphics.h"
#include "TrailblazerTypes.h"
#include "TrailblazerPQueue.h"
#include "random.h"

using namespace std;

struct Node{
	Loc coordinates;
	double weight;
	Color col;
	Node* parent;
};

/* Function: shortestPath
 * 
 * Finds the shortest path between the locations given by start and end in the
 * specified world.	 The cost of moving from one edge to the next is specified
 * by the given cost function.	The resulting path is then returned as a
 * Vector<Loc> containing the locations to visit in the order in which they
 * would be visited.	If no path is found, this function should report an
 * error.
 *
 * In Part Two of this assignment, you will need to add an additional parameter
 * to this function that represents the heuristic to use while performing the
 * search.  Make sure to update both this implementation prototype and the
 * function prototype in Trailblazer.h.
 */


int findset(int ind, Vector<int> & all, Vector<int> & size);
void unionset(int first, int second, Vector<int> & all, Vector<int> & size);
void makeGrey(Grid<Node*> & graph);
void dijkstra(Grid<double> & world, Grid<Node*> & graph, Loc start, Loc end, Vector<Loc> & result,
             double costFunction(Loc from, Loc to, Grid<double>& world), 
			 double heuristic(Loc start, Loc end, Grid<double>& world));

Vector<Loc> shortestPath(Loc start, Loc end, Grid<double>& world, 
						double costFunction(Loc one, Loc two, Grid<double>& world), 
						double heuristic(Loc start, Loc end, Grid<double>& world)){

	Vector<Loc> result;
	Grid<Node*> graph(world.numRows(), world.numCols());
	makeGrey(graph);
	graph[start.row][start.col]->col = YELLOW;
	colorCell(world, start, YELLOW);
	graph[start.row][start.col]->weight = 0;
	dijkstra(world, graph, start, end, result, costFunction, heuristic);
	for (int i=0; i<result.size()/2; i++){
		swap(result[i], result[result.size() - i - 1]);
	}

    return result;
}

Set<Edge> createMaze(int numRows, int numCols) {
	Set<Edge> allEdges, minTree;
	TrailblazerPQueue<Edge> pq;
	for (int i=0; i<numRows; i++){
		for (int j=0; j<numCols - 1; j++){
			Loc right = makeLoc(i, j+1);
			Loc down  = makeLoc(i+1, j);
			Loc start = makeLoc(i,   j);
			if (j+1 < numCols){ allEdges.add(makeEdge(start, right)); }
			if (i+1 < numRows){ allEdges.add(makeEdge(start, down));  }
		}
	}
	foreach (Edge e in allEdges){ pq.enqueue(e, randomReal(0, 1)); }
	Vector<int> all  (numRows * numCols);
	Vector<int> size (numRows * numCols);
	int count = numRows * numCols;
	for (int i=0; i<all.size(); i++){
		all[i] = i;
		size[i] = 0;
	}
	Edge e;
	while (count > 1){
		e = pq.dequeueMin();
		int start = e.start.row * numRows + e.start.col;
		int end = e.end.row * numRows + e.end.col;
		
		if (findset(start, all, size) != findset(end, all, size))
		{
			int first = findset(start, all, size);
			int second = findset(end, all, size);
			if (first != second){
				if (size[first] < size[second]) all[first] = second;
				else all[second] = first;
				size[first] += size[second];
				count--;
				minTree.add(e);
			}
		}
	}

    return minTree;
}

int findset(int ind, Vector<int> & all, Vector<int> & size)
{
	if (ind == all[ind]) return ind;
	all[ind] = findset(all[ind], all, size);
	return all[ind];
}

void unionset(int first, int second, Vector<int> & all, Vector<int> & size)
{
	first = findset(first, all, size);
	second = findset(second, all, size);
	if (first != second){
		if (size[first] < size[second])
			swap(first, second);
		all[second] = first;
		if (size[first] == size[second])
			all[first]++;
	}
}

void makeGrey(Grid<Node*> & graph)
{
	for (int i=0; i<graph.numRows(); i++){
		for (int j=0; j<graph.numCols(); j++){
			Node * newNode = new Node;
			newNode->col = GRAY;
			newNode->parent = NULL;
			Loc newloc;
			newloc.row = i;
			newloc.col = j;
			newNode->coordinates = newloc;
			graph[i][j] = newNode;
		}
	}
}

void dijkstra(Grid<double> & world, Grid<Node*> & graph, Loc start, Loc end, Vector<Loc> & result,
             double costFunction(Loc from, Loc to, Grid<double>& world), 
			 double heuristic(Loc start, Loc end, Grid<double>& world))
{
	Node* curr;
	TrailblazerPQueue<Node*> pq;
	graph[start.row][start.col]->weight = 0;
	pq.enqueue(graph[start.row][start.col], heuristic(start, end, world));
	while (!pq.isEmpty())
	{
		curr = pq.dequeueMin();
		colorCell(world, curr->coordinates, GREEN);
		curr->col = GREEN;
		if (curr->coordinates == end) break;
		for (int i=curr->coordinates.row - 1; i<curr->coordinates.row + 2; i++){
			for (int j=curr->coordinates.col - 1; j <curr->coordinates.col + 2; j++){
				if (i >= 0 && j >= 0 && i < graph.numRows() && j < graph.numCols()){
					double dis = curr->weight + costFunction(curr->coordinates, graph[i][j]->coordinates, world);
					if (graph[i][j]->col == GRAY){
						Loc help;
						help.row = i; help.col = j;
						colorCell(world, help, YELLOW);
						graph[i][j]->col = YELLOW;
						graph[i][j]->weight = dis;
						graph[i][j]->parent = curr;
						pq.enqueue(graph[i][j], graph[i][j]->weight + heuristic(graph[i][j]->coordinates, end, world));
					}else if (graph[i][j]->col == YELLOW && graph[i][j]->weight > dis){
						graph[i][j]->weight = dis;
						graph[i][j]->parent = curr;
						pq.decreaseKey(graph[i][j], dis + heuristic(graph[i][j]->coordinates, end, world));
					}
				}
			}
		}
	}

	while (curr != NULL)
	{
		result.push_back(curr->coordinates);
		curr = curr->parent;
	}
}