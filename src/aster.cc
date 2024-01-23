/*
Part of Aster Transcript Assembler
(c) 2024 by Xiaofei Carl Zang, Mingfu Shao, and The Pennsylvania State University.
See LICENSE for licensing.
*/

#include <algorithm>
#include <cmath>
#include "aster.h"
#include "basic_algo.h"

aster::aster(const splice_graph &g, const hyper_set &h)
	: origr(g), gr(g), hs(h)
{
	// prepare
    topological_sort_vertices();
	topological_sort_index_edges();
	make_stats();
	aggressive_purge_intersecting_edges();


	if(asterMode == aster_mode::STAT_ONLY) print_stats();
	assemble();
	get_transcripts();	
}

int aster::assemble()
{	
	if(asterMode == aster_mode::STAT_ONLY) return 0;
	if(gr.num_edges() == 0) return 0;
	if(gr.num_vertices() == 2) return 0;
	if(gr.num_vertices() > 1000) //FIXME:
	{
		cerr << "graph too big to process with D&C " << gr.num_vertices() << endl; 
		return 0;
	}

	assert(gr.num_vertices() > 2);
	//CLEAN: balance?
	if (true)
	{
		for(int i = 1; i < gr.num_vertices() - 1; i++) balance_vertex(i);
		for(int i = 1; i < gr.num_vertices() - 1; i++) balance_vertex(i);
	}
	divide_conquer();

	if (verbose >= 2)  for(int i = 0; i < paths.size(); i++) paths[i].print(i);
	return 0;
}

int aster::divide_conquer()
{
	if(verbose >= 3) gr.graphviz(gr.gid + ".dot", tp2v_to_string()); //CLEAN:
	assert(gr.num_vertices() > 2);
	assert(tp2v.size() == gr.num_vertices());
	int s = 0;									// tp2v index
	int t = gr.num_vertices() - 1;				// tp2v index
	aster_result res;
	divide_conquer(s, t, res);
	paths = res.subpaths;
	assert(paths.size() > 0);
	for(const path & p : paths)	assert(origr.valid_path(p.v));
	return 0;
}

// divide_conquer(i ,j) solves a subproblem between 
// i, j are tp2v indices
int aster::divide_conquer(int source, int target, aster_result& res)
{
	assert(res.subpaths.size() == 0);
	assert(res.dist == -1);
	assert(source <= target);
	assert(source < tp2v.size() && target < tp2v.size());
	int s = tp2v[source];
	int t = tp2v[target];
	assert(s <= t);
	assert(s < gr.num_vertices() && t < gr.num_vertices() && s >= 0 && t >= 0);

	// if(gr.degree(s) == 0 && s < t) divide_conquer(source + 1, target, res);
	// if(gr.degree(t) == 0 && s < t) divide_conquer(target - 1, source, res);

	comb_strat st = comb_strat::GREEDY_MIN; //TODO: which is better

	if (divide_conquer_single_vertex(source, target, res))			
	{
		dnc_counter_single ++;
		return 0;
	} 
	if (divide_conquer_unitig(source, target, res))					
	{
		dnc_counter_unitig ++;
		return 0;
	}
	if (divide_conquer_abutting(source, target, res))				
	{
		dnc_counter_abutting ++;
		return 0;
	}
	if (divide_conquer_nested_subgraphs(source, target, res)) 	
	{
		dnc_counter_nested ++;
		return 0;
	}
	if (divide_conquer_disjoint_at_pivot(source, target, res, st))		
	{
		dnc_counter_disjoint ++;
		return 0;
	}
	
	assert(0);
	return -1;
}

/* remove abutting edge, then call divide_conquer(source, target, res) again */
bool aster::divide_conquer_abutting(int source, int target, aster_result& res)
{	
	assert(source < tp2v.size() && target < tp2v.size());
	int s = tp2v[source];
	int t = tp2v[target];
	assert(s < gr.num_vertices() && t < gr.num_vertices() && s >= 0 && t >= 0);
	assert(s < t - 1);
	assert(source < target - 1);
	assert(gr.out_degree(s) >= 1 || gr.in_degree(t) >= 1);

	if(! gr.edge_exists(s,t)) return false;


	if(verbose >= 2) 
	{
		string msg = "aster D&C with a direct abutting edge between vertex [" + to_string(s) + ", " + to_string(t) + "]"; 
		msg += " (topoIndex [" + to_string(source) + "," + to_string(target) + "])";
		cout << msg << endl;
	}

	// remove abutting edge
	edge_descriptor e = gr.edge(s, t).first;
	assert(e != null_edge);
	double w = gr.get_edge_weight(e);
	char strand = gr.get_edge_info(e).strand;
	int eIdx = e2i.at(e);
	i2e[eIdx] = null_edge;
	e2i.erase(e);
	gr.remove_edge(e);
	assert(gr.check_path(s, t));

	divide_conquer(source, target, res);

	// push back abutting subpath
	int shortestPathSize = res.subpaths.front().v.size();
	for(const path& p: res.subpaths) 
	{
		if(p.v.size() < shortestPathSize) shortestPathSize = p.v.size();
	}
	int eventSize = shortestPathSize - 2;
	assert(eventSize >= 1);
	path p({s, t}, w);
	res.subpaths.push_back(p);
	res.dist += event_size_penalty(eventSize);

	// put back abutting edge
	edge_descriptor e_new = gr.add_edge(s, t);
	edge_info ei;
	ei.weight = w;
	ei.strand = strand;
	gr.set_edge_info(e_new, ei);
	gr.set_edge_weight(e_new, w);
	i2e[eIdx] = e_new;
	e2i.insert({e_new, eIdx});

	return true;
}	

bool aster::divide_conquer_disjoint_subgraphs(int source, int target, aster_result& res)
{
	assert(source < tp2v.size() && target < tp2v.size());
	int s = tp2v[source];
	int t = tp2v[target];
	assert(s < gr.num_vertices() && t < gr.num_vertices() && s >= 0 && t >= 0);
	assert(s < t - 1);
	assert(source < target - 1);
	assert(gr.out_degree(s) >= 1 || gr.in_degree(t) >= 1);
	if(gr.edge_exists(s,t)) return false;
	
	// examine if disjoint subgraphs; as the graph is DFS topo-sorted
	int disjointPoint = -1;
	for(int i = source; i < target; i++)
	{
		int ss = tp2v[i];
		int tt = tp2v[i + 1];
		if(gr.edge_exists(ss, tt)) continue;
		disjointPoint = i;
		break;
	}
	if(disjointPoint == -1) return false;

	if(verbose >= 2) 
	{
		string msg = "aster D&C with nested subgraphs, vertex [" + to_string(s) + "," + to_string(t) + "]"; 
		msg += " (topoIndex [" + to_string(source) + "," + to_string(target) + "])";
		cout << msg << endl;
	}

	for(int i = source + 1; i <= disjointPoint; i++)
	{
		for(int j = disjointPoint + 1; j <= target - 1; j++)
		{
			bool b = gr.check_path(tp2v[i], tp2v[j]);
			if (b) cout << "nested subgraph wrong " << i << " " << tp2v[i] << " " << j << " " << tp2v[j]  << endl;
			assert(!b);
		}
	}

	aster_result res1;
	divide_conquer(source, disjointPoint, res1);
	assert(res1.subpaths.size() > 0);
	for(path& p: res1.subpaths)
	{
		assert(p.v.back() != t);
		p.v.push_back(t);
		assert(origr.valid_path(p.v));
	}
	aster_result res2;
	divide_conquer(disjointPoint + 1, target, res2);
	assert(res2.subpaths.size() > 0);
	for(path& p: res2.subpaths)
	{
		assert(p.v.front() != s);
		p.v.insert(p.v.begin(), s);
		assert(origr.valid_path(p.v));
	}

	//TODO: assert res1 res2 subpaths are completely disjoint except s/t

	int shortestPathSize1 = find_shortest_path(res1);
	int shortestPathSize2 = find_shortest_path(res2);
	assert(shortestPathSize1 - 2 >= 1);
	assert(shortestPathSize2 - 2 >= 1);
	int eventSize = shortestPathSize1 - 2 + shortestPathSize2 - 2;
	assert(eventSize >= 1);
	res.dist = event_size_penalty(eventSize) + res1.dist + res2.dist;
	res.subpaths.clear();
	res.subpaths.insert(res.subpaths.begin(), res1.subpaths.begin(), res1.subpaths.end());
	res.subpaths.insert(res.subpaths.begin(), res2.subpaths.begin(), res2.subpaths.end());
	return true;
}

/* check if the graph can be divided to two disjoint graphs at a pivot k
** divide and conquer the problem to dnc[s, k], dnc [k, t]
** not disjoint: (s,t) edge exists, or multiple subgraphs can be found
*/
bool aster::divide_conquer_disjoint_at_pivot(int source, int target, aster_result& res, comb_strat st)
{
	assert(source < tp2v.size() && target < tp2v.size());
	int s = tp2v[source];
	int t = tp2v[target];
	assert(s < gr.num_vertices() && t < gr.num_vertices() && s >= 0 && t >= 0);
	assert(s < t - 1);
	assert(source < target - 1);
	assert(gr.out_degree(s) >= 1 || gr.in_degree(t) >= 1);
	if(gr.edge_exists(s,t)) return false;
	
	if(verbose >= 2)
	{
		string msg = "aster D&C with disjoint graphs at pivot, vertex [" + to_string(s) + ", " + to_string(t) + "]"; 
		msg += " (topoIndex [" + to_string(source) + "," + to_string(target) + "])";
		cout << msg << endl;
	}

	int pivot = divide_conquer_find_pivot(source, target);
	assert(pivot > source);
	assert(pivot < target);

	aster_result res1;
	divide_conquer(source, pivot, res1);
	assert(res1.subpaths.size() > 0);
	aster_result res2;
	divide_conquer(pivot, target, res2);
	assert(res2.subpaths.size() > 0);

	divide_conquer_combine(res1, res2, pivot, res, st);
	
	return true;
}

/*  combines subpaths of left and right sides of pivot k */
int aster::divide_conquer_combine(aster_result& res1,  aster_result& res2, int pivot, aster_result& comb, comb_strat st)
{
	int k = tp2v[pivot];
	assert(res1.subpaths.size() > 0);
	assert(res2.subpaths.size() > 0);
	for(const path& p: res1.subpaths)	
	{
		assert(p.v.size() >= 1); 
		assert(p.v.back() == k);
	}
	for(const path& p: res2.subpaths)	
	{
		assert(p.v.size() >= 1); 
		assert(p.v.front() == k);
	}

	int index1 = find_shortest_path(res1);
	int index2 = find_shortest_path(res2);
	assert(index1 >= 0);
	assert(index2 >= 0);
	assert(index1 < res1.subpaths.size());
	assert(index2 < res2.subpaths.size());
	
	comb.subpaths.clear();
	comb.dist = -1;	
	const path& rAnchor = res2.subpaths[index2];
	for(int i = 0; i < res1.subpaths.size(); i++)	
	{
		if(i == index1) continue;
		const path& p = res1.subpaths[i];
		double abd = p.abd < rAnchor.abd? p.abd: rAnchor.abd;
		vector<int> v;
		v.insert(v.end(), p.v.begin(), p.v.end());
		v.insert(v.end(), next(rAnchor.v.begin()), rAnchor.v.end());
		assert(origr.valid_path(v));
		comb.subpaths.push_back(path(v, abd));
	}
	const path& lAnchor = res1.subpaths[index1];
	for(int i = 0; i < res2.subpaths.size(); i++)	
	{
		if(i == index2) continue;
		const path& p = res2.subpaths[i];
		double abd = p.abd < lAnchor.abd? p.abd: lAnchor.abd;
		vector<int> v;
		v.insert(v.end(), lAnchor.v.begin(), lAnchor.v.end());
		v.insert(v.end(), next(p.v.begin()), p.v.end());
		assert(origr.valid_path(v));
		comb.subpaths.push_back(path(v, abd));
	}
	double abd = lAnchor.abd < rAnchor.abd? lAnchor.abd: rAnchor.abd;
	vector<int> v;
	v.insert(v.end(), lAnchor.v.begin(), lAnchor.v.end());
	v.insert(v.end(), next(rAnchor.v.begin()), rAnchor.v.end());
	assert(origr.valid_path(v));
	comb.subpaths.push_back(path(v, abd));

	assert(comb.subpaths.size() >= res1.subpaths.size() + res2.subpaths.size() - 1);
	return 0;
}

int aster::divide_conquer_find_pivot(int source, int target)
{
	assert(source < target - 1);
	int s = tp2v[source];
	int t = tp2v[target];
	PEEI peei = gr.out_edges(s);
	int k = -1;
	for(edge_iterator it1 = peei.first, it2 = peei.second; it1 != it2; it1++)
	{
		edge_descriptor e = *it1;
		assert(e != null_edge);
		int tt = e->target();
		assert(tt > s);
		if(tt >= t) continue;
		if(tt > k) k = tt;
	}
	assert(k >= 0);

	int pivot = v2tp.at(k);
	if(verbose >= 2) 
	{
		string msg = "\t inside ["+ to_string(s) + ", " + to_string(t) + "], " + " pivot = " + to_string(tp2v[pivot]);
		msg += " (topoIndex = " + to_string(pivot) + ")";
		cout << msg << endl;
	}
	assert(pivot > source);
	assert(pivot < target);

	//assertion
	splice_graph gr2(gr);
	gr2.clear_vertex(k);
	assert(! gr2.check_path(s, t));

	return pivot;
}

/* examine if dnc unitig between source to target; if true, populate res */
bool aster::divide_conquer_unitig(int source, int target, aster_result& res)
{
	assert(source < tp2v.size() && target < tp2v.size());
	int s = tp2v[source];
	int t = tp2v[target];
	assert(s < gr.num_vertices() && t < gr.num_vertices() && s >= 0 && t >= 0);
	assert(s < t);
	assert(gr.out_degree(s) >= 1 || gr.in_degree(t) >= 1);

	int n = gr.compute_num_paths(s, t, 2);
	assert(n >= 1);
	if(n > 1) return false;	

	vector<int> unitig;
	int    ss    = s;
	bool   _avg_ = false;        // average if true, geom mean if false
	double w     = _avg_? 0: 1;
	int    c     = 0;
	while(true)
	{
		unitig.push_back(ss);
		assert(s <= t);
		if (ss == t) break;
		if (gr.out_degree(ss) > 1) return false;
		edge_descriptor e = (*gr.out_edges(ss).first);
		if(_avg_) w += gr.get_edge_weight(e);
		else w  *= gr.get_edge_weight(e);
		c  += 1;
		ss  = e->target();
	}

	if(verbose >= 2) 
	{
		string msg = "aster D&C with unitig subgraph, vertex [" + to_string(s) + ", " + to_string(t) + "]"; 
		msg += " (topoIndex [" + to_string(source) + "," + to_string(target) + "])";
		cout << msg << endl;
		cout << "\t unitig path is: ";
		printv(unitig);
		cout << endl;
	}

	if(_avg_) w = w / double(c);
	else 	  w = pow(w, 1.0/c);
	res.subpaths.push_back(path(unitig, w));
	res.dist = 0;
	return true;

	/*
	int n = gr.compute_num_paths(s, t, 2);
	assert(n >= 1);
	if(n > 1) return false;	
	vector<edge_descriptor> edgePath;
	double abd = gr.compute_maximum_st_path_w(edgePath, s, t);
	assert(edgePath.size() > 0);
	vector<int> vertexPath;
	edge_path_to_vertex_path(edgePath, vertexPath);
	res.subpaths.push_back(path(vertexPath, abd));
	res.dist = 0;
	*/
	return true;
}

/* WON'T USE. CASE TAKEN BY divide_conquer_unitig
* examine if dnc single vertex; if true, populate res 
*/
bool aster::divide_conquer_single_vertex(int source, int target, aster_result& res)
{
	assert(source < tp2v.size() && target < tp2v.size());
	int s = tp2v[source];
	int t = tp2v[target];
	assert(s < gr.num_vertices() && t < gr.num_vertices() && s >= 0 && t >= 0);
	assert(s <= t);
	if(s != t) return false;
	assert(source == target);
	
	if(verbose >= 2) 
	{
		string msg = "aster D&C with single-vertex subgraph, vertex [" + to_string(s) + ", " + to_string(t) + "]"; 
		msg += " (topoIndex [" + to_string(source) + "," + to_string(target) + "])";
		cout << msg << endl;
	}

	res.subpaths.clear();
	res.dist = -1;
	if(gr.degree(s) == 0)  return true;
	double abd = gr.get_vertex_weight(s);
	res.subpaths.push_back(path({s}, abd)); 
	res.dist = 0;
	return true;
}

/*
* DFS based topological sort 
* This guarantees all disjoint subgraphs are gathered together
* assuming no two vertices have the same <lpos, rpos>
*/
int aster::topological_sort_vertices()
{
	assert(gr.num_vertices() >= 2);
	tp2v.clear();

	vector<bool> visited(gr.num_vertices(), false);
	for(int i = 0; i < gr.num_vertices(); i++)	
	{
		topological_sort_vertices_visit(i, visited);
	}
	reverse(tp2v.begin(), tp2v.end());
	v2tp.resize(tp2v.size());
	for(int i = 0; i < tp2v.size(); i++) v2tp[tp2v[i]] = i;

	// assertions topo sorted
	for(int i = 0; i < gr.num_vertices(); i++)	
	{
		if(i >= 1) assert(gr.get_vertex_info(i - 1).rpos <= gr.get_vertex_info(i).lpos);
		if (i < gr.num_vertices() - 1) 		// proceeding out edges
		{
			PEEI ei = gr.out_edges(i);
			edge_iterator it1 = ei.first; 
			edge_iterator it2 = ei.second;
			for(; it1 != it2; it1++) assert(gr.get_vertex_info(i).rpos <= gr.get_vertex_info((*it1)->target()).lpos);
		}
		if (i >= 1) 						// trailing in edges
		{
			PEEI ei = gr.in_edges(i);
			edge_iterator it1 = ei.first; 
			edge_iterator it2 = ei.second;
			for(; it1 != it2; it1++) assert(gr.get_vertex_info(i).lpos >= gr.get_vertex_info((*it1)->source()).rpos);
		}
	}
	// assertions disjoint sorted
	printv(tp2v);
	cout << endl;
	for(int i = 0; i < gr.num_vertices() - 2; i++)	
	{
		bool previouslyIsDisjoint = false;
		int ss = tp2v[i];
		int tt = tp2v[i + 1];
		bool hasEdge = gr.edge_exists(ss, tt);
		if(hasEdge) continue;
		previouslyIsDisjoint = true;

		for(int j = i + 2; j < gr.num_vertices() - 1; j++)
		{
			int tt = tp2v[j];
			cout << ss << " " << tt << endl;
			bool hasEdge = gr.edge_exists(ss, tt);
			if(previouslyIsDisjoint) assert(! hasEdge);
		}
	}
	assert(tp2v.front() == 0);
	assert(tp2v.back() == gr.num_vertices() - 1);
	assert(tp2v.size() == gr.num_vertices());
	
	return 0;
}

/* visit a node i in DFS, push i to tp2v */
int aster::topological_sort_vertices_visit(int i, vector<bool>& visited)
{
	if(visited[i]) return 0;
	visited[i] = true;

	if(gr.out_degree(i) == 0)
	{
		tp2v.push_back(i);
		return 0;
	}

	PEEI peei;
	edge_iterator it1, it2;
	for(peei = gr.out_edges(i), it1 = peei.first, it2 = peei.second; it1 != it2; it1++)
	{
		int j = (*it1)->target();
		topological_sort_vertices_visit(j , visited);
	}
	tp2v.push_back(i);
	return 0;
}

/* Edges are sorted by their <source, target> represented by tp2v
*  edges should be sorted & fetched independent of splice_graph implementation
*/ 
int aster::topological_sort_index_edges()
{
	map<pair<int, int>, edge_descriptor> sortedEdges;
	PEEI pei = gr.edges(); 
	for(edge_iterator it1 = pei.first, it2 = pei.second; it1 != it2; it1++)
	{
		edge_descriptor e = *it1;
		sortedEdges.insert({{v2tp[e->source()], v2tp[e->target()]}, e});	
	}
	for(auto it = next(sortedEdges.begin()); it != sortedEdges.end(); it++) 
	{
		auto pv = prev(it); 
		int p1 = it->first.first, p2 = it->first.second;
		int q1 = pv->first.first, q2 = pv->first.second;
		assert(p1 >= q1);
		assert(p1 > q1 || (p1 == q1 && p2 > q2));
	}

	i2e.clear();
	e2i.clear();
	int index = 0;
	for(const auto& iie: sortedEdges)
	{
		edge_descriptor e = iie.second;
		e2i.insert({e, index});
		i2e.push_back(e);
		index++;
	}
	assert(e2i.size() == index); 
	assert(i2e.size() == index);
	return 0;
}

/*
*  aggresively remove intersecting edges, whichever is topilocially smaller
*/
bool aster::aggressive_purge_intersecting_edges()
{
	if(gr.num_vertices() > 1000) //FIXME:
	{
		cerr << "graph too big to pruge, #vertex = " << gr.num_vertices() << endl;
		return 0;
	}

	int purgeCount = 0;
	for (int i = 0; i < i2e.size() - 1; i ++)
	{
		edge_descriptor edge1 = i2e[i];
		if (edge1 == null_edge) continue;
		for (int j = i + 1; j < i2e.size(); j ++)
		{
			edge_descriptor edge2 = i2e[j];
			if (edge2 == null_edge) continue;
			if (!gr.intersect(edge1, edge2)) continue;
			i2e[i] = null_edge;
			e2i.erase(edge1);
			gr.remove_edge(edge1);
			purgeCount ++;
			break;
		}
	}
	// if (verbose >= 2 && purgeCount <  1) cout << "graph does not have intersecting edges." << endl;
	if (verbose >= 2 && purgeCount >= 1) cout << "removed " << purgeCount << " intersecting edges from graph.\nSort graph again" << endl;
	gr.refine_splice_graph();
	assert(gr.check_nested());
	if (purgeCount >= 1)
	{
		topological_sort_vertices();
		topological_sort_index_edges();
	}
	return (purgeCount >= 1);
}

int aster::balance_vertex(int vertexIndex)
{
	int v = vertexIndex;
	if(gr.degree(v) <= 0) return 0;

	edge_iterator it1, it2;
	PEEI pei;
	double w1 = 0, w2 = 0;
	for(pei = gr.in_edges(v), it1 = pei.first, it2 = pei.second; it1 != it2; it1++)
	{
		double w = gr.get_edge_weight(*it1);
		w1 += w;
	}
	for(pei = gr.out_edges(v), it1 = pei.first, it2 = pei.second; it1 != it2; it1++)
	{
		double w = gr.get_edge_weight(*it1);
		w2 += w;
	}

	assert(w1 >= SMIN);
	assert(w2 >= SMIN);

	// use sqrt-meature
	double ww = sqrt(w1 * w2);

	double r1 = ww / w1;
	double r2 = ww / w2;

	double m1 = 0, m2 = 0;
	for(pei = gr.in_edges(v), it1 = pei.first, it2 = pei.second; it1 != it2; it1++)
	{
		double wx = gr.get_edge_weight(*it1);
		double wy = wx * r1;
		if(wy < min_guaranteed_edge_weight)
		{
			m1 += (min_guaranteed_edge_weight - wy);
			wy = min_guaranteed_edge_weight;
		}
		gr.set_edge_weight(*it1, wy);
	}
	for(pei = gr.out_edges(v), it1 = pei.first, it2 = pei.second; it1 != it2; it1++)
	{
		double wx = gr.get_edge_weight(*it1);
		double wy = wx * r2;
		if(wy < min_guaranteed_edge_weight)
		{
			m2 += min_guaranteed_edge_weight - wy;
			wy = min_guaranteed_edge_weight;
		}
		gr.set_edge_weight(*it1, wy);
	}

	if(m1 > m2)
	{
		edge_descriptor e = gr.max_out_edge(v);
		double w = gr.get_edge_weight(e);
		gr.set_edge_weight(e, w + m1 - m2);
	}
	else if(m1 < m2)
	{
		edge_descriptor e = gr.max_in_edge(v);
		double w = gr.get_edge_weight(e);
		gr.set_edge_weight(e, w + m2 - m1);
	}

	return 0;
}

int aster::edge_path_to_vertex_path(const VE& edgePath, VI& vertexPath)
{
	if (edgePath.size() == 0) return 0;
	vertexPath.clear();
	vertexPath.resize(edgePath.size() + 1);
	vertexPath[0] = edgePath.front()->source();
	int i = 1;
	for(const edge_descriptor& e: edgePath)
	{
		int s = e->source();
		int t = e->target();
		assert(s == vertexPath[i - 1]);
		vertexPath[i] = t;
		i ++;
	}
	assert(vertexPath.size() >= 2);
	assert(origr.valid_path(vertexPath));
	return 0;
}

// return index of longest path in res.subpaths
int aster::find_longest_path(const aster_result& res)
{
	int longestPathIndex = -1;
	int longestPathSize = -1;
	if (res.subpaths.size() == 0) return longestPathIndex;

	longestPathIndex = 0;
	longestPathSize = res.subpaths.front().v.size();
	for(int i = 0; i < res.subpaths.size(); i++) 
	{
		const path& p = res.subpaths[i];
		if(p.v.size() <= longestPathSize) continue;
		longestPathSize = p.v.size();
		longestPathIndex = i;
	}
	return longestPathIndex;
}

// return index of shortest path in res.subpaths
int aster::find_shortest_path(const aster_result& res)
{
	int shortestPathIndex = -1;
	int shortestPathSize = -1;
	if (res.subpaths.size() == 0) return shortestPathIndex;

	shortestPathIndex = 0;
	shortestPathSize = res.subpaths.front().v.size();
	for(int i = 0; i < res.subpaths.size(); i++) 
	{
		const path& p = res.subpaths[i];
		if(p.v.size() >= shortestPathSize) continue;
		shortestPathSize = p.v.size();
		shortestPathIndex = i;
	}
	return shortestPathIndex;
}


int aster::get_transcripts()
{
	if(asterMode == aster_mode::STAT_ONLY) return 0;
	if(gr.num_edges() == 0) return 0;
	if(gr.num_vertices() == 2) return 0;
	assert(paths.size() > 0);
	trsts.clear();	
	non_full_trsts.clear();
	origr.output_transcripts1(trsts, non_full_trsts, paths);
	if (verbose >= 2) cout << "Aster assembled bundle " << gr.gid.c_str() << endl;
	return 0;
}


int aster::make_stats()
{
	num_graph ++;
	for(int i = 1; i < gr.num_vertices() - 1; i++)
	{
		if(gr.degree(i) == 0) continue;
		num_exon ++;
	}
	for (int i = 0; i < i2e.size(); i ++) 
	{
		if(i2e[i]->source() == 0) continue; 
		if(i2e[i]->target() == gr.num_vertices() - 1) continue;
		num_intron++;
	}
	
	// num intersecting: introns, intron_pairs, graphs
	bool intersecting = false;
	for (int i = 0; i < i2e.size() - 1; i ++)
	{
		bool intersecting_edge = false;
		if(i2e[i] == null_edge) continue;
		for (int j = i + 1; j < i2e.size(); j ++)
		{
			if (i2e[j] == null_edge) continue;
			if (!gr.intersect(i2e[i], i2e[j])) continue;
			num_intersecting_intron_pair ++;
			intersecting_edge = true;
			intersecting = true;
		}
		
		if (!intersecting_edge) continue;
		num_intersecting_intron_count++;
		if(i == i2e.size() - 2) num_intersecting_intron_count ++;
	}

	if(intersecting) num_intersecting_graph ++;
	if(!gr.check_nested()) num_intersecting_graph2 ++;

	if(verbose >= 3 && num_graph % 100 == 0)	print_stats();
	return 0;
}

int aster::print_stats()
{
	cout << "aster print stats" << endl;
	cout << "\t num graph " << num_graph  << endl; 
	cout << "\t num_intersecting_graph " << num_intersecting_graph  << endl; 
	cout << "\t num_intersecting_graph2 " << num_intersecting_graph2  << endl; 
	cout << "\t num intron " << num_intron << endl;
	cout << "\t num_exon " << num_exon  << endl; 
	cout << "\t num_intersecting_intron_count " << num_intersecting_intron_count << endl;
	cout << "\t num_intersecting_intron_pair " << num_intersecting_intron_pair << endl;
	cout << "\t D&C counter: ";
	cout << "\t " << dnc_counter_single;
	cout << "\t " << dnc_counter_unitig;
	cout << "\t " << dnc_counter_abutting;
	cout << "\t " << dnc_counter_nested;
	cout << "\t " << dnc_counter_disjoint;
	cout << endl;
	cout << "aster printed stats" << endl;
	return 0;
}

string aster::tp2v_to_string()
{
	string tp2vString = gr.gid + " DFS TopoSorted vertex index vector:";
	for (int i = 0; i < tp2v.size(); i++)
	{
		if (i % 10 == 0)  tp2vString += "\n\t[" + to_string(i) + "]:";
		tp2vString += to_string(tp2v[i]);
		tp2vString += " ";
	}
	return tp2vString;
}

// exponential penalty guarantees to violate triangle inequality
int aster::event_size_penalty(int eventSize)
{
	assert(eventSize >= 0);
	return pow(2, eventSize) - 1;
}

/*
* return: 
*		positive int: number of edits
*		assertion error: size not positive
*/
int aster::path_distance(const path& p1, const path& p2)
{
	throw runtime_error("aster::path_distancen not implemented yet");
	const vector<int>& v1 = p1.v; 
	const vector<int>& v2 = p2.v;
	assert(v1.size() > 0);
	assert(v2.size() > 0);
	if(v1.back() > v2.front()) return -1;
	if(v2.back() > v1.front()) return -1;

	int edits = - basic_algo::ref_sw_query_nw(v1, v2, -1, -2, 0);

	return edits; //TODO: not finished
}

