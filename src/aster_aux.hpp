/*
Part of Aster Assembler
(c) 2024 by Xiaofei Carl Zang, Mingfu Shao, and The Pennsylvania State University.
See LICENSE for licensing.
*/

#ifndef __ASTER_AUX_H__
#define __ASTER_AUX_H__

#include <stdexcept>
#include "splice_graph.h"
#include "hyper_set.h"
#include "equation.h"
#include "router.h"
#include "path.h"
#include "config.h"

typedef map< edge_descriptor, vector<int> > MEV;
typedef pair< edge_descriptor, vector<int> > PEV;
typedef pair< vector<int>, vector<int> > PVV;
typedef pair<PEE, int> PPEEI;
typedef map<PEE, int> MPEEI;
typedef pair<int, int> PI;
typedef map<int, int> MI;
typedef vector<int> VI;

struct aster_index
{
private:
	vector<int> indices;	
	
public:
	inline aster_index() {};
	inline aster_index(vector<int> v): indices(v) {};
	inline const vector<int>& get_index() const {return indices;}
	inline int s() const {return indices.front();}
	inline int t() const {return indices.back();}
	inline int at(int i) const {return indices.at(i);}
	inline int size() const {return indices.size();}
	// split an index at index i, to left and right, both inclusive of i
	inline void split(int i, aster_index left, aster_index right) const
	{
		assert(i > 0 && i < size() - 1);
		assert(size() >= 3);
		left = aster_index(vector<int>(indices.begin(), indices.begin() + i + 1));
		right = aster_index(vector<int>(indices.begin() + i, indices.end()));
		assert(left.size() >= 2);
		assert(right.size() >= 2);
		assert(left.t() == right.s());
		assert(left.size() + right.size() == size() - 1);
	}
	inline void erase_itertor(vector<int>::iterator it) 
	{
		assert(std::distance(indices.begin(), it) >= 0);
		assert(std::distance(indices.end(), it)   <= 0);
		indices.erase(it);
	}
	inline void erase_index(int i) 
	{
		assert(i >= 0 && i <= size() - 1);
		indices.erase(indices.begin() + i);
	}
	inline void erase_element(int i) 
	{
		assert(i >= s() && i <= t());
		auto it = find(indices.begin(), indices.end(), i);
		assert(it != indices.end());
		indices.erase(it);
	}
	inline bool find_index(int i) const {return find(indices.begin(), indices.end(), i) != indices.end();}
	inline bool operator< (const aster_index& ai) const {return indices < ai.get_index();}
};

struct aster_result
{
	vector<path> subpaths;	// predicted paths, original v index, inclusive
	int dist = -1;
	inline aster_result() {};
	inline aster_result(const vector<int>& v, double w): subpaths({path(v, w)}), dist(0) {};
	inline void clear() {subpaths.clear(); dist = -1;}
};

class aster_error : public runtime_error 
{
public:
    aster_error(const char* message) : runtime_error(message) {}
};

#endif