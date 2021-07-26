#include <vector>
#include <queue>
#include "upcxx_comp.h"

#ifdef USE_TSL_ROBIN_MAP
#include <tsl/robin_map.h>
#else
#include "robin_hood.h"
#endif

#include "log.h"
#include "utils.h"

#ifdef USE_TSL_ROBIN_MAP
template<typename K, typename V>
using unordered_map = tsl::robin_map<K,V>;
#else
template<typename K, typename V>
using unordered_map = robin_hood::unordered_map<K,V>;
#endif

template<typename K, typename V>
using MapIterator = typename unordered_map<K,V>::const_iterator;

template<typename T> uint32_t consistent_hash(const T &v);
template<typename T> T empty();

template<> inline uint32_t consistent_hash(const std::string &str) {
	uint32_t h = 2166136261;
	for (size_t i = 0; i < str.size(); i++) {
		h = h ^ uint32_t(int8_t(str[i]));
		h = h * 16777619;
	}
	return h;
}
template<> inline uint32_t consistent_hash(const uint32_t &i) {
	return consistent_hash(std::to_string(i));
}


template<typename T>
struct NodeWeight {

	T id;
	float w;
};

template<typename T>
class NodeWeightComparator {
public:
	// Comparator function
	bool operator()(NodeWeight<T> &left, NodeWeight<T> &right) {

		//remove low weight first ,Accorading to std pop
		//Removes the element on top of the priority_queue, effectively reducing its size by one.
		//The element removed is the one with the highest value.
		return left.w > right.w;
	}
};

template<typename READID_TYPE>
using EdgeWeightQueue = std::priority_queue<NodeWeight<READID_TYPE>,
std::vector<NodeWeight<READID_TYPE>>, NodeWeightComparator<READID_TYPE>>;

class DistrGraph {

private:
	// store the local unordered map in a distributed object to access from RPCs
	using dobj_map_t = upcxx::dist_object<unordered_map<uint32_t,EdgeWeightQueue<uint32_t>> >;
	dobj_map_t local_map;
	// map the key to a target process
	uint32_t get_target_rank(const uint32_t &key) {
		return consistent_hash<uint32_t>(key) % upcxx::rank_n();
	}

	uint32_t k_neigbhor;
public:
	// initialize the local map
	DistrGraph(uint32_t k_neigbhor) :
			local_map(unordered_map<uint32_t, EdgeWeightQueue<uint32_t>>()), k_neigbhor(
					k_neigbhor) {
	}

	upcxx::future<> add_edges(
			std::vector<std::tuple<uint32_t, uint32_t, float>> &edges) {
		std::unordered_map<uint32_t, std::vector<uint32_t>> grouped_A;
		std::unordered_map<uint32_t, std::vector<uint32_t>> grouped_B;
		std::unordered_map<uint32_t, std::vector<float>> grouped_W;
		for (auto &x : edges) {
			auto rank = get_target_rank(std::get<0>(x));
			grouped_A[rank].push_back(std::get<0>(x));
			grouped_B[rank].push_back(std::get<1>(x));
			grouped_W[rank].push_back(std::get<2>(x));
		}
		for (auto &x : edges) {
			auto rank = get_target_rank(std::get<1>(x));
			grouped_A[rank].push_back(std::get<1>(x));
			grouped_B[rank].push_back(std::get<0>(x));
			grouped_W[rank].push_back(std::get<2>(x));
		}

		upcxx::future<> fut_all = upcxx::make_future();

		uint32_t i = 0;
		for (auto &kv : grouped_A) {
			uint32_t target_rank = kv.first;
			auto &va = kv.second;
			auto &vb = grouped_B.at(target_rank);
			auto &vw = grouped_W.at(target_rank);
			auto k = k_neigbhor;
			upcxx::future<> fut = upcxx::rpc(target_rank,
					[k](dobj_map_t &lmap, upcxx::view<uint32_t> va,
							upcxx::view<uint32_t> vb, upcxx::view<float> vw) {
						auto it1 = va.begin();
						auto it2 = vb.begin();
						auto it3 = vw.begin();
						for (; it1 != va.end(); it1++, it2++, it3++) {
							if (true) {
								auto &node = lmap->operator[](*it1);
								node.push( { *it2, *it3 });
								while (node.size() > k) {
									node.pop();
								}
							}
							if (true) {
								auto &node = lmap->operator[](*it2);
								node.push( { *it1, *it3 });
								while (node.size() > k) {
									node.pop();
								}
							}
						}
					}

					, local_map, upcxx::make_view(va.begin(), va.end()),
					upcxx::make_view(vb.begin(), vb.end()),
					upcxx::make_view(vw.begin(), vw.end()));
			fut_all = upcxx::when_all(fut_all, fut);
			if (i++ % 10 == 0) {
				upcxx::progress();
			}
		}
		return fut_all;
	}

	int dump_edges(const std::string &filepath, char sep) {
		int stat = 0;
		std::ostream *myfile_pointer = 0;
		if (sparc::endswith(filepath, ".gz")) {
			myfile_pointer = new ogzstream(filepath.c_str());

		} else {
			myfile_pointer = new std::ofstream(filepath.c_str());
		}
		std::ostream &myfile = *myfile_pointer;
		uint64_t n = 0;
		for (unordered_map<uint32_t, EdgeWeightQueue<uint32_t>>::iterator itor =
				local_map->begin(); itor != local_map->end(); itor++, n++) {

			auto &tmp_q = itor->second;
			while (!tmp_q.empty()) {
				auto q_element = tmp_q.top();
				myfile << itor->first << sep << q_element.id << " "
						<< q_element.w << "\n";
				tmp_q.pop();
			}

		}

		if (sparc::endswith(filepath, ".gz")) {
			((ogzstream&) myfile).close();
		} else {
			((std::ofstream&) myfile).close();
		}
		delete myfile_pointer;
		myinfo("Wrote %ld nodes", n);

		upcxx::barrier();
		uint64_t sum_n = upcxx_reduce_all(n, upcxx_op_add).wait();
		if (upcxx::rank_me() == 0) {
			myinfo("# of total nodes = %ld", sum_n);
		}
		return stat;

	}

};
