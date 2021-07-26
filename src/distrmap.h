#include <unordered_set>
#include <upcxx/upcxx.hpp>

#ifdef USE_TSL_ROBIN_MAP
#include <tsl/robin_map.h>
#else
#include "robin_hood.h"
#endif

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

template<typename K, typename V>
bool TRUE_FILTER(const K&, const V&) {
	return true;
}

template<typename K, typename V> class DistrMap {
private:
	// store the local unordered map in a distributed object to access from RPCs
	using dobj_map_t = upcxx::dist_object<unordered_map<K,V> >;
	dobj_map_t local_map;
	// map the key to a target process
	uint32_t get_target_rank(const K &key) {
		return consistent_hash<K>(key) % upcxx::rank_n();
	}
public:
	// initialize the local map
	DistrMap() :
			local_map(unordered_map<K, V>()) {
	}

	upcxx::future<> incr(const K &key, const V &val) {
		return upcxx::rpc(get_target_rank(key),
				[](dobj_map_t &lmap, K key, V val) {
					lmap->operator[](key) += val;
				}, local_map, key, val);
	}

	upcxx::future<> incr(const std::vector<K> &arr) {
		std::unordered_map<uint32_t, std::vector<K>> grouped;
		for (auto &a : arr) {
			grouped[get_target_rank(a)].push_back(a);
		}
		upcxx::future<> fut_all = upcxx::make_future();

		uint32_t i = 0;
		for (auto &kv : grouped) {
			uint32_t target_rank = kv.first;
			std::vector<K> &v = kv.second;
			upcxx::future<> fut = upcxx::rpc(target_rank,
					[](dobj_map_t &lmap, upcxx::view<K> val) {
						for (auto a : val) {
							lmap->operator[](a)++;}
						}
					, local_map, upcxx::make_view(v.begin(), v.end()));
			fut_all = upcxx::when_all(fut_all, fut);
		}
		return fut_all;
	}

	template<typename PairV> upcxx::future<> value_set_insert(
			const std::vector<std::pair<K, PairV>> &values) {
		unordered_map<uint32_t, std::vector<K> > grouped_K;
		unordered_map<uint32_t, std::vector<PairV> > grouped_V;
		for (const std::pair<K, PairV> &a : values) {
			auto key = get_target_rank(a.first);
			grouped_K[key].push_back(a.first);
			grouped_V[key].push_back(a.second);
		}
		upcxx::future<> fut_all = upcxx::make_future();
		uint32_t i = 0;
		for (auto &kv : grouped_K) {
			uint32_t target_rank = kv.first;
			std::vector<K> &k = kv.second;
			std::vector<PairV> &v = grouped_V.at(target_rank);
			upcxx::future<> fut = upcxx::rpc(target_rank,
					[](dobj_map_t &lmap, upcxx::view<K> keys,
							upcxx::view<PairV> vals) {
						auto it1 = keys.begin();
						auto it2 = vals.begin();
						for (; it1 != keys.end(); it1++, it2++) {
							lmap->operator[](*it1).insert(*it2);
						}
					}, local_map, upcxx::make_view(k.begin(), k.end()),
					upcxx::make_view(v.begin(), v.end()));
			fut_all = upcxx::when_all(fut_all, fut);
			if (i++ % 10 == 0) {
				upcxx::progress();
			}

		}
		return fut_all;

	}

	template<typename SetV> upcxx::future<> value_set_insert(const K &key,
			const SetV &val) {
		return upcxx::rpc(get_target_rank(key),
				[](dobj_map_t &lmap, K key, SetV val) {
					lmap->operator[](key).insert(val);
				}, local_map, key, val);
	}

	// insert a key-value pair into the hash table
	upcxx::future<> insert(const K &key, const V &val) {
		// the RPC returns an empty upcxx::future by default
		return upcxx::rpc(get_target_rank(key),
		// lambda to insert the key-value pair
				[](dobj_map_t &lmap, K key, V val) {
					// insert into the local map at the target
					lmap->insert( { key, val });
				}, local_map, key, val);
	}
	// find a key and return associated value in a future
	upcxx::future<V> find(const std::string &key) {
		return upcxx::rpc(get_target_rank(key),
		// lambda to find the key in the local map
				[](dobj_map_t &lmap, K key) -> std::string {
					auto elem = lmap->find(key);
					// no key found
					if (elem == lmap->end())
						return empty<V>();
					// the key was found, return the value
					return elem->second;
				}, local_map, key);
	}

	template<typename VALFUN> int dump(const std::string &filepath, char sep,
			VALFUN fun) {
		return dump(filepath, sep, fun, TRUE_FILTER<K, V>);
	}

	template<typename VALFUN, typename FILTER> int dump(
			const std::string &filepath, char sep, VALFUN fun, FILTER filter) {
		int stat = 0;
		std::ostream *myfile_pointer = 0;
		if (sparc::endswith(filepath, ".gz")) {
			myfile_pointer = new ogzstream(filepath.c_str());

		} else {
			myfile_pointer = new std::ofstream(filepath.c_str());
		}
		std::ostream &myfile = *myfile_pointer;

		MapIterator<K, V> it = local_map->begin();
		uint64_t n = 0;
		for (; it != local_map->end(); it++) {
			if (filter(it->first, it->second)) {
				myfile << it->first << sep << fun(it->second) << std::endl;
				++n;
			}
		}
		if (sparc::endswith(filepath, ".gz")) {
			((ogzstream&) myfile).close();
		} else {
			((std::ofstream&) myfile).close();
		}
		delete myfile_pointer;
		myinfo("Wrote %ld records", n);
		return stat;

	}

};
