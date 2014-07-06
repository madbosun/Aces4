/*
 * lru_array_policy.h
 *
 *  Created on: May 15, 2014
 *      Author: njindal
 */

#ifndef LRU_ARRAY_POLICY_H_
#define LRU_ARRAY_POLICY_H_

#include "sip_interface.h"
#include "id_block_map.h"
#include "block_id.h"

#include <list>

namespace sip {

class BlockId;
template <class BLOCK_TYPE>

/**
 * Policy for block replacement on servers.
 * This policy maintains an LRU data structure on least recently used arrays.
 * From the array, an arbitrary block is offered up for removal.
 */
class LRUArrayPolicy {
public:
	LRUArrayPolicy(IdBlockMap<BLOCK_TYPE>& block_map) : block_map_(block_map) {};
    
    /*! Mark a block as recently used  */
	void touch(const BlockId& bid) {
		 // Search for array_id; if present, place at front of queue

		int array_id = bid.array_id();
		std::list<int>::iterator it;
		for (it = lru_list_.begin(); it != lru_list_.end(); ++it) {
			if (array_id == *it)
				break;
		}

		if (it != lru_list_.end())
			lru_list_.erase(it);

		lru_list_.push_front(array_id);

	}

    /*! Removes all blocks with given array id */
    void remove_all_blocks_for_array(int array_id){
    	lru_list_.remove(array_id);
    }

    /*! Get the next block for removal */
	BlockId get_next_block_for_removal(){
		/** Return an arbitrary block from the least recently used array
		 */
		sip::check(!lru_list_.empty(),
				"No blocks have been touched, yet block requested for flushing");
		while (!lru_list_.empty()) {
			int to_remove_array = lru_list_.back();
			typename IdBlockMap<BLOCK_TYPE>::PerArrayMap* array_map = block_map_.per_array_map(to_remove_array);
			typename IdBlockMap<BLOCK_TYPE>::PerArrayMap::iterator it = array_map->begin();
			if (it == array_map->end())
				lru_list_.pop_back();
			else
				return it->first;
		}
		sip::fail("No blocks to remove !", current_line());
	}

private:
	std::list<int> lru_list_; /*! List of least recently ysed arrays  */
    IdBlockMap<BLOCK_TYPE>& block_map_;

};

} /* namespace sip */

#endif /* LRU_ARRAY_POLICY_H_ */
