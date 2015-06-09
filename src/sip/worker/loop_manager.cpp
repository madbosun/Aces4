/*
 * loop_manager.cpp
 *
 *  Created on: Aug 23, 2013
 *      Author: basbas
 */

#include "loop_manager.h"
#include <sstream>
#include "interpreter.h"

namespace sip {
LoopManager::LoopManager() :
		to_exit_(false) {
}
LoopManager::~LoopManager() {
}
void LoopManager::do_set_to_exit() {
	to_exit_ = true;
}
std::ostream& operator<<(std::ostream& os, const LoopManager &obj) {
	os << obj.to_string();
	return os;
}

std::string LoopManager::to_string() const {
	return std::string("LoopManager");
}


//+++++++++++++++++++++++++++++++++++++++++

DoLoop::DoLoop(int index_id, DataManager & data_manager,
		const SipTables & sip_tables) :
		data_manager_(data_manager), sip_tables_(sip_tables), index_id_(
				index_id), first_time_(true) {
	lower_seg_ = sip_tables_.lower_seg(index_id);
	upper_bound_ = lower_seg_ + sip_tables_.num_segments(index_id);
//	sip::check_and_warn(lower_seg_ < upper_bound_,
//			std::string("doloop has empty range"),
//			Interpreter::global_interpreter->line_number());
}

DoLoop::~DoLoop() {
}

bool DoLoop::do_update() {
//	std::cout << "DEBUG DoLoop:update:46 \n" << this->to_string() << std::endl;
	if (to_exit_)
		return false;
	int current_value;
	if (first_time_) {  //initialize index to lower value
		first_time_ = false;
		sip::check(
				data_manager_.index_value(index_id_)
						== DataManager::undefined_index_value,
				"SIAL or SIP error, index " + sip_tables_.index_name(index_id_)
						+ " already has value before loop",
				Interpreter::global_interpreter->line_number());
		current_value = lower_seg_;
	} else { //not the first time through loop.  Get the current value and try to increment it
		current_value = data_manager_.index_value(index_id_);
		++current_value;
	}
	if (current_value < upper_bound_) {
		data_manager_.set_index_value(index_id_, current_value);
		return true;
	}
	//If here on first time through, the range is empty. This leaves index undefined, which is required behavior.
	return false;
}
void DoLoop::do_finalize() {
	data_manager_.set_index_undefined(index_id_);
}

std::string DoLoop::to_string() const {
	std::stringstream ss;
	ss << "index_id_=" << sip_tables_.index_name(index_id_);
	ss << ", lower_seg_=" << lower_seg_;
	ss << ", upper_bound_=" << upper_bound_;
	ss << ", current= " << data_manager_.index_value_to_string(index_id_);
	return ss.str();
}

std::ostream& operator<<(std::ostream& os, const DoLoop &obj) {
	os << obj.to_string();
	return os;
}

//+++++++++++++++++++++++++++++++++++

SubindexDoLoop::SubindexDoLoop(int subindex_id, DataManager & data_manager,
		const SipTables & sip_tables) :
		DoLoop(subindex_id, data_manager, sip_tables) {
	sip::check(sip_tables_.is_subindex(subindex_id),
			"Attempting subindex do loop with non-subindex loop variable");
	parent_id_ = sip_tables_.parent_index(subindex_id);
	parent_value_ = data_manager_.index_value(parent_id_);
	lower_seg_ = 1;  //subindices always start at 1
	upper_bound_ = lower_seg_
			+ sip_tables_.num_subsegments(subindex_id, parent_value_);
	sip::check_and_warn(lower_seg_ < upper_bound_,
			std::string("SubindexDoLoop has empty range"),
			Interpreter::global_interpreter->line_number());
}

std::string SubindexDoLoop::to_string() const {
	std::stringstream ss;
	ss << DoLoop::to_string();
	ss << ", parent_id_=" << parent_id_;
	ss << ", parent_value_=" << parent_value_;
	return ss.str();
}

std::ostream& operator<<(std::ostream& os, const SubindexDoLoop &obj) {
	os << obj.to_string();
	return os;
}

SubindexDoLoop::~SubindexDoLoop() {
}


//++++++++++++++++++++++++++++++++++++++++

//note that the max number of indices allowed by the implementation is MAX_RANK.  This limitation is
// due to the structure of the pardo instruction inherited from aces3
SequentialPardoLoop::SequentialPardoLoop(int num_indices,
		const int (&index_id)[MAX_RANK], DataManager & data_manager,
		const SipTables & sip_tables) :
		data_manager_(data_manager), sip_tables_(sip_tables), num_indices_(
				num_indices), first_time_(true) {
	std::copy(index_id + 0, index_id + MAX_RANK, index_id_ + 0);
	for (int i = 0; i < num_indices; ++i) {
		lower_seg_[i] = sip_tables_.lower_seg(index_id_[i]);
		upper_bound_[i] = lower_seg_[i]
				+ sip_tables_.num_segments(index_id_[i]);
		sip::check(lower_seg_[i] < upper_bound_[i],
				"Pardo loop index " + sip_tables_.index_name(index_id_[i])
						+ " has empty range",
				Interpreter::global_interpreter->line_number());
	}
//	std::cout << "SequentialPardoLoop::SequentialPardoLoop at line " << Interpreter::global_interpreter->line_number()  << std::endl;
//
//		Interpreter::global_interpreter->set_index_value(index_ids_[i],
//				first_segments_[i]);
}

SequentialPardoLoop::~SequentialPardoLoop() {
}

bool SequentialPardoLoop::do_update() {
//	std::cout << "DEBUG: SequentialPardoLoop:112 \n" << this->to_string()
//			<< std::endl;
	if (to_exit_)
		return false;
	int current_value;
	if (first_time_) {
		first_time_ = false;
		//initialize values of all indices
		for (int i = 0; i < num_indices_; ++i) {
			if (lower_seg_[i] >= upper_bound_[i])
				return false; //this loop has an empty range in at least one dimension.
			sip::check(
					data_manager_.index_value(index_id_[i])
							== DataManager::undefined_index_value,
					"SIAL or SIP error, index "
							+ sip_tables_.index_name(index_id_[i])
							+ " already has value before loop",
					Interpreter::global_interpreter->line_number());
			data_manager_.set_index_value(index_id_[i], lower_seg_[i]);
		}
		return true;
	} else {
		for (int i = 0; i < num_indices_; ++i) {
			current_value = data_manager_.index_value(index_id_[i]);
			++current_value;
			if (current_value < upper_bound_[i]) { //increment current index and return
				data_manager_.set_index_value(index_id_[i], current_value);
				return true;
			} else { //wrap around and handle next index
				data_manager_.set_index_value(index_id_[i], lower_seg_[i]);
			}
		} //if here, then all indices are at their max value
	}
	return false;
}

void SequentialPardoLoop::do_finalize() {
	for (int i = 0; i < num_indices_; ++i) {
		data_manager_.set_index_undefined(index_id_[i]);
	}

}

std::string SequentialPardoLoop::to_string() const {
	std::stringstream ss;
	ss << "Sequential Pardo Loop:  num_indices=" << num_indices_ << std::endl;
	ss << "index_ids_=[";
	for (int i = 0; i < num_indices_; ++i) {
		ss << (i == 0 ? "" : ",") << sip_tables_.index_name(index_id_[i]);
	}
	ss << "] lower_seg_=[";
	for (int i = 0; i < num_indices_; ++i) {
		ss << (i == 0 ? "" : ",") << lower_seg_[i];
	}
	ss << "] upper_bound_=[";
	for (int i = 0; i < num_indices_; ++i) {
		ss << (i == 0 ? "" : ",") << upper_bound_[i];
	}
	ss << "] current= [";
	for (int i = 0; i < num_indices_; ++i) {
		ss << (i == 0 ? "" : ",")
				<< data_manager_.index_value_to_string(index_id_[i]);
	}
	ss << "]";
	return ss.str();
}
;

#ifdef HAVE_MPI

//++++++++++++++++++++++++++++++++++++

StaticTaskAllocParallelPardoLoop::StaticTaskAllocParallelPardoLoop(
		int num_indices, const int (&index_id)[MAX_RANK],
		DataManager & data_manager, const SipTables & sip_tables,
		SIPMPIAttr & sip_mpi_attr) :
		data_manager_(data_manager), sip_tables_(sip_tables), num_indices_(
				num_indices), first_time_(true), iteration_(0), sip_mpi_attr_(
				sip_mpi_attr) {

	std::copy(index_id + 0, index_id + MAX_RANK, index_id_ + 0);
	for (int i = 0; i < num_indices; ++i) {
		lower_seg_[i] = sip_tables_.lower_seg(index_id_[i]);
		upper_bound_[i] = lower_seg_[i]
				+ sip_tables_.num_segments(index_id_[i]);
		sip::check(lower_seg_[i] < upper_bound_[i],
				"Pardo loop index " + sip_tables_.index_name(index_id_[i])
						+ " has empty range",
				Interpreter::global_interpreter->line_number());
	}
}

StaticTaskAllocParallelPardoLoop::~StaticTaskAllocParallelPardoLoop() {
}

inline bool StaticTaskAllocParallelPardoLoop::increment_indices() {
	bool more = false; 	// More iterations?
	int current_value;
	for (int i = 0; i < num_indices_; ++i) {
		current_value = data_manager_.index_value(index_id_[i]);
		++current_value;
		if (current_value < upper_bound_[i]) {
			//increment current index and return
			data_manager_.set_index_value(index_id_[i], current_value);
			more = true;
			break;
		} else {
			//wrap around and handle next index
			data_manager_.set_index_value(index_id_[i], lower_seg_[i]);
		}
	} //if here, then all indices are at their max value
	return more;
}

inline bool StaticTaskAllocParallelPardoLoop::initialize_indices() {
	//initialize values of all indices
	bool more_iterations = true;
	for (int i = 0; i < num_indices_; ++i) {
		if (lower_seg_[i] >= upper_bound_[i]) {
			more_iterations = false; //this loop has an empty range in at least one dimension.
			return more_iterations;
		}

		sip::check(
				data_manager_.index_value(index_id_[i])
						== DataManager::undefined_index_value,
				"SIAL or SIP error, index "
						+ sip_tables_.index_name(index_id_[i])
						+ " already has value before loop",
				Interpreter::global_interpreter->line_number());
		data_manager_.set_index_value(index_id_[i], lower_seg_[i]);
	}
	more_iterations = true;
	return more_iterations;
}

bool StaticTaskAllocParallelPardoLoop::do_update() {

	if (to_exit_)
		return false;

	int company_rank = sip_mpi_attr_.company_rank();
	int num_workers = sip_mpi_attr_.num_workers();

	if (first_time_) {
		first_time_ = false;
		bool more_iters = initialize_indices();
		while (more_iters && iteration_ % num_workers != company_rank) {
			more_iters = increment_indices();
			iteration_++;
		}
		return more_iters;
	} else {
		iteration_++;
		bool more_iters = increment_indices();
		while (more_iters && iteration_ % num_workers != company_rank) {
			more_iters = increment_indices();
			iteration_++;
		}
		return more_iters;
	}
}

void StaticTaskAllocParallelPardoLoop::do_finalize() {
	for (int i = 0; i < num_indices_; ++i) {
		data_manager_.set_index_undefined(index_id_[i]);
	}
}

std::string StaticTaskAllocParallelPardoLoop::to_string() const {
	std::stringstream ss;
	ss << "Static Task Allocation Parallel Pardo Loop:  num_indices="
			<< num_indices_ << std::endl;
	ss << "index_ids_=[";
	for (int i = 0; i < num_indices_; ++i) {
		ss << (i == 0 ? "" : ",") << sip_tables_.index_name(index_id_[i]);
	}
	ss << "] lower_seg_=[";
	for (int i = 0; i < num_indices_; ++i) {
		ss << (i == 0 ? "" : ",") << lower_seg_[i];
	}
	ss << "] upper_bound_=[";
	for (int i = 0; i < num_indices_; ++i) {
		ss << (i == 0 ? "" : ",") << upper_bound_[i];
	}
	ss << "] current= [";
	for (int i = 0; i < num_indices_; ++i) {
		ss << (i == 0 ? "" : ",")
				<< data_manager_.index_value_to_string(index_id_[i]);
	}
	ss << "]";
	return ss.str();
}

std::ostream& operator<<(std::ostream& os,
		const StaticTaskAllocParallelPardoLoop &obj) {
	os << obj.to_string();
	return os;
}

//++++++++++++++++++++++++++++++++++++++++++++

BalancedTaskAllocParallelPardoLoop::BalancedTaskAllocParallelPardoLoop(
		int num_indices, const int (&index_id)[MAX_RANK],
		DataManager & data_manager, const SipTables & sip_tables,
		SIPMPIAttr & sip_mpi_attr, int num_where_clauses,
		Interpreter* interpreter, long& iteration) :
		data_manager_(data_manager), sip_tables_(sip_tables), num_indices_(
				num_indices), first_time_(true), iteration_(iteration), sip_mpi_attr_(
				sip_mpi_attr), num_where_clauses_(num_where_clauses), company_rank_(
				sip_mpi_attr.company_rank()), num_workers_(
				sip_mpi_attr_.num_workers()), interpreter_(interpreter) {

	std::copy(index_id + 0, index_id + MAX_RANK, index_id_ + 0);
	for (int i = 0; i < num_indices; ++i) {
		lower_seg_[i] = sip_tables_.lower_seg(index_id_[i]);
		upper_bound_[i] = lower_seg_[i]
				+ sip_tables_.num_segments(index_id_[i]);
		sip::check(lower_seg_[i] < upper_bound_[i],
				"Pardo loop index " + sip_tables_.index_name(index_id_[i])
						+ " has empty range",
				Interpreter::global_interpreter->line_number());
	}
}

BalancedTaskAllocParallelPardoLoop::~BalancedTaskAllocParallelPardoLoop() {
}

inline bool BalancedTaskAllocParallelPardoLoop::increment_indices() {
	bool more = false; 	// More iterations?
	int current_value;
	for (int i = 0; i < num_indices_; ++i) {
		current_value = data_manager_.index_value(index_id_[i]);
		++current_value;
		if (current_value < upper_bound_[i]) {
			//increment current index and return
			data_manager_.set_index_value(index_id_[i], current_value);
			more = true;
			break;
		} else {
			//wrap around and handle next index
			data_manager_.set_index_value(index_id_[i], lower_seg_[i]);
		}
	} //if here, then all indices are at their max value
	return more;
}

inline bool BalancedTaskAllocParallelPardoLoop::initialize_indices() {
	//initialize values of all indices
	bool more_iterations = true;
	for (int i = 0; i < num_indices_; ++i) {
		if (lower_seg_[i] >= upper_bound_[i]) {
			more_iterations = false; //this loop has an empty range in at least one dimension.
			return more_iterations;
		}
		sip::check(
				data_manager_.index_value(index_id_[i])
						== DataManager::undefined_index_value,
				"SIAL or SIP error, index "
						+ sip_tables_.index_name(index_id_[i])
						+ " already has value before loop",
				Interpreter::global_interpreter->line_number());
		data_manager_.set_index_value(index_id_[i], lower_seg_[i]);
	}
	more_iterations = true;
	return more_iterations;
}

bool BalancedTaskAllocParallelPardoLoop::do_update() {

	if (to_exit_)
		return false;
	bool more_iters;
	if (first_time_) {
		first_time_ = false;
		more_iters = initialize_indices();
	} else {
		more_iters = increment_indices();
	}

	while(more_iters){
		bool where_clauses_value = interpreter_->interpret_where(num_where_clauses_);
		//if true, the pc will be after the last where clause
		//otherwise it is undefined
		if(where_clauses_value){
			iteration_++;
			if ((iteration_-1) % num_workers_ == company_rank_){
//				std::cout << "rank " << company_rank_ << " executing iteration";
//				std::cout <<  " [";
//				for (int i = 0; i < num_indices_; ++i) {
//					std::cout << data_manager_.index_value(index_id_[i]) << ",";
//				}
//				std::cout << "]" << std::endl << std::flush;
				return true;
			}
		}
		more_iters = increment_indices();
	}
	return more_iters; //this should be false here
}

void BalancedTaskAllocParallelPardoLoop::do_finalize() {
	for (int i = 0; i < num_indices_; ++i) {
		data_manager_.set_index_undefined(index_id_[i]);
	}
}

std::string BalancedTaskAllocParallelPardoLoop::to_string() const {
	std::stringstream ss;
	ss << "Balanced Task Allocation Parallel Pardo Loop:  num_indices="
			<< num_indices_ << std::endl;
	ss << "index_ids_=[";
	for (int i = 0; i < num_indices_; ++i) {
		ss << (i == 0 ? "" : ",") << sip_tables_.index_name(index_id_[i]);
	}
	ss << "] lower_seg_=[";
	for (int i = 0; i < num_indices_; ++i) {
		ss << (i == 0 ? "" : ",") << lower_seg_[i];
	}
	ss << "] upper_bound_=[";
	for (int i = 0; i < num_indices_; ++i) {
		ss << (i == 0 ? "" : ",") << upper_bound_[i];
	}
	ss << "] current= [";
	for (int i = 0; i < num_indices_; ++i) {
		ss << (i == 0 ? "" : ",")
				<< data_manager_.index_value_to_string(index_id_[i]);
	}
	ss << "]";
	return ss.str();
}

std::ostream& operator<<(std::ostream& os,
		const BalancedTaskAllocParallelPardoLoop &obj) {
	os << obj.to_string();
	return os;
}

/*

 naming scheme: _<simple index fragment indices>_<occ (o), virt (v), ao (a) in ifrag>_<o,v,a in jfrag>

*/

/*
-------------------------------------------
            _ij_aa__
-------------------------------------------

      PARDO ifrag, jfrag, mu, nu #GETLINE: Fragment_ij_aa__
      WHERE jfrag != ifrag
      where (int)elst_dist[ifrag,jfrag] == ifrag
      where (int)SwAO_frag[(index)mu] == ifrag
      where (int)SwAO_frag[(index)nu] == ifrag

*/

Fragment_ij_aa__PardoLoopManager::Fragment_ij_aa__PardoLoopManager(
		int num_indices, const int (&index_id)[MAX_RANK],
		DataManager & data_manager, const SipTables & sip_tables,
		SIPMPIAttr & sip_mpi_attr, int num_where_clauses,
		Interpreter* interpreter, long& iteration) :
		data_manager_(data_manager), sip_tables_(sip_tables), num_indices_(
				num_indices), first_time_(true), iteration_(iteration), sip_mpi_attr_(
				sip_mpi_attr), num_where_clauses_(num_where_clauses), company_rank_(
				sip_mpi_attr.company_rank()), num_workers_(
				sip_mpi_attr_.num_workers()), interpreter_(interpreter) {

	std::copy(index_id + 0, index_id + MAX_RANK, index_id_ + 0);
	for (int i = 0; i < num_indices; ++i) {
		lower_seg_[i] = sip_tables_.lower_seg(index_id_[i]);
		upper_bound_[i] = lower_seg_[i]
				+ sip_tables_.num_segments(index_id_[i]);
		sip::check(lower_seg_[i] < upper_bound_[i],
				"Pardo loop index " + sip_tables_.index_name(index_id_[i])
						+ " has empty range",
				Interpreter::global_interpreter->line_number());
	}
}

Fragment_ij_aa__PardoLoopManager::~Fragment_ij_aa__PardoLoopManager() {}


inline bool Fragment_ij_aa__PardoLoopManager::initialize_indices() {
	//initialize values of all indices
	bool more_iterations = true;
	for (int i = 0; i < num_indices_; ++i) {
		if (lower_seg_[i] >= upper_bound_[i]) {
			more_iterations = false; //this loop has an empty range in at least one dimension.
			return more_iterations;
		}
		sip::check(
				data_manager_.index_value(index_id_[i])
						== DataManager::undefined_index_value,
				"SIAL or SIP error, index "
						+ sip_tables_.index_name(index_id_[i])
						+ " already has value before loop",
				Interpreter::global_interpreter->line_number());
		data_manager_.set_index_value(index_id_[i], lower_seg_[i]);
		index_values_[i] = lower_seg_[i];
	}
	more_iterations = true;
	return more_iterations;
}

double Fragment_ij_aa__PardoLoopManager::return_val_elst_dist(int index1,int index2) {
	    int elst_dist_array_slot = sip_tables_.array_slot(std::string("elst_dist"));
	    const index_value_array_t indices_elst_dist = {index1, index2, unused_index_value, unused_index_value, unused_index_value, unused_index_value};
	    BlockId bid_elst_dist(elst_dist_array_slot, indices_elst_dist);
	    Block::BlockPtr bptr_elst_dist = data_manager_.block_manager_.get_block_for_reading(bid_elst_dist);
   	    double val_elst_dist = (int)(bptr_elst_dist->get_data()[0]);
	    return val_elst_dist;
}

double Fragment_ij_aa__PardoLoopManager::return_val_swao_frag(int index1) {
	    int swao_frag_array_slot = sip_tables_.array_slot(std::string("swao_frag"));
	    const index_value_array_t indices_swao_frag = {index1, unused_index_value, unused_index_value, unused_index_value, unused_index_value, unused_index_value};
	    BlockId bid_swao_frag(swao_frag_array_slot, indices_swao_frag);
	    Block::BlockPtr bptr_swao_frag = data_manager_.block_manager_.get_block_for_reading(bid_swao_frag);
	    double val_swao_frag = (int)(bptr_swao_frag->get_data()[0]);
	    return val_swao_frag;
}

double Fragment_ij_aa__PardoLoopManager::return_val_swocca_frag(int index1) {
	    int swocca_frag_array_slot = sip_tables_.array_slot(std::string("swocca_frag"));
	    const index_value_array_t indices_swocca_frag = {index1, unused_index_value, unused_index_value, unused_index_value, unused_index_value, unused_index_value};
	    BlockId bid_swocca_frag(swocca_frag_array_slot, indices_swocca_frag);
	    Block::BlockPtr bptr_swocca_frag = data_manager_.block_manager_.get_block_for_reading(bid_swocca_frag);
	    double val_swocca_frag = (int)(bptr_swocca_frag->get_data()[0]);
	    return val_swocca_frag;
}

double Fragment_ij_aa__PardoLoopManager::return_val_swvirta_frag(int index1) {
	    int swvirta_frag_array_slot = sip_tables_.array_slot(std::string("swvirta_frag"));
	    const index_value_array_t indices_swvirta_frag = {index1, unused_index_value, unused_index_value, unused_index_value, unused_index_value, unused_index_value};
	    BlockId bid_swvirta_frag(swvirta_frag_array_slot, indices_swvirta_frag);
	    Block::BlockPtr bptr_swvirta_frag = data_manager_.block_manager_.get_block_for_reading(bid_swvirta_frag);
	    double val_swvirta_frag = (int)(bptr_swvirta_frag->get_data()[0]);
	    return val_swvirta_frag;
}

bool Fragment_ij_aa__PardoLoopManager::increment_special(){
	// Go over indices in this order
	// ifrag, jfrag
	bool more_ifrag_jfrag = false;
	bool more = false;

restart_ij_aa___ifrag_jfrag:
	{
		more_ifrag_jfrag = false; 	// More iterations?
		for (int i = 0; i < 2; ++i) {
			int current_value = data_manager_.index_value(index_id_[i]);
			++current_value;
			if (current_value < upper_bound_[i]) {
				//increment current index and return
				data_manager_.set_index_value(index_id_[i], current_value);
				index_values_[i] = current_value;
				more_ifrag_jfrag = true;
				break;
			} else {
				//wrap around and handle next index
				data_manager_.set_index_value(index_id_[i], lower_seg_[i]);
				index_values_[i] = lower_seg_[i];
			}
		}

		double val_elst_dist = return_val_elst_dist(index_values_[0],index_values_[1]);

		// where ifrag != jfrag
		// where elst_dist[ifrag,jfrag] == ifrag
		if (!(index_values_[0] != index_values_[1] && val_elst_dist == index_values_[0])){
			if (more_ifrag_jfrag){
				goto restart_ij_aa___ifrag_jfrag;
			} else {
				goto return_increment_special_ij_aa__;
			}
		}
	}

restart_ij_aa___1: // ao in ifrag
	{
	    more = false;
	    int frag = 0;
	    int i = 2;
	    int current_value = data_manager_.index_value(index_id_[i]);
	    ++current_value;
	    if (current_value < upper_bound_[i]) {
		    //increment current index and return
		    data_manager_.set_index_value(index_id_[i], current_value);
		    index_values_[i] = current_value;
		    more = true;
	    } else {
		    //wrap around and handle next index
		    data_manager_.set_index_value(index_id_[i], lower_seg_[i]);
		    index_values_[i] = lower_seg_[i];
	    }

	    double val_swao_frag = return_val_swao_frag(index_values_[i]);
	    // if SwAO_frag[(index)] != ifrag break out of loop.
	    if (val_swao_frag != index_values_[frag]) { 
		    if (more) {
			    goto restart_ij_aa___1;
		    } else {
			    goto restart_ij_aa___ifrag_jfrag;
		    }
	    }
	}

restart_ij_aa___2: // ao in ifrag
	{
	    more = false;
	    int frag = 0;
	    int i = 3;
	    int current_value = data_manager_.index_value(index_id_[i]);
	    ++current_value;
	    if (current_value < upper_bound_[i]) {
		    //increment current index and return
		    data_manager_.set_index_value(index_id_[i], current_value);
		    index_values_[i] = current_value;
		    more = true;
	    } else {
		    //wrap around and handle next index
		    data_manager_.set_index_value(index_id_[i], lower_seg_[i]);
		    index_values_[i] = lower_seg_[i];
	    }

	    double val_swao_frag = return_val_swao_frag(index_values_[i]);
	    // if SwAO_frag[(index)] != ifrag break out of loop.
	    if (val_swao_frag != index_values_[frag]) { 
		    if (more) {
			    goto restart_ij_aa___2;
		    } else {
			    goto restart_ij_aa___ifrag_jfrag;
		    }
	    }
	}


return_increment_special_ij_aa__:
	return more_ifrag_jfrag;
}

bool Fragment_ij_aa__PardoLoopManager::do_update() {

	if (to_exit_)
		return false;
	bool more_iters = false;
	if (first_time_) {
		first_time_ = false;
		bool first_iteration = initialize_indices();
		if (first_iteration){
			bool do_first_iteration = interpreter_->interpret_where(num_where_clauses_);
			if (do_first_iteration)
				return true;
			else
				more_iters = increment_special();
		}
	} else {
		more_iters = increment_special();
	}

	while(more_iters){
		interpreter_->skip_where_clauses(num_where_clauses_);
		iteration_++;
		if ((iteration_-1) % num_workers_ == company_rank_){
			return true;
		}
		more_iters = increment_special();
	}
	return more_iters; //this should be false here
}

void Fragment_ij_aa__PardoLoopManager::do_finalize() {
	for (int i = 0; i < num_indices_; ++i) {
		data_manager_.set_index_undefined(index_id_[i]);
	}
}

std::string Fragment_ij_aa__PardoLoopManager::to_string() const {
	std::stringstream ss;
	ss << "Fragment_ij_aa__PardoLoopManager:  num_indices="
			<< num_indices_ << std::endl;
	ss << "index_ids_=[";
	for (int i = 0; i < num_indices_; ++i) {
		ss << (i == 0 ? "" : ",") << sip_tables_.index_name(index_id_[i]);
	}
	ss << "] lower_seg_=[";
	for (int i = 0; i < num_indices_; ++i) {
		ss << (i == 0 ? "" : ",") << lower_seg_[i];
	}
	ss << "] upper_bound_=[";
	for (int i = 0; i < num_indices_; ++i) {
		ss << (i == 0 ? "" : ",") << upper_bound_[i];
	}
	ss << "] current= [";
	for (int i = 0; i < num_indices_; ++i) {
		ss << (i == 0 ? "" : ",")
				<< data_manager_.index_value_to_string(index_id_[i]);
	}
	ss << "]";
	return ss.str();
}

std::ostream& operator<<(std::ostream& os,
		const Fragment_ij_aa__PardoLoopManager &obj) {
	os << obj.to_string();
	return os;
}

/*
-------------------------------------------
            _ij_ao_vo_
-------------------------------------------

      PARDO ifrag, jfrag, mu, i, b, j #GETLINE: Fragment_ij_ao_vo_
      where (int)elst_dist[ifrag,jfrag] == ifrag
      where (int)SwAO_frag[(index)mu] == ifrag
      where (int)SwOccA_frag[(index)i] == ifrag
      where (int)SwVirtA_frag[(index)b] == jfrag
      where (int)SwOccA_frag[(index)j] == jfrag

*/


Fragment_ij_ao_vo_PardoLoopManager::Fragment_ij_ao_vo_PardoLoopManager(
		int num_indices, const int (&index_id)[MAX_RANK],
		DataManager & data_manager, const SipTables & sip_tables,
		SIPMPIAttr & sip_mpi_attr, int num_where_clauses,
		Interpreter* interpreter, long& iteration) :
		data_manager_(data_manager), sip_tables_(sip_tables), num_indices_(
				num_indices), first_time_(true), iteration_(iteration), sip_mpi_attr_(
				sip_mpi_attr), num_where_clauses_(num_where_clauses), company_rank_(
				sip_mpi_attr.company_rank()), num_workers_(
				sip_mpi_attr_.num_workers()), interpreter_(interpreter) {

	std::copy(index_id + 0, index_id + MAX_RANK, index_id_ + 0);
	for (int i = 0; i < num_indices; ++i) {
		lower_seg_[i] = sip_tables_.lower_seg(index_id_[i]);
		upper_bound_[i] = lower_seg_[i]
				+ sip_tables_.num_segments(index_id_[i]);
		sip::check(lower_seg_[i] < upper_bound_[i],
				"Pardo loop index " + sip_tables_.index_name(index_id_[i])
						+ " has empty range",
				Interpreter::global_interpreter->line_number());
	}
}

Fragment_ij_ao_vo_PardoLoopManager::~Fragment_ij_ao_vo_PardoLoopManager() {}


inline bool Fragment_ij_ao_vo_PardoLoopManager::initialize_indices() {
	//initialize values of all indices
	bool more_iterations = true;
	for (int i = 0; i < num_indices_; ++i) {
		if (lower_seg_[i] >= upper_bound_[i]) {
			more_iterations = false; //this loop has an empty range in at least one dimension.
			return more_iterations;
		}
		sip::check(
				data_manager_.index_value(index_id_[i])
						== DataManager::undefined_index_value,
				"SIAL or SIP error, index "
						+ sip_tables_.index_name(index_id_[i])
						+ " already has value before loop",
				Interpreter::global_interpreter->line_number());
		data_manager_.set_index_value(index_id_[i], lower_seg_[i]);
		index_values_[i] = lower_seg_[i];
	}
	more_iterations = true;
	return more_iterations;
}

double Fragment_ij_ao_vo_PardoLoopManager::return_val_elst_dist(int index1,int index2) {
	    int elst_dist_array_slot = sip_tables_.array_slot(std::string("elst_dist"));
	    const index_value_array_t indices_elst_dist = {index1, index2, unused_index_value, unused_index_value, unused_index_value, unused_index_value};
	    BlockId bid_elst_dist(elst_dist_array_slot, indices_elst_dist);
	    Block::BlockPtr bptr_elst_dist = data_manager_.block_manager_.get_block_for_reading(bid_elst_dist);
   	    double val_elst_dist = (int)(bptr_elst_dist->get_data()[0]);
	    return val_elst_dist;
}

double Fragment_ij_ao_vo_PardoLoopManager::return_val_swao_frag(int index1) {
	    int swao_frag_array_slot = sip_tables_.array_slot(std::string("swao_frag"));
	    const index_value_array_t indices_swao_frag = {index1, unused_index_value, unused_index_value, unused_index_value, unused_index_value, unused_index_value};
	    BlockId bid_swao_frag(swao_frag_array_slot, indices_swao_frag);
	    Block::BlockPtr bptr_swao_frag = data_manager_.block_manager_.get_block_for_reading(bid_swao_frag);
	    double val_swao_frag = (int)(bptr_swao_frag->get_data()[0]);
	    return val_swao_frag;
}

double Fragment_ij_ao_vo_PardoLoopManager::return_val_swocca_frag(int index1) {
	    int swocca_frag_array_slot = sip_tables_.array_slot(std::string("swocca_frag"));
	    const index_value_array_t indices_swocca_frag = {index1, unused_index_value, unused_index_value, unused_index_value, unused_index_value, unused_index_value};
	    BlockId bid_swocca_frag(swocca_frag_array_slot, indices_swocca_frag);
	    Block::BlockPtr bptr_swocca_frag = data_manager_.block_manager_.get_block_for_reading(bid_swocca_frag);
	    double val_swocca_frag = (int)(bptr_swocca_frag->get_data()[0]);
	    return val_swocca_frag;
}

double Fragment_ij_ao_vo_PardoLoopManager::return_val_swvirta_frag(int index1) {
	    int swvirta_frag_array_slot = sip_tables_.array_slot(std::string("swvirta_frag"));
	    const index_value_array_t indices_swvirta_frag = {index1, unused_index_value, unused_index_value, unused_index_value, unused_index_value, unused_index_value};
	    BlockId bid_swvirta_frag(swvirta_frag_array_slot, indices_swvirta_frag);
	    Block::BlockPtr bptr_swvirta_frag = data_manager_.block_manager_.get_block_for_reading(bid_swvirta_frag);
	    double val_swvirta_frag = (int)(bptr_swvirta_frag->get_data()[0]);
	    return val_swvirta_frag;
}

bool Fragment_ij_ao_vo_PardoLoopManager::increment_special(){
	// Go over indices in this order
	// ifrag, jfrag
	bool more_ifrag_jfrag = false;
	bool more = false;

restart_ij_ao_vo__ifrag_jfrag:
	{
		more_ifrag_jfrag = false; 	// More iterations?
		for (int i = 0; i < 2; ++i) {
			int current_value = data_manager_.index_value(index_id_[i]);
			++current_value;
			if (current_value < upper_bound_[i]) {
				//increment current index and return
				data_manager_.set_index_value(index_id_[i], current_value);
				index_values_[i] = current_value;
				more_ifrag_jfrag = true;
				break;
			} else {
				//wrap around and handle next index
				data_manager_.set_index_value(index_id_[i], lower_seg_[i]);
				index_values_[i] = lower_seg_[i];
			}
		}

		double val_elst_dist = return_val_elst_dist(index_values_[0],index_values_[1]);

		// elst_dist[ifrag,jfrag] != ifrag
		if (val_elst_dist != index_values_[0]){
			if (more_ifrag_jfrag){
				goto restart_ij_ao_vo__ifrag_jfrag;
			} else {
				goto return_increment_special_ij_ao_vo_;
			}
		}
	}

restart_ij_ao_vo__1: // ao in ifrag
	{
	    more = false;
	    int frag = 0;
	    int i = 2;
	    int current_value = data_manager_.index_value(index_id_[i]);
	    ++current_value;
	    if (current_value < upper_bound_[i]) {
		    //increment current index and return
		    data_manager_.set_index_value(index_id_[i], current_value);
		    index_values_[i] = current_value;
		    more = true;
	    } else {
		    //wrap around and handle next index
		    data_manager_.set_index_value(index_id_[i], lower_seg_[i]);
		    index_values_[i] = lower_seg_[i];
	    }

	    double val_swao_frag = return_val_swao_frag(index_values_[i]);
	    // if SwAO_frag[(index)] != ifrag break out of loop.
	    if (val_swao_frag != index_values_[frag]) { 
		    if (more) {
			    goto restart_ij_ao_vo__1;
		    } else {
			    goto restart_ij_ao_vo__ifrag_jfrag;
		    }
	    }
	}

restart_ij_ao_vo__2: // occ in ifrag
	{
	    more = false;
	    int frag = 0;
	    int i = 3;
	    int current_value = data_manager_.index_value(index_id_[i]);
	    ++current_value;
	    if (current_value < upper_bound_[i]) {
		    //increment current index and return
		    data_manager_.set_index_value(index_id_[i], current_value);
		    index_values_[i] = current_value;
		    more = true;
	    } else {
		    //wrap around and handle next index
		    data_manager_.set_index_value(index_id_[i], lower_seg_[i]);
		    index_values_[i] = lower_seg_[i];
	    }

	    double val_swocca_frag = return_val_swocca_frag(index_values_[i]);
	    // if SwOccA_frag[(index)] != ifrag break out of loop.
	    if (val_swocca_frag != index_values_[frag]) { 
		    if (more) {
			    goto restart_ij_ao_vo__2;
		    } else {
			    goto restart_ij_ao_vo__1;
		    }
	    }
	}

restart_ij_ao_vo__3: // virt in jfrag
	{
	    more = false;
	    int frag = 1;
	    int i = 4;
	    int current_value = data_manager_.index_value(index_id_[i]);
	    ++current_value;
	    if (current_value < upper_bound_[i]) {
		    //increment current index and return
		    data_manager_.set_index_value(index_id_[i], current_value);
		    index_values_[i] = current_value;
		    more = true;
	    } else {
		    //wrap around and handle next index
		    data_manager_.set_index_value(index_id_[i], lower_seg_[i]);
		    index_values_[i] = lower_seg_[i];
	    }

	    double val_swvirta_frag = return_val_swvirta_frag(index_values_[i]);
	    // if SwVirtA_frag[(index)] != jfrag break out of loop.
	    if (val_swvirta_frag != index_values_[frag]) { 
		    if (more) {
			    goto restart_ij_ao_vo__3;
		    } else {
			    goto restart_ij_ao_vo__2;
		    }
	    }
	}

restart_ij_ao_vo__4: // occ in jfrag
	{
	    more = false;
	    int frag = 1;
	    int i = 5;
	    int current_value = data_manager_.index_value(index_id_[i]);
	    ++current_value;
	    if (current_value < upper_bound_[i]) {
		    //increment current index and return
		    data_manager_.set_index_value(index_id_[i], current_value);
		    index_values_[i] = current_value;
		    more = true;
	    } else {
		    //wrap around and handle next index
		    data_manager_.set_index_value(index_id_[i], lower_seg_[i]);
		    index_values_[i] = lower_seg_[i];
	    }

	    double val_swocca_frag = return_val_swocca_frag(index_values_[i]);
	    // if SwOccA_frag[(index)] != jfrag break out of loop.
	    if (val_swocca_frag != index_values_[frag]) { 
		    if (more) {
			    goto restart_ij_ao_vo__4;
		    } else {
			    goto restart_ij_ao_vo__3;
		    }
	    }
	}


return_increment_special_ij_ao_vo_:
	return more_ifrag_jfrag;
}

bool Fragment_ij_ao_vo_PardoLoopManager::do_update() {

	if (to_exit_)
		return false;
	bool more_iters = false;
	if (first_time_) {
		first_time_ = false;
		bool first_iteration = initialize_indices();
		if (first_iteration){
			bool do_first_iteration = interpreter_->interpret_where(num_where_clauses_);
			if (do_first_iteration)
				return true;
			else
				more_iters = increment_special();
		}
	} else {
		more_iters = increment_special();
	}

	while(more_iters){
		interpreter_->skip_where_clauses(num_where_clauses_);
		iteration_++;
		if ((iteration_-1) % num_workers_ == company_rank_){
			return true;
		}
		more_iters = increment_special();
	}
	return more_iters; //this should be false here
}

void Fragment_ij_ao_vo_PardoLoopManager::do_finalize() {
	for (int i = 0; i < num_indices_; ++i) {
		data_manager_.set_index_undefined(index_id_[i]);
	}
}

std::string Fragment_ij_ao_vo_PardoLoopManager::to_string() const {
	std::stringstream ss;
	ss << "Fragment_ij_ao_vo_PardoLoopManager:  num_indices="
			<< num_indices_ << std::endl;
	ss << "index_ids_=[";
	for (int i = 0; i < num_indices_; ++i) {
		ss << (i == 0 ? "" : ",") << sip_tables_.index_name(index_id_[i]);
	}
	ss << "] lower_seg_=[";
	for (int i = 0; i < num_indices_; ++i) {
		ss << (i == 0 ? "" : ",") << lower_seg_[i];
	}
	ss << "] upper_bound_=[";
	for (int i = 0; i < num_indices_; ++i) {
		ss << (i == 0 ? "" : ",") << upper_bound_[i];
	}
	ss << "] current= [";
	for (int i = 0; i < num_indices_; ++i) {
		ss << (i == 0 ? "" : ",")
				<< data_manager_.index_value_to_string(index_id_[i]);
	}
	ss << "]";
	return ss.str();
}

std::ostream& operator<<(std::ostream& os,
		const Fragment_ij_ao_vo_PardoLoopManager &obj) {
	os << obj.to_string();
	return os;
}

#endif

} /* namespace sip */
