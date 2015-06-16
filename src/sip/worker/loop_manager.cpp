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
// ----------------------------------------------------------------------------------------
/*!
 Fragment_ special naming scheme: _<simple index fragment indices>_<occ (o), virt (v), ao (a) in ifrag>_<o,v,a in jfrag>
*/

/*!
 -------------------------------------------
 _ij_aa__
 -------------------------------------------
 
 PARDO ifrag, jfrag, mu, nu #GETLINE: Fragment_ij_aa__
 WHERE jfrag != ifrag
 where (int)elst_dist[ifrag,jfrag] == ifrag
 where (int)SwAO_frag[(index)mu] == ifrag
 where (int)SwAO_frag[(index)nu] == ifrag
 
 */
    
    /*
     for each new special fragment where clause pattern, this should be the only thing realy changed.
     see comment above for fragment_special_where_clause syntax
     */
    bool Fragment_ij_aa__PardoLoopManager::where_clause(int index) {
        bool where_;
        int ifrag = 0;
        int jfrag = 1;
        int ao = 1;
        int occ = 2;
        int virt = 3;
        int elst = 4;
        int rcut = 5;
        int NE = 0;

        switch (index) {
            case 0:
            case 1:
                where_ = fragment_special_where_clause(elst,jfrag,ifrag) && fragment_special_where_clause(NE,jfrag,ifrag);
                break;
            case 2:
                where_ = fragment_special_where_clause(ao,index,ifrag);
                break;
            case 3:
                where_ = fragment_special_where_clause(ao,index,ifrag);
                break;
            default:
                where_ = false;
        }
        return where_;
    }
    
    bool Fragment_ij_aa__PardoLoopManager::do_update() {
        if (to_exit_)
            return false;
        bool more_iters;
        bool where_clauses_value;
        
        interpreter_->skip_where_clauses(num_where_clauses_);
        
        if (first_time_) {
            first_time_ = false;
            more_iters = initialize_indices();
            bool tmp = increment_single_index(1);
            
            where_clauses_value = true;
            for (int i = 1; i < num_indices_; ++i) {
                where_clauses_value = where_clauses_value && where_clause(i);
            }
            if(where_clauses_value){
                iteration_++;
                if ((iteration_-1) % num_workers_ == company_rank_){
                    return true;
                }
            }
        }
        
        int index_restart[MAX_RANK];
        
        for (int i = 0; i < num_indices_; ++i) {
            index_restart[i] = index_values_[i];
        }
        
        int loop_count = iteration_;
        for (int index_i = index_restart[0]; index_i < upper_bound_[0]; ++index_i) {
            index_values_[0] = index_i;
            data_manager_.set_index_value(index_id_[0], index_i);
            
        for (int index_j = index_restart[1]; index_j < upper_bound_[1]; ++index_j) {
            index_values_[1] = index_j;
            data_manager_.set_index_value(index_id_[1], index_j);
            
            if (where_clause(1)) {
                for (int index_2 = index_restart[2]; index_2 < upper_bound_[2]; ++index_2) {
                    index_values_[2] = index_2;
                    data_manager_.set_index_value(index_id_[2], index_2);
                        
            if (where_clause(2)) {
                for (int index_3 = index_restart[3]; index_3 < upper_bound_[3]; ++index_3) {
                    index_values_[3] = index_3;
                    data_manager_.set_index_value(index_id_[3], index_3);
                    
            if (where_clause(3)) {
                /*
                for (int i = 0; i < num_indices_; ++i) {
                    std::cout << index_values_[i] << ", ";
                }
                std::cout << loop_count << " " << iteration_ << std::endl;
                */

            if (loop_count > iteration_) {
                iteration_++;
                if ((iteration_-1) % num_workers_ == company_rank_){
                    return true;
                }
            }
            ++loop_count;
        } // where 3
            } // index_3
            for (int i = 3; i < num_indices_; ++i) { index_restart[i] = lower_seg_[i]; }
        } // where 2
            } // index_2
            for (int i = 2; i < num_indices_; ++i) { index_restart[i] = lower_seg_[i]; }
        } // where 1
        } // index_j
            for (int i = 1; i < num_indices_; ++i) { index_restart[i] = lower_seg_[i]; }
        } // index_i
        
        return false; //this should be false here
    }


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
    
    form_elst_dist();
    //form_rcut_dist();
    form_swao_frag();
    //form_swocca_frag();
    //form_swvirta_frag();
}

Fragment_ij_aa__PardoLoopManager::~Fragment_ij_aa__PardoLoopManager() {}

bool Fragment_ij_aa__PardoLoopManager::increment_simple_pair(){
    bool more = false;      // More iterations?
    int current_value;
    for (int i = 0; i < 2; ++i) {
        more = increment_single_index(i);
        if (more) {break;}
    } //if here, then all indices are at their max value
    return more;
}

bool Fragment_ij_aa__PardoLoopManager::increment_single_index(int index){
    bool more = false; 	// More iterations?
    int current_value;
    current_value = data_manager_.index_value(index_id_[index]);
    ++current_value;
    if (current_value < upper_bound_[index]) {
        //increment current index and return
        index_values_[index] = current_value;
        data_manager_.set_index_value(index_id_[index], current_value);
        more = true;
    } else {
        //wrap around and handle next index
        index_values_[index] = lower_seg_[index];
        data_manager_.set_index_value(index_id_[index], lower_seg_[index]);
    }
    return more;
}
    
    
    inline bool Fragment_ij_aa__PardoLoopManager::increment_all() {
        bool more = false;      // More iterations?
        int current_value;
        for (int i = 0; i < num_indices_; ++i) {
            more = increment_single_index(i);
            if (more) {break;}
        } //if here, then all indices are at their max value
        return more;
    }
    
    void Fragment_ij_aa__PardoLoopManager::form_elst_dist() {
        int elst_dist_array_slot = sip_tables_.array_slot(std::string("elst_dist"));
        Block::BlockPtr bptr_elst_dist = data_manager_.contiguous_array_manager().get_array(elst_dist_array_slot);
        double *val_elst_dist = bptr_elst_dist->get_data();
        int elst_dist_size = bptr_elst_dist->size();
        int num_frag = upper_bound_[0] - lower_seg_[0];
        elst_dist.resize(num_frag);
        for (int i = 0; i < num_frag; ++i) {
            //std::cout << i << std::endl;
            elst_dist[i].resize(num_frag);
        }
        int ij = 0;
        for (int i = 0; i < num_frag; ++i) {
            for (int j = 0; j < num_frag; ++j) {
                elst_dist[j][i] = (int)val_elst_dist[ij];
                //std::cout << i << " " << j << " " << ij << " " << elst_dist[j][i] << " " << val_elst_dist[ij]<< std::endl;
                ++ij;
            }
        }
        return;
    }

    void Fragment_ij_aa__PardoLoopManager::form_rcut_dist() {
        int rcut_dist_array_slot = sip_tables_.array_slot(std::string("rcut_dist"));
        Block::BlockPtr bptr_rcut_dist = data_manager_.contiguous_array_manager().get_array(rcut_dist_array_slot);
        double *val_rcut_dist = bptr_rcut_dist->get_data();
        int rcut_dist_size = bptr_rcut_dist->size();
        int num_frag = upper_bound_[0] - lower_seg_[0];
        rcut_dist.resize(num_frag);
        for (int i = 0; i < num_frag; ++i) {
            //std::cout << i << std::endl;
            rcut_dist[i].resize(num_frag);
        }
        int ij = 0;
        for (int i = 0; i < num_frag; ++i) {
            for (int j = 0; j < num_frag; ++j) {
                rcut_dist[j][i] = (int)val_rcut_dist[ij];
                //std::cout << i << " " << j << " " << ij << " " << rcut_dist[j][i] << " " << val_rcut_dist[ij]<< std::endl;
                ++ij;
            }
        }
        return;
    }
    
    void Fragment_ij_aa__PardoLoopManager::form_swao_frag() {
        int swao_frag_array_slot = sip_tables_.array_slot(std::string("swao_frag"));
        Block::BlockPtr bptr_swao_frag = data_manager_.contiguous_array_manager().get_array(swao_frag_array_slot);
        double *val_swao_frag = bptr_swao_frag->get_data();
        int swao_frag_size = bptr_swao_frag->size();
        swao_frag.resize(swao_frag_size);
        for (int i = 0; i < swao_frag_size; ++i) {
            swao_frag[i] = (int)val_swao_frag[i];
            //std::cout << "vec " << swao_frag[i] << " data " << val_swao_frag[i] << std::endl;
        }
        return;
    }
    
    void Fragment_ij_aa__PardoLoopManager::form_swocca_frag() {
        int swocca_frag_array_slot = sip_tables_.array_slot(std::string("swocca_frag"));
        Block::BlockPtr bptr_swocca_frag = data_manager_.contiguous_array_manager().get_array(swocca_frag_array_slot);
        double *val_swocca_frag = bptr_swocca_frag->get_data();
        int swocca_frag_size = bptr_swocca_frag->size();
        swocca_frag.resize(swocca_frag_size);
        for (int i = 0; i < swocca_frag_size; ++i) {
            swocca_frag[i] = (int)val_swocca_frag[i];
            //std::cout << "vec " << swocca_frag[i] << " data " << val_swocca_frag[i] << std::endl;
        }
        return;
    }
    
    void Fragment_ij_aa__PardoLoopManager::form_swvirta_frag() {
        int swvirta_frag_array_slot = sip_tables_.array_slot(std::string("swvirta_frag"));
        Block::BlockPtr bptr_swvirta_frag = data_manager_.contiguous_array_manager().get_array(swvirta_frag_array_slot);
        double *val_swvirta_frag = bptr_swvirta_frag->get_data();
        int swvirta_frag_size = bptr_swvirta_frag->size();
        swvirta_frag.resize(swvirta_frag_size);
        for (int i = 0; i < swvirta_frag_size; ++i) {
            swvirta_frag[i] = (int)val_swvirta_frag[i];
            //std::cout << "vec " << swvirta_frag[i] << " data " << val_swvirta_frag[i] << std::endl;
        }
        return;
    }
    
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
            index_values_[i] = lower_seg_[i];
            data_manager_.set_index_value(index_id_[i], lower_seg_[i]);
        }
        more_iterations = true;
        return more_iterations;
    }

/*!
 the generic where clauses are:
 typ = 1: map AO to fragment
 typ = 2: map Occ to fragment
 typ = 3: map Virt to fragment
 typ = 4: map elst_dist(ifrag,jfrag) == ifrag
 typ = 5: map rcut_dist(ifrag,jfrag) == ifrag
 typ = 0: ifrag != jfrag
 */
bool Fragment_ij_aa__PardoLoopManager::fragment_special_where_clause(int typ, int index, int frag) {
    bool where_clause;
    int ij = 0;
    switch (typ) {
            
        case 1: // check against ao index
            //where_clause = return_val_swao_frag(index_values_[index])     == index_values_[frag];
            where_clause = swao_frag[index_values_[index] - lower_seg_[index]]     == index_values_[frag];
            break;
            
        case 2: // check against occ index
            //where_clause = return_val_swocca_frag(index_values_[index])   == index_values_[frag];
            where_clause = swocca_frag[index_values_[index] - lower_seg_[index]]   == index_values_[frag];
            break;
            
        case 3: // check against virt index
            //where_clause = return_val_swvirta_frag(index_values_[index])  == index_values_[frag];
            where_clause = swvirta_frag[index_values_[index] - lower_seg_[index]]  == index_values_[frag];
            break;
            
        case 4: // check elst_dist[ifrag,jfrag] == ifrag
            //where_clause = return_val_elst_dist(index_values_[frag], index_values_[index]) == index_values_[frag];
            //ij = index_values_[frag] - 1 + (upper_bound_[0])*(index_values_[index] - 1);
            //std::cout << index_values_[index] << " " << index_values_[frag] << " " << ij << " " << elst_dist[ij] << std::endl;
            //where_clause = elst_dist[ij] == index_values_[frag];
            where_clause = elst_dist[index_values_[frag] - lower_seg_[frag]][index_values_[index] - lower_seg_[index]] == index_values_[frag];
            break;
            
        case 5: // check rcut_dist[ifrag,jfrag] == ifrag
            //where_clause = return_val_rcut_dist(index_values_[frag], index_values_[index]) == index_values_[frag];
            //ij = index_values_[index] - 1 + upper_bound_[0]*(index_values_[frag] - 1);
            //where_clause = rcut_dist[ij] == index_values_[frag];
            //break;
            where_clause = rcut_dist[index_values_[frag] - lower_seg_[frag]][index_values_[index] - lower_seg_[index]] == index_values_[frag];
            break;
            
        case 0:
            where_clause = index_values_[frag] != index_values_[index];
            break;
        default:
            where_clause = false;
    }
    //std::cout << "where_clause " << typ << " " << index_values_[index] << " " << index_values_[frag] << " " << where_clause << std::endl;
    return where_clause;
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


/*!
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

    
    /*
     for each new special fragment where clause pattern, this should be the only thing realy changed.
     see comment above for fragment_special_where_clause syntax
     */
    bool Fragment_ij_ao_vo_PardoLoopManager::where_clause(int index) {
        bool where_;
        int ifrag = 0;
        int jfrag = 1;
        int ao = 1;
        int occ = 2;
        int virt = 3;
        int elst = 4;
        int rcut = 5;

        
        switch (index) {
            case 0:
            case 1:
                where_ = fragment_special_where_clause(elst,jfrag,ifrag);
                break;
            case 2:
                where_ = fragment_special_where_clause(ao,index,ifrag);
                break;
            case 3:
                where_ = fragment_special_where_clause(occ,index,ifrag);
                break;
            case 4:
                where_ = fragment_special_where_clause(virt,index,jfrag);
                break;
            case 5:
                where_ = fragment_special_where_clause(occ,index,jfrag);
                break;
            default:
                where_ = false;
        }
        return where_;
    }
    
    bool Fragment_ij_ao_vo_PardoLoopManager::do_update() {
        if (to_exit_)
            return false;
        bool more_iters;
        bool where_clauses_value;
        
        interpreter_->skip_where_clauses(num_where_clauses_);
        
        if (first_time_) {
            first_time_ = false;
            more_iters = initialize_indices();
            
            where_clauses_value = true;
            for (int i = 1; i < num_indices_; ++i) {
                where_clauses_value = where_clauses_value && where_clause(i);
            }
            if(where_clauses_value){
                iteration_++;
                if ((iteration_-1) % num_workers_ == company_rank_){
                    return true;
                }
            }
        }
        
        int index_restart[MAX_RANK];
        
        for (int i = 0; i < num_indices_; ++i) {
            index_restart[i] = index_values_[i];
        }
        
        int loop_count = iteration_;
        for (int index_i = index_restart[0]; index_i < upper_bound_[0]; ++index_i) {
            index_values_[0] = index_i;
            data_manager_.set_index_value(index_id_[0], index_i);
            
        for (int index_j = index_restart[1]; index_j < upper_bound_[1]; ++index_j) {
            index_values_[1] = index_j;
            data_manager_.set_index_value(index_id_[1], index_j);
                
            if (where_clause(1)) {
                for (int index_2 = index_restart[2]; index_2 < upper_bound_[2]; ++index_2) {
                    index_values_[2] = index_2;
                    data_manager_.set_index_value(index_id_[2], index_2);
                        
            if (where_clause(2)) {
                for (int index_3 = index_restart[3]; index_3 < upper_bound_[3]; ++index_3) {
                    index_values_[3] = index_3;
                    data_manager_.set_index_value(index_id_[3], index_3);
                    
            if (where_clause(3)) {
                for (int index_4 = index_restart[4]; index_4 < upper_bound_[4]; ++index_4) {
                    index_values_[4] = index_4;
                    data_manager_.set_index_value(index_id_[4], index_4);
                    
            if (where_clause(4)) {
                for (int index_5 = index_restart[5]; index_5 < upper_bound_[5]; ++index_5) {
                    index_values_[5] = index_5;
                    data_manager_.set_index_value(index_id_[5], index_5);
                    
            if (where_clause(5)) {
            if (loop_count > iteration_) {
                iteration_++;
                if ((iteration_-1) % num_workers_ == company_rank_){
                    return true;
                }
            }
            ++loop_count;
        } // where 5
            } // index_5
            for (int i = 5; i < num_indices_; ++i) { index_restart[i] = lower_seg_[i]; }
        } // where 4
            } // index_4
            for (int i = 4; i < num_indices_; ++i) { index_restart[i] = lower_seg_[i]; }
        } // where 3
            } // index_3
            for (int i = 3; i < num_indices_; ++i) { index_restart[i] = lower_seg_[i]; }
        } // where 2
            } // index_2
            for (int i = 2; i < num_indices_; ++i) { index_restart[i] = lower_seg_[i]; }
        } // where 1
            } // index_j
            for (int i = 1; i < num_indices_; ++i) { index_restart[i] = lower_seg_[i]; }
        } // index_i
        
        return false; //this should be false here
    }


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
    
    form_elst_dist();
    //form_rcut_dist();
    form_swao_frag();
    form_swocca_frag();
    form_swvirta_frag();
}

Fragment_ij_ao_vo_PardoLoopManager::~Fragment_ij_ao_vo_PardoLoopManager() {}

bool Fragment_ij_ao_vo_PardoLoopManager::increment_simple_pair(){
    bool more = false;      // More iterations?
    int current_value;
    for (int i = 0; i < 2; ++i) {
        more = increment_single_index(i);
        if (more) {break;}
    } //if here, then all indices are at their max value
    return more;
}

bool Fragment_ij_ao_vo_PardoLoopManager::increment_single_index(int index){
    bool more = false; 	// More iterations?
    int current_value;
    current_value = data_manager_.index_value(index_id_[index]);
    ++current_value;
    if (current_value < upper_bound_[index]) {
        //increment current index and return
        index_values_[index] = current_value;
        data_manager_.set_index_value(index_id_[index], current_value);
        more = true;
    } else {
        //wrap around and handle next index
        index_values_[index] = lower_seg_[index];
        data_manager_.set_index_value(index_id_[index], lower_seg_[index]);
    }
    return more;
}
    
    
    inline bool Fragment_ij_ao_vo_PardoLoopManager::increment_all() {
        bool more = false;      // More iterations?
        int current_value;
        for (int i = 0; i < num_indices_; ++i) {
            more = increment_single_index(i);
            if (more) {break;}
        } //if here, then all indices are at their max value
        return more;
    }
    
    void Fragment_ij_ao_vo_PardoLoopManager::form_elst_dist() {
        int elst_dist_array_slot = sip_tables_.array_slot(std::string("elst_dist"));
        Block::BlockPtr bptr_elst_dist = data_manager_.contiguous_array_manager().get_array(elst_dist_array_slot);
        double *val_elst_dist = bptr_elst_dist->get_data();
        int elst_dist_size = bptr_elst_dist->size();
        int num_frag = upper_bound_[0] - lower_seg_[0];
        elst_dist.resize(num_frag);
        for (int i = 0; i < num_frag; ++i) {
            //std::cout << i << std::endl;
            elst_dist[i].resize(num_frag);
        }
        int ij = 0;
        for (int i = 0; i < num_frag; ++i) {
            for (int j = 0; j < num_frag; ++j) {
                elst_dist[j][i] = (int)val_elst_dist[ij];
                //std::cout << i << " " << j << " " << ij << " " << elst_dist[j][i] << " " << val_elst_dist[ij]<< std::endl;
                ++ij;
            }
        }
        return;
    }

    void Fragment_ij_ao_vo_PardoLoopManager::form_rcut_dist() {
        int rcut_dist_array_slot = sip_tables_.array_slot(std::string("rcut_dist"));
        Block::BlockPtr bptr_rcut_dist = data_manager_.contiguous_array_manager().get_array(rcut_dist_array_slot);
        double *val_rcut_dist = bptr_rcut_dist->get_data();
        int rcut_dist_size = bptr_rcut_dist->size();
        int num_frag = upper_bound_[0] - lower_seg_[0];
        rcut_dist.resize(num_frag);
        for (int i = 0; i < num_frag; ++i) {
            //std::cout << i << std::endl;
            rcut_dist[i].resize(num_frag);
        }
        int ij = 0;
        for (int i = 0; i < num_frag; ++i) {
            for (int j = 0; j < num_frag; ++j) {
                rcut_dist[j][i] = (int)val_rcut_dist[ij];
                //std::cout << i << " " << j << " " << ij << " " << rcut_dist[j][i] << " " << val_rcut_dist[ij]<< std::endl;
                ++ij;
            }
        }
        return;
    }
    
    void Fragment_ij_ao_vo_PardoLoopManager::form_swao_frag() {
        int swao_frag_array_slot = sip_tables_.array_slot(std::string("swao_frag"));
        Block::BlockPtr bptr_swao_frag = data_manager_.contiguous_array_manager().get_array(swao_frag_array_slot);
        double *val_swao_frag = bptr_swao_frag->get_data();
        int swao_frag_size = bptr_swao_frag->size();
        swao_frag.resize(swao_frag_size);
        for (int i = 0; i < swao_frag_size; ++i) {
            swao_frag[i] = (int)val_swao_frag[i];
            //std::cout << "vec " << swao_frag[i] << " data " << val_swao_frag[i] << std::endl;
        }
        return;
    }
    
    void Fragment_ij_ao_vo_PardoLoopManager::form_swocca_frag() {
        int swocca_frag_array_slot = sip_tables_.array_slot(std::string("swocca_frag"));
        Block::BlockPtr bptr_swocca_frag = data_manager_.contiguous_array_manager().get_array(swocca_frag_array_slot);
        double *val_swocca_frag = bptr_swocca_frag->get_data();
        int swocca_frag_size = bptr_swocca_frag->size();
        swocca_frag.resize(swocca_frag_size);
        for (int i = 0; i < swocca_frag_size; ++i) {
            swocca_frag[i] = (int)val_swocca_frag[i];
            //std::cout << "vec " << swocca_frag[i] << " data " << val_swocca_frag[i] << std::endl;
        }
        return;
    }
    
    void Fragment_ij_ao_vo_PardoLoopManager::form_swvirta_frag() {
        int swvirta_frag_array_slot = sip_tables_.array_slot(std::string("swvirta_frag"));
        Block::BlockPtr bptr_swvirta_frag = data_manager_.contiguous_array_manager().get_array(swvirta_frag_array_slot);
        double *val_swvirta_frag = bptr_swvirta_frag->get_data();
        int swvirta_frag_size = bptr_swvirta_frag->size();
        swvirta_frag.resize(swvirta_frag_size);
        for (int i = 0; i < swvirta_frag_size; ++i) {
            swvirta_frag[i] = (int)val_swvirta_frag[i];
            //std::cout << "vec " << swvirta_frag[i] << " data " << val_swvirta_frag[i] << std::endl;
        }
        return;
    }
    
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
            index_values_[i] = lower_seg_[i];
            data_manager_.set_index_value(index_id_[i], lower_seg_[i]);
        }
        more_iterations = true;
        return more_iterations;
    }

/*
 the generic where clauses are:
 typ = 1: map AO to fragment
 typ = 2: map Occ to fragment
 typ = 3: map Virt to fragment
 typ = 4: map elst_dist(ifrag,jfrag) == ifrag
 typ = 5: map rcut_dist(ifrag,jfrag) == ifrag
 typ = 0: ifrag != jfrag
 */
bool Fragment_ij_ao_vo_PardoLoopManager::fragment_special_where_clause(int typ, int index, int frag) {
    bool where_clause;
    int ij = 0;
    switch (typ) {
            
        case 1: // check against ao index
            //where_clause = return_val_swao_frag(index_values_[index])     == index_values_[frag];
            where_clause = swao_frag[index_values_[index] - lower_seg_[index]]     == index_values_[frag];
            break;
            
        case 2: // check against occ index
            //where_clause = return_val_swocca_frag(index_values_[index])   == index_values_[frag];
            where_clause = swocca_frag[index_values_[index] - lower_seg_[index]]   == index_values_[frag];
            break;
            
        case 3: // check against virt index
            //where_clause = return_val_swvirta_frag(index_values_[index])  == index_values_[frag];
            where_clause = swvirta_frag[index_values_[index] - lower_seg_[index]]  == index_values_[frag];
            break;
            
        case 4: // check elst_dist[ifrag,jfrag] == ifrag
            //where_clause = return_val_elst_dist(index_values_[frag], index_values_[index]) == index_values_[frag];
            //ij = index_values_[frag] - 1 + (upper_bound_[0])*(index_values_[index] - 1);
            //std::cout << index_values_[index] << " " << index_values_[frag] << " " << ij << " " << elst_dist[ij] << std::endl;
            //where_clause = elst_dist[ij] == index_values_[frag];
            where_clause = elst_dist[index_values_[frag] - lower_seg_[frag]][index_values_[index] - lower_seg_[index]] == index_values_[frag];
            break;
            
        case 5: // check rcut_dist[ifrag,jfrag] == ifrag
            //where_clause = return_val_rcut_dist(index_values_[frag], index_values_[index]) == index_values_[frag];
            //ij = index_values_[index] - 1 + upper_bound_[0]*(index_values_[frag] - 1);
            //where_clause = rcut_dist[ij] == index_values_[frag];
            //break;
            where_clause = rcut_dist[index_values_[frag] - lower_seg_[frag]][index_values_[index] - lower_seg_[index]] == index_values_[frag];
            break;
            
        case 0:
            where_clause = index_values_[frag] != index_values_[index];
            break;
        default:
            where_clause = false;
    }
    //std::cout << "where_clause " << typ << " " << index_values_[index] << " " << index_values_[frag] << " " << where_clause << std::endl;
    return where_clause;
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
    
    /*!
     -------------------------------------------
     _i_aa__
     -------------------------------------------
     
     PARDO ifrag, mu,....#GETLINE: Fragment_i_aa__
     where (int)SwAO_frag[(index)mu] == ifrag
     .
     .
     .
     */
    
    
    /*
     for each new special fragment where clause pattern, this should be the only thing realy changed.
     see comment above for fragment_special_where_clause syntax
     */
    bool Fragment_i_aa__PardoLoopManager::where_clause(int index) {
        bool where_;
        int ifrag = 0;
        int jfrag = 1;
        int ao = 1;
        int occ = 2;
        int virt = 3;
        int elst = 4;
        int rcut = 5;
        
        
        switch (index) {
            case 0:
                where_ = true;
            case 1:
            case 2:
            case 3:
            case 4:
            case 5:
                where_ = fragment_special_where_clause(ao,index,ifrag);
                break;
            default:
                where_ = false;
        }
        return where_;
    }
    
    bool Fragment_i_aa__PardoLoopManager::do_update() {
        if (to_exit_)
            return false;
        bool more_iters;
        bool where_clauses_value;
        
        interpreter_->skip_where_clauses(num_where_clauses_);
        
        if (first_time_) {
            first_time_ = false;
            more_iters = initialize_indices();
            
            where_clauses_value = true;
            for (int i = 1; i < num_indices_; ++i) {
                where_clauses_value = where_clauses_value && where_clause(i);
            }
            if(where_clauses_value){
                iteration_++;
                if ((iteration_-1) % num_workers_ == company_rank_){
                    return true;
                }
            }
        }
        
        int index_restart[MAX_RANK];
        
        for (int i = 0; i < num_indices_; ++i) {
            index_restart[i] = index_values_[i];
        }
        
        int loop_count = iteration_;
        for (int index_i = index_restart[0]; index_i < upper_bound_[0]; ++index_i) {
            index_values_[0] = index_i;
            data_manager_.set_index_value(index_id_[0], index_i);
                
                for (int index_1 = index_restart[1]; index_1 < upper_bound_[1]; ++index_1) {
                    index_values_[1] = index_1;
                    data_manager_.set_index_value(index_id_[1], index_1);
                        
            if (where_clause(1)) {
                for (int index_2 = index_restart[2]; index_2 < upper_bound_[2]; ++index_2) {
                    index_values_[2] = index_2;
                    data_manager_.set_index_value(index_id_[2], index_2);
                    
            if (where_clause(2)) {
            if (loop_count > iteration_) {
                iteration_++;
                if ((iteration_-1) % num_workers_ == company_rank_){
                    return true;
                }
            }
            ++loop_count;
        } // where 2
            } // index_2
            for (int i = 2; i < num_indices_; ++i) { index_restart[i] = lower_seg_[i]; }
        } // where 1
            } // index_1
            for (int i = 1; i < num_indices_; ++i) { index_restart[i] = lower_seg_[i]; }
        } // index_i
        
        return false; //this should be false here
    }
    
    Fragment_i_aa__PardoLoopManager::Fragment_i_aa__PardoLoopManager(
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
        
        //form_elst_dist();
        //form_rcut_dist();
        form_swao_frag();
        //form_swocca_frag();
        //form_swvirta_frag();
    }
    
    Fragment_i_aa__PardoLoopManager::~Fragment_i_aa__PardoLoopManager() {}
    
    bool Fragment_i_aa__PardoLoopManager::increment_simple_pair(){
        bool more = false;      // More iterations?
        int current_value;
        for (int i = 0; i < 2; ++i) {
            more = increment_single_index(i);
            if (more) {break;}
        } //if here, then all indices are at their max value
        return more;
    }
    
    bool Fragment_i_aa__PardoLoopManager::increment_single_index(int index){
        bool more = false; 	// More iterations?
        int current_value;
        current_value = data_manager_.index_value(index_id_[index]);
        ++current_value;
        if (current_value < upper_bound_[index]) {
            //increment current index and return
            index_values_[index] = current_value;
            data_manager_.set_index_value(index_id_[index], current_value);
            more = true;
        } else {
            //wrap around and handle next index
            index_values_[index] = lower_seg_[index];
            data_manager_.set_index_value(index_id_[index], lower_seg_[index]);
        }
        return more;
    }
    
    
    inline bool Fragment_i_aa__PardoLoopManager::increment_all() {
        bool more = false;      // More iterations?
        int current_value;
        for (int i = 0; i < num_indices_; ++i) {
            more = increment_single_index(i);
            if (more) {break;}
        } //if here, then all indices are at their max value
        return more;
    }
    
    void Fragment_i_aa__PardoLoopManager::form_elst_dist() {
        int elst_dist_array_slot = sip_tables_.array_slot(std::string("elst_dist"));
        Block::BlockPtr bptr_elst_dist = data_manager_.contiguous_array_manager().get_array(elst_dist_array_slot);
        double *val_elst_dist = bptr_elst_dist->get_data();
        int elst_dist_size = bptr_elst_dist->size();
        int num_frag = upper_bound_[0] - lower_seg_[0];
        elst_dist.resize(num_frag);
        for (int i = 0; i < num_frag; ++i) {
            //std::cout << i << std::endl;
            elst_dist[i].resize(num_frag);
        }
        int ij = 0;
        for (int i = 0; i < num_frag; ++i) {
            for (int j = 0; j < num_frag; ++j) {
                elst_dist[j][i] = (int)val_elst_dist[ij];
                //std::cout << i << " " << j << " " << ij << " " << elst_dist[j][i] << " " << val_elst_dist[ij]<< std::endl;
                ++ij;
            }
        }
        return;
    }
    
    void Fragment_i_aa__PardoLoopManager::form_rcut_dist() {
        int rcut_dist_array_slot = sip_tables_.array_slot(std::string("rcut_dist"));
        Block::BlockPtr bptr_rcut_dist = data_manager_.contiguous_array_manager().get_array(rcut_dist_array_slot);
        double *val_rcut_dist = bptr_rcut_dist->get_data();
        int rcut_dist_size = bptr_rcut_dist->size();
        int num_frag = upper_bound_[0] - lower_seg_[0];
        rcut_dist.resize(num_frag);
        for (int i = 0; i < num_frag; ++i) {
            //std::cout << i << std::endl;
            rcut_dist[i].resize(num_frag);
        }
        int ij = 0;
        for (int i = 0; i < num_frag; ++i) {
            for (int j = 0; j < num_frag; ++j) {
                rcut_dist[j][i] = (int)val_rcut_dist[ij];
                //std::cout << i << " " << j << " " << ij << " " << rcut_dist[j][i] << " " << val_rcut_dist[ij]<< std::endl;
                ++ij;
            }
        }
        return;
    }
    
    void Fragment_i_aa__PardoLoopManager::form_swao_frag() {
        int swao_frag_array_slot = sip_tables_.array_slot(std::string("swao_frag"));
        Block::BlockPtr bptr_swao_frag = data_manager_.contiguous_array_manager().get_array(swao_frag_array_slot);
        double *val_swao_frag = bptr_swao_frag->get_data();
        int swao_frag_size = bptr_swao_frag->size();
        swao_frag.resize(swao_frag_size);
        for (int i = 0; i < swao_frag_size; ++i) {
            swao_frag[i] = (int)val_swao_frag[i];
            //std::cout << "vec " << swao_frag[i] << " data " << val_swao_frag[i] << std::endl;
        }
        return;
    }
    
    void Fragment_i_aa__PardoLoopManager::form_swocca_frag() {
        int swocca_frag_array_slot = sip_tables_.array_slot(std::string("swocca_frag"));
        Block::BlockPtr bptr_swocca_frag = data_manager_.contiguous_array_manager().get_array(swocca_frag_array_slot);
        double *val_swocca_frag = bptr_swocca_frag->get_data();
        int swocca_frag_size = bptr_swocca_frag->size();
        swocca_frag.resize(swocca_frag_size);
        for (int i = 0; i < swocca_frag_size; ++i) {
            swocca_frag[i] = (int)val_swocca_frag[i];
            //std::cout << "vec " << swocca_frag[i] << " data " << val_swocca_frag[i] << std::endl;
        }
        return;
    }
    
    void Fragment_i_aa__PardoLoopManager::form_swvirta_frag() {
        int swvirta_frag_array_slot = sip_tables_.array_slot(std::string("swvirta_frag"));
        Block::BlockPtr bptr_swvirta_frag = data_manager_.contiguous_array_manager().get_array(swvirta_frag_array_slot);
        double *val_swvirta_frag = bptr_swvirta_frag->get_data();
        int swvirta_frag_size = bptr_swvirta_frag->size();
        swvirta_frag.resize(swvirta_frag_size);
        for (int i = 0; i < swvirta_frag_size; ++i) {
            swvirta_frag[i] = (int)val_swvirta_frag[i];
            //std::cout << "vec " << swvirta_frag[i] << " data " << val_swvirta_frag[i] << std::endl;
        }
        return;
    }
    
    inline bool Fragment_i_aa__PardoLoopManager::initialize_indices() {
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
            index_values_[i] = lower_seg_[i];
            data_manager_.set_index_value(index_id_[i], lower_seg_[i]);
        }
        more_iterations = true;
        return more_iterations;
    }
    
    /*
     the generic where clauses are:
     typ = 1: map AO to fragment
     typ = 2: map Occ to fragment
     typ = 3: map Virt to fragment
     typ = 4: map elst_dist(ifrag,jfrag) == ifrag
     typ = 5: map rcut_dist(ifrag,jfrag) == ifrag
     typ = 0: ifrag != jfrag
     */
    bool Fragment_i_aa__PardoLoopManager::fragment_special_where_clause(int typ, int index, int frag) {
        bool where_clause;
        int ij = 0;
        switch (typ) {
                
            case 1: // check against ao index
                //where_clause = return_val_swao_frag(index_values_[index])     == index_values_[frag];
                where_clause = swao_frag[index_values_[index] - lower_seg_[index]]     == index_values_[frag];
                break;
                
            case 2: // check against occ index
                //where_clause = return_val_swocca_frag(index_values_[index])   == index_values_[frag];
                where_clause = swocca_frag[index_values_[index] - lower_seg_[index]]   == index_values_[frag];
                break;
                
            case 3: // check against virt index
                //where_clause = return_val_swvirta_frag(index_values_[index])  == index_values_[frag];
                where_clause = swvirta_frag[index_values_[index] - lower_seg_[index]]  == index_values_[frag];
                break;
                
            case 4: // check elst_dist[ifrag,jfrag] == ifrag
                //where_clause = return_val_elst_dist(index_values_[frag], index_values_[index]) == index_values_[frag];
                //ij = index_values_[frag] - 1 + (upper_bound_[0])*(index_values_[index] - 1);
                //std::cout << index_values_[index] << " " << index_values_[frag] << " " << ij << " " << elst_dist[ij] << std::endl;
                //where_clause = elst_dist[ij] == index_values_[frag];
                where_clause = elst_dist[index_values_[frag] - lower_seg_[frag]][index_values_[index] - lower_seg_[index]] == index_values_[frag];
                break;
                
            case 5: // check rcut_dist[ifrag,jfrag] == ifrag
                //where_clause = return_val_rcut_dist(index_values_[frag], index_values_[index]) == index_values_[frag];
                //ij = index_values_[index] - 1 + upper_bound_[0]*(index_values_[frag] - 1);
                //where_clause = rcut_dist[ij] == index_values_[frag];
                //break;
                where_clause = rcut_dist[index_values_[frag] - lower_seg_[frag]][index_values_[index] - lower_seg_[index]] == index_values_[frag];
                break;
                
            case 0:
                where_clause = index_values_[frag] != index_values_[index];
                break;
            default:
                where_clause = false;
        }
        //std::cout << "where_clause " << typ << " " << index_values_[index] << " " << index_values_[frag] << " " << where_clause << std::endl;
        return where_clause;
    }
    
    
    void Fragment_i_aa__PardoLoopManager::do_finalize() {
        for (int i = 0; i < num_indices_; ++i) {
            data_manager_.set_index_undefined(index_id_[i]);
        }
    }
    
    std::string Fragment_i_aa__PardoLoopManager::to_string() const {
        std::stringstream ss;
        ss << "Fragment_i_aa__PardoLoopManager:  num_indices="
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
                             const Fragment_i_aa__PardoLoopManager &obj) {
        os << obj.to_string();
        return os;
    }

    /*!
     -------------------------------------------
     _i_aavo__
     -------------------------------------------
     */
    
    
    /*
     for each new special fragment where clause pattern, this should be the only thing realy changed.
     see comment above for fragment_special_where_clause syntax
     */
    bool Fragment_i_aavo__PardoLoopManager::where_clause(int index) {
        bool where_;
        int ifrag = 0;
        int jfrag = 1;
        int ao = 1;
        int occ = 2;
        int virt = 3;
        int elst = 4;
        int rcut = 5;
        
        
        switch (index) {
            case 0:
                where_ = true;
            case 1:
            case 2:
            case 3:
            case 4:
            case 5:
                where_ = fragment_special_where_clause(ao,index,ifrag);
                break;
            default:
                where_ = false;
        }
        return where_;
    }
    
    bool Fragment_i_aavo__PardoLoopManager::do_update() {
        if (to_exit_)
            return false;
        bool more_iters;
        bool where_clauses_value;
        
        interpreter_->skip_where_clauses(num_where_clauses_);
        
        if (first_time_) {
            first_time_ = false;
            more_iters = initialize_indices();
            
            where_clauses_value = true;
            for (int i = 1; i < num_indices_; ++i) {
                where_clauses_value = where_clauses_value && where_clause(i);
            }
            if(where_clauses_value){
                iteration_++;
                if ((iteration_-1) % num_workers_ == company_rank_){
                    return true;
                }
            }
        }
        
        int index_restart[MAX_RANK];
        
        for (int i = 0; i < num_indices_; ++i) {
            index_restart[i] = index_values_[i];
        }
        
        int loop_count = iteration_;
        for (int index_i = index_restart[0]; index_i < upper_bound_[0]; ++index_i) {
            index_values_[0] = index_i;
            data_manager_.set_index_value(index_id_[0], index_i);
                
                for (int index_1 = index_restart[1]; index_1 < upper_bound_[1]; ++index_1) {
                    index_values_[1] = index_1;
                    data_manager_.set_index_value(index_id_[1], index_1);
                        
            if (where_clause(1)) {
                for (int index_2 = index_restart[2]; index_2 < upper_bound_[2]; ++index_2) {
                    index_values_[2] = index_2;
                    data_manager_.set_index_value(index_id_[2], index_2);
                    
            if (where_clause(2)) {
                for (int index_3 = index_restart[3]; index_3 < upper_bound_[3]; ++index_3) {
                    index_values_[3] = index_3;
                    data_manager_.set_index_value(index_id_[3], index_3);
                    
            if (where_clause(3)) {
                for (int index_4 = index_restart[4]; index_4 < upper_bound_[4]; ++index_4) {
                    index_values_[4] = index_4;
                    data_manager_.set_index_value(index_id_[4], index_4);
                    
            if (where_clause(4)) {
            if (loop_count > iteration_) {
                iteration_++;
                if ((iteration_-1) % num_workers_ == company_rank_){
                    return true;
                }
            }
            ++loop_count;
        } // where 4
            } // index_4
            for (int i = 4; i < num_indices_; ++i) { index_restart[i] = lower_seg_[i]; }
        } // where 3
            } // index_3
            for (int i = 3; i < num_indices_; ++i) { index_restart[i] = lower_seg_[i]; }
        } // where 2
            } // index_2
            for (int i = 2; i < num_indices_; ++i) { index_restart[i] = lower_seg_[i]; }
        } // where 1
            } // index_1
            for (int i = 1; i < num_indices_; ++i) { index_restart[i] = lower_seg_[i]; }
        } // index_i
        
        return false; //this should be false here
    }
    
    Fragment_i_aavo__PardoLoopManager::Fragment_i_aavo__PardoLoopManager(
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
        
        //form_elst_dist();
        //form_rcut_dist();
        form_swao_frag();
        //form_swocca_frag();
        //form_swvirta_frag();
    }
    
    Fragment_i_aavo__PardoLoopManager::~Fragment_i_aavo__PardoLoopManager() {}
    
    bool Fragment_i_aavo__PardoLoopManager::increment_simple_pair(){
        bool more = false;      // More iterations?
        int current_value;
        for (int i = 0; i < 2; ++i) {
            more = increment_single_index(i);
            if (more) {break;}
        } //if here, then all indices are at their max value
        return more;
    }
    
    bool Fragment_i_aavo__PardoLoopManager::increment_single_index(int index){
        bool more = false; 	// More iterations?
        int current_value;
        current_value = data_manager_.index_value(index_id_[index]);
        ++current_value;
        if (current_value < upper_bound_[index]) {
            //increment current index and return
            index_values_[index] = current_value;
            data_manager_.set_index_value(index_id_[index], current_value);
            more = true;
        } else {
            //wrap around and handle next index
            index_values_[index] = lower_seg_[index];
            data_manager_.set_index_value(index_id_[index], lower_seg_[index]);
        }
        return more;
    }
    
    
    inline bool Fragment_i_aavo__PardoLoopManager::increment_all() {
        bool more = false;      // More iterations?
        int current_value;
        for (int i = 0; i < num_indices_; ++i) {
            more = increment_single_index(i);
            if (more) {break;}
        } //if here, then all indices are at their max value
        return more;
    }
    
    void Fragment_i_aavo__PardoLoopManager::form_elst_dist() {
        int elst_dist_array_slot = sip_tables_.array_slot(std::string("elst_dist"));
        Block::BlockPtr bptr_elst_dist = data_manager_.contiguous_array_manager().get_array(elst_dist_array_slot);
        double *val_elst_dist = bptr_elst_dist->get_data();
        int elst_dist_size = bptr_elst_dist->size();
        int num_frag = upper_bound_[0] - lower_seg_[0];
        elst_dist.resize(num_frag);
        for (int i = 0; i < num_frag; ++i) {
            //std::cout << i << std::endl;
            elst_dist[i].resize(num_frag);
        }
        int ij = 0;
        for (int i = 0; i < num_frag; ++i) {
            for (int j = 0; j < num_frag; ++j) {
                elst_dist[j][i] = (int)val_elst_dist[ij];
                //std::cout << i << " " << j << " " << ij << " " << elst_dist[j][i] << " " << val_elst_dist[ij]<< std::endl;
                ++ij;
            }
        }
        return;
    }
    
    void Fragment_i_aavo__PardoLoopManager::form_rcut_dist() {
        int rcut_dist_array_slot = sip_tables_.array_slot(std::string("rcut_dist"));
        Block::BlockPtr bptr_rcut_dist = data_manager_.contiguous_array_manager().get_array(rcut_dist_array_slot);
        double *val_rcut_dist = bptr_rcut_dist->get_data();
        int rcut_dist_size = bptr_rcut_dist->size();
        int num_frag = upper_bound_[0] - lower_seg_[0];
        rcut_dist.resize(num_frag);
        for (int i = 0; i < num_frag; ++i) {
            //std::cout << i << std::endl;
            rcut_dist[i].resize(num_frag);
        }
        int ij = 0;
        for (int i = 0; i < num_frag; ++i) {
            for (int j = 0; j < num_frag; ++j) {
                rcut_dist[j][i] = (int)val_rcut_dist[ij];
                //std::cout << i << " " << j << " " << ij << " " << rcut_dist[j][i] << " " << val_rcut_dist[ij]<< std::endl;
                ++ij;
            }
        }
        return;
    }
    
    void Fragment_i_aavo__PardoLoopManager::form_swao_frag() {
        int swao_frag_array_slot = sip_tables_.array_slot(std::string("swao_frag"));
        Block::BlockPtr bptr_swao_frag = data_manager_.contiguous_array_manager().get_array(swao_frag_array_slot);
        double *val_swao_frag = bptr_swao_frag->get_data();
        int swao_frag_size = bptr_swao_frag->size();
        swao_frag.resize(swao_frag_size);
        for (int i = 0; i < swao_frag_size; ++i) {
            swao_frag[i] = (int)val_swao_frag[i];
            //std::cout << "vec " << swao_frag[i] << " data " << val_swao_frag[i] << std::endl;
        }
        return;
    }
    
    void Fragment_i_aavo__PardoLoopManager::form_swocca_frag() {
        int swocca_frag_array_slot = sip_tables_.array_slot(std::string("swocca_frag"));
        Block::BlockPtr bptr_swocca_frag = data_manager_.contiguous_array_manager().get_array(swocca_frag_array_slot);
        double *val_swocca_frag = bptr_swocca_frag->get_data();
        int swocca_frag_size = bptr_swocca_frag->size();
        swocca_frag.resize(swocca_frag_size);
        for (int i = 0; i < swocca_frag_size; ++i) {
            swocca_frag[i] = (int)val_swocca_frag[i];
            //std::cout << "vec " << swocca_frag[i] << " data " << val_swocca_frag[i] << std::endl;
        }
        return;
    }
    
    void Fragment_i_aavo__PardoLoopManager::form_swvirta_frag() {
        int swvirta_frag_array_slot = sip_tables_.array_slot(std::string("swvirta_frag"));
        Block::BlockPtr bptr_swvirta_frag = data_manager_.contiguous_array_manager().get_array(swvirta_frag_array_slot);
        double *val_swvirta_frag = bptr_swvirta_frag->get_data();
        int swvirta_frag_size = bptr_swvirta_frag->size();
        swvirta_frag.resize(swvirta_frag_size);
        for (int i = 0; i < swvirta_frag_size; ++i) {
            swvirta_frag[i] = (int)val_swvirta_frag[i];
            //std::cout << "vec " << swvirta_frag[i] << " data " << val_swvirta_frag[i] << std::endl;
        }
        return;
    }
    
    inline bool Fragment_i_aavo__PardoLoopManager::initialize_indices() {
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
            index_values_[i] = lower_seg_[i];
            data_manager_.set_index_value(index_id_[i], lower_seg_[i]);
        }
        more_iterations = true;
        return more_iterations;
    }
    
    /*
     the generic where clauses are:
     typ = 1: map AO to fragment
     typ = 2: map Occ to fragment
     typ = 3: map Virt to fragment
     typ = 4: map elst_dist(ifrag,jfrag) == ifrag
     typ = 5: map rcut_dist(ifrag,jfrag) == ifrag
     typ = 0: ifrag != jfrag
     */
    bool Fragment_i_aavo__PardoLoopManager::fragment_special_where_clause(int typ, int index, int frag) {
        bool where_clause;
        int ij = 0;
        switch (typ) {
                
            case 1: // check against ao index
                //where_clause = return_val_swao_frag(index_values_[index])     == index_values_[frag];
                where_clause = swao_frag[index_values_[index] - lower_seg_[index]]     == index_values_[frag];
                break;
                
            case 2: // check against occ index
                //where_clause = return_val_swocca_frag(index_values_[index])   == index_values_[frag];
                where_clause = swocca_frag[index_values_[index] - lower_seg_[index]]   == index_values_[frag];
                break;
                
            case 3: // check against virt index
                //where_clause = return_val_swvirta_frag(index_values_[index])  == index_values_[frag];
                where_clause = swvirta_frag[index_values_[index] - lower_seg_[index]]  == index_values_[frag];
                break;
                
            case 4: // check elst_dist[ifrag,jfrag] == ifrag
                //where_clause = return_val_elst_dist(index_values_[frag], index_values_[index]) == index_values_[frag];
                //ij = index_values_[frag] - 1 + (upper_bound_[0])*(index_values_[index] - 1);
                //std::cout << index_values_[index] << " " << index_values_[frag] << " " << ij << " " << elst_dist[ij] << std::endl;
                //where_clause = elst_dist[ij] == index_values_[frag];
                where_clause = elst_dist[index_values_[frag] - lower_seg_[frag]][index_values_[index] - lower_seg_[index]] == index_values_[frag];
                break;
                
            case 5: // check rcut_dist[ifrag,jfrag] == ifrag
                //where_clause = return_val_rcut_dist(index_values_[frag], index_values_[index]) == index_values_[frag];
                //ij = index_values_[index] - 1 + upper_bound_[0]*(index_values_[frag] - 1);
                //where_clause = rcut_dist[ij] == index_values_[frag];
                //break;
                where_clause = rcut_dist[index_values_[frag] - lower_seg_[frag]][index_values_[index] - lower_seg_[index]] == index_values_[frag];
                break;
                
            case 0:
                where_clause = index_values_[frag] != index_values_[index];
                break;
            default:
                where_clause = false;
        }
        //std::cout << "where_clause " << typ << " " << index_values_[index] << " " << index_values_[frag] << " " << where_clause << std::endl;
        return where_clause;
    }
    
    
    void Fragment_i_aavo__PardoLoopManager::do_finalize() {
        for (int i = 0; i < num_indices_; ++i) {
            data_manager_.set_index_undefined(index_id_[i]);
        }
    }
    
    std::string Fragment_i_aavo__PardoLoopManager::to_string() const {
        std::stringstream ss;
        ss << "Fragment_i_aavo__PardoLoopManager:  num_indices="
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
                             const Fragment_i_aavo__PardoLoopManager &obj) {
        os << obj.to_string();
        return os;
    }

    /*!
     -------------------------------------------
     _i_aavv__
     -------------------------------------------
     */
    
    
    /*
     for each new special fragment where clause pattern, this should be the only thing realy changed.
     see comment above for fragment_special_where_clause syntax
     */
    bool Fragment_i_aavv__PardoLoopManager::where_clause(int index) {
        bool where_;
        int ifrag = 0;
        int jfrag = 1;
        int ao = 1;
        int occ = 2;
        int virt = 3;
        int elst = 4;
        int rcut = 5;
        
        
        switch (index) {
            case 0:
                where_ = true;
            case 1:
            case 2:
            case 3:
            case 4:
            case 5:
                where_ = fragment_special_where_clause(ao,index,ifrag);
                break;
            default:
                where_ = false;
        }
        return where_;
    }
    
    bool Fragment_i_aavv__PardoLoopManager::do_update() {
        if (to_exit_)
            return false;
        bool more_iters;
        bool where_clauses_value;
        
        interpreter_->skip_where_clauses(num_where_clauses_);
        
        if (first_time_) {
            first_time_ = false;
            more_iters = initialize_indices();
            
            where_clauses_value = true;
            for (int i = 1; i < num_indices_; ++i) {
                where_clauses_value = where_clauses_value && where_clause(i);
            }
            if(where_clauses_value){
                iteration_++;
                if ((iteration_-1) % num_workers_ == company_rank_){
                    return true;
                }
            }
        }
        
        int index_restart[MAX_RANK];
        
        for (int i = 0; i < num_indices_; ++i) {
            index_restart[i] = index_values_[i];
        }
        
        int loop_count = iteration_;
        for (int index_i = index_restart[0]; index_i < upper_bound_[0]; ++index_i) {
            index_values_[0] = index_i;
            data_manager_.set_index_value(index_id_[0], index_i);
                
                for (int index_1 = index_restart[1]; index_1 < upper_bound_[1]; ++index_1) {
                    index_values_[1] = index_1;
                    data_manager_.set_index_value(index_id_[1], index_1);
                        
            if (where_clause(1)) {
                for (int index_2 = index_restart[2]; index_2 < upper_bound_[2]; ++index_2) {
                    index_values_[2] = index_2;
                    data_manager_.set_index_value(index_id_[2], index_2);
                    
            if (where_clause(2)) {
                for (int index_3 = index_restart[3]; index_3 < upper_bound_[3]; ++index_3) {
                    index_values_[3] = index_3;
                    data_manager_.set_index_value(index_id_[3], index_3);
                    
            if (where_clause(3)) {
                for (int index_4 = index_restart[4]; index_4 < upper_bound_[4]; ++index_4) {
                    index_values_[4] = index_4;
                    data_manager_.set_index_value(index_id_[4], index_4);
                    
            if (where_clause(4)) {
            if (loop_count > iteration_) {
                iteration_++;
                if ((iteration_-1) % num_workers_ == company_rank_){
                    return true;
                }
            }
            ++loop_count;
        } // where 4
            } // index_4
            for (int i = 4; i < num_indices_; ++i) { index_restart[i] = lower_seg_[i]; }
        } // where 3
            } // index_3
            for (int i = 3; i < num_indices_; ++i) { index_restart[i] = lower_seg_[i]; }
        } // where 2
            } // index_2
            for (int i = 2; i < num_indices_; ++i) { index_restart[i] = lower_seg_[i]; }
        } // where 1
            } // index_1
            for (int i = 1; i < num_indices_; ++i) { index_restart[i] = lower_seg_[i]; }
        } // index_i
        
        return false; //this should be false here
    }
    
    Fragment_i_aavv__PardoLoopManager::Fragment_i_aavv__PardoLoopManager(
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
        
        //form_elst_dist();
        //form_rcut_dist();
        form_swao_frag();
        //form_swocca_frag();
        //form_swvirta_frag();
    }
    
    Fragment_i_aavv__PardoLoopManager::~Fragment_i_aavv__PardoLoopManager() {}
    
    bool Fragment_i_aavv__PardoLoopManager::increment_simple_pair(){
        bool more = false;      // More iterations?
        int current_value;
        for (int i = 0; i < 2; ++i) {
            more = increment_single_index(i);
            if (more) {break;}
        } //if here, then all indices are at their max value
        return more;
    }
    
    bool Fragment_i_aavv__PardoLoopManager::increment_single_index(int index){
        bool more = false; 	// More iterations?
        int current_value;
        current_value = data_manager_.index_value(index_id_[index]);
        ++current_value;
        if (current_value < upper_bound_[index]) {
            //increment current index and return
            index_values_[index] = current_value;
            data_manager_.set_index_value(index_id_[index], current_value);
            more = true;
        } else {
            //wrap around and handle next index
            index_values_[index] = lower_seg_[index];
            data_manager_.set_index_value(index_id_[index], lower_seg_[index]);
        }
        return more;
    }
    
    
    inline bool Fragment_i_aavv__PardoLoopManager::increment_all() {
        bool more = false;      // More iterations?
        int current_value;
        for (int i = 0; i < num_indices_; ++i) {
            more = increment_single_index(i);
            if (more) {break;}
        } //if here, then all indices are at their max value
        return more;
    }
    
    void Fragment_i_aavv__PardoLoopManager::form_elst_dist() {
        int elst_dist_array_slot = sip_tables_.array_slot(std::string("elst_dist"));
        Block::BlockPtr bptr_elst_dist = data_manager_.contiguous_array_manager().get_array(elst_dist_array_slot);
        double *val_elst_dist = bptr_elst_dist->get_data();
        int elst_dist_size = bptr_elst_dist->size();
        int num_frag = upper_bound_[0] - lower_seg_[0];
        elst_dist.resize(num_frag);
        for (int i = 0; i < num_frag; ++i) {
            //std::cout << i << std::endl;
            elst_dist[i].resize(num_frag);
        }
        int ij = 0;
        for (int i = 0; i < num_frag; ++i) {
            for (int j = 0; j < num_frag; ++j) {
                elst_dist[j][i] = (int)val_elst_dist[ij];
                //std::cout << i << " " << j << " " << ij << " " << elst_dist[j][i] << " " << val_elst_dist[ij]<< std::endl;
                ++ij;
            }
        }
        return;
    }
    
    void Fragment_i_aavv__PardoLoopManager::form_rcut_dist() {
        int rcut_dist_array_slot = sip_tables_.array_slot(std::string("rcut_dist"));
        Block::BlockPtr bptr_rcut_dist = data_manager_.contiguous_array_manager().get_array(rcut_dist_array_slot);
        double *val_rcut_dist = bptr_rcut_dist->get_data();
        int rcut_dist_size = bptr_rcut_dist->size();
        int num_frag = upper_bound_[0] - lower_seg_[0];
        rcut_dist.resize(num_frag);
        for (int i = 0; i < num_frag; ++i) {
            //std::cout << i << std::endl;
            rcut_dist[i].resize(num_frag);
        }
        int ij = 0;
        for (int i = 0; i < num_frag; ++i) {
            for (int j = 0; j < num_frag; ++j) {
                rcut_dist[j][i] = (int)val_rcut_dist[ij];
                //std::cout << i << " " << j << " " << ij << " " << rcut_dist[j][i] << " " << val_rcut_dist[ij]<< std::endl;
                ++ij;
            }
        }
        return;
    }
    
    void Fragment_i_aavv__PardoLoopManager::form_swao_frag() {
        int swao_frag_array_slot = sip_tables_.array_slot(std::string("swao_frag"));
        Block::BlockPtr bptr_swao_frag = data_manager_.contiguous_array_manager().get_array(swao_frag_array_slot);
        double *val_swao_frag = bptr_swao_frag->get_data();
        int swao_frag_size = bptr_swao_frag->size();
        swao_frag.resize(swao_frag_size);
        for (int i = 0; i < swao_frag_size; ++i) {
            swao_frag[i] = (int)val_swao_frag[i];
            //std::cout << "vec " << swao_frag[i] << " data " << val_swao_frag[i] << std::endl;
        }
        return;
    }
    
    void Fragment_i_aavv__PardoLoopManager::form_swocca_frag() {
        int swocca_frag_array_slot = sip_tables_.array_slot(std::string("swocca_frag"));
        Block::BlockPtr bptr_swocca_frag = data_manager_.contiguous_array_manager().get_array(swocca_frag_array_slot);
        double *val_swocca_frag = bptr_swocca_frag->get_data();
        int swocca_frag_size = bptr_swocca_frag->size();
        swocca_frag.resize(swocca_frag_size);
        for (int i = 0; i < swocca_frag_size; ++i) {
            swocca_frag[i] = (int)val_swocca_frag[i];
            //std::cout << "vec " << swocca_frag[i] << " data " << val_swocca_frag[i] << std::endl;
        }
        return;
    }
    
    void Fragment_i_aavv__PardoLoopManager::form_swvirta_frag() {
        int swvirta_frag_array_slot = sip_tables_.array_slot(std::string("swvirta_frag"));
        Block::BlockPtr bptr_swvirta_frag = data_manager_.contiguous_array_manager().get_array(swvirta_frag_array_slot);
        double *val_swvirta_frag = bptr_swvirta_frag->get_data();
        int swvirta_frag_size = bptr_swvirta_frag->size();
        swvirta_frag.resize(swvirta_frag_size);
        for (int i = 0; i < swvirta_frag_size; ++i) {
            swvirta_frag[i] = (int)val_swvirta_frag[i];
            //std::cout << "vec " << swvirta_frag[i] << " data " << val_swvirta_frag[i] << std::endl;
        }
        return;
    }
    
    inline bool Fragment_i_aavv__PardoLoopManager::initialize_indices() {
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
            index_values_[i] = lower_seg_[i];
            data_manager_.set_index_value(index_id_[i], lower_seg_[i]);
        }
        more_iterations = true;
        return more_iterations;
    }
    
    /*
     the generic where clauses are:
     typ = 1: map AO to fragment
     typ = 2: map Occ to fragment
     typ = 3: map Virt to fragment
     typ = 4: map elst_dist(ifrag,jfrag) == ifrag
     typ = 5: map rcut_dist(ifrag,jfrag) == ifrag
     typ = 0: ifrag != jfrag
     */
    bool Fragment_i_aavv__PardoLoopManager::fragment_special_where_clause(int typ, int index, int frag) {
        bool where_clause;
        int ij = 0;
        switch (typ) {
                
            case 1: // check against ao index
                //where_clause = return_val_swao_frag(index_values_[index])     == index_values_[frag];
                where_clause = swao_frag[index_values_[index] - lower_seg_[index]]     == index_values_[frag];
                break;
                
            case 2: // check against occ index
                //where_clause = return_val_swocca_frag(index_values_[index])   == index_values_[frag];
                where_clause = swocca_frag[index_values_[index] - lower_seg_[index]]   == index_values_[frag];
                break;
                
            case 3: // check against virt index
                //where_clause = return_val_swvirta_frag(index_values_[index])  == index_values_[frag];
                where_clause = swvirta_frag[index_values_[index] - lower_seg_[index]]  == index_values_[frag];
                break;
                
            case 4: // check elst_dist[ifrag,jfrag] == ifrag
                //where_clause = return_val_elst_dist(index_values_[frag], index_values_[index]) == index_values_[frag];
                //ij = index_values_[frag] - 1 + (upper_bound_[0])*(index_values_[index] - 1);
                //std::cout << index_values_[index] << " " << index_values_[frag] << " " << ij << " " << elst_dist[ij] << std::endl;
                //where_clause = elst_dist[ij] == index_values_[frag];
                where_clause = elst_dist[index_values_[frag] - lower_seg_[frag]][index_values_[index] - lower_seg_[index]] == index_values_[frag];
                break;
                
            case 5: // check rcut_dist[ifrag,jfrag] == ifrag
                //where_clause = return_val_rcut_dist(index_values_[frag], index_values_[index]) == index_values_[frag];
                //ij = index_values_[index] - 1 + upper_bound_[0]*(index_values_[frag] - 1);
                //where_clause = rcut_dist[ij] == index_values_[frag];
                //break;
                where_clause = rcut_dist[index_values_[frag] - lower_seg_[frag]][index_values_[index] - lower_seg_[index]] == index_values_[frag];
                break;
                
            case 0:
                where_clause = index_values_[frag] != index_values_[index];
                break;
            default:
                where_clause = false;
        }
        //std::cout << "where_clause " << typ << " " << index_values_[index] << " " << index_values_[frag] << " " << where_clause << std::endl;
        return where_clause;
    }
    
    
    void Fragment_i_aavv__PardoLoopManager::do_finalize() {
        for (int i = 0; i < num_indices_; ++i) {
            data_manager_.set_index_undefined(index_id_[i]);
        }
    }
    
    std::string Fragment_i_aavv__PardoLoopManager::to_string() const {
        std::stringstream ss;
        ss << "Fragment_i_aavv__PardoLoopManager:  num_indices="
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
                             const Fragment_i_aavv__PardoLoopManager &obj) {
        os << obj.to_string();
        return os;
    }

    /*!
     -------------------------------------------
     _i_aaoo__
     -------------------------------------------
     */
    
    
    /*
     for each new special fragment where clause pattern, this should be the only thing realy changed.
     see comment above for fragment_special_where_clause syntax
     */
    bool Fragment_i_aaoo__PardoLoopManager::where_clause(int index) {
        bool where_;
        int ifrag = 0;
        int jfrag = 1;
        int ao = 1;
        int occ = 2;
        int virt = 3;
        int elst = 4;
        int rcut = 5;
        
        
        switch (index) {
            case 0:
                where_ = true;
            case 1:
            case 2:
            case 3:
            case 4:
            case 5:
                where_ = fragment_special_where_clause(ao,index,ifrag);
                break;
            default:
                where_ = false;
        }
        return where_;
    }
    
    bool Fragment_i_aaoo__PardoLoopManager::do_update() {
        if (to_exit_)
            return false;
        bool more_iters;
        bool where_clauses_value;
        
        interpreter_->skip_where_clauses(num_where_clauses_);
        
        if (first_time_) {
            first_time_ = false;
            more_iters = initialize_indices();
            
            where_clauses_value = true;
            for (int i = 1; i < num_indices_; ++i) {
                where_clauses_value = where_clauses_value && where_clause(i);
            }
            if(where_clauses_value){
                iteration_++;
                if ((iteration_-1) % num_workers_ == company_rank_){
                    return true;
                }
            }
        }
        
        int index_restart[MAX_RANK];
        
        for (int i = 0; i < num_indices_; ++i) {
            index_restart[i] = index_values_[i];
        }
        
        int loop_count = iteration_;
        for (int index_i = index_restart[0]; index_i < upper_bound_[0]; ++index_i) {
            index_values_[0] = index_i;
            data_manager_.set_index_value(index_id_[0], index_i);
                
                for (int index_1 = index_restart[1]; index_1 < upper_bound_[1]; ++index_1) {
                    index_values_[1] = index_1;
                    data_manager_.set_index_value(index_id_[1], index_1);
                        
            if (where_clause(1)) {
                for (int index_2 = index_restart[2]; index_2 < upper_bound_[2]; ++index_2) {
                    index_values_[2] = index_2;
                    data_manager_.set_index_value(index_id_[2], index_2);
                    
            if (where_clause(2)) {
                for (int index_3 = index_restart[3]; index_3 < upper_bound_[3]; ++index_3) {
                    index_values_[3] = index_3;
                    data_manager_.set_index_value(index_id_[3], index_3);
                    
            if (where_clause(3)) {
                for (int index_4 = index_restart[4]; index_4 < upper_bound_[4]; ++index_4) {
                    index_values_[4] = index_4;
                    data_manager_.set_index_value(index_id_[4], index_4);
                    
            if (where_clause(4)) {
            if (loop_count > iteration_) {
                iteration_++;
                if ((iteration_-1) % num_workers_ == company_rank_){
                    return true;
                }
            }
            ++loop_count;
        } // where 4
            } // index_4
            for (int i = 4; i < num_indices_; ++i) { index_restart[i] = lower_seg_[i]; }
        } // where 3
            } // index_3
            for (int i = 3; i < num_indices_; ++i) { index_restart[i] = lower_seg_[i]; }
        } // where 2
            } // index_2
            for (int i = 2; i < num_indices_; ++i) { index_restart[i] = lower_seg_[i]; }
        } // where 1
            } // index_1
            for (int i = 1; i < num_indices_; ++i) { index_restart[i] = lower_seg_[i]; }
        } // index_i
        
        return false; //this should be false here
    }
    
    Fragment_i_aaoo__PardoLoopManager::Fragment_i_aaoo__PardoLoopManager(
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
        
        //form_elst_dist();
        //form_rcut_dist();
        form_swao_frag();
        //form_swocca_frag();
        //form_swvirta_frag();
    }
    
    Fragment_i_aaoo__PardoLoopManager::~Fragment_i_aaoo__PardoLoopManager() {}
    
    bool Fragment_i_aaoo__PardoLoopManager::increment_simple_pair(){
        bool more = false;      // More iterations?
        int current_value;
        for (int i = 0; i < 2; ++i) {
            more = increment_single_index(i);
            if (more) {break;}
        } //if here, then all indices are at their max value
        return more;
    }
    
    bool Fragment_i_aaoo__PardoLoopManager::increment_single_index(int index){
        bool more = false; 	// More iterations?
        int current_value;
        current_value = data_manager_.index_value(index_id_[index]);
        ++current_value;
        if (current_value < upper_bound_[index]) {
            //increment current index and return
            index_values_[index] = current_value;
            data_manager_.set_index_value(index_id_[index], current_value);
            more = true;
        } else {
            //wrap around and handle next index
            index_values_[index] = lower_seg_[index];
            data_manager_.set_index_value(index_id_[index], lower_seg_[index]);
        }
        return more;
    }
    
    
    inline bool Fragment_i_aaoo__PardoLoopManager::increment_all() {
        bool more = false;      // More iterations?
        int current_value;
        for (int i = 0; i < num_indices_; ++i) {
            more = increment_single_index(i);
            if (more) {break;}
        } //if here, then all indices are at their max value
        return more;
    }
    
    void Fragment_i_aaoo__PardoLoopManager::form_elst_dist() {
        int elst_dist_array_slot = sip_tables_.array_slot(std::string("elst_dist"));
        Block::BlockPtr bptr_elst_dist = data_manager_.contiguous_array_manager().get_array(elst_dist_array_slot);
        double *val_elst_dist = bptr_elst_dist->get_data();
        int elst_dist_size = bptr_elst_dist->size();
        int num_frag = upper_bound_[0] - lower_seg_[0];
        elst_dist.resize(num_frag);
        for (int i = 0; i < num_frag; ++i) {
            //std::cout << i << std::endl;
            elst_dist[i].resize(num_frag);
        }
        int ij = 0;
        for (int i = 0; i < num_frag; ++i) {
            for (int j = 0; j < num_frag; ++j) {
                elst_dist[j][i] = (int)val_elst_dist[ij];
                //std::cout << i << " " << j << " " << ij << " " << elst_dist[j][i] << " " << val_elst_dist[ij]<< std::endl;
                ++ij;
            }
        }
        return;
    }
    
    void Fragment_i_aaoo__PardoLoopManager::form_rcut_dist() {
        int rcut_dist_array_slot = sip_tables_.array_slot(std::string("rcut_dist"));
        Block::BlockPtr bptr_rcut_dist = data_manager_.contiguous_array_manager().get_array(rcut_dist_array_slot);
        double *val_rcut_dist = bptr_rcut_dist->get_data();
        int rcut_dist_size = bptr_rcut_dist->size();
        int num_frag = upper_bound_[0] - lower_seg_[0];
        rcut_dist.resize(num_frag);
        for (int i = 0; i < num_frag; ++i) {
            //std::cout << i << std::endl;
            rcut_dist[i].resize(num_frag);
        }
        int ij = 0;
        for (int i = 0; i < num_frag; ++i) {
            for (int j = 0; j < num_frag; ++j) {
                rcut_dist[j][i] = (int)val_rcut_dist[ij];
                //std::cout << i << " " << j << " " << ij << " " << rcut_dist[j][i] << " " << val_rcut_dist[ij]<< std::endl;
                ++ij;
            }
        }
        return;
    }
    
    void Fragment_i_aaoo__PardoLoopManager::form_swao_frag() {
        int swao_frag_array_slot = sip_tables_.array_slot(std::string("swao_frag"));
        Block::BlockPtr bptr_swao_frag = data_manager_.contiguous_array_manager().get_array(swao_frag_array_slot);
        double *val_swao_frag = bptr_swao_frag->get_data();
        int swao_frag_size = bptr_swao_frag->size();
        swao_frag.resize(swao_frag_size);
        for (int i = 0; i < swao_frag_size; ++i) {
            swao_frag[i] = (int)val_swao_frag[i];
            //std::cout << "vec " << swao_frag[i] << " data " << val_swao_frag[i] << std::endl;
        }
        return;
    }
    
    void Fragment_i_aaoo__PardoLoopManager::form_swocca_frag() {
        int swocca_frag_array_slot = sip_tables_.array_slot(std::string("swocca_frag"));
        Block::BlockPtr bptr_swocca_frag = data_manager_.contiguous_array_manager().get_array(swocca_frag_array_slot);
        double *val_swocca_frag = bptr_swocca_frag->get_data();
        int swocca_frag_size = bptr_swocca_frag->size();
        swocca_frag.resize(swocca_frag_size);
        for (int i = 0; i < swocca_frag_size; ++i) {
            swocca_frag[i] = (int)val_swocca_frag[i];
            //std::cout << "vec " << swocca_frag[i] << " data " << val_swocca_frag[i] << std::endl;
        }
        return;
    }
    
    void Fragment_i_aaoo__PardoLoopManager::form_swvirta_frag() {
        int swvirta_frag_array_slot = sip_tables_.array_slot(std::string("swvirta_frag"));
        Block::BlockPtr bptr_swvirta_frag = data_manager_.contiguous_array_manager().get_array(swvirta_frag_array_slot);
        double *val_swvirta_frag = bptr_swvirta_frag->get_data();
        int swvirta_frag_size = bptr_swvirta_frag->size();
        swvirta_frag.resize(swvirta_frag_size);
        for (int i = 0; i < swvirta_frag_size; ++i) {
            swvirta_frag[i] = (int)val_swvirta_frag[i];
            //std::cout << "vec " << swvirta_frag[i] << " data " << val_swvirta_frag[i] << std::endl;
        }
        return;
    }
    
    inline bool Fragment_i_aaoo__PardoLoopManager::initialize_indices() {
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
            index_values_[i] = lower_seg_[i];
            data_manager_.set_index_value(index_id_[i], lower_seg_[i]);
        }
        more_iterations = true;
        return more_iterations;
    }
    
    /*
     the generic where clauses are:
     typ = 1: map AO to fragment
     typ = 2: map Occ to fragment
     typ = 3: map Virt to fragment
     typ = 4: map elst_dist(ifrag,jfrag) == ifrag
     typ = 5: map rcut_dist(ifrag,jfrag) == ifrag
     typ = 0: ifrag != jfrag
     */
    bool Fragment_i_aaoo__PardoLoopManager::fragment_special_where_clause(int typ, int index, int frag) {
        bool where_clause;
        int ij = 0;
        switch (typ) {
                
            case 1: // check against ao index
                //where_clause = return_val_swao_frag(index_values_[index])     == index_values_[frag];
                where_clause = swao_frag[index_values_[index] - lower_seg_[index]]     == index_values_[frag];
                break;
                
            case 2: // check against occ index
                //where_clause = return_val_swocca_frag(index_values_[index])   == index_values_[frag];
                where_clause = swocca_frag[index_values_[index] - lower_seg_[index]]   == index_values_[frag];
                break;
                
            case 3: // check against virt index
                //where_clause = return_val_swvirta_frag(index_values_[index])  == index_values_[frag];
                where_clause = swvirta_frag[index_values_[index] - lower_seg_[index]]  == index_values_[frag];
                break;
                
            case 4: // check elst_dist[ifrag,jfrag] == ifrag
                //where_clause = return_val_elst_dist(index_values_[frag], index_values_[index]) == index_values_[frag];
                //ij = index_values_[frag] - 1 + (upper_bound_[0])*(index_values_[index] - 1);
                //std::cout << index_values_[index] << " " << index_values_[frag] << " " << ij << " " << elst_dist[ij] << std::endl;
                //where_clause = elst_dist[ij] == index_values_[frag];
                where_clause = elst_dist[index_values_[frag] - lower_seg_[frag]][index_values_[index] - lower_seg_[index]] == index_values_[frag];
                break;
                
            case 5: // check rcut_dist[ifrag,jfrag] == ifrag
                //where_clause = return_val_rcut_dist(index_values_[frag], index_values_[index]) == index_values_[frag];
                //ij = index_values_[index] - 1 + upper_bound_[0]*(index_values_[frag] - 1);
                //where_clause = rcut_dist[ij] == index_values_[frag];
                //break;
                where_clause = rcut_dist[index_values_[frag] - lower_seg_[frag]][index_values_[index] - lower_seg_[index]] == index_values_[frag];
                break;
                
            case 0:
                where_clause = index_values_[frag] != index_values_[index];
                break;
            default:
                where_clause = false;
        }
        //std::cout << "where_clause " << typ << " " << index_values_[index] << " " << index_values_[frag] << " " << where_clause << std::endl;
        return where_clause;
    }
    
    
    void Fragment_i_aaoo__PardoLoopManager::do_finalize() {
        for (int i = 0; i < num_indices_; ++i) {
            data_manager_.set_index_undefined(index_id_[i]);
        }
    }
    
    std::string Fragment_i_aaoo__PardoLoopManager::to_string() const {
        std::stringstream ss;
        ss << "Fragment_i_aaoo__PardoLoopManager:  num_indices="
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
                             const Fragment_i_aaoo__PardoLoopManager &obj) {
        os << obj.to_string();
        return os;
    }

/*!
 -------------------------------------------
 _ij_aaa__
 -------------------------------------------

 PARDO ifrag, jfrag, mu, i, b, j #GETLINE: Fragment_ij_aaa__
 where (int)elst_dist[ifrag,jfrag] == ifrag
 where (int)SwAO_frag[(index)mu] == ifrag
 */

    
    /*
     for each new special fragment where clause pattern, this should be the only thing realy changed.
     see comment above for fragment_special_where_clause syntax
     */
    bool Fragment_ij_aaa__PardoLoopManager::where_clause(int index) {
        bool where_;
        int ifrag = 0;
        int jfrag = 1;
        int ao = 1;
        int occ = 2;
        int virt = 3;
        int elst = 4;
        int rcut = 5;

        
        switch (index) {
            case 0:
            case 1:
                where_ = fragment_special_where_clause(elst,jfrag,ifrag);
                break;
            case 2:
            case 3:
            case 4:
            case 5:
                where_ = fragment_special_where_clause(ao,index,ifrag);
                break;
            default:
                where_ = false;
        }
        return where_;
    }
    
    bool Fragment_ij_aaa__PardoLoopManager::do_update() {
        if (to_exit_)
            return false;
        bool more_iters;
        bool where_clauses_value;
        
        interpreter_->skip_where_clauses(num_where_clauses_);
        
        if (first_time_) {
            first_time_ = false;
            more_iters = initialize_indices();
            
            where_clauses_value = true;
            for (int i = 1; i < num_indices_; ++i) {
                where_clauses_value = where_clauses_value && where_clause(i);
            }
            if(where_clauses_value){
                iteration_++;
                if ((iteration_-1) % num_workers_ == company_rank_){
                    return true;
                }
            }
        }
        
        int index_restart[MAX_RANK];
        
        for (int i = 0; i < num_indices_; ++i) {
            index_restart[i] = index_values_[i];
        }
        
        int loop_count = iteration_;
        for (int index_i = index_restart[0]; index_i < upper_bound_[0]; ++index_i) {
            index_values_[0] = index_i;
            data_manager_.set_index_value(index_id_[0], index_i);
            
        for (int index_j = index_restart[1]; index_j < upper_bound_[1]; ++index_j) {
            index_values_[1] = index_j;
            data_manager_.set_index_value(index_id_[1], index_j);
                
            if (where_clause(1)) {
                for (int index_2 = index_restart[2]; index_2 < upper_bound_[2]; ++index_2) {
                    index_values_[2] = index_2;
                    data_manager_.set_index_value(index_id_[2], index_2);
                        
            if (where_clause(2)) {
                for (int index_3 = index_restart[3]; index_3 < upper_bound_[3]; ++index_3) {
                    index_values_[3] = index_3;
                    data_manager_.set_index_value(index_id_[3], index_3);
                    
            if (where_clause(3)) {
                for (int index_4 = index_restart[4]; index_4 < upper_bound_[4]; ++index_4) {
                    index_values_[4] = index_4;
                    data_manager_.set_index_value(index_id_[4], index_4);
                    
            if (where_clause(4)) {
            if (loop_count > iteration_) {
                iteration_++;
                if ((iteration_-1) % num_workers_ == company_rank_){
                    return true;
                }
            }
            ++loop_count;
        } // where 4
            } // index_4
            for (int i = 4; i < num_indices_; ++i) { index_restart[i] = lower_seg_[i]; }
        } // where 3
            } // index_3
            for (int i = 3; i < num_indices_; ++i) { index_restart[i] = lower_seg_[i]; }
        } // where 2
            } // index_2
            for (int i = 2; i < num_indices_; ++i) { index_restart[i] = lower_seg_[i]; }
        } // where 1
            } // index_j
            for (int i = 1; i < num_indices_; ++i) { index_restart[i] = lower_seg_[i]; }
        } // index_i
        
        return false; //this should be false here
    }

Fragment_ij_aaa__PardoLoopManager::Fragment_ij_aaa__PardoLoopManager(
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
    
    form_elst_dist();
    //form_rcut_dist();
    form_swao_frag();
    //form_swocca_frag();
    //form_swvirta_frag();
}

Fragment_ij_aaa__PardoLoopManager::~Fragment_ij_aaa__PardoLoopManager() {}

bool Fragment_ij_aaa__PardoLoopManager::increment_simple_pair(){
    bool more = false;      // More iterations?
    int current_value;
    for (int i = 0; i < 2; ++i) {
        more = increment_single_index(i);
        if (more) {break;}
    } //if here, then all indices are at their max value
    return more;
}

bool Fragment_ij_aaa__PardoLoopManager::increment_single_index(int index){
    bool more = false; 	// More iterations?
    int current_value;
    current_value = data_manager_.index_value(index_id_[index]);
    ++current_value;
    if (current_value < upper_bound_[index]) {
        //increment current index and return
        index_values_[index] = current_value;
        data_manager_.set_index_value(index_id_[index], current_value);
        more = true;
    } else {
        //wrap around and handle next index
        index_values_[index] = lower_seg_[index];
        data_manager_.set_index_value(index_id_[index], lower_seg_[index]);
    }
    return more;
}
    
    
    inline bool Fragment_ij_aaa__PardoLoopManager::increment_all() {
        bool more = false;      // More iterations?
        int current_value;
        for (int i = 0; i < num_indices_; ++i) {
            more = increment_single_index(i);
            if (more) {break;}
        } //if here, then all indices are at their max value
        return more;
    }
    
    void Fragment_ij_aaa__PardoLoopManager::form_elst_dist() {
        int elst_dist_array_slot = sip_tables_.array_slot(std::string("elst_dist"));
        Block::BlockPtr bptr_elst_dist = data_manager_.contiguous_array_manager().get_array(elst_dist_array_slot);
        double *val_elst_dist = bptr_elst_dist->get_data();
        int elst_dist_size = bptr_elst_dist->size();
        int num_frag = upper_bound_[0] - lower_seg_[0];
        elst_dist.resize(num_frag);
        for (int i = 0; i < num_frag; ++i) {
            //std::cout << i << std::endl;
            elst_dist[i].resize(num_frag);
        }
        int ij = 0;
        for (int i = 0; i < num_frag; ++i) {
            for (int j = 0; j < num_frag; ++j) {
                elst_dist[j][i] = (int)val_elst_dist[ij];
                //std::cout << i << " " << j << " " << ij << " " << elst_dist[j][i] << " " << val_elst_dist[ij]<< std::endl;
                ++ij;
            }
        }
        return;
    }

    void Fragment_ij_aaa__PardoLoopManager::form_rcut_dist() {
        int rcut_dist_array_slot = sip_tables_.array_slot(std::string("rcut_dist"));
        Block::BlockPtr bptr_rcut_dist = data_manager_.contiguous_array_manager().get_array(rcut_dist_array_slot);
        double *val_rcut_dist = bptr_rcut_dist->get_data();
        int rcut_dist_size = bptr_rcut_dist->size();
        int num_frag = upper_bound_[0] - lower_seg_[0];
        rcut_dist.resize(num_frag);
        for (int i = 0; i < num_frag; ++i) {
            //std::cout << i << std::endl;
            rcut_dist[i].resize(num_frag);
        }
        int ij = 0;
        for (int i = 0; i < num_frag; ++i) {
            for (int j = 0; j < num_frag; ++j) {
                rcut_dist[j][i] = (int)val_rcut_dist[ij];
                //std::cout << i << " " << j << " " << ij << " " << rcut_dist[j][i] << " " << val_rcut_dist[ij]<< std::endl;
                ++ij;
            }
        }
        return;
    }
    
    void Fragment_ij_aaa__PardoLoopManager::form_swao_frag() {
        int swao_frag_array_slot = sip_tables_.array_slot(std::string("swao_frag"));
        Block::BlockPtr bptr_swao_frag = data_manager_.contiguous_array_manager().get_array(swao_frag_array_slot);
        double *val_swao_frag = bptr_swao_frag->get_data();
        int swao_frag_size = bptr_swao_frag->size();
        swao_frag.resize(swao_frag_size);
        for (int i = 0; i < swao_frag_size; ++i) {
            swao_frag[i] = (int)val_swao_frag[i];
            //std::cout << "vec " << swao_frag[i] << " data " << val_swao_frag[i] << std::endl;
        }
        return;
    }
    
    void Fragment_ij_aaa__PardoLoopManager::form_swocca_frag() {
        int swocca_frag_array_slot = sip_tables_.array_slot(std::string("swocca_frag"));
        Block::BlockPtr bptr_swocca_frag = data_manager_.contiguous_array_manager().get_array(swocca_frag_array_slot);
        double *val_swocca_frag = bptr_swocca_frag->get_data();
        int swocca_frag_size = bptr_swocca_frag->size();
        swocca_frag.resize(swocca_frag_size);
        for (int i = 0; i < swocca_frag_size; ++i) {
            swocca_frag[i] = (int)val_swocca_frag[i];
            //std::cout << "vec " << swocca_frag[i] << " data " << val_swocca_frag[i] << std::endl;
        }
        return;
    }
    
    void Fragment_ij_aaa__PardoLoopManager::form_swvirta_frag() {
        int swvirta_frag_array_slot = sip_tables_.array_slot(std::string("swvirta_frag"));
        Block::BlockPtr bptr_swvirta_frag = data_manager_.contiguous_array_manager().get_array(swvirta_frag_array_slot);
        double *val_swvirta_frag = bptr_swvirta_frag->get_data();
        int swvirta_frag_size = bptr_swvirta_frag->size();
        swvirta_frag.resize(swvirta_frag_size);
        for (int i = 0; i < swvirta_frag_size; ++i) {
            swvirta_frag[i] = (int)val_swvirta_frag[i];
            //std::cout << "vec " << swvirta_frag[i] << " data " << val_swvirta_frag[i] << std::endl;
        }
        return;
    }
    
    inline bool Fragment_ij_aaa__PardoLoopManager::initialize_indices() {
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
            index_values_[i] = lower_seg_[i];
            data_manager_.set_index_value(index_id_[i], lower_seg_[i]);
        }
        more_iterations = true;
        return more_iterations;
    }

/*
 the generic where clauses are:
 typ = 1: map AO to fragment
 typ = 2: map Occ to fragment
 typ = 3: map Virt to fragment
 typ = 4: map elst_dist(ifrag,jfrag) == ifrag
 typ = 5: map rcut_dist(ifrag,jfrag) == ifrag
 typ = 0: ifrag != jfrag
 */
bool Fragment_ij_aaa__PardoLoopManager::fragment_special_where_clause(int typ, int index, int frag) {
    bool where_clause;
    int ij = 0;
    switch (typ) {
            
        case 1: // check against ao index
            //where_clause = return_val_swao_frag(index_values_[index])     == index_values_[frag];
            where_clause = swao_frag[index_values_[index] - lower_seg_[index]]     == index_values_[frag];
            break;
            
        case 2: // check against occ index
            //where_clause = return_val_swocca_frag(index_values_[index])   == index_values_[frag];
            where_clause = swocca_frag[index_values_[index] - lower_seg_[index]]   == index_values_[frag];
            break;
            
        case 3: // check against virt index
            //where_clause = return_val_swvirta_frag(index_values_[index])  == index_values_[frag];
            where_clause = swvirta_frag[index_values_[index] - lower_seg_[index]]  == index_values_[frag];
            break;
            
        case 4: // check elst_dist[ifrag,jfrag] == ifrag
            //where_clause = return_val_elst_dist(index_values_[frag], index_values_[index]) == index_values_[frag];
            //ij = index_values_[frag] - 1 + (upper_bound_[0])*(index_values_[index] - 1);
            //std::cout << index_values_[index] << " " << index_values_[frag] << " " << ij << " " << elst_dist[ij] << std::endl;
            //where_clause = elst_dist[ij] == index_values_[frag];
            where_clause = elst_dist[index_values_[frag] - lower_seg_[frag]][index_values_[index] - lower_seg_[index]] == index_values_[frag];
            break;
            
        case 5: // check rcut_dist[ifrag,jfrag] == ifrag
            //where_clause = return_val_rcut_dist(index_values_[frag], index_values_[index]) == index_values_[frag];
            //ij = index_values_[index] - 1 + upper_bound_[0]*(index_values_[frag] - 1);
            //where_clause = rcut_dist[ij] == index_values_[frag];
            //break;
            where_clause = rcut_dist[index_values_[frag] - lower_seg_[frag]][index_values_[index] - lower_seg_[index]] == index_values_[frag];
            break;
            
        case 0:
            where_clause = index_values_[frag] != index_values_[index];
            break;
        default:
            where_clause = false;
    }
    //std::cout << "where_clause " << typ << " " << index_values_[index] << " " << index_values_[frag] << " " << where_clause << std::endl;
    return where_clause;
}


void Fragment_ij_aaa__PardoLoopManager::do_finalize() {
    for (int i = 0; i < num_indices_; ++i) {
        data_manager_.set_index_undefined(index_id_[i]);
    }
}

std::string Fragment_ij_aaa__PardoLoopManager::to_string() const {
    std::stringstream ss;
    ss << "Fragment_ij_aaa__PardoLoopManager:  num_indices="
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
                         const Fragment_ij_aaa__PardoLoopManager &obj) {
    os << obj.to_string();
    return os;
}

/*!
 -------------------------------------------
 _ij_aa_a_
 -------------------------------------------
 */

    
    /*
     for each new special fragment where clause pattern, this should be the only thing realy changed.
     see comment above for fragment_special_where_clause syntax
     */
    bool Fragment_ij_aa_a_PardoLoopManager::where_clause(int index) {
        bool where_;
        int ifrag = 0;
        int jfrag = 1;
        int ao = 1;
        int occ = 2;
        int virt = 3;
        int elst = 4;
        int rcut = 5;

        
        switch (index) {
            case 0:
            case 1:
                where_ = fragment_special_where_clause(elst,jfrag,ifrag);
                break;
            case 2:
            case 3:
                where_ = fragment_special_where_clause(ao,index,ifrag);
                break;
            case 4:
            case 5:
                where_ = fragment_special_where_clause(ao,index,jfrag);
                break;
            default:
                where_ = false;
        }
        return where_;
    }
    
    bool Fragment_ij_aa_a_PardoLoopManager::do_update() {
        if (to_exit_)
            return false;
        bool more_iters;
        bool where_clauses_value;
        
        interpreter_->skip_where_clauses(num_where_clauses_);
        
        if (first_time_) {
            first_time_ = false;
            more_iters = initialize_indices();
            
            where_clauses_value = true;
            for (int i = 1; i < num_indices_; ++i) {
                where_clauses_value = where_clauses_value && where_clause(i);
            }
            if(where_clauses_value){
                iteration_++;
                if ((iteration_-1) % num_workers_ == company_rank_){
                    return true;
                }
            }
        }
        
        int index_restart[MAX_RANK];
        
        for (int i = 0; i < num_indices_; ++i) {
            index_restart[i] = index_values_[i];
        }
        
        int loop_count = iteration_;
        for (int index_i = index_restart[0]; index_i < upper_bound_[0]; ++index_i) {
            index_values_[0] = index_i;
            data_manager_.set_index_value(index_id_[0], index_i);
            
        for (int index_j = index_restart[1]; index_j < upper_bound_[1]; ++index_j) {
            index_values_[1] = index_j;
            data_manager_.set_index_value(index_id_[1], index_j);
                
            if (where_clause(1)) {
                for (int index_2 = index_restart[2]; index_2 < upper_bound_[2]; ++index_2) {
                    index_values_[2] = index_2;
                    data_manager_.set_index_value(index_id_[2], index_2);
                        
            if (where_clause(2)) {
                for (int index_3 = index_restart[3]; index_3 < upper_bound_[3]; ++index_3) {
                    index_values_[3] = index_3;
                    data_manager_.set_index_value(index_id_[3], index_3);
                    
            if (where_clause(3)) {
                for (int index_4 = index_restart[4]; index_4 < upper_bound_[4]; ++index_4) {
                    index_values_[4] = index_4;
                    data_manager_.set_index_value(index_id_[4], index_4);
                    
            if (where_clause(4)) {
            if (loop_count > iteration_) {
                iteration_++;
                if ((iteration_-1) % num_workers_ == company_rank_){
                    return true;
                }
            }
            ++loop_count;
        } // where 4
            } // index_4
            for (int i = 4; i < num_indices_; ++i) { index_restart[i] = lower_seg_[i]; }
        } // where 3
            } // index_3
            for (int i = 3; i < num_indices_; ++i) { index_restart[i] = lower_seg_[i]; }
        } // where 2
            } // index_2
            for (int i = 2; i < num_indices_; ++i) { index_restart[i] = lower_seg_[i]; }
        } // where 1
            } // index_j
            for (int i = 1; i < num_indices_; ++i) { index_restart[i] = lower_seg_[i]; }
        } // index_i
        
        return false; //this should be false here
    }

Fragment_ij_aa_a_PardoLoopManager::Fragment_ij_aa_a_PardoLoopManager(
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
    
    form_elst_dist();
    //form_rcut_dist();
    form_swao_frag();
    //form_swocca_frag();
    //form_swvirta_frag();
}

Fragment_ij_aa_a_PardoLoopManager::~Fragment_ij_aa_a_PardoLoopManager() {}

bool Fragment_ij_aa_a_PardoLoopManager::increment_simple_pair(){
    bool more = false;      // More iterations?
    int current_value;
    for (int i = 0; i < 2; ++i) {
        more = increment_single_index(i);
        if (more) {break;}
    } //if here, then all indices are at their max value
    return more;
}

bool Fragment_ij_aa_a_PardoLoopManager::increment_single_index(int index){
    bool more = false; 	// More iterations?
    int current_value;
    current_value = data_manager_.index_value(index_id_[index]);
    ++current_value;
    if (current_value < upper_bound_[index]) {
        //increment current index and return
        index_values_[index] = current_value;
        data_manager_.set_index_value(index_id_[index], current_value);
        more = true;
    } else {
        //wrap around and handle next index
        index_values_[index] = lower_seg_[index];
        data_manager_.set_index_value(index_id_[index], lower_seg_[index]);
    }
    return more;
}
    
    
    inline bool Fragment_ij_aa_a_PardoLoopManager::increment_all() {
        bool more = false;      // More iterations?
        int current_value;
        for (int i = 0; i < num_indices_; ++i) {
            more = increment_single_index(i);
            if (more) {break;}
        } //if here, then all indices are at their max value
        return more;
    }
    
    void Fragment_ij_aa_a_PardoLoopManager::form_elst_dist() {
        int elst_dist_array_slot = sip_tables_.array_slot(std::string("elst_dist"));
        Block::BlockPtr bptr_elst_dist = data_manager_.contiguous_array_manager().get_array(elst_dist_array_slot);
        double *val_elst_dist = bptr_elst_dist->get_data();
        int elst_dist_size = bptr_elst_dist->size();
        int num_frag = upper_bound_[0] - lower_seg_[0];
        elst_dist.resize(num_frag);
        for (int i = 0; i < num_frag; ++i) {
            //std::cout << i << std::endl;
            elst_dist[i].resize(num_frag);
        }
        int ij = 0;
        for (int i = 0; i < num_frag; ++i) {
            for (int j = 0; j < num_frag; ++j) {
                elst_dist[j][i] = (int)val_elst_dist[ij];
                //std::cout << i << " " << j << " " << ij << " " << elst_dist[j][i] << " " << val_elst_dist[ij]<< std::endl;
                ++ij;
            }
        }
        return;
    }

    void Fragment_ij_aa_a_PardoLoopManager::form_rcut_dist() {
        int rcut_dist_array_slot = sip_tables_.array_slot(std::string("rcut_dist"));
        Block::BlockPtr bptr_rcut_dist = data_manager_.contiguous_array_manager().get_array(rcut_dist_array_slot);
        double *val_rcut_dist = bptr_rcut_dist->get_data();
        int rcut_dist_size = bptr_rcut_dist->size();
        int num_frag = upper_bound_[0] - lower_seg_[0];
        rcut_dist.resize(num_frag);
        for (int i = 0; i < num_frag; ++i) {
            //std::cout << i << std::endl;
            rcut_dist[i].resize(num_frag);
        }
        int ij = 0;
        for (int i = 0; i < num_frag; ++i) {
            for (int j = 0; j < num_frag; ++j) {
                rcut_dist[j][i] = (int)val_rcut_dist[ij];
                //std::cout << i << " " << j << " " << ij << " " << rcut_dist[j][i] << " " << val_rcut_dist[ij]<< std::endl;
                ++ij;
            }
        }
        return;
    }
    
    void Fragment_ij_aa_a_PardoLoopManager::form_swao_frag() {
        int swao_frag_array_slot = sip_tables_.array_slot(std::string("swao_frag"));
        Block::BlockPtr bptr_swao_frag = data_manager_.contiguous_array_manager().get_array(swao_frag_array_slot);
        double *val_swao_frag = bptr_swao_frag->get_data();
        int swao_frag_size = bptr_swao_frag->size();
        swao_frag.resize(swao_frag_size);
        for (int i = 0; i < swao_frag_size; ++i) {
            swao_frag[i] = (int)val_swao_frag[i];
            //std::cout << "vec " << swao_frag[i] << " data " << val_swao_frag[i] << std::endl;
        }
        return;
    }
    
    void Fragment_ij_aa_a_PardoLoopManager::form_swocca_frag() {
        int swocca_frag_array_slot = sip_tables_.array_slot(std::string("swocca_frag"));
        Block::BlockPtr bptr_swocca_frag = data_manager_.contiguous_array_manager().get_array(swocca_frag_array_slot);
        double *val_swocca_frag = bptr_swocca_frag->get_data();
        int swocca_frag_size = bptr_swocca_frag->size();
        swocca_frag.resize(swocca_frag_size);
        for (int i = 0; i < swocca_frag_size; ++i) {
            swocca_frag[i] = (int)val_swocca_frag[i];
            //std::cout << "vec " << swocca_frag[i] << " data " << val_swocca_frag[i] << std::endl;
        }
        return;
    }
    
    void Fragment_ij_aa_a_PardoLoopManager::form_swvirta_frag() {
        int swvirta_frag_array_slot = sip_tables_.array_slot(std::string("swvirta_frag"));
        Block::BlockPtr bptr_swvirta_frag = data_manager_.contiguous_array_manager().get_array(swvirta_frag_array_slot);
        double *val_swvirta_frag = bptr_swvirta_frag->get_data();
        int swvirta_frag_size = bptr_swvirta_frag->size();
        swvirta_frag.resize(swvirta_frag_size);
        for (int i = 0; i < swvirta_frag_size; ++i) {
            swvirta_frag[i] = (int)val_swvirta_frag[i];
            //std::cout << "vec " << swvirta_frag[i] << " data " << val_swvirta_frag[i] << std::endl;
        }
        return;
    }
    
    inline bool Fragment_ij_aa_a_PardoLoopManager::initialize_indices() {
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
            index_values_[i] = lower_seg_[i];
            data_manager_.set_index_value(index_id_[i], lower_seg_[i]);
        }
        more_iterations = true;
        return more_iterations;
    }

/*
 the generic where clauses are:
 typ = 1: map AO to fragment
 typ = 2: map Occ to fragment
 typ = 3: map Virt to fragment
 typ = 4: map elst_dist(ifrag,jfrag) == ifrag
 typ = 5: map rcut_dist(ifrag,jfrag) == ifrag
 typ = 0: ifrag != jfrag
 */
bool Fragment_ij_aa_a_PardoLoopManager::fragment_special_where_clause(int typ, int index, int frag) {
    bool where_clause;
    int ij = 0;
    switch (typ) {
            
        case 1: // check against ao index
            //where_clause = return_val_swao_frag(index_values_[index])     == index_values_[frag];
            where_clause = swao_frag[index_values_[index] - lower_seg_[index]]     == index_values_[frag];
            break;
            
        case 2: // check against occ index
            //where_clause = return_val_swocca_frag(index_values_[index])   == index_values_[frag];
            where_clause = swocca_frag[index_values_[index] - lower_seg_[index]]   == index_values_[frag];
            break;
            
        case 3: // check against virt index
            //where_clause = return_val_swvirta_frag(index_values_[index])  == index_values_[frag];
            where_clause = swvirta_frag[index_values_[index] - lower_seg_[index]]  == index_values_[frag];
            break;
            
        case 4: // check elst_dist[ifrag,jfrag] == ifrag
            //where_clause = return_val_elst_dist(index_values_[frag], index_values_[index]) == index_values_[frag];
            //ij = index_values_[frag] - 1 + (upper_bound_[0])*(index_values_[index] - 1);
            //std::cout << index_values_[index] << " " << index_values_[frag] << " " << ij << " " << elst_dist[ij] << std::endl;
            //where_clause = elst_dist[ij] == index_values_[frag];
            where_clause = elst_dist[index_values_[frag] - lower_seg_[frag]][index_values_[index] - lower_seg_[index]] == index_values_[frag];
            break;
            
        case 5: // check rcut_dist[ifrag,jfrag] == ifrag
            //where_clause = return_val_rcut_dist(index_values_[frag], index_values_[index]) == index_values_[frag];
            //ij = index_values_[index] - 1 + upper_bound_[0]*(index_values_[frag] - 1);
            //where_clause = rcut_dist[ij] == index_values_[frag];
            //break;
            where_clause = rcut_dist[index_values_[frag] - lower_seg_[frag]][index_values_[index] - lower_seg_[index]] == index_values_[frag];
            break;
            
        case 0:
            where_clause = index_values_[frag] != index_values_[index];
            break;
        default:
            where_clause = false;
    }
    //std::cout << "where_clause " << typ << " " << index_values_[index] << " " << index_values_[frag] << " " << where_clause << std::endl;
    return where_clause;
}


void Fragment_ij_aa_a_PardoLoopManager::do_finalize() {
    for (int i = 0; i < num_indices_; ++i) {
        data_manager_.set_index_undefined(index_id_[i]);
    }
}

std::string Fragment_ij_aa_a_PardoLoopManager::to_string() const {
    std::stringstream ss;
    ss << "Fragment_ij_aa_a_PardoLoopManager:  num_indices="
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
                         const Fragment_ij_aa_a_PardoLoopManager &obj) {
    os << obj.to_string();
    return os;
}

/*!
 -------------------------------------------
 _ij_ao_ao_
 -------------------------------------------
 */

    
    /*
     for each new special fragment where clause pattern, this should be the only thing realy changed.
     see comment above for fragment_special_where_clause syntax
     */
    bool Fragment_ij_ao_ao_PardoLoopManager::where_clause(int index) {
        bool where_;
        int ifrag = 0;
        int jfrag = 1;
        int ao = 1;
        int occ = 2;
        int virt = 3;
        int elst = 4;
        int rcut = 5;

        
        switch (index) {
            case 0:
            case 1:
                where_ = fragment_special_where_clause(elst,jfrag,ifrag);
                break;
            case 2:
                where_ = fragment_special_where_clause(ao,index,ifrag);
                break;
            case 3:
                where_ = fragment_special_where_clause(occ,index,ifrag);
                break;
            case 4:
                where_ = fragment_special_where_clause(ao,index,jfrag);
                break;
            case 5:
                where_ = fragment_special_where_clause(occ,index,jfrag);
                break;
            default:
                where_ = false;
        }
        return where_;
    }
    
    bool Fragment_ij_ao_ao_PardoLoopManager::do_update() {
        if (to_exit_)
            return false;
        bool more_iters;
        bool where_clauses_value;
        
        interpreter_->skip_where_clauses(num_where_clauses_);
        
        if (first_time_) {
            first_time_ = false;
            more_iters = initialize_indices();
            
            where_clauses_value = true;
            for (int i = 1; i < num_indices_; ++i) {
                where_clauses_value = where_clauses_value && where_clause(i);
            }
            if(where_clauses_value){
                iteration_++;
                if ((iteration_-1) % num_workers_ == company_rank_){
                    return true;
                }
            }
        }
        
        int index_restart[MAX_RANK];
        
        for (int i = 0; i < num_indices_; ++i) {
            index_restart[i] = index_values_[i];
        }
        
        int loop_count = iteration_;
        for (int index_i = index_restart[0]; index_i < upper_bound_[0]; ++index_i) {
            index_values_[0] = index_i;
            data_manager_.set_index_value(index_id_[0], index_i);
            
        for (int index_j = index_restart[1]; index_j < upper_bound_[1]; ++index_j) {
            index_values_[1] = index_j;
            data_manager_.set_index_value(index_id_[1], index_j);
                
            if (where_clause(1)) {
                for (int index_2 = index_restart[2]; index_2 < upper_bound_[2]; ++index_2) {
                    index_values_[2] = index_2;
                    data_manager_.set_index_value(index_id_[2], index_2);
                        
            if (where_clause(2)) {
                for (int index_3 = index_restart[3]; index_3 < upper_bound_[3]; ++index_3) {
                    index_values_[3] = index_3;
                    data_manager_.set_index_value(index_id_[3], index_3);
                    
            if (where_clause(3)) {
                for (int index_4 = index_restart[4]; index_4 < upper_bound_[4]; ++index_4) {
                    index_values_[4] = index_4;
                    data_manager_.set_index_value(index_id_[4], index_4);
                    
            if (where_clause(4)) {
                for (int index_5 = index_restart[5]; index_5 < upper_bound_[5]; ++index_5) {
                    index_values_[5] = index_5;
                    data_manager_.set_index_value(index_id_[5], index_5);
                    
            if (where_clause(5)) {
            if (loop_count > iteration_) {
                iteration_++;
                if ((iteration_-1) % num_workers_ == company_rank_){
                    return true;
                }
            }
            ++loop_count;
        } // where 5
            } // index_5
            for (int i = 5; i < num_indices_; ++i) { index_restart[i] = lower_seg_[i]; }
        } // where 4
            } // index_4
            for (int i = 4; i < num_indices_; ++i) { index_restart[i] = lower_seg_[i]; }
        } // where 3
            } // index_3
            for (int i = 3; i < num_indices_; ++i) { index_restart[i] = lower_seg_[i]; }
        } // where 2
            } // index_2
            for (int i = 2; i < num_indices_; ++i) { index_restart[i] = lower_seg_[i]; }
        } // where 1
            } // index_j
            for (int i = 1; i < num_indices_; ++i) { index_restart[i] = lower_seg_[i]; }
        } // index_i
        
        return false; //this should be false here
    }
    
Fragment_ij_ao_ao_PardoLoopManager::Fragment_ij_ao_ao_PardoLoopManager(
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
    
    form_elst_dist();
    //form_rcut_dist();
    form_swao_frag();
    form_swocca_frag();
    //form_swvirta_frag();
}

Fragment_ij_ao_ao_PardoLoopManager::~Fragment_ij_ao_ao_PardoLoopManager() {}

bool Fragment_ij_ao_ao_PardoLoopManager::increment_simple_pair(){
    bool more = false;      // More iterations?
    int current_value;
    for (int i = 0; i < 2; ++i) {
        more = increment_single_index(i);
        if (more) {break;}
    } //if here, then all indices are at their max value
    return more;
}

bool Fragment_ij_ao_ao_PardoLoopManager::increment_single_index(int index){
    bool more = false; 	// More iterations?
    int current_value;
    current_value = data_manager_.index_value(index_id_[index]);
    ++current_value;
    if (current_value < upper_bound_[index]) {
        //increment current index and return
        index_values_[index] = current_value;
        data_manager_.set_index_value(index_id_[index], current_value);
        more = true;
    } else {
        //wrap around and handle next index
        index_values_[index] = lower_seg_[index];
        data_manager_.set_index_value(index_id_[index], lower_seg_[index]);
    }
    return more;
}
    
    
    inline bool Fragment_ij_ao_ao_PardoLoopManager::increment_all() {
        bool more = false;      // More iterations?
        int current_value;
        for (int i = 0; i < num_indices_; ++i) {
            more = increment_single_index(i);
            if (more) {break;}
        } //if here, then all indices are at their max value
        return more;
    }
    
    void Fragment_ij_ao_ao_PardoLoopManager::form_elst_dist() {
        int elst_dist_array_slot = sip_tables_.array_slot(std::string("elst_dist"));
        Block::BlockPtr bptr_elst_dist = data_manager_.contiguous_array_manager().get_array(elst_dist_array_slot);
        double *val_elst_dist = bptr_elst_dist->get_data();
        int elst_dist_size = bptr_elst_dist->size();
        int num_frag = upper_bound_[0] - lower_seg_[0];
        elst_dist.resize(num_frag);
        for (int i = 0; i < num_frag; ++i) {
            //std::cout << i << std::endl;
            elst_dist[i].resize(num_frag);
        }
        int ij = 0;
        for (int i = 0; i < num_frag; ++i) {
            for (int j = 0; j < num_frag; ++j) {
                elst_dist[j][i] = (int)val_elst_dist[ij];
                //std::cout << i << " " << j << " " << ij << " " << elst_dist[j][i] << " " << val_elst_dist[ij]<< std::endl;
                ++ij;
            }
        }
        return;
    }

    void Fragment_ij_ao_ao_PardoLoopManager::form_rcut_dist() {
        int rcut_dist_array_slot = sip_tables_.array_slot(std::string("rcut_dist"));
        Block::BlockPtr bptr_rcut_dist = data_manager_.contiguous_array_manager().get_array(rcut_dist_array_slot);
        double *val_rcut_dist = bptr_rcut_dist->get_data();
        int rcut_dist_size = bptr_rcut_dist->size();
        int num_frag = upper_bound_[0] - lower_seg_[0];
        rcut_dist.resize(num_frag);
        for (int i = 0; i < num_frag; ++i) {
            //std::cout << i << std::endl;
            rcut_dist[i].resize(num_frag);
        }
        int ij = 0;
        for (int i = 0; i < num_frag; ++i) {
            for (int j = 0; j < num_frag; ++j) {
                rcut_dist[j][i] = (int)val_rcut_dist[ij];
                //std::cout << i << " " << j << " " << ij << " " << rcut_dist[j][i] << " " << val_rcut_dist[ij]<< std::endl;
                ++ij;
            }
        }
        return;
    }
    
    void Fragment_ij_ao_ao_PardoLoopManager::form_swao_frag() {
        int swao_frag_array_slot = sip_tables_.array_slot(std::string("swao_frag"));
        Block::BlockPtr bptr_swao_frag = data_manager_.contiguous_array_manager().get_array(swao_frag_array_slot);
        double *val_swao_frag = bptr_swao_frag->get_data();
        int swao_frag_size = bptr_swao_frag->size();
        swao_frag.resize(swao_frag_size);
        for (int i = 0; i < swao_frag_size; ++i) {
            swao_frag[i] = (int)val_swao_frag[i];
            //std::cout << "vec " << swao_frag[i] << " data " << val_swao_frag[i] << std::endl;
        }
        return;
    }
    
    void Fragment_ij_ao_ao_PardoLoopManager::form_swocca_frag() {
        int swocca_frag_array_slot = sip_tables_.array_slot(std::string("swocca_frag"));
        Block::BlockPtr bptr_swocca_frag = data_manager_.contiguous_array_manager().get_array(swocca_frag_array_slot);
        double *val_swocca_frag = bptr_swocca_frag->get_data();
        int swocca_frag_size = bptr_swocca_frag->size();
        swocca_frag.resize(swocca_frag_size);
        for (int i = 0; i < swocca_frag_size; ++i) {
            swocca_frag[i] = (int)val_swocca_frag[i];
            //std::cout << "vec " << swocca_frag[i] << " data " << val_swocca_frag[i] << std::endl;
        }
        return;
    }
    
    void Fragment_ij_ao_ao_PardoLoopManager::form_swvirta_frag() {
        int swvirta_frag_array_slot = sip_tables_.array_slot(std::string("swvirta_frag"));
        Block::BlockPtr bptr_swvirta_frag = data_manager_.contiguous_array_manager().get_array(swvirta_frag_array_slot);
        double *val_swvirta_frag = bptr_swvirta_frag->get_data();
        int swvirta_frag_size = bptr_swvirta_frag->size();
        swvirta_frag.resize(swvirta_frag_size);
        for (int i = 0; i < swvirta_frag_size; ++i) {
            swvirta_frag[i] = (int)val_swvirta_frag[i];
            //std::cout << "vec " << swvirta_frag[i] << " data " << val_swvirta_frag[i] << std::endl;
        }
        return;
    }
    
    inline bool Fragment_ij_ao_ao_PardoLoopManager::initialize_indices() {
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
            index_values_[i] = lower_seg_[i];
            data_manager_.set_index_value(index_id_[i], lower_seg_[i]);
        }
        more_iterations = true;
        return more_iterations;
    }

/*
 the generic where clauses are:
 typ = 1: map AO to fragment
 typ = 2: map Occ to fragment
 typ = 3: map Virt to fragment
 typ = 4: map elst_dist(ifrag,jfrag) == ifrag
 typ = 5: map rcut_dist(ifrag,jfrag) == ifrag
 typ = 0: ifrag != jfrag
 */
bool Fragment_ij_ao_ao_PardoLoopManager::fragment_special_where_clause(int typ, int index, int frag) {
    bool where_clause;
    int ij = 0;
    switch (typ) {
            
        case 1: // check against ao index
            //where_clause = return_val_swao_frag(index_values_[index])     == index_values_[frag];
            where_clause = swao_frag[index_values_[index] - lower_seg_[index]]     == index_values_[frag];
            break;
            
        case 2: // check against occ index
            //where_clause = return_val_swocca_frag(index_values_[index])   == index_values_[frag];
            where_clause = swocca_frag[index_values_[index] - lower_seg_[index]]   == index_values_[frag];
            break;
            
        case 3: // check against virt index
            //where_clause = return_val_swvirta_frag(index_values_[index])  == index_values_[frag];
            where_clause = swvirta_frag[index_values_[index] - lower_seg_[index]]  == index_values_[frag];
            break;
            
        case 4: // check elst_dist[ifrag,jfrag] == ifrag
            //where_clause = return_val_elst_dist(index_values_[frag], index_values_[index]) == index_values_[frag];
            //ij = index_values_[frag] - 1 + (upper_bound_[0])*(index_values_[index] - 1);
            //std::cout << index_values_[index] << " " << index_values_[frag] << " " << ij << " " << elst_dist[ij] << std::endl;
            //where_clause = elst_dist[ij] == index_values_[frag];
            where_clause = elst_dist[index_values_[frag] - lower_seg_[frag]][index_values_[index] - lower_seg_[index]] == index_values_[frag];
            break;
            
        case 5: // check rcut_dist[ifrag,jfrag] == ifrag
            //where_clause = return_val_rcut_dist(index_values_[frag], index_values_[index]) == index_values_[frag];
            //ij = index_values_[index] - 1 + upper_bound_[0]*(index_values_[frag] - 1);
            //where_clause = rcut_dist[ij] == index_values_[frag];
            //break;
            where_clause = rcut_dist[index_values_[frag] - lower_seg_[frag]][index_values_[index] - lower_seg_[index]] == index_values_[frag];
            break;
            
        case 0:
            where_clause = index_values_[frag] != index_values_[index];
            break;
        default:
            where_clause = false;
    }
    //std::cout << "where_clause " << typ << " " << index_values_[index] << " " << index_values_[frag] << " " << where_clause << std::endl;
    return where_clause;
}


void Fragment_ij_ao_ao_PardoLoopManager::do_finalize() {
    for (int i = 0; i < num_indices_; ++i) {
        data_manager_.set_index_undefined(index_id_[i]);
    }
}

std::string Fragment_ij_ao_ao_PardoLoopManager::to_string() const {
    std::stringstream ss;
    ss << "Fragment_ij_ao_ao_PardoLoopManager:  num_indices="
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
                         const Fragment_ij_ao_ao_PardoLoopManager &obj) {
    os << obj.to_string();
    return os;
}

/*!
 -------------------------------------------
 _ij_aa_oo_
 -------------------------------------------
 */

    
    /*
     for each new special fragment where clause pattern, this should be the only thing realy changed.
     see comment above for fragment_special_where_clause syntax
     */
    bool Fragment_ij_aa_oo_PardoLoopManager::where_clause(int index) {
        bool where_;
        int ifrag = 0;
        int jfrag = 1;
        int ao = 1;
        int occ = 2;
        int virt = 3;
        int elst = 4;
        int rcut = 5;

        
        switch (index) {
            case 0:
            case 1:
                where_ = fragment_special_where_clause(elst,jfrag,ifrag);
                break;
            case 2:
                where_ = fragment_special_where_clause(ao,index,ifrag);
                break;
            case 3:
                where_ = fragment_special_where_clause(ao,index,ifrag);
                break;
            case 4:
                where_ = fragment_special_where_clause(occ,index,jfrag);
                break;
            case 5:
                where_ = fragment_special_where_clause(occ,index,jfrag);
                break;
            default:
                where_ = false;
        }
        return where_;
    }
    
    bool Fragment_ij_aa_oo_PardoLoopManager::do_update() {
        if (to_exit_)
            return false;
        bool more_iters;
        bool where_clauses_value;
        
        interpreter_->skip_where_clauses(num_where_clauses_);
        
        if (first_time_) {
            first_time_ = false;
            more_iters = initialize_indices();
            
            where_clauses_value = true;
            for (int i = 1; i < num_indices_; ++i) {
                where_clauses_value = where_clauses_value && where_clause(i);
            }
            if(where_clauses_value){
                iteration_++;
                if ((iteration_-1) % num_workers_ == company_rank_){
                    return true;
                }
            }
        }
        
        int index_restart[MAX_RANK];
        
        for (int i = 0; i < num_indices_; ++i) {
            index_restart[i] = index_values_[i];
        }
        
        int loop_count = iteration_;
        for (int index_i = index_restart[0]; index_i < upper_bound_[0]; ++index_i) {
            index_values_[0] = index_i;
            data_manager_.set_index_value(index_id_[0], index_i);
            
        for (int index_j = index_restart[1]; index_j < upper_bound_[1]; ++index_j) {
            index_values_[1] = index_j;
            data_manager_.set_index_value(index_id_[1], index_j);
                
            if (where_clause(1)) {
                for (int index_2 = index_restart[2]; index_2 < upper_bound_[2]; ++index_2) {
                    index_values_[2] = index_2;
                    data_manager_.set_index_value(index_id_[2], index_2);
                        
            if (where_clause(2)) {
                for (int index_3 = index_restart[3]; index_3 < upper_bound_[3]; ++index_3) {
                    index_values_[3] = index_3;
                    data_manager_.set_index_value(index_id_[3], index_3);
                    
            if (where_clause(3)) {
                for (int index_4 = index_restart[4]; index_4 < upper_bound_[4]; ++index_4) {
                    index_values_[4] = index_4;
                    data_manager_.set_index_value(index_id_[4], index_4);
                    
            if (where_clause(4)) {
                for (int index_5 = index_restart[5]; index_5 < upper_bound_[5]; ++index_5) {
                    index_values_[5] = index_5;
                    data_manager_.set_index_value(index_id_[5], index_5);
                    
            if (where_clause(5)) {
            if (loop_count > iteration_) {
                iteration_++;
                if ((iteration_-1) % num_workers_ == company_rank_){
                    return true;
                }
            }
            ++loop_count;
        } // where 5
            } // index_5
            for (int i = 5; i < num_indices_; ++i) { index_restart[i] = lower_seg_[i]; }
        } // where 4
            } // index_4
            for (int i = 4; i < num_indices_; ++i) { index_restart[i] = lower_seg_[i]; }
        } // where 3
            } // index_3
            for (int i = 3; i < num_indices_; ++i) { index_restart[i] = lower_seg_[i]; }
        } // where 2
            } // index_2
            for (int i = 2; i < num_indices_; ++i) { index_restart[i] = lower_seg_[i]; }
        } // where 1
            } // index_j
            for (int i = 1; i < num_indices_; ++i) { index_restart[i] = lower_seg_[i]; }
        } // index_i
        
        return false; //this should be false here
    }

Fragment_ij_aa_oo_PardoLoopManager::Fragment_ij_aa_oo_PardoLoopManager(
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
    
    form_elst_dist();
    //form_rcut_dist();
    form_swao_frag();
    form_swocca_frag();
    //form_swvirta_frag();
}

Fragment_ij_aa_oo_PardoLoopManager::~Fragment_ij_aa_oo_PardoLoopManager() {}

bool Fragment_ij_aa_oo_PardoLoopManager::increment_simple_pair(){
    bool more = false;      // More iterations?
    int current_value;
    for (int i = 0; i < 2; ++i) {
        more = increment_single_index(i);
        if (more) {break;}
    } //if here, then all indices are at their max value
    return more;
}

bool Fragment_ij_aa_oo_PardoLoopManager::increment_single_index(int index){
    bool more = false; 	// More iterations?
    int current_value;
    current_value = data_manager_.index_value(index_id_[index]);
    ++current_value;
    if (current_value < upper_bound_[index]) {
        //increment current index and return
        index_values_[index] = current_value;
        data_manager_.set_index_value(index_id_[index], current_value);
        more = true;
    } else {
        //wrap around and handle next index
        index_values_[index] = lower_seg_[index];
        data_manager_.set_index_value(index_id_[index], lower_seg_[index]);
    }
    return more;
}
    
    
    inline bool Fragment_ij_aa_oo_PardoLoopManager::increment_all() {
        bool more = false;      // More iterations?
        int current_value;
        for (int i = 0; i < num_indices_; ++i) {
            more = increment_single_index(i);
            if (more) {break;}
        } //if here, then all indices are at their max value
        return more;
    }
    
    void Fragment_ij_aa_oo_PardoLoopManager::form_elst_dist() {
        int elst_dist_array_slot = sip_tables_.array_slot(std::string("elst_dist"));
        Block::BlockPtr bptr_elst_dist = data_manager_.contiguous_array_manager().get_array(elst_dist_array_slot);
        double *val_elst_dist = bptr_elst_dist->get_data();
        int elst_dist_size = bptr_elst_dist->size();
        int num_frag = upper_bound_[0] - lower_seg_[0];
        elst_dist.resize(num_frag);
        for (int i = 0; i < num_frag; ++i) {
            //std::cout << i << std::endl;
            elst_dist[i].resize(num_frag);
        }
        int ij = 0;
        for (int i = 0; i < num_frag; ++i) {
            for (int j = 0; j < num_frag; ++j) {
                elst_dist[j][i] = (int)val_elst_dist[ij];
                //std::cout << i << " " << j << " " << ij << " " << elst_dist[j][i] << " " << val_elst_dist[ij]<< std::endl;
                ++ij;
            }
        }
        return;
    }

    void Fragment_ij_aa_oo_PardoLoopManager::form_rcut_dist() {
        int rcut_dist_array_slot = sip_tables_.array_slot(std::string("rcut_dist"));
        Block::BlockPtr bptr_rcut_dist = data_manager_.contiguous_array_manager().get_array(rcut_dist_array_slot);
        double *val_rcut_dist = bptr_rcut_dist->get_data();
        int rcut_dist_size = bptr_rcut_dist->size();
        int num_frag = upper_bound_[0] - lower_seg_[0];
        rcut_dist.resize(num_frag);
        for (int i = 0; i < num_frag; ++i) {
            //std::cout << i << std::endl;
            rcut_dist[i].resize(num_frag);
        }
        int ij = 0;
        for (int i = 0; i < num_frag; ++i) {
            for (int j = 0; j < num_frag; ++j) {
                rcut_dist[j][i] = (int)val_rcut_dist[ij];
                //std::cout << i << " " << j << " " << ij << " " << rcut_dist[j][i] << " " << val_rcut_dist[ij]<< std::endl;
                ++ij;
            }
        }
        return;
    }
    
    void Fragment_ij_aa_oo_PardoLoopManager::form_swao_frag() {
        int swao_frag_array_slot = sip_tables_.array_slot(std::string("swao_frag"));
        Block::BlockPtr bptr_swao_frag = data_manager_.contiguous_array_manager().get_array(swao_frag_array_slot);
        double *val_swao_frag = bptr_swao_frag->get_data();
        int swao_frag_size = bptr_swao_frag->size();
        swao_frag.resize(swao_frag_size);
        for (int i = 0; i < swao_frag_size; ++i) {
            swao_frag[i] = (int)val_swao_frag[i];
            //std::cout << "vec " << swao_frag[i] << " data " << val_swao_frag[i] << std::endl;
        }
        return;
    }
    
    void Fragment_ij_aa_oo_PardoLoopManager::form_swocca_frag() {
        int swocca_frag_array_slot = sip_tables_.array_slot(std::string("swocca_frag"));
        Block::BlockPtr bptr_swocca_frag = data_manager_.contiguous_array_manager().get_array(swocca_frag_array_slot);
        double *val_swocca_frag = bptr_swocca_frag->get_data();
        int swocca_frag_size = bptr_swocca_frag->size();
        swocca_frag.resize(swocca_frag_size);
        for (int i = 0; i < swocca_frag_size; ++i) {
            swocca_frag[i] = (int)val_swocca_frag[i];
            //std::cout << "vec " << swocca_frag[i] << " data " << val_swocca_frag[i] << std::endl;
        }
        return;
    }
    
    void Fragment_ij_aa_oo_PardoLoopManager::form_swvirta_frag() {
        int swvirta_frag_array_slot = sip_tables_.array_slot(std::string("swvirta_frag"));
        Block::BlockPtr bptr_swvirta_frag = data_manager_.contiguous_array_manager().get_array(swvirta_frag_array_slot);
        double *val_swvirta_frag = bptr_swvirta_frag->get_data();
        int swvirta_frag_size = bptr_swvirta_frag->size();
        swvirta_frag.resize(swvirta_frag_size);
        for (int i = 0; i < swvirta_frag_size; ++i) {
            swvirta_frag[i] = (int)val_swvirta_frag[i];
            //std::cout << "vec " << swvirta_frag[i] << " data " << val_swvirta_frag[i] << std::endl;
        }
        return;
    }
    
    inline bool Fragment_ij_aa_oo_PardoLoopManager::initialize_indices() {
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
            index_values_[i] = lower_seg_[i];
            data_manager_.set_index_value(index_id_[i], lower_seg_[i]);
        }
        more_iterations = true;
        return more_iterations;
    }

/*
 the generic where clauses are:
 typ = 1: map AO to fragment
 typ = 2: map Occ to fragment
 typ = 3: map Virt to fragment
 typ = 4: map elst_dist(ifrag,jfrag) == ifrag
 typ = 5: map rcut_dist(ifrag,jfrag) == ifrag
 typ = 0: ifrag != jfrag
 */
bool Fragment_ij_aa_oo_PardoLoopManager::fragment_special_where_clause(int typ, int index, int frag) {
    bool where_clause;
    int ij = 0;
    switch (typ) {
            
        case 1: // check against ao index
            //where_clause = return_val_swao_frag(index_values_[index])     == index_values_[frag];
            where_clause = swao_frag[index_values_[index] - lower_seg_[index]]     == index_values_[frag];
            break;
            
        case 2: // check against occ index
            //where_clause = return_val_swocca_frag(index_values_[index])   == index_values_[frag];
            where_clause = swocca_frag[index_values_[index] - lower_seg_[index]]   == index_values_[frag];
            break;
            
        case 3: // check against virt index
            //where_clause = return_val_swvirta_frag(index_values_[index])  == index_values_[frag];
            where_clause = swvirta_frag[index_values_[index] - lower_seg_[index]]  == index_values_[frag];
            break;
            
        case 4: // check elst_dist[ifrag,jfrag] == ifrag
            //where_clause = return_val_elst_dist(index_values_[frag], index_values_[index]) == index_values_[frag];
            //ij = index_values_[frag] - 1 + (upper_bound_[0])*(index_values_[index] - 1);
            //std::cout << index_values_[index] << " " << index_values_[frag] << " " << ij << " " << elst_dist[ij] << std::endl;
            //where_clause = elst_dist[ij] == index_values_[frag];
            where_clause = elst_dist[index_values_[frag] - lower_seg_[frag]][index_values_[index] - lower_seg_[index]] == index_values_[frag];
            break;
            
        case 5: // check rcut_dist[ifrag,jfrag] == ifrag
            //where_clause = return_val_rcut_dist(index_values_[frag], index_values_[index]) == index_values_[frag];
            //ij = index_values_[index] - 1 + upper_bound_[0]*(index_values_[frag] - 1);
            //where_clause = rcut_dist[ij] == index_values_[frag];
            //break;
            where_clause = rcut_dist[index_values_[frag] - lower_seg_[frag]][index_values_[index] - lower_seg_[index]] == index_values_[frag];
            break;
            
        case 0:
            where_clause = index_values_[frag] != index_values_[index];
            break;
        default:
            where_clause = false;
    }
    //std::cout << "where_clause " << typ << " " << index_values_[index] << " " << index_values_[frag] << " " << where_clause << std::endl;
    return where_clause;
}


void Fragment_ij_aa_oo_PardoLoopManager::do_finalize() {
    for (int i = 0; i < num_indices_; ++i) {
        data_manager_.set_index_undefined(index_id_[i]);
    }
}

std::string Fragment_ij_aa_oo_PardoLoopManager::to_string() const {
    std::stringstream ss;
    ss << "Fragment_ij_aa_oo_PardoLoopManager:  num_indices="
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
                         const Fragment_ij_aa_oo_PardoLoopManager &obj) {
    os << obj.to_string();
    return os;
}

/*!
 -------------------------------------------
 _ij_aoa_o_
 -------------------------------------------
 */

    
    /*
     for each new special fragment where clause pattern, this should be the only thing realy changed.
     see comment above for fragment_special_where_clause syntax
     */
    bool Fragment_ij_aoa_o_PardoLoopManager::where_clause(int index) {
        bool where_;
        int ifrag = 0;
        int jfrag = 1;
        int ao = 1;
        int occ = 2;
        int virt = 3;
        int elst = 4;
        int rcut = 5;

        
        switch (index) {
            case 0:
            case 1:
                where_ = fragment_special_where_clause(elst,jfrag,ifrag);
                break;
            case 2:
                where_ = fragment_special_where_clause(ao,index,ifrag);
                break;
            case 3:
                where_ = fragment_special_where_clause(occ,index,ifrag);
                break;
            case 4:
                where_ = fragment_special_where_clause(ao,index,ifrag);
                break;
            case 5:
                where_ = fragment_special_where_clause(occ,index,jfrag);
                break;
            default:
                where_ = false;
        }
        return where_;
    }
    
    bool Fragment_ij_aoa_o_PardoLoopManager::do_update() {
        if (to_exit_)
            return false;
        bool more_iters;
        bool where_clauses_value;
        
        interpreter_->skip_where_clauses(num_where_clauses_);
        
        if (first_time_) {
            first_time_ = false;
            more_iters = initialize_indices();
            
            where_clauses_value = true;
            for (int i = 1; i < num_indices_; ++i) {
                where_clauses_value = where_clauses_value && where_clause(i);
            }
            if(where_clauses_value){
                iteration_++;
                if ((iteration_-1) % num_workers_ == company_rank_){
                    return true;
                }
            }
        }
        
        int index_restart[MAX_RANK];
        
        for (int i = 0; i < num_indices_; ++i) {
            index_restart[i] = index_values_[i];
        }
        
        int loop_count = iteration_;
        for (int index_i = index_restart[0]; index_i < upper_bound_[0]; ++index_i) {
            index_values_[0] = index_i;
            data_manager_.set_index_value(index_id_[0], index_i);
            
        for (int index_j = index_restart[1]; index_j < upper_bound_[1]; ++index_j) {
            index_values_[1] = index_j;
            data_manager_.set_index_value(index_id_[1], index_j);
                
            if (where_clause(1)) {
                for (int index_2 = index_restart[2]; index_2 < upper_bound_[2]; ++index_2) {
                    index_values_[2] = index_2;
                    data_manager_.set_index_value(index_id_[2], index_2);
                        
            if (where_clause(2)) {
                for (int index_3 = index_restart[3]; index_3 < upper_bound_[3]; ++index_3) {
                    index_values_[3] = index_3;
                    data_manager_.set_index_value(index_id_[3], index_3);
                    
            if (where_clause(3)) {
                for (int index_4 = index_restart[4]; index_4 < upper_bound_[4]; ++index_4) {
                    index_values_[4] = index_4;
                    data_manager_.set_index_value(index_id_[4], index_4);
                    
            if (where_clause(4)) {
                for (int index_5 = index_restart[5]; index_5 < upper_bound_[5]; ++index_5) {
                    index_values_[5] = index_5;
                    data_manager_.set_index_value(index_id_[5], index_5);
                    
            if (where_clause(5)) {
            if (loop_count > iteration_) {
                iteration_++;
                if ((iteration_-1) % num_workers_ == company_rank_){
                    return true;
                }
            }
            ++loop_count;
        } // where 5
            } // index_5
            for (int i = 5; i < num_indices_; ++i) { index_restart[i] = lower_seg_[i]; }
        } // where 4
            } // index_4
            for (int i = 4; i < num_indices_; ++i) { index_restart[i] = lower_seg_[i]; }
        } // where 3
            } // index_3
            for (int i = 3; i < num_indices_; ++i) { index_restart[i] = lower_seg_[i]; }
        } // where 2
            } // index_2
            for (int i = 2; i < num_indices_; ++i) { index_restart[i] = lower_seg_[i]; }
        } // where 1
            } // index_j
            for (int i = 1; i < num_indices_; ++i) { index_restart[i] = lower_seg_[i]; }
        } // index_i
        
        return false; //this should be false here
    }

Fragment_ij_aoa_o_PardoLoopManager::Fragment_ij_aoa_o_PardoLoopManager(
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
    
    form_elst_dist();
    //form_rcut_dist();
    form_swao_frag();
    form_swocca_frag();
    //form_swvirta_frag();
}

Fragment_ij_aoa_o_PardoLoopManager::~Fragment_ij_aoa_o_PardoLoopManager() {}

bool Fragment_ij_aoa_o_PardoLoopManager::increment_simple_pair(){
    bool more = false;      // More iterations?
    int current_value;
    for (int i = 0; i < 2; ++i) {
        more = increment_single_index(i);
        if (more) {break;}
    } //if here, then all indices are at their max value
    return more;
}

bool Fragment_ij_aoa_o_PardoLoopManager::increment_single_index(int index){
    bool more = false; 	// More iterations?
    int current_value;
    current_value = data_manager_.index_value(index_id_[index]);
    ++current_value;
    if (current_value < upper_bound_[index]) {
        //increment current index and return
        index_values_[index] = current_value;
        data_manager_.set_index_value(index_id_[index], current_value);
        more = true;
    } else {
        //wrap around and handle next index
        index_values_[index] = lower_seg_[index];
        data_manager_.set_index_value(index_id_[index], lower_seg_[index]);
    }
    return more;
}
    
    
    inline bool Fragment_ij_aoa_o_PardoLoopManager::increment_all() {
        bool more = false;      // More iterations?
        int current_value;
        for (int i = 0; i < num_indices_; ++i) {
            more = increment_single_index(i);
            if (more) {break;}
        } //if here, then all indices are at their max value
        return more;
    }
    
    void Fragment_ij_aoa_o_PardoLoopManager::form_elst_dist() {
        int elst_dist_array_slot = sip_tables_.array_slot(std::string("elst_dist"));
        Block::BlockPtr bptr_elst_dist = data_manager_.contiguous_array_manager().get_array(elst_dist_array_slot);
        double *val_elst_dist = bptr_elst_dist->get_data();
        int elst_dist_size = bptr_elst_dist->size();
        int num_frag = upper_bound_[0] - lower_seg_[0];
        elst_dist.resize(num_frag);
        for (int i = 0; i < num_frag; ++i) {
            //std::cout << i << std::endl;
            elst_dist[i].resize(num_frag);
        }
        int ij = 0;
        for (int i = 0; i < num_frag; ++i) {
            for (int j = 0; j < num_frag; ++j) {
                elst_dist[j][i] = (int)val_elst_dist[ij];
                //std::cout << i << " " << j << " " << ij << " " << elst_dist[j][i] << " " << val_elst_dist[ij]<< std::endl;
                ++ij;
            }
        }
        return;
    }

    void Fragment_ij_aoa_o_PardoLoopManager::form_rcut_dist() {
        int rcut_dist_array_slot = sip_tables_.array_slot(std::string("rcut_dist"));
        Block::BlockPtr bptr_rcut_dist = data_manager_.contiguous_array_manager().get_array(rcut_dist_array_slot);
        double *val_rcut_dist = bptr_rcut_dist->get_data();
        int rcut_dist_size = bptr_rcut_dist->size();
        int num_frag = upper_bound_[0] - lower_seg_[0];
        rcut_dist.resize(num_frag);
        for (int i = 0; i < num_frag; ++i) {
            //std::cout << i << std::endl;
            rcut_dist[i].resize(num_frag);
        }
        int ij = 0;
        for (int i = 0; i < num_frag; ++i) {
            for (int j = 0; j < num_frag; ++j) {
                rcut_dist[j][i] = (int)val_rcut_dist[ij];
                //std::cout << i << " " << j << " " << ij << " " << rcut_dist[j][i] << " " << val_rcut_dist[ij]<< std::endl;
                ++ij;
            }
        }
        return;
    }
    
    void Fragment_ij_aoa_o_PardoLoopManager::form_swao_frag() {
        int swao_frag_array_slot = sip_tables_.array_slot(std::string("swao_frag"));
        Block::BlockPtr bptr_swao_frag = data_manager_.contiguous_array_manager().get_array(swao_frag_array_slot);
        double *val_swao_frag = bptr_swao_frag->get_data();
        int swao_frag_size = bptr_swao_frag->size();
        swao_frag.resize(swao_frag_size);
        for (int i = 0; i < swao_frag_size; ++i) {
            swao_frag[i] = (int)val_swao_frag[i];
            //std::cout << "vec " << swao_frag[i] << " data " << val_swao_frag[i] << std::endl;
        }
        return;
    }
    
    void Fragment_ij_aoa_o_PardoLoopManager::form_swocca_frag() {
        int swocca_frag_array_slot = sip_tables_.array_slot(std::string("swocca_frag"));
        Block::BlockPtr bptr_swocca_frag = data_manager_.contiguous_array_manager().get_array(swocca_frag_array_slot);
        double *val_swocca_frag = bptr_swocca_frag->get_data();
        int swocca_frag_size = bptr_swocca_frag->size();
        swocca_frag.resize(swocca_frag_size);
        for (int i = 0; i < swocca_frag_size; ++i) {
            swocca_frag[i] = (int)val_swocca_frag[i];
            //std::cout << "vec " << swocca_frag[i] << " data " << val_swocca_frag[i] << std::endl;
        }
        return;
    }
    
    void Fragment_ij_aoa_o_PardoLoopManager::form_swvirta_frag() {
        int swvirta_frag_array_slot = sip_tables_.array_slot(std::string("swvirta_frag"));
        Block::BlockPtr bptr_swvirta_frag = data_manager_.contiguous_array_manager().get_array(swvirta_frag_array_slot);
        double *val_swvirta_frag = bptr_swvirta_frag->get_data();
        int swvirta_frag_size = bptr_swvirta_frag->size();
        swvirta_frag.resize(swvirta_frag_size);
        for (int i = 0; i < swvirta_frag_size; ++i) {
            swvirta_frag[i] = (int)val_swvirta_frag[i];
            //std::cout << "vec " << swvirta_frag[i] << " data " << val_swvirta_frag[i] << std::endl;
        }
        return;
    }
    
    inline bool Fragment_ij_aoa_o_PardoLoopManager::initialize_indices() {
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
            index_values_[i] = lower_seg_[i];
            data_manager_.set_index_value(index_id_[i], lower_seg_[i]);
        }
        more_iterations = true;
        return more_iterations;
    }

/*
 the generic where clauses are:
 typ = 1: map AO to fragment
 typ = 2: map Occ to fragment
 typ = 3: map Virt to fragment
 typ = 4: map elst_dist(ifrag,jfrag) == ifrag
 typ = 5: map rcut_dist(ifrag,jfrag) == ifrag
 typ = 0: ifrag != jfrag
 */
bool Fragment_ij_aoa_o_PardoLoopManager::fragment_special_where_clause(int typ, int index, int frag) {
    bool where_clause;
    int ij = 0;
    switch (typ) {
            
        case 1: // check against ao index
            //where_clause = return_val_swao_frag(index_values_[index])     == index_values_[frag];
            where_clause = swao_frag[index_values_[index] - lower_seg_[index]]     == index_values_[frag];
            break;
            
        case 2: // check against occ index
            //where_clause = return_val_swocca_frag(index_values_[index])   == index_values_[frag];
            where_clause = swocca_frag[index_values_[index] - lower_seg_[index]]   == index_values_[frag];
            break;
            
        case 3: // check against virt index
            //where_clause = return_val_swvirta_frag(index_values_[index])  == index_values_[frag];
            where_clause = swvirta_frag[index_values_[index] - lower_seg_[index]]  == index_values_[frag];
            break;
            
        case 4: // check elst_dist[ifrag,jfrag] == ifrag
            //where_clause = return_val_elst_dist(index_values_[frag], index_values_[index]) == index_values_[frag];
            //ij = index_values_[frag] - 1 + (upper_bound_[0])*(index_values_[index] - 1);
            //std::cout << index_values_[index] << " " << index_values_[frag] << " " << ij << " " << elst_dist[ij] << std::endl;
            //where_clause = elst_dist[ij] == index_values_[frag];
            where_clause = elst_dist[index_values_[frag] - lower_seg_[frag]][index_values_[index] - lower_seg_[index]] == index_values_[frag];
            break;
            
        case 5: // check rcut_dist[ifrag,jfrag] == ifrag
            //where_clause = return_val_rcut_dist(index_values_[frag], index_values_[index]) == index_values_[frag];
            //ij = index_values_[index] - 1 + upper_bound_[0]*(index_values_[frag] - 1);
            //where_clause = rcut_dist[ij] == index_values_[frag];
            //break;
            where_clause = rcut_dist[index_values_[frag] - lower_seg_[frag]][index_values_[index] - lower_seg_[index]] == index_values_[frag];
            break;
            
        case 0:
            where_clause = index_values_[frag] != index_values_[index];
            break;
        default:
            where_clause = false;
    }
    //std::cout << "where_clause " << typ << " " << index_values_[index] << " " << index_values_[frag] << " " << where_clause << std::endl;
    return where_clause;
}


void Fragment_ij_aoa_o_PardoLoopManager::do_finalize() {
    for (int i = 0; i < num_indices_; ++i) {
        data_manager_.set_index_undefined(index_id_[i]);
    }
}

std::string Fragment_ij_aoa_o_PardoLoopManager::to_string() const {
    std::stringstream ss;
    ss << "Fragment_ij_aoa_o_PardoLoopManager:  num_indices="
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
                         const Fragment_ij_aoa_o_PardoLoopManager &obj) {
    os << obj.to_string();
    return os;
}

/*!
 -------------------------------------------
 _ij_av_oo_
 -------------------------------------------
 */
    
    /*
     for each new special fragment where clause pattern, this should be the only thing realy changed.
     see comment above for fragment_special_where_clause syntax
     */
    bool Fragment_ij_av_oo_PardoLoopManager::where_clause(int index) {
        bool where_;
        int ifrag = 0;
        int jfrag = 1;
        int ao = 1;
        int occ = 2;
        int virt = 3;
        int elst = 4;
        int rcut = 5;

        
        switch (index) {
            case 0:
            case 1:
                where_ = fragment_special_where_clause(elst,jfrag,ifrag);
                break;
            case 2:
                where_ = fragment_special_where_clause(ao,index,ifrag);
                break;
            case 3:
                where_ = fragment_special_where_clause(virt,index,ifrag);
                break;
            case 4:
                where_ = fragment_special_where_clause(occ,index,jfrag);
                break;
            case 5:
                where_ = fragment_special_where_clause(occ,index,jfrag);
                break;
            default:
                where_ = false;
        }
        return where_;
    }
    
    bool Fragment_ij_av_oo_PardoLoopManager::do_update() {
        if (to_exit_)
            return false;
        bool more_iters;
        bool where_clauses_value;
        
        interpreter_->skip_where_clauses(num_where_clauses_);
        
        if (first_time_) {
            first_time_ = false;
            more_iters = initialize_indices();
            
            where_clauses_value = true;
            for (int i = 1; i < num_indices_; ++i) {
                where_clauses_value = where_clauses_value && where_clause(i);
            }
            if(where_clauses_value){
                iteration_++;
                if ((iteration_-1) % num_workers_ == company_rank_){
                    return true;
                }
            }
        }
        
        int index_restart[MAX_RANK];
        
        for (int i = 0; i < num_indices_; ++i) {
            index_restart[i] = index_values_[i];
        }
        
        int loop_count = iteration_;
        for (int index_i = index_restart[0]; index_i < upper_bound_[0]; ++index_i) {
            index_values_[0] = index_i;
            data_manager_.set_index_value(index_id_[0], index_i);
            
        for (int index_j = index_restart[1]; index_j < upper_bound_[1]; ++index_j) {
            index_values_[1] = index_j;
            data_manager_.set_index_value(index_id_[1], index_j);
                
            if (where_clause(1)) {
                for (int index_2 = index_restart[2]; index_2 < upper_bound_[2]; ++index_2) {
                    index_values_[2] = index_2;
                    data_manager_.set_index_value(index_id_[2], index_2);
                        
            if (where_clause(2)) {
                for (int index_3 = index_restart[3]; index_3 < upper_bound_[3]; ++index_3) {
                    index_values_[3] = index_3;
                    data_manager_.set_index_value(index_id_[3], index_3);
                    
            if (where_clause(3)) {
                for (int index_4 = index_restart[4]; index_4 < upper_bound_[4]; ++index_4) {
                    index_values_[4] = index_4;
                    data_manager_.set_index_value(index_id_[4], index_4);
                    
            if (where_clause(4)) {
                for (int index_5 = index_restart[5]; index_5 < upper_bound_[5]; ++index_5) {
                    index_values_[5] = index_5;
                    data_manager_.set_index_value(index_id_[5], index_5);
                    
            if (where_clause(5)) {
            if (loop_count > iteration_) {
                iteration_++;
                if ((iteration_-1) % num_workers_ == company_rank_){
                    return true;
                }
            }
            ++loop_count;
        } // where 5
            } // index_5
            for (int i = 5; i < num_indices_; ++i) { index_restart[i] = lower_seg_[i]; }
        } // where 4
            } // index_4
            for (int i = 4; i < num_indices_; ++i) { index_restart[i] = lower_seg_[i]; }
        } // where 3
            } // index_3
            for (int i = 3; i < num_indices_; ++i) { index_restart[i] = lower_seg_[i]; }
        } // where 2
            } // index_2
            for (int i = 2; i < num_indices_; ++i) { index_restart[i] = lower_seg_[i]; }
        } // where 1
            } // index_j
            for (int i = 1; i < num_indices_; ++i) { index_restart[i] = lower_seg_[i]; }
        } // index_i
        
        return false; //this should be false here
    }

Fragment_ij_av_oo_PardoLoopManager::Fragment_ij_av_oo_PardoLoopManager(
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
    
    form_elst_dist();
    //form_rcut_dist();
    form_swao_frag();
    form_swocca_frag();
    form_swvirta_frag();
}

Fragment_ij_av_oo_PardoLoopManager::~Fragment_ij_av_oo_PardoLoopManager() {}

bool Fragment_ij_av_oo_PardoLoopManager::increment_simple_pair(){
    bool more = false;      // More iterations?
    int current_value;
    for (int i = 0; i < 2; ++i) {
        more = increment_single_index(i);
        if (more) {break;}
    } //if here, then all indices are at their max value
    return more;
}

bool Fragment_ij_av_oo_PardoLoopManager::increment_single_index(int index){
    bool more = false; 	// More iterations?
    int current_value;
    current_value = data_manager_.index_value(index_id_[index]);
    ++current_value;
    if (current_value < upper_bound_[index]) {
        //increment current index and return
        index_values_[index] = current_value;
        data_manager_.set_index_value(index_id_[index], current_value);
        more = true;
    } else {
        //wrap around and handle next index
        index_values_[index] = lower_seg_[index];
        data_manager_.set_index_value(index_id_[index], lower_seg_[index]);
    }
    return more;
}
    
    
    inline bool Fragment_ij_av_oo_PardoLoopManager::increment_all() {
        bool more = false;      // More iterations?
        int current_value;
        for (int i = 0; i < num_indices_; ++i) {
            more = increment_single_index(i);
            if (more) {break;}
        } //if here, then all indices are at their max value
        return more;
    }
    
    void Fragment_ij_av_oo_PardoLoopManager::form_elst_dist() {
        int elst_dist_array_slot = sip_tables_.array_slot(std::string("elst_dist"));
        Block::BlockPtr bptr_elst_dist = data_manager_.contiguous_array_manager().get_array(elst_dist_array_slot);
        double *val_elst_dist = bptr_elst_dist->get_data();
        int elst_dist_size = bptr_elst_dist->size();
        int num_frag = upper_bound_[0] - lower_seg_[0];
        elst_dist.resize(num_frag);
        for (int i = 0; i < num_frag; ++i) {
            //std::cout << i << std::endl;
            elst_dist[i].resize(num_frag);
        }
        int ij = 0;
        for (int i = 0; i < num_frag; ++i) {
            for (int j = 0; j < num_frag; ++j) {
                elst_dist[j][i] = (int)val_elst_dist[ij];
                //std::cout << i << " " << j << " " << ij << " " << elst_dist[j][i] << " " << val_elst_dist[ij]<< std::endl;
                ++ij;
            }
        }
        return;
    }

    void Fragment_ij_av_oo_PardoLoopManager::form_rcut_dist() {
        int rcut_dist_array_slot = sip_tables_.array_slot(std::string("rcut_dist"));
        Block::BlockPtr bptr_rcut_dist = data_manager_.contiguous_array_manager().get_array(rcut_dist_array_slot);
        double *val_rcut_dist = bptr_rcut_dist->get_data();
        int rcut_dist_size = bptr_rcut_dist->size();
        int num_frag = upper_bound_[0] - lower_seg_[0];
        rcut_dist.resize(num_frag);
        for (int i = 0; i < num_frag; ++i) {
            //std::cout << i << std::endl;
            rcut_dist[i].resize(num_frag);
        }
        int ij = 0;
        for (int i = 0; i < num_frag; ++i) {
            for (int j = 0; j < num_frag; ++j) {
                rcut_dist[j][i] = (int)val_rcut_dist[ij];
                //std::cout << i << " " << j << " " << ij << " " << rcut_dist[j][i] << " " << val_rcut_dist[ij]<< std::endl;
                ++ij;
            }
        }
        return;
    }
    
    void Fragment_ij_av_oo_PardoLoopManager::form_swao_frag() {
        int swao_frag_array_slot = sip_tables_.array_slot(std::string("swao_frag"));
        Block::BlockPtr bptr_swao_frag = data_manager_.contiguous_array_manager().get_array(swao_frag_array_slot);
        double *val_swao_frag = bptr_swao_frag->get_data();
        int swao_frag_size = bptr_swao_frag->size();
        swao_frag.resize(swao_frag_size);
        for (int i = 0; i < swao_frag_size; ++i) {
            swao_frag[i] = (int)val_swao_frag[i];
            //std::cout << "vec " << swao_frag[i] << " data " << val_swao_frag[i] << std::endl;
        }
        return;
    }
    
    void Fragment_ij_av_oo_PardoLoopManager::form_swocca_frag() {
        int swocca_frag_array_slot = sip_tables_.array_slot(std::string("swocca_frag"));
        Block::BlockPtr bptr_swocca_frag = data_manager_.contiguous_array_manager().get_array(swocca_frag_array_slot);
        double *val_swocca_frag = bptr_swocca_frag->get_data();
        int swocca_frag_size = bptr_swocca_frag->size();
        swocca_frag.resize(swocca_frag_size);
        for (int i = 0; i < swocca_frag_size; ++i) {
            swocca_frag[i] = (int)val_swocca_frag[i];
            //std::cout << "vec " << swocca_frag[i] << " data " << val_swocca_frag[i] << std::endl;
        }
        return;
    }
    
    void Fragment_ij_av_oo_PardoLoopManager::form_swvirta_frag() {
        int swvirta_frag_array_slot = sip_tables_.array_slot(std::string("swvirta_frag"));
        Block::BlockPtr bptr_swvirta_frag = data_manager_.contiguous_array_manager().get_array(swvirta_frag_array_slot);
        double *val_swvirta_frag = bptr_swvirta_frag->get_data();
        int swvirta_frag_size = bptr_swvirta_frag->size();
        swvirta_frag.resize(swvirta_frag_size);
        for (int i = 0; i < swvirta_frag_size; ++i) {
            swvirta_frag[i] = (int)val_swvirta_frag[i];
            //std::cout << "vec " << swvirta_frag[i] << " data " << val_swvirta_frag[i] << std::endl;
        }
        return;
    }
    
    inline bool Fragment_ij_av_oo_PardoLoopManager::initialize_indices() {
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
            index_values_[i] = lower_seg_[i];
            data_manager_.set_index_value(index_id_[i], lower_seg_[i]);
        }
        more_iterations = true;
        return more_iterations;
    }

/*
 the generic where clauses are:
 typ = 1: map AO to fragment
 typ = 2: map Occ to fragment
 typ = 3: map Virt to fragment
 typ = 4: map elst_dist(ifrag,jfrag) == ifrag
 typ = 5: map rcut_dist(ifrag,jfrag) == ifrag
 typ = 0: ifrag != jfrag
 */
bool Fragment_ij_av_oo_PardoLoopManager::fragment_special_where_clause(int typ, int index, int frag) {
    bool where_clause;
    int ij = 0;
    switch (typ) {
            
        case 1: // check against ao index
            //where_clause = return_val_swao_frag(index_values_[index])     == index_values_[frag];
            where_clause = swao_frag[index_values_[index] - lower_seg_[index]]     == index_values_[frag];
            break;
            
        case 2: // check against occ index
            //where_clause = return_val_swocca_frag(index_values_[index])   == index_values_[frag];
            where_clause = swocca_frag[index_values_[index] - lower_seg_[index]]   == index_values_[frag];
            break;
            
        case 3: // check against virt index
            //where_clause = return_val_swvirta_frag(index_values_[index])  == index_values_[frag];
            where_clause = swvirta_frag[index_values_[index] - lower_seg_[index]]  == index_values_[frag];
            break;
            
        case 4: // check elst_dist[ifrag,jfrag] == ifrag
            //where_clause = return_val_elst_dist(index_values_[frag], index_values_[index]) == index_values_[frag];
            //ij = index_values_[frag] - 1 + (upper_bound_[0])*(index_values_[index] - 1);
            //std::cout << index_values_[index] << " " << index_values_[frag] << " " << ij << " " << elst_dist[ij] << std::endl;
            //where_clause = elst_dist[ij] == index_values_[frag];
            where_clause = elst_dist[index_values_[frag] - lower_seg_[frag]][index_values_[index] - lower_seg_[index]] == index_values_[frag];
            break;
            
        case 5: // check rcut_dist[ifrag,jfrag] == ifrag
            //where_clause = return_val_rcut_dist(index_values_[frag], index_values_[index]) == index_values_[frag];
            //ij = index_values_[index] - 1 + upper_bound_[0]*(index_values_[frag] - 1);
            //where_clause = rcut_dist[ij] == index_values_[frag];
            //break;
            where_clause = rcut_dist[index_values_[frag] - lower_seg_[frag]][index_values_[index] - lower_seg_[index]] == index_values_[frag];
            break;
            
        case 0:
            where_clause = index_values_[frag] != index_values_[index];
            break;
        default:
            where_clause = false;
    }
    //std::cout << "where_clause " << typ << " " << index_values_[index] << " " << index_values_[frag] << " " << where_clause << std::endl;
    return where_clause;
}


void Fragment_ij_av_oo_PardoLoopManager::do_finalize() {
    for (int i = 0; i < num_indices_; ++i) {
        data_manager_.set_index_undefined(index_id_[i]);
    }
}

std::string Fragment_ij_av_oo_PardoLoopManager::to_string() const {
    std::stringstream ss;
    ss << "Fragment_ij_av_oo_PardoLoopManager:  num_indices="
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
                         const Fragment_ij_av_oo_PardoLoopManager &obj) {
    os << obj.to_string();
    return os;
}
    

/*!
 -------------------------------------------
 _ij_ao_oo_
 -------------------------------------------
 */
    
    /*
     for each new special fragment where clause pattern, this should be the only thing realy changed.
     see comment above for fragment_special_where_clause syntax
     */
    bool Fragment_ij_ao_oo_PardoLoopManager::where_clause(int index) {
        bool where_;
        int ifrag = 0;
        int jfrag = 1;
        int ao = 1;
        int occ = 2;
        int virt = 3;
        int elst = 4;
        int rcut = 5;

        
        switch (index) {
            case 0:
            case 1:
                where_ = fragment_special_where_clause(elst,jfrag,ifrag);
                break;
            case 2:
                where_ = fragment_special_where_clause(ao,index,ifrag);
                break;
            case 3:
                where_ = fragment_special_where_clause(occ,index,ifrag);
                break;
            case 4:
                where_ = fragment_special_where_clause(occ,index,jfrag);
                break;
            case 5:
                where_ = fragment_special_where_clause(occ,index,jfrag);
                break;
            default:
                where_ = false;
        }
        return where_;
    }
    
    bool Fragment_ij_ao_oo_PardoLoopManager::do_update() {
        if (to_exit_)
            return false;
        bool more_iters;
        bool where_clauses_value;
        
        interpreter_->skip_where_clauses(num_where_clauses_);
        
        if (first_time_) {
            first_time_ = false;
            more_iters = initialize_indices();
            
            where_clauses_value = true;
            for (int i = 1; i < num_indices_; ++i) {
                where_clauses_value = where_clauses_value && where_clause(i);
            }
            if(where_clauses_value){
                iteration_++;
                if ((iteration_-1) % num_workers_ == company_rank_){
                    return true;
                }
            }
        }
        
        int index_restart[MAX_RANK];
        
        for (int i = 0; i < num_indices_; ++i) {
            index_restart[i] = index_values_[i];
        }
        
        int loop_count = iteration_;
        for (int index_i = index_restart[0]; index_i < upper_bound_[0]; ++index_i) {
            index_values_[0] = index_i;
            data_manager_.set_index_value(index_id_[0], index_i);
            
        for (int index_j = index_restart[1]; index_j < upper_bound_[1]; ++index_j) {
            index_values_[1] = index_j;
            data_manager_.set_index_value(index_id_[1], index_j);
                
            if (where_clause(1)) {
                for (int index_2 = index_restart[2]; index_2 < upper_bound_[2]; ++index_2) {
                    index_values_[2] = index_2;
                    data_manager_.set_index_value(index_id_[2], index_2);
                        
            if (where_clause(2)) {
                for (int index_3 = index_restart[3]; index_3 < upper_bound_[3]; ++index_3) {
                    index_values_[3] = index_3;
                    data_manager_.set_index_value(index_id_[3], index_3);
                    
            if (where_clause(3)) {
                for (int index_4 = index_restart[4]; index_4 < upper_bound_[4]; ++index_4) {
                    index_values_[4] = index_4;
                    data_manager_.set_index_value(index_id_[4], index_4);
                    
            if (where_clause(4)) {
                for (int index_5 = index_restart[5]; index_5 < upper_bound_[5]; ++index_5) {
                    index_values_[5] = index_5;
                    data_manager_.set_index_value(index_id_[5], index_5);
                    
            if (where_clause(5)) {
            if (loop_count > iteration_) {
                iteration_++;
                if ((iteration_-1) % num_workers_ == company_rank_){
                    return true;
                }
            }
            ++loop_count;
        } // where 5
            } // index_5
            for (int i = 5; i < num_indices_; ++i) { index_restart[i] = lower_seg_[i]; }
        } // where 4
            } // index_4
            for (int i = 4; i < num_indices_; ++i) { index_restart[i] = lower_seg_[i]; }
        } // where 3
            } // index_3
            for (int i = 3; i < num_indices_; ++i) { index_restart[i] = lower_seg_[i]; }
        } // where 2
            } // index_2
            for (int i = 2; i < num_indices_; ++i) { index_restart[i] = lower_seg_[i]; }
        } // where 1
            } // index_j
            for (int i = 1; i < num_indices_; ++i) { index_restart[i] = lower_seg_[i]; }
        } // index_i
        
        return false; //this should be false here
    }

Fragment_ij_ao_oo_PardoLoopManager::Fragment_ij_ao_oo_PardoLoopManager(
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
    
    form_elst_dist();
    //form_rcut_dist();
    form_swao_frag();
    form_swocca_frag();
    //form_swvirta_frag();
}

Fragment_ij_ao_oo_PardoLoopManager::~Fragment_ij_ao_oo_PardoLoopManager() {}

bool Fragment_ij_ao_oo_PardoLoopManager::increment_simple_pair(){
    bool more = false;      // More iterations?
    int current_value;
    for (int i = 0; i < 2; ++i) {
        more = increment_single_index(i);
        if (more) {break;}
    } //if here, then all indices are at their max value
    return more;
}

bool Fragment_ij_ao_oo_PardoLoopManager::increment_single_index(int index){
    bool more = false; 	// More iterations?
    int current_value;
    current_value = data_manager_.index_value(index_id_[index]);
    ++current_value;
    if (current_value < upper_bound_[index]) {
        //increment current index and return
        index_values_[index] = current_value;
        data_manager_.set_index_value(index_id_[index], current_value);
        more = true;
    } else {
        //wrap around and handle next index
        index_values_[index] = lower_seg_[index];
        data_manager_.set_index_value(index_id_[index], lower_seg_[index]);
    }
    return more;
}
    
    
    inline bool Fragment_ij_ao_oo_PardoLoopManager::increment_all() {
        bool more = false;      // More iterations?
        int current_value;
        for (int i = 0; i < num_indices_; ++i) {
            more = increment_single_index(i);
            if (more) {break;}
        } //if here, then all indices are at their max value
        return more;
    }
    
    void Fragment_ij_ao_oo_PardoLoopManager::form_elst_dist() {
        int elst_dist_array_slot = sip_tables_.array_slot(std::string("elst_dist"));
        Block::BlockPtr bptr_elst_dist = data_manager_.contiguous_array_manager().get_array(elst_dist_array_slot);
        double *val_elst_dist = bptr_elst_dist->get_data();
        int elst_dist_size = bptr_elst_dist->size();
        int num_frag = upper_bound_[0] - lower_seg_[0];
        elst_dist.resize(num_frag);
        for (int i = 0; i < num_frag; ++i) {
            //std::cout << i << std::endl;
            elst_dist[i].resize(num_frag);
        }
        int ij = 0;
        for (int i = 0; i < num_frag; ++i) {
            for (int j = 0; j < num_frag; ++j) {
                elst_dist[j][i] = (int)val_elst_dist[ij];
                //std::cout << i << " " << j << " " << ij << " " << elst_dist[j][i] << " " << val_elst_dist[ij]<< std::endl;
                ++ij;
            }
        }
        return;
    }

    void Fragment_ij_ao_oo_PardoLoopManager::form_rcut_dist() {
        int rcut_dist_array_slot = sip_tables_.array_slot(std::string("rcut_dist"));
        Block::BlockPtr bptr_rcut_dist = data_manager_.contiguous_array_manager().get_array(rcut_dist_array_slot);
        double *val_rcut_dist = bptr_rcut_dist->get_data();
        int rcut_dist_size = bptr_rcut_dist->size();
        int num_frag = upper_bound_[0] - lower_seg_[0];
        rcut_dist.resize(num_frag);
        for (int i = 0; i < num_frag; ++i) {
            //std::cout << i << std::endl;
            rcut_dist[i].resize(num_frag);
        }
        int ij = 0;
        for (int i = 0; i < num_frag; ++i) {
            for (int j = 0; j < num_frag; ++j) {
                rcut_dist[j][i] = (int)val_rcut_dist[ij];
                //std::cout << i << " " << j << " " << ij << " " << rcut_dist[j][i] << " " << val_rcut_dist[ij]<< std::endl;
                ++ij;
            }
        }
        return;
    }
    
    void Fragment_ij_ao_oo_PardoLoopManager::form_swao_frag() {
        int swao_frag_array_slot = sip_tables_.array_slot(std::string("swao_frag"));
        Block::BlockPtr bptr_swao_frag = data_manager_.contiguous_array_manager().get_array(swao_frag_array_slot);
        double *val_swao_frag = bptr_swao_frag->get_data();
        int swao_frag_size = bptr_swao_frag->size();
        swao_frag.resize(swao_frag_size);
        for (int i = 0; i < swao_frag_size; ++i) {
            swao_frag[i] = (int)val_swao_frag[i];
            //std::cout << "vec " << swao_frag[i] << " data " << val_swao_frag[i] << std::endl;
        }
        return;
    }
    
    void Fragment_ij_ao_oo_PardoLoopManager::form_swocca_frag() {
        int swocca_frag_array_slot = sip_tables_.array_slot(std::string("swocca_frag"));
        Block::BlockPtr bptr_swocca_frag = data_manager_.contiguous_array_manager().get_array(swocca_frag_array_slot);
        double *val_swocca_frag = bptr_swocca_frag->get_data();
        int swocca_frag_size = bptr_swocca_frag->size();
        swocca_frag.resize(swocca_frag_size);
        for (int i = 0; i < swocca_frag_size; ++i) {
            swocca_frag[i] = (int)val_swocca_frag[i];
            //std::cout << "vec " << swocca_frag[i] << " data " << val_swocca_frag[i] << std::endl;
        }
        return;
    }
    
    void Fragment_ij_ao_oo_PardoLoopManager::form_swvirta_frag() {
        int swvirta_frag_array_slot = sip_tables_.array_slot(std::string("swvirta_frag"));
        Block::BlockPtr bptr_swvirta_frag = data_manager_.contiguous_array_manager().get_array(swvirta_frag_array_slot);
        double *val_swvirta_frag = bptr_swvirta_frag->get_data();
        int swvirta_frag_size = bptr_swvirta_frag->size();
        swvirta_frag.resize(swvirta_frag_size);
        for (int i = 0; i < swvirta_frag_size; ++i) {
            swvirta_frag[i] = (int)val_swvirta_frag[i];
            //std::cout << "vec " << swvirta_frag[i] << " data " << val_swvirta_frag[i] << std::endl;
        }
        return;
    }
    
    inline bool Fragment_ij_ao_oo_PardoLoopManager::initialize_indices() {
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
            index_values_[i] = lower_seg_[i];
            data_manager_.set_index_value(index_id_[i], lower_seg_[i]);
        }
        more_iterations = true;
        return more_iterations;
    }

/*
 the generic where clauses are:
 typ = 1: map AO to fragment
 typ = 2: map Occ to fragment
 typ = 3: map Virt to fragment
 typ = 4: map elst_dist(ifrag,jfrag) == ifrag
 typ = 5: map rcut_dist(ifrag,jfrag) == ifrag
 typ = 0: ifrag != jfrag
 */
bool Fragment_ij_ao_oo_PardoLoopManager::fragment_special_where_clause(int typ, int index, int frag) {
    bool where_clause;
    int ij = 0;
    switch (typ) {
            
        case 1: // check against ao index
            //where_clause = return_val_swao_frag(index_values_[index])     == index_values_[frag];
            where_clause = swao_frag[index_values_[index] - lower_seg_[index]]     == index_values_[frag];
            break;
            
        case 2: // check against occ index
            //where_clause = return_val_swocca_frag(index_values_[index])   == index_values_[frag];
            where_clause = swocca_frag[index_values_[index] - lower_seg_[index]]   == index_values_[frag];
            break;
            
        case 3: // check against virt index
            //where_clause = return_val_swvirta_frag(index_values_[index])  == index_values_[frag];
            where_clause = swvirta_frag[index_values_[index] - lower_seg_[index]]  == index_values_[frag];
            break;
            
        case 4: // check elst_dist[ifrag,jfrag] == ifrag
            //where_clause = return_val_elst_dist(index_values_[frag], index_values_[index]) == index_values_[frag];
            //ij = index_values_[frag] - 1 + (upper_bound_[0])*(index_values_[index] - 1);
            //std::cout << index_values_[index] << " " << index_values_[frag] << " " << ij << " " << elst_dist[ij] << std::endl;
            //where_clause = elst_dist[ij] == index_values_[frag];
            where_clause = elst_dist[index_values_[frag] - lower_seg_[frag]][index_values_[index] - lower_seg_[index]] == index_values_[frag];
            break;
            
        case 5: // check rcut_dist[ifrag,jfrag] == ifrag
            //where_clause = return_val_rcut_dist(index_values_[frag], index_values_[index]) == index_values_[frag];
            //ij = index_values_[index] - 1 + upper_bound_[0]*(index_values_[frag] - 1);
            //where_clause = rcut_dist[ij] == index_values_[frag];
            //break;
            where_clause = rcut_dist[index_values_[frag] - lower_seg_[frag]][index_values_[index] - lower_seg_[index]] == index_values_[frag];
            break;
            
        case 0:
            where_clause = index_values_[frag] != index_values_[index];
            break;
        default:
            where_clause = false;
    }
    //std::cout << "where_clause " << typ << " " << index_values_[index] << " " << index_values_[frag] << " " << where_clause << std::endl;
    return where_clause;
}


void Fragment_ij_ao_oo_PardoLoopManager::do_finalize() {
    for (int i = 0; i < num_indices_; ++i) {
        data_manager_.set_index_undefined(index_id_[i]);
    }
}

std::string Fragment_ij_ao_oo_PardoLoopManager::to_string() const {
    std::stringstream ss;
    ss << "Fragment_ij_ao_oo_PardoLoopManager:  num_indices="
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
                         const Fragment_ij_ao_oo_PardoLoopManager &obj) {
    os << obj.to_string();
    return os;
}
    
/*!
 -------------------------------------------
 _ij_oo_ao_
 -------------------------------------------
 */
    
    /*
     for each new special fragment where clause pattern, this should be the only thing realy changed.
     see comment above for fragment_special_where_clause syntax
     */
    bool Fragment_ij_oo_ao_PardoLoopManager::where_clause(int index) {
        bool where_;
        int ifrag = 0;
        int jfrag = 1;
        int ao = 1;
        int occ = 2;
        int virt = 3;
        int elst = 4;
        int rcut = 5;

        
        switch (index) {
            case 0:
            case 1:
                where_ = fragment_special_where_clause(elst,jfrag,ifrag);
                break;
            case 2:
                where_ = fragment_special_where_clause(occ,index,ifrag);
                break;
            case 3:
                where_ = fragment_special_where_clause(occ,index,ifrag);
                break;
            case 4:
                where_ = fragment_special_where_clause(ao,index,jfrag);
                break;
            case 5:
                where_ = fragment_special_where_clause(occ,index,jfrag);
                break;
            default:
                where_ = false;
        }
        return where_;
    }
    
    bool Fragment_ij_oo_ao_PardoLoopManager::do_update() {
        if (to_exit_)
            return false;
        bool more_iters;
        bool where_clauses_value;
        
        interpreter_->skip_where_clauses(num_where_clauses_);
        
        if (first_time_) {
            first_time_ = false;
            more_iters = initialize_indices();
            
            where_clauses_value = true;
            for (int i = 1; i < num_indices_; ++i) {
                where_clauses_value = where_clauses_value && where_clause(i);
            }
            if(where_clauses_value){
                iteration_++;
                if ((iteration_-1) % num_workers_ == company_rank_){
                    return true;
                }
            }
        }
        
        int index_restart[MAX_RANK];
        
        for (int i = 0; i < num_indices_; ++i) {
            index_restart[i] = index_values_[i];
        }
        
        int loop_count = iteration_;
        for (int index_i = index_restart[0]; index_i < upper_bound_[0]; ++index_i) {
            index_values_[0] = index_i;
            data_manager_.set_index_value(index_id_[0], index_i);
            
        for (int index_j = index_restart[1]; index_j < upper_bound_[1]; ++index_j) {
            index_values_[1] = index_j;
            data_manager_.set_index_value(index_id_[1], index_j);
                
            if (where_clause(1)) {
                for (int index_2 = index_restart[2]; index_2 < upper_bound_[2]; ++index_2) {
                    index_values_[2] = index_2;
                    data_manager_.set_index_value(index_id_[2], index_2);
                        
            if (where_clause(2)) {
                for (int index_3 = index_restart[3]; index_3 < upper_bound_[3]; ++index_3) {
                    index_values_[3] = index_3;
                    data_manager_.set_index_value(index_id_[3], index_3);
                    
            if (where_clause(3)) {
                for (int index_4 = index_restart[4]; index_4 < upper_bound_[4]; ++index_4) {
                    index_values_[4] = index_4;
                    data_manager_.set_index_value(index_id_[4], index_4);
                    
            if (where_clause(4)) {
                for (int index_5 = index_restart[5]; index_5 < upper_bound_[5]; ++index_5) {
                    index_values_[5] = index_5;
                    data_manager_.set_index_value(index_id_[5], index_5);
                    
            if (where_clause(5)) {
            if (loop_count > iteration_) {
                iteration_++;
                if ((iteration_-1) % num_workers_ == company_rank_){
                    return true;
                }
            }
            ++loop_count;
        } // where 5
            } // index_5
            for (int i = 5; i < num_indices_; ++i) { index_restart[i] = lower_seg_[i]; }
        } // where 4
            } // index_4
            for (int i = 4; i < num_indices_; ++i) { index_restart[i] = lower_seg_[i]; }
        } // where 3
            } // index_3
            for (int i = 3; i < num_indices_; ++i) { index_restart[i] = lower_seg_[i]; }
        } // where 2
            } // index_2
            for (int i = 2; i < num_indices_; ++i) { index_restart[i] = lower_seg_[i]; }
        } // where 1
            } // index_j
            for (int i = 1; i < num_indices_; ++i) { index_restart[i] = lower_seg_[i]; }
        } // index_i
        
        return false; //this should be false here
    }

Fragment_ij_oo_ao_PardoLoopManager::Fragment_ij_oo_ao_PardoLoopManager(
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
    
    form_elst_dist();
    //form_rcut_dist();
    form_swao_frag();
    form_swocca_frag();
    //form_swvirta_frag();
}

Fragment_ij_oo_ao_PardoLoopManager::~Fragment_ij_oo_ao_PardoLoopManager() {}

bool Fragment_ij_oo_ao_PardoLoopManager::increment_simple_pair(){
    bool more = false;      // More iterations?
    int current_value;
    for (int i = 0; i < 2; ++i) {
        more = increment_single_index(i);
        if (more) {break;}
    } //if here, then all indices are at their max value
    return more;
}

bool Fragment_ij_oo_ao_PardoLoopManager::increment_single_index(int index){
    bool more = false; 	// More iterations?
    int current_value;
    current_value = data_manager_.index_value(index_id_[index]);
    ++current_value;
    if (current_value < upper_bound_[index]) {
        //increment current index and return
        index_values_[index] = current_value;
        data_manager_.set_index_value(index_id_[index], current_value);
        more = true;
    } else {
        //wrap around and handle next index
        index_values_[index] = lower_seg_[index];
        data_manager_.set_index_value(index_id_[index], lower_seg_[index]);
    }
    return more;
}
    
    
    inline bool Fragment_ij_oo_ao_PardoLoopManager::increment_all() {
        bool more = false;      // More iterations?
        int current_value;
        for (int i = 0; i < num_indices_; ++i) {
            more = increment_single_index(i);
            if (more) {break;}
        } //if here, then all indices are at their max value
        return more;
    }
    
    void Fragment_ij_oo_ao_PardoLoopManager::form_elst_dist() {
        int elst_dist_array_slot = sip_tables_.array_slot(std::string("elst_dist"));
        Block::BlockPtr bptr_elst_dist = data_manager_.contiguous_array_manager().get_array(elst_dist_array_slot);
        double *val_elst_dist = bptr_elst_dist->get_data();
        int elst_dist_size = bptr_elst_dist->size();
        int num_frag = upper_bound_[0] - lower_seg_[0];
        elst_dist.resize(num_frag);
        for (int i = 0; i < num_frag; ++i) {
            //std::cout << i << std::endl;
            elst_dist[i].resize(num_frag);
        }
        int ij = 0;
        for (int i = 0; i < num_frag; ++i) {
            for (int j = 0; j < num_frag; ++j) {
                elst_dist[j][i] = (int)val_elst_dist[ij];
                //std::cout << i << " " << j << " " << ij << " " << elst_dist[j][i] << " " << val_elst_dist[ij]<< std::endl;
                ++ij;
            }
        }
        return;
    }

    void Fragment_ij_oo_ao_PardoLoopManager::form_rcut_dist() {
        int rcut_dist_array_slot = sip_tables_.array_slot(std::string("rcut_dist"));
        Block::BlockPtr bptr_rcut_dist = data_manager_.contiguous_array_manager().get_array(rcut_dist_array_slot);
        double *val_rcut_dist = bptr_rcut_dist->get_data();
        int rcut_dist_size = bptr_rcut_dist->size();
        int num_frag = upper_bound_[0] - lower_seg_[0];
        rcut_dist.resize(num_frag);
        for (int i = 0; i < num_frag; ++i) {
            //std::cout << i << std::endl;
            rcut_dist[i].resize(num_frag);
        }
        int ij = 0;
        for (int i = 0; i < num_frag; ++i) {
            for (int j = 0; j < num_frag; ++j) {
                rcut_dist[j][i] = (int)val_rcut_dist[ij];
                //std::cout << i << " " << j << " " << ij << " " << rcut_dist[j][i] << " " << val_rcut_dist[ij]<< std::endl;
                ++ij;
            }
        }
        return;
    }
    
    void Fragment_ij_oo_ao_PardoLoopManager::form_swao_frag() {
        int swao_frag_array_slot = sip_tables_.array_slot(std::string("swao_frag"));
        Block::BlockPtr bptr_swao_frag = data_manager_.contiguous_array_manager().get_array(swao_frag_array_slot);
        double *val_swao_frag = bptr_swao_frag->get_data();
        int swao_frag_size = bptr_swao_frag->size();
        swao_frag.resize(swao_frag_size);
        for (int i = 0; i < swao_frag_size; ++i) {
            swao_frag[i] = (int)val_swao_frag[i];
            //std::cout << "vec " << swao_frag[i] << " data " << val_swao_frag[i] << std::endl;
        }
        return;
    }
    
    void Fragment_ij_oo_ao_PardoLoopManager::form_swocca_frag() {
        int swocca_frag_array_slot = sip_tables_.array_slot(std::string("swocca_frag"));
        Block::BlockPtr bptr_swocca_frag = data_manager_.contiguous_array_manager().get_array(swocca_frag_array_slot);
        double *val_swocca_frag = bptr_swocca_frag->get_data();
        int swocca_frag_size = bptr_swocca_frag->size();
        swocca_frag.resize(swocca_frag_size);
        for (int i = 0; i < swocca_frag_size; ++i) {
            swocca_frag[i] = (int)val_swocca_frag[i];
            //std::cout << "vec " << swocca_frag[i] << " data " << val_swocca_frag[i] << std::endl;
        }
        return;
    }
    
    void Fragment_ij_oo_ao_PardoLoopManager::form_swvirta_frag() {
        int swvirta_frag_array_slot = sip_tables_.array_slot(std::string("swvirta_frag"));
        Block::BlockPtr bptr_swvirta_frag = data_manager_.contiguous_array_manager().get_array(swvirta_frag_array_slot);
        double *val_swvirta_frag = bptr_swvirta_frag->get_data();
        int swvirta_frag_size = bptr_swvirta_frag->size();
        swvirta_frag.resize(swvirta_frag_size);
        for (int i = 0; i < swvirta_frag_size; ++i) {
            swvirta_frag[i] = (int)val_swvirta_frag[i];
            //std::cout << "vec " << swvirta_frag[i] << " data " << val_swvirta_frag[i] << std::endl;
        }
        return;
    }
    
    inline bool Fragment_ij_oo_ao_PardoLoopManager::initialize_indices() {
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
            index_values_[i] = lower_seg_[i];
            data_manager_.set_index_value(index_id_[i], lower_seg_[i]);
        }
        more_iterations = true;
        return more_iterations;
    }

/*
 the generic where clauses are:
 typ = 1: map AO to fragment
 typ = 2: map Occ to fragment
 typ = 3: map Virt to fragment
 typ = 4: map elst_dist(ifrag,jfrag) == ifrag
 typ = 5: map rcut_dist(ifrag,jfrag) == ifrag
 typ = 0: ifrag != jfrag
 */
bool Fragment_ij_oo_ao_PardoLoopManager::fragment_special_where_clause(int typ, int index, int frag) {
    bool where_clause;
    int ij = 0;
    switch (typ) {
            
        case 1: // check against ao index
            //where_clause = return_val_swao_frag(index_values_[index])     == index_values_[frag];
            where_clause = swao_frag[index_values_[index] - lower_seg_[index]]     == index_values_[frag];
            break;
            
        case 2: // check against occ index
            //where_clause = return_val_swocca_frag(index_values_[index])   == index_values_[frag];
            where_clause = swocca_frag[index_values_[index] - lower_seg_[index]]   == index_values_[frag];
            break;
            
        case 3: // check against virt index
            //where_clause = return_val_swvirta_frag(index_values_[index])  == index_values_[frag];
            where_clause = swvirta_frag[index_values_[index] - lower_seg_[index]]  == index_values_[frag];
            break;
            
        case 4: // check elst_dist[ifrag,jfrag] == ifrag
            //where_clause = return_val_elst_dist(index_values_[frag], index_values_[index]) == index_values_[frag];
            //ij = index_values_[frag] - 1 + (upper_bound_[0])*(index_values_[index] - 1);
            //std::cout << index_values_[index] << " " << index_values_[frag] << " " << ij << " " << elst_dist[ij] << std::endl;
            //where_clause = elst_dist[ij] == index_values_[frag];
            where_clause = elst_dist[index_values_[frag] - lower_seg_[frag]][index_values_[index] - lower_seg_[index]] == index_values_[frag];
            break;
            
        case 5: // check rcut_dist[ifrag,jfrag] == ifrag
            //where_clause = return_val_rcut_dist(index_values_[frag], index_values_[index]) == index_values_[frag];
            //ij = index_values_[index] - 1 + upper_bound_[0]*(index_values_[frag] - 1);
            //where_clause = rcut_dist[ij] == index_values_[frag];
            //break;
            where_clause = rcut_dist[index_values_[frag] - lower_seg_[frag]][index_values_[index] - lower_seg_[index]] == index_values_[frag];
            break;
            
        case 0:
            where_clause = index_values_[frag] != index_values_[index];
            break;
        default:
            where_clause = false;
    }
    //std::cout << "where_clause " << typ << " " << index_values_[index] << " " << index_values_[frag] << " " << where_clause << std::endl;
    return where_clause;
}


void Fragment_ij_oo_ao_PardoLoopManager::do_finalize() {
    for (int i = 0; i < num_indices_; ++i) {
        data_manager_.set_index_undefined(index_id_[i]);
    }
}

std::string Fragment_ij_oo_ao_PardoLoopManager::to_string() const {
    std::stringstream ss;
    ss << "Fragment_ij_oo_ao_PardoLoopManager:  num_indices="
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
                         const Fragment_ij_oo_ao_PardoLoopManager &obj) {
    os << obj.to_string();
    return os;
}

/*!
 -------------------------------------------
 _ij_aoo_o_
 -------------------------------------------
 */
    
    /*
     for each new special fragment where clause pattern, this should be the only thing realy changed.
     see comment above for fragment_special_where_clause syntax
     */
    bool Fragment_ij_aoo_o_PardoLoopManager::where_clause(int index) {
        bool where_;
        int ifrag = 0;
        int jfrag = 1;
        int ao = 1;
        int occ = 2;
        int virt = 3;
        int elst = 4;
        int rcut = 5;

        
        switch (index) {
            case 0:
            case 1:
                where_ = fragment_special_where_clause(elst,jfrag,ifrag);
                break;
            case 2:
                where_ = fragment_special_where_clause(ao,index,ifrag);
                break;
            case 3:
                where_ = fragment_special_where_clause(occ,index,ifrag);
                break;
            case 4:
                where_ = fragment_special_where_clause(occ,index,ifrag);
                break;
            case 5:
                where_ = fragment_special_where_clause(occ,index,jfrag);
                break;
            default:
                where_ = false;
        }
        return where_;
    }
    
    bool Fragment_ij_aoo_o_PardoLoopManager::do_update() {
        if (to_exit_)
            return false;
        bool more_iters;
        bool where_clauses_value;
        
        interpreter_->skip_where_clauses(num_where_clauses_);
        
        if (first_time_) {
            first_time_ = false;
            more_iters = initialize_indices();
            
            where_clauses_value = true;
            for (int i = 1; i < num_indices_; ++i) {
                where_clauses_value = where_clauses_value && where_clause(i);
            }
            if(where_clauses_value){
                iteration_++;
                if ((iteration_-1) % num_workers_ == company_rank_){
                    return true;
                }
            }
        }
        
        int index_restart[MAX_RANK];
        
        for (int i = 0; i < num_indices_; ++i) {
            index_restart[i] = index_values_[i];
        }
        
        int loop_count = iteration_;
        for (int index_i = index_restart[0]; index_i < upper_bound_[0]; ++index_i) {
            index_values_[0] = index_i;
            data_manager_.set_index_value(index_id_[0], index_i);
            
        for (int index_j = index_restart[1]; index_j < upper_bound_[1]; ++index_j) {
            index_values_[1] = index_j;
            data_manager_.set_index_value(index_id_[1], index_j);
                
            if (where_clause(1)) {
                for (int index_2 = index_restart[2]; index_2 < upper_bound_[2]; ++index_2) {
                    index_values_[2] = index_2;
                    data_manager_.set_index_value(index_id_[2], index_2);
                        
            if (where_clause(2)) {
                for (int index_3 = index_restart[3]; index_3 < upper_bound_[3]; ++index_3) {
                    index_values_[3] = index_3;
                    data_manager_.set_index_value(index_id_[3], index_3);
                    
            if (where_clause(3)) {
                for (int index_4 = index_restart[4]; index_4 < upper_bound_[4]; ++index_4) {
                    index_values_[4] = index_4;
                    data_manager_.set_index_value(index_id_[4], index_4);
                    
            if (where_clause(4)) {
                for (int index_5 = index_restart[5]; index_5 < upper_bound_[5]; ++index_5) {
                    index_values_[5] = index_5;
                    data_manager_.set_index_value(index_id_[5], index_5);
                    
            if (where_clause(5)) {
            if (loop_count > iteration_) {
                iteration_++;
                if ((iteration_-1) % num_workers_ == company_rank_){
                    return true;
                }
            }
            ++loop_count;
        } // where 5
            } // index_5
            for (int i = 5; i < num_indices_; ++i) { index_restart[i] = lower_seg_[i]; }
        } // where 4
            } // index_4
            for (int i = 4; i < num_indices_; ++i) { index_restart[i] = lower_seg_[i]; }
        } // where 3
            } // index_3
            for (int i = 3; i < num_indices_; ++i) { index_restart[i] = lower_seg_[i]; }
        } // where 2
            } // index_2
            for (int i = 2; i < num_indices_; ++i) { index_restart[i] = lower_seg_[i]; }
        } // where 1
            } // index_j
            for (int i = 1; i < num_indices_; ++i) { index_restart[i] = lower_seg_[i]; }
        } // index_i
        
        return false; //this should be false here
    }

Fragment_ij_aoo_o_PardoLoopManager::Fragment_ij_aoo_o_PardoLoopManager(
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
    
    form_elst_dist();
    //form_rcut_dist();
    form_swao_frag();
    form_swocca_frag();
    //form_swvirta_frag();
}

Fragment_ij_aoo_o_PardoLoopManager::~Fragment_ij_aoo_o_PardoLoopManager() {}

bool Fragment_ij_aoo_o_PardoLoopManager::increment_simple_pair(){
    bool more = false;      // More iterations?
    int current_value;
    for (int i = 0; i < 2; ++i) {
        more = increment_single_index(i);
        if (more) {break;}
    } //if here, then all indices are at their max value
    return more;
}

bool Fragment_ij_aoo_o_PardoLoopManager::increment_single_index(int index){
    bool more = false; 	// More iterations?
    int current_value;
    current_value = data_manager_.index_value(index_id_[index]);
    ++current_value;
    if (current_value < upper_bound_[index]) {
        //increment current index and return
        index_values_[index] = current_value;
        data_manager_.set_index_value(index_id_[index], current_value);
        more = true;
    } else {
        //wrap around and handle next index
        index_values_[index] = lower_seg_[index];
        data_manager_.set_index_value(index_id_[index], lower_seg_[index]);
    }
    return more;
}
    
    
    inline bool Fragment_ij_aoo_o_PardoLoopManager::increment_all() {
        bool more = false;      // More iterations?
        int current_value;
        for (int i = 0; i < num_indices_; ++i) {
            more = increment_single_index(i);
            if (more) {break;}
        } //if here, then all indices are at their max value
        return more;
    }
    
    void Fragment_ij_aoo_o_PardoLoopManager::form_elst_dist() {
        int elst_dist_array_slot = sip_tables_.array_slot(std::string("elst_dist"));
        Block::BlockPtr bptr_elst_dist = data_manager_.contiguous_array_manager().get_array(elst_dist_array_slot);
        double *val_elst_dist = bptr_elst_dist->get_data();
        int elst_dist_size = bptr_elst_dist->size();
        int num_frag = upper_bound_[0] - lower_seg_[0];
        elst_dist.resize(num_frag);
        for (int i = 0; i < num_frag; ++i) {
            //std::cout << i << std::endl;
            elst_dist[i].resize(num_frag);
        }
        int ij = 0;
        for (int i = 0; i < num_frag; ++i) {
            for (int j = 0; j < num_frag; ++j) {
                elst_dist[j][i] = (int)val_elst_dist[ij];
                //std::cout << i << " " << j << " " << ij << " " << elst_dist[j][i] << " " << val_elst_dist[ij]<< std::endl;
                ++ij;
            }
        }
        return;
    }

    void Fragment_ij_aoo_o_PardoLoopManager::form_rcut_dist() {
        int rcut_dist_array_slot = sip_tables_.array_slot(std::string("rcut_dist"));
        Block::BlockPtr bptr_rcut_dist = data_manager_.contiguous_array_manager().get_array(rcut_dist_array_slot);
        double *val_rcut_dist = bptr_rcut_dist->get_data();
        int rcut_dist_size = bptr_rcut_dist->size();
        int num_frag = upper_bound_[0] - lower_seg_[0];
        rcut_dist.resize(num_frag);
        for (int i = 0; i < num_frag; ++i) {
            //std::cout << i << std::endl;
            rcut_dist[i].resize(num_frag);
        }
        int ij = 0;
        for (int i = 0; i < num_frag; ++i) {
            for (int j = 0; j < num_frag; ++j) {
                rcut_dist[j][i] = (int)val_rcut_dist[ij];
                //std::cout << i << " " << j << " " << ij << " " << rcut_dist[j][i] << " " << val_rcut_dist[ij]<< std::endl;
                ++ij;
            }
        }
        return;
    }
    
    void Fragment_ij_aoo_o_PardoLoopManager::form_swao_frag() {
        int swao_frag_array_slot = sip_tables_.array_slot(std::string("swao_frag"));
        Block::BlockPtr bptr_swao_frag = data_manager_.contiguous_array_manager().get_array(swao_frag_array_slot);
        double *val_swao_frag = bptr_swao_frag->get_data();
        int swao_frag_size = bptr_swao_frag->size();
        swao_frag.resize(swao_frag_size);
        for (int i = 0; i < swao_frag_size; ++i) {
            swao_frag[i] = (int)val_swao_frag[i];
            //std::cout << "vec " << swao_frag[i] << " data " << val_swao_frag[i] << std::endl;
        }
        return;
    }
    
    void Fragment_ij_aoo_o_PardoLoopManager::form_swocca_frag() {
        int swocca_frag_array_slot = sip_tables_.array_slot(std::string("swocca_frag"));
        Block::BlockPtr bptr_swocca_frag = data_manager_.contiguous_array_manager().get_array(swocca_frag_array_slot);
        double *val_swocca_frag = bptr_swocca_frag->get_data();
        int swocca_frag_size = bptr_swocca_frag->size();
        swocca_frag.resize(swocca_frag_size);
        for (int i = 0; i < swocca_frag_size; ++i) {
            swocca_frag[i] = (int)val_swocca_frag[i];
            //std::cout << "vec " << swocca_frag[i] << " data " << val_swocca_frag[i] << std::endl;
        }
        return;
    }
    
    void Fragment_ij_aoo_o_PardoLoopManager::form_swvirta_frag() {
        int swvirta_frag_array_slot = sip_tables_.array_slot(std::string("swvirta_frag"));
        Block::BlockPtr bptr_swvirta_frag = data_manager_.contiguous_array_manager().get_array(swvirta_frag_array_slot);
        double *val_swvirta_frag = bptr_swvirta_frag->get_data();
        int swvirta_frag_size = bptr_swvirta_frag->size();
        swvirta_frag.resize(swvirta_frag_size);
        for (int i = 0; i < swvirta_frag_size; ++i) {
            swvirta_frag[i] = (int)val_swvirta_frag[i];
            //std::cout << "vec " << swvirta_frag[i] << " data " << val_swvirta_frag[i] << std::endl;
        }
        return;
    }
    
    inline bool Fragment_ij_aoo_o_PardoLoopManager::initialize_indices() {
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
            index_values_[i] = lower_seg_[i];
            data_manager_.set_index_value(index_id_[i], lower_seg_[i]);
        }
        more_iterations = true;
        return more_iterations;
    }

/*
 the generic where clauses are:
 typ = 1: map AO to fragment
 typ = 2: map Occ to fragment
 typ = 3: map Virt to fragment
 typ = 4: map elst_dist(ifrag,jfrag) == ifrag
 typ = 5: map rcut_dist(ifrag,jfrag) == ifrag
 typ = 0: ifrag != jfrag
 */
bool Fragment_ij_aoo_o_PardoLoopManager::fragment_special_where_clause(int typ, int index, int frag) {
    bool where_clause;
    int ij = 0;
    switch (typ) {
            
        case 1: // check against ao index
            //where_clause = return_val_swao_frag(index_values_[index])     == index_values_[frag];
            where_clause = swao_frag[index_values_[index] - lower_seg_[index]]     == index_values_[frag];
            break;
            
        case 2: // check against occ index
            //where_clause = return_val_swocca_frag(index_values_[index])   == index_values_[frag];
            where_clause = swocca_frag[index_values_[index] - lower_seg_[index]]   == index_values_[frag];
            break;
            
        case 3: // check against virt index
            //where_clause = return_val_swvirta_frag(index_values_[index])  == index_values_[frag];
            where_clause = swvirta_frag[index_values_[index] - lower_seg_[index]]  == index_values_[frag];
            break;
            
        case 4: // check elst_dist[ifrag,jfrag] == ifrag
            //where_clause = return_val_elst_dist(index_values_[frag], index_values_[index]) == index_values_[frag];
            //ij = index_values_[frag] - 1 + (upper_bound_[0])*(index_values_[index] - 1);
            //std::cout << index_values_[index] << " " << index_values_[frag] << " " << ij << " " << elst_dist[ij] << std::endl;
            //where_clause = elst_dist[ij] == index_values_[frag];
            where_clause = elst_dist[index_values_[frag] - lower_seg_[frag]][index_values_[index] - lower_seg_[index]] == index_values_[frag];
            break;
            
        case 5: // check rcut_dist[ifrag,jfrag] == ifrag
            //where_clause = return_val_rcut_dist(index_values_[frag], index_values_[index]) == index_values_[frag];
            //ij = index_values_[index] - 1 + upper_bound_[0]*(index_values_[frag] - 1);
            //where_clause = rcut_dist[ij] == index_values_[frag];
            //break;
            where_clause = rcut_dist[index_values_[frag] - lower_seg_[frag]][index_values_[index] - lower_seg_[index]] == index_values_[frag];
            break;
            
        case 0:
            where_clause = index_values_[frag] != index_values_[index];
            break;
        default:
            where_clause = false;
    }
    //std::cout << "where_clause " << typ << " " << index_values_[index] << " " << index_values_[frag] << " " << where_clause << std::endl;
    return where_clause;
}


void Fragment_ij_aoo_o_PardoLoopManager::do_finalize() {
    for (int i = 0; i < num_indices_; ++i) {
        data_manager_.set_index_undefined(index_id_[i]);
    }
}

std::string Fragment_ij_aoo_o_PardoLoopManager::to_string() const {
    std::stringstream ss;
    ss << "Fragment_ij_aoo_o_PardoLoopManager:  num_indices="
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
                         const Fragment_ij_aoo_o_PardoLoopManager &obj) {
    os << obj.to_string();
    return os;
}
    
    
    /*
     // this bock is commented out because we want to load the arrays at once.  These will break now also because the included values are statics instead of locals.  keep around just in case.
     
     double Fragment_ij_ao_vo_PardoLoopManager::return_val_rcut_dist(int index1,int index2) {
     int rcut_dist_array_slot = sip_tables_.array_slot(std::string("rcut_dist"));
     const index_value_array_t indices_rcut_dist = {index1, index2, unused_index_value, unused_index_value, unused_index_value, unused_index_value};
     BlockId bid_rcut_dist(rcut_dist_array_slot, indices_rcut_dist);
     Block::BlockPtr bptr_rcut_dist = data_manager_.block_manager_.get_block_for_reading(bid_rcut_dist);
     double val_rcut_dist = (int)(bptr_rcut_dist->get_data()[0]);
     return val_rcut_dist;
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
     */

#endif

} /* namespace sip */
