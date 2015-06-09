/*
 * loop_manager.h
 *
 *  Created on: Aug 23, 2013
 *      Author: basbas
 */

#ifndef LOOP_MANAGER_H_
#define LOOP_MANAGER_H_

#include "config.h"
#include "aces_defs.h"
#include "sip.h"
#include "array_constants.h"
#include "data_manager.h"

#include "sip_mpi_attr.h"

namespace sip {

class Interpreter;

/*! Base class for loop managers. */
class LoopManager {
public:
	LoopManager();
	virtual ~LoopManager();
	/* sets the value of the managed indices to the next value and returns true
	 * If there are no more values, the loop is done, return false;
	 */
	bool update() {
		return do_update();
	}
	/** reset the value of all managed indices to their undefined value */
	void finalize() {
		do_finalize();
	}
	/** set flag to immediately exit loop */
	void set_to_exit() {
		do_set_to_exit();
	}

	friend std::ostream& operator<<(std::ostream&, const LoopManager &);
protected:
	bool to_exit_;
	virtual std::string to_string() const;
private:
	virtual bool do_update() = 0;
	virtual void do_finalize() = 0;
	virtual void do_set_to_exit();
	DISALLOW_COPY_AND_ASSIGN(LoopManager);
};

class DoLoop: public LoopManager {
public:
	DoLoop(int index_id, DataManager & data_manager,
			const SipTables & sip_tables);
	virtual ~DoLoop();
	friend std::ostream& operator<<(std::ostream&, const DoLoop &);
protected:
	bool first_time_;
	int index_id_;
	int lower_seg_;
	int upper_bound_;  // first_segment_ <= current_value_ < upper_bound

	DataManager & data_manager_;
	const SipTables & sip_tables_;

	virtual std::string to_string() const;
private:
	virtual bool do_update();
	virtual void do_finalize();

	DISALLOW_COPY_AND_ASSIGN(DoLoop);
};

class SubindexDoLoop: public DoLoop {
public:
	SubindexDoLoop(int subindex_id, DataManager & data_manager,
			const SipTables & sip_tables);
	virtual ~SubindexDoLoop();
	friend std::ostream& operator<<(std::ostream&, const SubindexDoLoop &);
private:
	int parent_id_;
	int parent_value_;
	virtual std::string to_string() const;DISALLOW_COPY_AND_ASSIGN(SubindexDoLoop);
};
class SequentialPardoLoop: public LoopManager {
public:
	SequentialPardoLoop(int num_indices, const int (&index_ids)[MAX_RANK],
			DataManager & data_manager, const SipTables & sip_tables);
	virtual ~SequentialPardoLoop();
	friend std::ostream& operator<<(std::ostream&, const SequentialPardoLoop &);
private:
	virtual std::string to_string() const;
	virtual bool do_update();
	virtual void do_finalize();
	bool first_time_;
	int num_indices_;
	sip::index_selector_t index_id_;
	sip::index_value_array_t lower_seg_;
	sip::index_value_array_t upper_bound_;

	DataManager & data_manager_;
	const SipTables & sip_tables_;

	DISALLOW_COPY_AND_ASSIGN(SequentialPardoLoop);
};

#ifdef HAVE_MPI
class StaticTaskAllocParallelPardoLoop: public LoopManager {
public:
	StaticTaskAllocParallelPardoLoop(int num_indices,
			const int (&index_ids)[MAX_RANK], DataManager & data_manager,
			const SipTables & sip_tables, SIPMPIAttr& sip_mpi_attr);
	virtual ~StaticTaskAllocParallelPardoLoop();
	friend std::ostream& operator<<(std::ostream&,
			const StaticTaskAllocParallelPardoLoop &);
private:
	virtual std::string to_string() const;
	virtual bool do_update();
	virtual void do_finalize();
	bool first_time_;
	int num_indices_;
	long iteration_;
	sip::index_selector_t index_id_;
	sip::index_value_array_t lower_seg_;
	sip::index_value_array_t upper_bound_;

	DataManager & data_manager_;
	const SipTables & sip_tables_;
	SIPMPIAttr & sip_mpi_attr_;


	bool increment_indices();
	bool initialize_indices();

	DISALLOW_COPY_AND_ASSIGN(StaticTaskAllocParallelPardoLoop);

};




class BalancedTaskAllocParallelPardoLoop: public LoopManager {
public:
	BalancedTaskAllocParallelPardoLoop(int num_indices,
			const int (&index_ids)[MAX_RANK], DataManager & data_manager,
			const SipTables & sip_tables, SIPMPIAttr& sip_mpi_attr,
			int num_where_clauses, Interpreter* interpreter, long& iteration);
	virtual ~BalancedTaskAllocParallelPardoLoop();
	friend std::ostream& operator<<(std::ostream&,
			const BalancedTaskAllocParallelPardoLoop &);
private:
	virtual std::string to_string() const;
	virtual bool do_update();
	virtual void do_finalize();
	bool first_time_;
	int num_indices_;
	long& iteration_;
	index_selector_t index_id_;
	index_value_array_t lower_seg_;
	index_value_array_t upper_bound_;
	int num_where_clauses_;

	DataManager & data_manager_;
	const SipTables & sip_tables_;
	SIPMPIAttr & sip_mpi_attr_;
	int company_rank_;
	int num_workers_;
	Interpreter* interpreter_;

	bool increment_indices();
	bool initialize_indices();

	DISALLOW_COPY_AND_ASSIGN(BalancedTaskAllocParallelPardoLoop);

};

/**
 * Special Loop manager for Fragment code mcpt2_corr_lowmem.
 * Optimizes where clause lookups that involves looking up blocks.
 *
 *  naming scheme: _<simple index fragment indices>_<occ (o), virt (v), ao (a)
 *  in ifrag>_<o,v,a in jfrag>
 *
 *  */

class Fragment_ij_aa__PardoLoopManager: public LoopManager{
public:
	Fragment_ij_aa__PardoLoopManager(int num_indices,
			const int (&index_ids)[MAX_RANK], DataManager & data_manager,
			const SipTables & sip_tables, SIPMPIAttr& sip_mpi_attr,
			int num_where_clauses, Interpreter* interpreter, long& iteration);
	virtual ~Fragment_ij_aa__PardoLoopManager();
	friend std::ostream& operator<<(std::ostream&,
				const Fragment_ij_aa__PardoLoopManager &);
private:
	virtual std::string to_string() const;
	virtual bool do_update();
	virtual void do_finalize();
	bool first_time_;
	int num_indices_;
	long& iteration_;
	index_selector_t index_id_;
	index_value_array_t lower_seg_;
	index_value_array_t upper_bound_;
	int num_where_clauses_;

	DataManager & data_manager_;
	const SipTables & sip_tables_;
	SIPMPIAttr & sip_mpi_attr_;
	int company_rank_;
	int num_workers_;
	Interpreter* interpreter_;

	int index_values_[MAX_RANK];

	bool increment_special();
	bool initialize_indices();

        int index1;
        int index2;
	double return_val_elst_dist(int index1, int index2);
	double return_val_swao_frag(int index1);
	double return_val_swocca_frag(int index1);
	double return_val_swvirta_frag(int index1);
};

class Fragment_ij_ao_vo_PardoLoopManager: public LoopManager{
public:
	Fragment_ij_ao_vo_PardoLoopManager(int num_indices,
			const int (&index_ids)[MAX_RANK], DataManager & data_manager,
			const SipTables & sip_tables, SIPMPIAttr& sip_mpi_attr,
			int num_where_clauses, Interpreter* interpreter, long& iteration);
	virtual ~Fragment_ij_ao_vo_PardoLoopManager();
	friend std::ostream& operator<<(std::ostream&,
				const Fragment_ij_ao_vo_PardoLoopManager &);
private:
	virtual std::string to_string() const;
	virtual bool do_update();
	virtual void do_finalize();
	bool first_time_;
	int num_indices_;
	long& iteration_;
	index_selector_t index_id_;
	index_value_array_t lower_seg_;
	index_value_array_t upper_bound_;
	int num_where_clauses_;

	DataManager & data_manager_;
	const SipTables & sip_tables_;
	SIPMPIAttr & sip_mpi_attr_;
	int company_rank_;
	int num_workers_;
	Interpreter* interpreter_;

	int index_values_[MAX_RANK];

	bool increment_special();
	bool initialize_indices();

        int index1;
        int index2;
	double return_val_elst_dist(int index1, int index2);
	double return_val_swao_frag(int index1);
	double return_val_swocca_frag(int index1);
	double return_val_swvirta_frag(int index1);
};

#endif

} /* namespace sip */
#endif /* LOOP_MANAGER_H_ */
