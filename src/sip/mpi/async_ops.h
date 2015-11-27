/*
 * async_state.h
 *
 *  Created on: Jun 22, 2015
 *      Author: basbas
 */

#ifndef ASYNC_STATE_H_
#define ASYNC_STATE_H_

#include <list>
#include <mpi.h>
#include "sip.h"
#include "sip_mpi_constants.h"
#include "barrier_support.h"
#include "sip_mpi_utils.h"


namespace sip {

class ServerBlock;

/**
 * Base class for pending asynchronous operations on a block
 *
 * These objects are owned by a data structure belonging to the block and constructed when
 * inserted and destructed when removed.  These objects are not copyable--there should be
 * a single instance per operation.
 *
 * Semantic requirements and coding standards related to asynchronous operations:
 *
 * The constructor of each class should perform the asynchronous send
 * or post the asynchronous receive that derived objects of this class are handling.
 *
 * The do_handle method should perform tasks necessary to complete the operation
 * after it has been enabled.  This may require sending a message, removing
 * items from a data structure, or doing nothing (in the case the op is an
 * asynchronous mpi send).
 *
 * Closure: Subclasses should contain enough information for their methods
 * to be callable from any context.
 *
 * Nonblocking:  An enabled operation should be able to be executed to completion without
 * blocking.
 *
 *
 * This class is a wrapper for these operations which ensures that calls are idempotent.
 *
 *
 */
class AsyncBase {
public:
	enum AsyncStatus {
		WAITING = 1, READY = 2, DONE = 3,
	};

	AsyncBase(int pc);

	virtual ~AsyncBase();

	/** Tests the status of this operation
	 *
	 * @return  true if the operation is enabled or completed, and false otherwise
	 */
	bool test();

	/**
	 * Waits for this operation to be enabled or completed.
	 * On return, async_status_ is READY or DONE.
	 */
	void wait();

	bool is_done();

	/**  Attempts to handle this operation.
	 *
	 * If the operation cannot be performed, then async_status remains WAITING and the method returns false
	 * If the operation can be performed, it is, asynch_status is set to DONE, and the method returns true
	 *
	 * @return  true if the operation has been completed and async_status_ set to DONE
	 */
	bool try_handle();

	/** returns true if the asynchronous operation represented by this object
	 * may modify the block
	 * @return
	 */
	bool is_write();

	/** returns true if the asynchronous operation represented by this object is a get operation
	 * at the worker.
	 *
	 * @return
	 */
	bool is_get();

	friend std::ostream& operator<<(std::ostream& os, const AsyncBase &obj);

protected:
	//subclasses should invoke AsyncStatus::to_string() in overriding methods
	//to handle private variables in this base class.
	virtual std::string to_string() const;

private:
	AsyncStatus async_status_; //initially WAITING, must be DONE before destruction
	int pc_; //index into optable of the sialx instruction that generated this operation

	//precondition:  async_state == WAITNG
	//postcondition: returned value == "operation is enabled"
	virtual bool do_test()=0;

	//precondition:   async_state == WAITING
	//postcondition:  "operation is enabled"
	virtual void do_wait()=0;

	//precondition:  async_state == READY (or WAITING and test has returned true)
	//postcondition:  op completed
	virtual void do_handle()=0;

	//returns true if this operation may modify the block
	virtual bool do_is_write()=0;

	//returns true if this operation is a worker get
	virtual bool do_is_get(){return false;}

	DISALLOW_COPY_AND_ASSIGN (AsyncBase);
};


/** This class manages the pending asynchronous
 * operations on a single block.  It contains
 * a list of pending operations.
 *
 * At the server, get_block_for_reading should wait for all pending write  ops,
 * and and get_block_for_writing stop and wait for all pending operations
 * to be completed before returning. get_block_for_accumulate,
 * on the other hand does not wait, but adds its async_op to the list
 * and returns.
 *
 * All asynchronous operations are handled in the order they are submitted to this list.
 *
 * Methods in this class are responsible for appropriately maintaining
 * the invariant that num_pending_writes_ is the number of ops in the
 * pending list which may modify the block (i.e. for which is_write()
 * will return true.
 *
 */
class BlockAsyncManager {
public:
	BlockAsyncManager():
		num_pending_writes_(0),
		num_pending_gets_(0){
	}
	~BlockAsyncManager() {
		check(pending_.empty(), "deleting block with pending async ops");
	}

	void add_async(AsyncBase* async);

	/**
	 *
	 * @return true if this block has pending operations
	 */
	bool has_pending() {
		return !pending_.empty();
	}

	/**
	 * Attempts to handle all pending items in the list, in order.  If an item is found that
	 * is not enabled, false is returned. Otherwise, the list will be empty and true is returned.
	 * On a return value of true, if the caller is keeping track of blocks with pending ops,
	 * this block should be deleted from that set.
	 *
	 * @return  false if the pending list for the block is not on return empty.
	 *
	 */
	bool try_handle_all_test_none_pending() {
		std::list<AsyncBase*>::iterator it = pending_.begin();
		while (it != pending_.end()) {
			bool res = (*it)->try_handle();
			if (res) {
				if ((*it)->is_write()) {
					num_pending_writes_--;
				}
				delete *it;
				it = pending_.erase(it);
			} else
				return false;
		}
		return true;
	}

	/**
	 * Attempts to handle the first item in the list.  If the first item is not enabled, the routine
	 * returns immediately.  Otherwise, the item will be handled and removed from the list.  If the
	 * first item is actually done, it will be removed from the list, and the next item tried.
	 * If the list is empty, or becomes empty after removing done operations, true is returned,
	 * otherwise false.
	 *
	 * @return  false if the pending list for the block is not empty after attempting to handle one op.
	 *
	 * From the point of view of the caller, try_handle_all and try_handle look the
	 * same, a return value of true means that all pending ops for this block have
	 * been handled.  The difference is how long the methods keep trying.
	 *
	 */
	bool try_handle_test_none_pending() {
		std::list<AsyncBase*>::iterator it = pending_.begin();
		//remove already done ops, if any, from list
		while (it != pending_.end() && (*it)->is_done()) {
			if ((*it)->is_write()) {
				num_pending_writes_--;
			}
			delete *it;
			it = pending_.erase(it);
		}
		if (it == pending_.end()) {
			//nothing left to handle
			return true;
		}
		//try to handle first op
		bool res = (*it)->try_handle();
		if (res) {
			if ((*it)->is_write()) {
				num_pending_writes_--;
			}
			delete *it;
			it = pending_.erase(it);
		}
		return pending_.empty();

	}

	/** Handles all pending asyncs on this block, waiting for each one in
	 * turn to be enabled.
	 *
	 * Postcondition:  there are no pending ops on this block
	 * */

	void wait_all() {
		std::list<AsyncBase*>::iterator it = pending_.begin();
		while (it != pending_.end()) {
			(*it)->wait();
			bool res = (*it)->try_handle();
			check(res, "in BlockAsyncManager::wait--handle ready op failed");
			if ((*it)->is_write()) {
				num_pending_writes_--;
			}
			delete *it;
			it = pending_.erase(it);
		}
	}

	/** Handles all pending writes on this block, waiting for each one in
	 * turn to be enabled.  Since ops must be handled in order, any pending
	 * reads preceding a write in the list are also awaited
	 *
	 * Postcondition:  there are no pending write ops on this block
	 * */
	void wait_for_writes() {
		std::list<AsyncBase*>::iterator it = pending_.begin();
		while (it != pending_.end() && num_pending_writes_>0) {
			(*it)->wait();
			bool res = (*it)->try_handle();
			check(res, "in BlockAsyncManager::wait--handle ready op failed");
			if ((*it)->is_write()) {
				num_pending_writes_--;
			}
			delete *it;
			it = pending_.erase(it);
		}
	}

private:
	std::list<AsyncBase*> pending_;
	int num_pending_writes_;
	int num_pending_gets_; //only useful on worker to avoid duplicate gets.
	DISALLOW_COPY_AND_ASSIGN (BlockAsyncManager);
};

} /* namespace sip */

#endif /* ASYNC_STATE_H_ */
