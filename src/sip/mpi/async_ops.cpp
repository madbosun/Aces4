/*
 * async_state.cpp
 *
 *  Created on: Jun 22, 2015
 *      Author: basbas
 */

#include "async_ops.h"
#include "server_block.h"
#include "disk_backed_block_map.h"
#include "sip_server.h"

namespace sip {


AsyncBase::AsyncBase(int pc) :
			async_status_(WAITING), pc_(pc) {
	}

AsyncBase::~AsyncBase() {
		check(async_status_ == DONE,
				"attempting to deconstruct unhandled request");
	}

	/** Tests the status of this operation
	 *
	 * @return  true if the operation is enabled or completed, and false otherwise
	 */
	bool AsyncBase::test() {
		if (async_status_ != WAITING)
			return true;
		if (do_test()) {
			async_status_ = READY;
			return true;
		}
		return false;
	}

	/**
	 * Waits for this operation to be enabled or completed.
	 * On return, async_status_ is READY or DONE.
	 */
	void AsyncBase::wait() {
		if (async_status_ != WAITING)
			return;
		do_wait();
		async_status_ = READY;
	}

	bool AsyncBase::is_done() {
		return async_status_ == DONE;
	}

	/**  Attempts to handle this operation.
	 *
	 * If the operation cannot be performed, then async_status remains WAITING and the method returns false
	 * If the operation can be performed, it is, asynch_status is set to DONE, and the method returns true
	 *
	 * @return  true if the operation has been completed and async_status_ set to DONE
	 */
	bool AsyncBase::try_handle() {
		if (async_status_ == READY || ((async_status_ == WAITING) && test())) {
			do_handle();  //this should always succeed if called when enabled
			async_status_ = DONE;
			return true;
		}
		return async_status_ == DONE;
	}

	/** returns true if the asynchronous operation represented by this object
	 * may modify the block
	 * @return
	 */
	bool AsyncBase::is_write() {
		return do_is_write();
	}

	bool AsyncBase::is_get(){
		return do_is_get();
	}

	std::ostream& operator<<(std::ostream& os, const AsyncBase &obj) {
		os << obj.to_string();
		return os;
	}


	//subclasses should invoke AsyncStatus::to_string() in overriding methods
	std::string AsyncBase::to_string() const {
		std::stringstream ss;
		ss << "async_status_=" << async_status_;
		ss << " pc_=" << pc_ << " ";
		return ss.str();
	}





/*************ServerBlockAsyncManager*************/


void BlockAsyncManager::add_async(AsyncBase* async){
	pending_.push_back(async);
	if (async->is_write()) num_pending_writes_++;
	if (async->is_get()) num_pending_gets_++;
}



} /* namespace sip */
