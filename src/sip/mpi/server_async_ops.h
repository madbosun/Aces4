/*
 * server_async_ops.h
 *
 *  Created on: Nov 27, 2015
 *      Author: basbas
 */

#ifndef SERVER_ASYNC_OPS_H_
#define SERVER_ASYNC_OPS_H_

#include "async_ops.h"
#include <mpi.h>

namespace sip {
class ServerBlock;

/** Represents asynchronous send of block by server in response to get operation
 * Once the send operation is complete, no additional "handling" needs to be done.
 *
 * The block used in the operation should have been obtained by a call to
 * get_block_for reading, which should wait for any pending writes to complete.
 */
class ServerGetAsync: public AsyncBase {
public:
	//asynchronous send with response performed in constructor.
	ServerGetAsync(int mpi_source, int get_tag, ServerBlock* block, int pc);
	virtual ~ServerGetAsync() {
	}

private:
	MPI_Request mpi_request_;

	bool do_test() {
		int flag = 0;
		MPI_Test(&mpi_request_, &flag, MPI_STATUS_IGNORE);
		return flag;
	}

	void do_wait() {
		MPI_Wait(&mpi_request_, MPI_STATUS_IGNORE);
	}

	void do_handle() {
		//nothing to do here since we have either returned from wait or test returned true
	}

	bool do_is_write() {
		return false;
	}

	std::string to_string() const {
		std::stringstream ss;
		ss << "GetAsync";
		ss << AsyncBase::to_string();
		return ss.str();
	}
	DISALLOW_COPY_AND_ASSIGN (ServerGetAsync);
};

/** Represents asynchronous put_accumulate operation at server
 *
 * Instance will be created in response to put_accumulate message to handle
 * the second message containing the data and to perform the accumulate operation.
 *
 * Instances of this class own the temporary buffer used to receive the data.
 *
 * Constructor:
 *             allocates temporary buffer,
 *             posts Irecv for block data.
 *
 * do_handle: sends ack for data message to source
 *           performs the accumulate operation of temp data into block
 *
 * Destructor: deletes the temp data buffer.
 *
 *The constructor does not send the "reply" message to the source to let it know that
 *the data message can be sent.  This should be done AFTER this object is constructed
 *in order to ensure that the Irecv has been posted, and thus ensuring that the
 *large message can bypass the MPI buffers
 *
 *TODO:  reuse the temp buffers?
 *
 */
class ServerPutAccumulateDataAsync: public AsyncBase {
public:
	ServerPutAccumulateDataAsync(int mpi_source, int put_accumulate_data_tag,
			ServerBlock* block, int pc);
	virtual ~ServerPutAccumulateDataAsync();
private:
	ServerBlock* block_;
	int mpi_source_;
	int tag_;
	MPI_Request mpi_request_;
	double* temp_; //buffer to receive data.  created in constructor, deleted in destructor

	bool do_test();
	void do_wait();
	void do_handle();
	bool do_is_write() {
		return true;
	}
	virtual std::string to_string() const;
	DISALLOW_COPY_AND_ASSIGN (ServerPutAccumulateDataAsync);
};

/** Represents asynchronous put data operation at server.
 *
 * Will be created in response to put message.
 *
 * This class assumes that there is at most one unhandled put message per block at any time.
 * This property is guaranteed because the block should have been obtained by a call
 * to get_block_for_writing,  This method should call the block's
 * CommunicationState.wait method which waits for all pending ops to complete.
 *
 * Constructor: sends reply to source,
 *             posts Irecv for block data
 *             The receive buffer is the block's data array, which should exist already
 *
 * do_handle: sends ack for data message to source
 *
 *The constructor does not send the "reply" message to the source to let it know that
 *the data message can be sent.  This should be done AFTER this object is constructed
 *order to ensure that the Irecv has been posted (and thus ensuring that the message
 *bypasses MPI buffers
 *
 */
class ServerPutDataAsync: public AsyncBase {

public:
	ServerPutDataAsync(int mpi_source, int put_data_tag, ServerBlock* block, int pc);
	virtual ~ServerPutDataAsync() {
	}

private:
	ServerBlock* block_;
	int mpi_source_;
	int tag_;
	MPI_Request mpi_request_;

	bool do_test();
	void do_wait();
	void do_handle();
	bool do_is_write() {
		return true;
	}
	std::string to_string() const;
	DISALLOW_COPY_AND_ASSIGN (ServerPutDataAsync);
};



} /* namespace sip */

#endif /* SERVER_ASYNC_OPS_H_ */
