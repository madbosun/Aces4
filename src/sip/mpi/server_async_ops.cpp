/*
 * server_sync_ops.cpp
 *
 *  Created on: Nov 27, 2015
 *      Author: basbas
 */

#include "server_async_ops.h"
#include "server_block.h"
#include "sip_server.h"

namespace sip {

/************* GetAsync ************/

ServerGetAsync::ServerGetAsync(int mpi_source, int get_tag, ServerBlock* block, int pc) :
		AsyncBase(pc), mpi_request_() {
	//send block
	SIPMPIUtils::check_err(
			MPI_Isend(block->get_data(), block->size(), MPI_DOUBLE, mpi_source,
					get_tag, MPI_COMM_WORLD, &mpi_request_), __LINE__,
			__FILE__);
}



/**  PutAccumulateAsync *******************/

ServerPutAccumulateDataAsync::ServerPutAccumulateDataAsync(int mpi_source,
		int put_accumulate_data_tag, ServerBlock* block, int pc) :
		AsyncBase(pc), mpi_source_(mpi_source), tag_(put_accumulate_data_tag),  block_(
				block), mpi_request_(MPI_REQUEST_NULL) {
	//allocate temp buffer
    size_t size = block_->size();
	temp_ = SIPServer::global_sipserver->disk_backed_block_map_.allocate_data(size, false);
	//post receive
	SIPMPIUtils::check_err(
			MPI_Irecv(temp_, block->size(), MPI_DOUBLE, mpi_source,
					put_accumulate_data_tag, MPI_COMM_WORLD, &mpi_request_), __LINE__,
			__FILE__);
	//reply should be sent in calling routing AFTER this object has been constructed,
	//so the Irecv will already be posted.
}

ServerPutAccumulateDataAsync::~ServerPutAccumulateDataAsync() {
	check(temp_ != NULL,
			"destructor of PutAccumulateDataAsync invoked with null data_");
	//delete[] temp_;
	SIPServer::global_sipserver->disk_backed_block_map_.free_data(temp_, block_->size());
}

bool ServerPutAccumulateDataAsync::do_test() {
	int flag = 0;
	MPI_Status status;
	MPI_Test(&mpi_request_, &flag, &status);
	if (flag) {
		//check that received message was expected size
		int count;
		MPI_Get_count(&status, MPI_DOUBLE, &count);
		check(count == block_->size(), "count != block_->size()");
	}
	return flag;
}

void ServerPutAccumulateDataAsync::do_handle() {
	//send ack to worker
	SIPMPIUtils::check_err(
			MPI_Send(0, 0, MPI_INT, mpi_source_, tag_, MPI_COMM_WORLD),
			__LINE__, __FILE__);
	//accumulate received data into block
	block_->accumulate_data(block_->size(), temp_);
}

void ServerPutAccumulateDataAsync::do_wait() {
	MPI_Status status;
	MPI_Wait(&mpi_request_, &status);
	//check that received message was expected size
	int count;
	MPI_Get_count(&status, MPI_DOUBLE, &count);
	check(count == block_->size(), "count != block_->size()");
}

std::string ServerPutAccumulateDataAsync::to_string() const {
	std::stringstream ss;
	ss << "PutAccumulateDataAsync";
	ss << AsyncBase::to_string();
	ss << " source=" << mpi_source_ << " tag=" << tag_;
	int transaction;
	SIPMPIConstants::MessageType_t message_type;
	BarrierSupport::decode_tag(tag_, message_type, transaction);
	ss << " message type="
			<< SIPMPIConstants::messageTypeToName(
					SIPMPIConstants::intToMessageType(message_type))
			<< " transaction=" << transaction << std::cout;
	return ss.str();
}


/******************* PutDataAsync **********/

ServerPutDataAsync::ServerPutDataAsync(int mpi_source, int put_data_tag, ServerBlock* block,
		int pc) :
		AsyncBase(pc), mpi_source_(mpi_source), tag_(put_data_tag), block_(
				block), mpi_request_(MPI_REQUEST_NULL) {
	//post receive
	SIPMPIUtils::check_err(
			MPI_Irecv(block->get_data(), block->size(), MPI_DOUBLE, mpi_source,
					put_data_tag, MPI_COMM_WORLD, &mpi_request_), __LINE__,
			__FILE__);
}

bool ServerPutDataAsync::do_test() {
	int flag = 0;
	MPI_Status status;
	MPI_Test(&mpi_request_, &flag, &status);
	if (flag) {
		int count;
		MPI_Get_count(&status, MPI_DOUBLE, &count);
		check(count == block_->size(), "count != block->size()");
	}
	return flag;
}

void ServerPutDataAsync::do_wait() {
	MPI_Status status;
	MPI_Wait(&mpi_request_, &status);
	//check that received message was expected size
	int count;
	MPI_Get_count(&status, MPI_DOUBLE, &count);
	check(count == block_->size(), "count != block_->size()");
}

void ServerPutDataAsync::do_handle() {
	//send ack to worker
	SIPMPIUtils::check_err(
			MPI_Send(0, 0, MPI_INT, mpi_source_, tag_, MPI_COMM_WORLD),
			__LINE__, __FILE__);
}

std::string ServerPutDataAsync::to_string() const {
	std::stringstream ss;
	ss << "PutDataAsync";
	ss << AsyncBase::to_string();
	ss << " source=" << mpi_source_ << " tag=" << tag_;
	int transaction;
	SIPMPIConstants::MessageType_t message_type;
	BarrierSupport::decode_tag(tag_, message_type, transaction);
	ss << " message type="
			<< SIPMPIConstants::messageTypeToName(
					SIPMPIConstants::intToMessageType(message_type))
			<< " transaction=" << transaction << std::cout;
	return ss.str();
}



} /* namespace sip */
