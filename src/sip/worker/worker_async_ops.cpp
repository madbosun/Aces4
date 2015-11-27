/*
 * worker_async_ops..cpp
 *
 *  Created on: Nov 27, 2015
 *      Author: basbas
 */

#include "worker_async_ops.h"

namespace sip {


WorkerGetAsync::WorkerGetAsync(int pc):
	AsyncBase(pc), mpi_request_(NULL){
	//a pointer to the mpi_request will be passed to the MPI_IRecv call, which will initialize it.
}

WorkerGetAsync::~WorkerGetAsync() {

}

std::string WorkerGetAsync::to_string() const {
	std::stringstream ss;
	ss << "WorkerGetAsync";
	ss << AsyncBase::to_string();
	ss << "mpi_request_=" << mpi_request_==NULL ? "NULL" : " not null";
	ss << std::endl;
	return ss.str();
}

WorkerPutAsync::WorkerPutAsync(int pc):
	AsyncBase(pc), mpi_request_(NULL){
	//a pointer to the mpi_request will be passed to the MPI_IRecv call, which will initialize it.
}

WorkerPutAsync::~WorkerPutAsync() {

}

std::string WorkerPutAsync::to_string() const {
	std::stringstream ss;
	ss << "WorkerPutAsync";
	ss << AsyncBase::to_string();
	ss << "mpi_request_=" << mpi_request_==NULL ? "NULL" : " not null";
	ss << std::endl;
	return ss.str();
}




} /* namespace sip */
