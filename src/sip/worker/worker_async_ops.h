/*
 * worker_async_ops.h
 *
 *  Created on: Nov 27, 2015
 *      Author: basbas
 */

#ifndef WORKER_ASYNC_OPS_H_
#define WORKER_ASYNC_OPS_H_


#include "async_ops.h"

#ifdef HAVE_MPI
#include <mpi.h>
#endif

namespace sip {


/** Represents an in progress get operation at a worker
 */
class WorkerGetAsync: public AsyncBase {
public:
	explicit WorkerGetAsync(int pc);
	virtual ~WorkerGetAsync();

	MPI_Request* mpi_request();

private:
	MPI_Request mpi_request_;

	bool do_test();
	void do_wait();
	void do_handle();
	bool do_is_write();
	bool do_is_get();
	std::string to_string() const;

	DISALLOW_COPY_AND_ASSIGN (WorkerGetAsync);
};


/**
 * Async Op for both Put and PutAccumulate requests at workers
 */
class WorkerPutAsync: public AsyncBase {
public:
	explicit WorkerPutAsync(int pc);
	virtual ~WorkerPutAsync();

	MPI_Request* mpi_request();

private:
	MPI_Request mpi_request_;

	bool do_test();
	void do_wait();
	void do_handle();
	bool do_is_write();
	std::string to_string() const;

	DISALLOW_COPY_AND_ASSIGN (WorkerPutAsync);
};







/**Inlined functions for WorkerGetAsync**********/

inline MPI_Request* WorkerGetAsync::mpi_request(){return &mpi_request_;}

inline bool WorkerGetAsync::do_test() {
	int flag = 0;
	MPI_Test(&mpi_request_, &flag, MPI_STATUS_IGNORE);
	return flag;
}

inline void WorkerGetAsync::do_wait() {
	MPI_Wait(&mpi_request_, MPI_STATUS_IGNORE);
}

inline void WorkerGetAsync::do_handle() {}

//This may seem counterintuitive but the buffer belonging to the operation will
//be written with data from the server.
//
inline bool WorkerGetAsync::do_is_write() {
	return true;
}

inline bool WorkerGetAsync::do_is_get(){
	return true;
}


/** Inlined functins for WorkerPutAsync */

inline MPI_Request* WorkerPutAsync::mpi_request(){return &mpi_request_;}

inline bool WorkerPutAsync::do_test() {
	int flag = 0;
	MPI_Test(&mpi_request_, &flag, MPI_STATUS_IGNORE);
	return flag;
}

inline void WorkerPutAsync::do_wait() {
	MPI_Wait(&mpi_request_, MPI_STATUS_IGNORE);
}

inline void WorkerPutAsync::do_handle() {}

//This may seem counterintuitive but the buffer belonging to the operation will
//be written with data from the server.
//
inline bool WorkerPutAsync::do_is_write() {
	return false;
}


} /* namespace sip */

#endif /* WORKER_ASYNC_OPS_H_ */
