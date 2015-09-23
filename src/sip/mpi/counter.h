/*
 * counter.h
 *
 *  Created on: Jan 19, 2015
 *      Author: basbas
 */

#ifndef COUNTER_H_
#define COUNTER_H_

#include <ctime>
#include <ostream>
#include <sstream>
#include <vector>
#include <iomanip>
#include <limits>
#include "sip.h"
#include "sip_tables.h"
#ifdef HAVE_MPI
#include "sip_mpi_attr.h"
#endif
namespace sip {
template<typename T, typename R> class Timer;
template<typename T, typename R> std::ostream& operator<<(std::ostream&,
		const Timer<T, R>&);

/**Basic monotonic counter */
template<typename T>
class Counter {
public:
	Counter() :
			value_(0) {
	}
	~Counter() {
	}
	;
	void inc(size_t delta = 1) {
		value_ += delta;
	}
	size_t get_value() {
		return value_;
	}
	void reset() {
		value_ = 0;
	}
	void gather() {
		static_cast<T*>(this)->gather_impl();
	}
	void reduce() {
		static_cast<T*>(this)->reduce_impl();
	}
	friend std::ostream& operator<<<>(std::ostream& os, const Counter<T>& obj);
protected:
	size_t value_;
	std::ostream& stream_out(std::ostream& os) const {
		os << value_;
		return os;
	}
	DISALLOW_COPY_AND_ASSIGN(Counter);
};

template<typename T>
std::ostream& operator<<(std::ostream& os, const Counter<T>& obj) {
	return static_cast<const T&>(obj).stream_out(os);
}

template<typename T>
class NoopCounter{
public:
	NoopCounter(const MPI_Comm&){}
	void inc(T delta = 0){}
	void dec(T delta = 0){}
	T get_value(){}
	void reset(){}
	void gather(){}
	void reduce(){}
	void set(T val){}
	std::ostream& stream_out(std::ostream& os) const{ return os;}
	friend std::ostream& operator<<<>(std::ostream& os, const NoopCounter<T>& obj);
};

template<typename T>
std::ostream& operator<<(std::ostream& os, const NoopCounter<T>& obj) {
	return os;
}


class SingleNodeCounter: public Counter<SingleNodeCounter> {
public:
	void gather_impl() {
	}
	void reduce_impl() {
	}
	template<typename T> friend class Counter;
	template<typename T> friend std::ostream& operator<<(std::ostream& os,
			const Counter<T>& obj);
protected:
	std::ostream& stream_out(std::ostream& os) const {
		os << "value_=" << "value";
		return os;
	}
	DISALLOW_COPY_AND_ASSIGN(SingleNodeCounter);
};




class MPICounter: public Counter<MPICounter> {
public:
	explicit MPICounter(const MPI_Comm& comm) :
			comm_(comm), gathered_vals_(0), reduced_val_(0), gather_done_(
					false), reduce_done_(false) {
	}
	template<typename T> friend class Counter;
	template<typename T> friend std::ostream& operator<<(std::ostream& os,
			const Counter<T>& obj);
protected:
	const MPI_Comm& comm_;
	std::vector<size_t> gathered_vals_;
	size_t reduced_val_;
	bool gather_done_;
	bool reduce_done_;

	std::ostream& stream_out(std::ostream& os) const {
		if (!gather_done_ && !reduce_done_) {
			//print own data
			int rank;
			MPI_Comm_rank(comm_, &rank);
			os << "company rank," << rank << ", value_," << value_ << std::endl;
			return os;
		}
		int comm_size;
		MPI_Comm_size(comm_, &comm_size);

		if (gather_done_) {
			//output gathered values in comma separated list
			os << "value";
			std::vector<size_t>::const_iterator it = gathered_vals_.begin();
			for (int j = 0; j < comm_size; ++j) {
				os << ',' << *it;
				++it;
			}
		}
		if (reduce_done_) {
			os << std::endl << "reduced mean=";
			os
					<< static_cast<double>(reduced_val_)
							/ static_cast<double>(comm_size);
		}
		os << std::endl;
		return os;
	}
	void gather_impl() {
//		check(sizeof(MPI_UNSIGNED_LONG) == sizeof(size_t),
//				"mismatch in  mpi and c++ type ");
		int rank;
		int comm_size;
		MPI_Comm_rank(comm_, &rank);
		MPI_Comm_size(comm_, &comm_size);
		if (rank == 0) {
			gathered_vals_.resize(comm_size, 0);
			gather_done_ = true;
		}
		if (comm_size > 1) {
			MPI_Gather(&value_, 1, MPI_UNSIGNED_LONG, gathered_vals_.data(), 1,
					MPI_UNSIGNED_LONG, 0, comm_);
		} else {
			gathered_vals_.at(0) = value_;
		}
	}
	void reduce_impl() {
		int rank;
		int comm_size;
		MPI_Comm_rank(comm_, &rank);
		MPI_Comm_size(comm_, &comm_size);

		if (comm_size > 1) {
			MPI_Reduce(&value_, &reduced_val_, 1, MPI_UNSIGNED_LONG, MPI_SUM, 0,
					comm_);
		} else {
			reduced_val_ = value_;
		}
		if (rank == 0) {
			reduce_done_ = true;
		}
	}

private:
	DISALLOW_COPY_AND_ASSIGN(MPICounter);
};

template<typename T>
class CounterList {
public:
	explicit CounterList(size_t size, bool filter = true) :
			size_(size), list_(size,0), filter_(filter) {
	}
	~CounterList() {
	}
	void inc(int index, size_t delta = 1) {
		list_[index] += delta;
	}
	size_t get_value(int index) {
		return list_[index];
	}
	void reset() {
		std::fill(list_.begin(), list_.end(), 0);
	}
	void gather() {
		static_cast<T*>(this)->gather_impl();
	}
	void reduce() {
		static_cast<T*>(this)->reduce_impl();
	}
	friend std::ostream& operator<<<>(std::ostream& os,
			const CounterList<T>& obj);
protected:
	std::vector<size_t> list_;
	size_t size_;
	bool filter_;
	DISALLOW_COPY_AND_ASSIGN(CounterList);
};

template<typename T>
std::ostream& operator<<(std::ostream& os, const CounterList<T>& obj) {
	return static_cast<const T&>(obj).stream_out(os);
}

class MPICounterList: public CounterList<MPICounterList> {
public:
	MPICounterList(const MPI_Comm& comm, size_t size) :
			CounterList<MPICounterList>(size), comm_(comm), gathered_vals_(0), reduced_vals_(
					0), reduce_done_(false), gather_done_(false) {
	}
	template<typename T> friend class CounterList;
	template<typename T> friend std::ostream& operator<<(std::ostream& os,
			const CounterList<T>& obj);
protected:
	const MPI_Comm& comm_;
	std::vector<size_t> gathered_vals_;
	std::vector<size_t> reduced_vals_;
	bool gather_done_;
	bool reduce_done_;
	void gather_impl() {
		int rank;
		int comm_size;
		MPI_Comm_rank(comm_, &rank);
		MPI_Comm_size(comm_, &comm_size);
		if (rank == 0) {
			gathered_vals_.resize(comm_size * size_, 0);
			gather_done_ = true;
		}
		if (comm_size > 1) {
			MPI_Gather(list_.data(), size_, MPI_UNSIGNED_LONG,
					gathered_vals_.data(), size_, MPI_UNSIGNED_LONG, 0, comm_);
		} else {
			gathered_vals_ = list_;
		}

	}
	void reduce_impl() {
		int rank;
		int comm_size;
		MPI_Comm_rank(comm_, &rank);
		MPI_Comm_size(comm_, &comm_size);
		if (rank == 0) {
			reduced_vals_.resize(size_, 0);
			reduce_done_ = true;
		}
		if (comm_size > 1) {
			MPI_Reduce(list_.data(), reduced_vals_.data(), size_,
					MPI_UNSIGNED_LONG, MPI_SUM, 0, comm_);
		} else {
			reduced_vals_ = list_;
		}
	}
	std::ostream& stream_out(std::ostream& os) const{
		if (!gather_done_ && !reduce_done_) {
			//output own data
			int i = 0;
			std::vector<size_t>::const_iterator it = list_.begin();
			while (i < size_) {
				os << i << ',' << *it << std::endl;
				++i;
				++it;
			}
			return os;
		}
		int comm_size;
		MPI_Comm_size(comm_, &comm_size);
		if (gather_done_) {
			std::vector<size_t>::const_iterator it;
			int i = 0;
			while (i < size_) {
				os << i;
				it = gathered_vals_.begin() + i;
				for (int j = 0; j < comm_size; ++j) {
					os << ',' << *it;
					it += size_;
				}
				os << std::endl;
				++i;
			}
			os << std::endl;
		}
		if (reduce_done_) {
			double dcomm_size = static_cast<double>(comm_size);
			int i = 0;
			std::vector<size_t>::const_iterator it = reduced_vals_.begin();
			while (i < size_) {
				if (filter_ && *it != 0){
				os << i << ',' << static_cast<double>(*it) / dcomm_size
						<< std::endl;
				}
				++i;
				++it;
			}
			os << std::endl;
		}

		return os;
	}
private:
	DISALLOW_COPY_AND_ASSIGN(MPICounterList);
};



/** Enhanced "counter" that can be decremented (and may have a value < 0) and remembers its maximum */
template<typename T, typename D>
class MaxCounter {
public:
	explicit MaxCounter(D init = 0) :
			value_(init), max_(init) {
	}

	~MaxCounter() {
	}
	;
	void inc(D delta = 1) {
		if (std::numeric_limits<D>::max() - value_  <= delta){
			std::cerr << "value, delta = " << value_ << ',' << delta;
		}
		value_ += delta;
		if (value_ > max_){
			max_ = value_;
		}
	}
	void dec(long delta = 1){
		value_ -= delta;
		if (value_ > max_){
			max_ = value_;
		}
	}
	D get_value() {
		return value_;
	}
	void set(D val){
		value_ = val;
		if (value_>max_) {
			max_ = value_;
		}
	}
	long get_max(){
		return max_;
	}
	void reset() {
		value_ = 0;
		max_=0;
	}
	void gather() {
		static_cast<T*>(this)->gather_impl();
	}
	void reduce() {
		static_cast<T*>(this)->reduce_impl();
	}
	template<typename T, typename D> friend class MaxCounter;
	template<typename T, typename D> friend std::ostream& operator<<(std::ostream& os,
			const MaxCounter<T,D>& obj);
protected:
	D value_;
	D max_;
	std::ostream& stream_out(std::ostream& os) const {
		os << value_;
		return os;
	}
	DISALLOW_COPY_AND_ASSIGN(MaxCounter);
};

template<typename T, typename D>
std::ostream& operator<<(std::ostream& os, const MaxCounter<T,D>& obj) {
	return static_cast<const T&>(obj).stream_out(os);
}

//class MaxCounter {
//public:
//	MaxCounter() :
//			value_(0), max_(0) {
//	}
//	virtual ~MaxCounter() {
//	}
//	;
//	void inc(long delta = 1) {
//		value_ += delta;
//		max_ = (value_ > max_) ? value_ : max_;
//	}
//	void set(long value) {
//		value_ = value;
//		max_ = (value_ > max_) ? value_ : max_;
//	}
//	long get_value() {
//		return value_;
//	}
//	long get_max() {
//		return max_;
//	}
//	void reset() {
//		value_ = 0;
//		max_ = 0;
//	}
//	friend std::ostream& operator<<(std::ostream& os, const MaxCounter& obj) {
//		return obj.stream_out(os);
//	}
//protected:
//	long value_;
//	long max_;
//	virtual std::ostream& stream_out(std::ostream& os) const {
//		os << value_ << ',' << max_;
//		return os;
//	}
//	DISALLOW_COPY_AND_ASSIGN(MaxCounter);
//};


class MPIMaxCounter: public MaxCounter<MPIMaxCounter,long long> {
public:
	explicit MPIMaxCounter(const MPI_Comm& comm, long long init=0) : MaxCounter(init),
			comm_(comm), gathered_value_(0), reduced_value_(0), gathered_max_(0), reduced_max_(0)
	, gather_done_(false), reduce_done_(false){
	}
	~MPIMaxCounter() {
	}

	template<typename T, typename D> friend class MaxCounter;
	template<typename T, typename D> friend std::ostream& operator<<(std::ostream& os,
			const MaxCounter<T,D>& obj);
protected:

	const MPI_Comm& comm_;
	std::vector<long long> gathered_value_;
	long long reduced_value_;
	std::vector<long long> gathered_max_;
	long long reduced_max_;
	bool gather_done_;
	bool reduce_done_;

	void gather_impl(){
		int rank;
		int comm_size;
		MPI_Comm_rank(comm_, &rank);
		MPI_Comm_size(comm_, &comm_size);
		if (rank == 0){
			gathered_value_.resize(comm_size,0);
			gathered_max_.resize(comm_size,0);
			gather_done_=true;
		}
		if(comm_size>1){
			MPI_Gather(&value_, 1, MPI_LONG_LONG, gathered_value_.data(), 1, MPI_LONG_LONG, 0, comm_);
			MPI_Gather(&max_, 1, MPI_LONG_LONG, gathered_max_.data(), 1, MPI_LONG_LONG, 0, comm_);
		}
		else{
			gathered_value_[0] = value_;
			gathered_max_[0] = max_;
		}
	}
	void reduce_impl(){
		int rank;
		int comm_size;
		MPI_Comm_rank(comm_, &rank);
		MPI_Comm_size(comm_, &comm_size);
		if (rank == 0){
			reduce_done_=true;
		}
		if(comm_size>1){
			MPI_Reduce(&value_, &reduced_value_, 1, MPI_LONG_LONG, MPI_SUM, 0, comm_);
			MPI_Reduce(&max_, &reduced_max_, 1, MPI_LONG_LONG, MPI_MAX, 0, comm_);
		}
		else{
			reduced_value_ = value_;
			reduced_max_ = max_;
		}
	}


	std::ostream& stream_out(std::ostream& os) const{
		if (!gather_done_ && !reduce_done_) {
			//print own data
			int rank;
			MPI_Comm_rank(comm_, &rank);
			os << "company rank," << rank << ", value_," << value_ << " max_," << max_ << std::endl;
			return os;
		}
		int comm_size;
		MPI_Comm_size(comm_, &comm_size);

		if (gather_done_) {
			//output gathered values in comma separated list
			os << "value_";
			std::vector<long long>::const_iterator it = gathered_value_.begin();
			for (int j = 0; j < comm_size; ++j) {
				os << ',' << *it;
				++it;
			}
			os << std::endl;
			os << "max_";
			std::vector<long long>::const_iterator itm = gathered_max_.begin();
			for (int j = 0; j < comm_size; ++j) {
				os << ',' << *itm;
				++itm;
			}
			os << std::endl;
		}
		if (reduce_done_) {
			os << std::endl << "reduced mean";
			os << static_cast<double>(reduced_value_)
							/ static_cast<double>(comm_size);
			os << "reduced max";
			os << reduced_max_;
			os << std::endl;
		}
		return os;
	}
private:
	DISALLOW_COPY_AND_ASSIGN(MPIMaxCounter);
};

//class MaxCounterList {
//public:
//	explicit MaxCounterList(size_t size) :
//			list_(size) {
//	}
//	virtual ~MaxCounterList() {
//	}
//	void inc(int index, size_t delta = 1) {
//		list_[index].inc(delta);
//	}
//	void set(int index, long value) {
//		list_[index].set(value);
//	}
//	size_t get_value(int index) {
//		return list_[index].get_value();
//	}
//	void reset() {
//		for (std::vector<MaxCounter>::iterator it = list_.begin();
//				it != list_.end(); ++it) {
//			it->reset();
//		}
//	}
//	friend std::ostream& operator<<(std::ostream& os,
//			const MaxCounterList& obj) {
//		return obj.stream_out(os);
//	}
//protected:
//	std::vector<MaxCounter> list_;
//	virtual std::ostream& stream_out(std::ostream& os) const;
//private:
//	DISALLOW_COPY_AND_ASSIGN(MaxCounterList);
//};
//
//class PMaxCounterList: public MaxCounterList {
//public:
//	PMaxCounterList(const MPI_Comm& comm, size_t size) :
//			MaxCounterList(size), comm_(comm), gathered_vals_(0), reduced_vals_(
//					0) {
//	}
//	virtual ~PMaxCounterList() {
//	}
//	void gather();
//	void reduce();
//protected:
//	const MPI_Comm& comm_;
//	std::vector<long> gathered_vals_;
//	std::vector<long> reduced_vals_;
//	virtual std::ostream& stream_out(std::ostream& os) const;
//private:
//	DISALLOW_COPY_AND_ASSIGN(PMaxCounterList);
//};

template<typename T, typename R> class Timer;

template<typename T, typename R>
class Timer {
public:
	Timer() :
			total_(R()), start_(R()), max_(R()), num_epochs_(0), on_(false) {
	}
	void start() {
		start_ = get_time();
		check(!on_, "starting timer already on");
		on_ = true;
	}
	R get_time() {
		return static_cast<T*>(this)->get_time_impl();
	}

	double diff(R start, R end) {
		return static_cast<T*>(this)->diff_impl(start, end);
	}
	void pause() {
		double elapsed = diff(start_, get_time());
		total_ += elapsed;
		if (elapsed > max_) {
			max_ = elapsed;
		}
		check(on_, "pausing time that is not on");
		on_ = false;
	}
	void gather() {
		static_cast<T*>(this)->gather_impl();
	}
	void reduce() {
		static_cast<T*>(this)->reduce_impl();
	}

	double get_mean(){
		return static_cast<T*>(this)->get_mean_impl();
	}

	double get_total() {
		return total_;
	}
	double get_max() {
		return max_;
	}
	unsigned long get_num_epochs() {
		return num_epochs_;
	}

	friend std::ostream& operator<<<>(std::ostream& os, const Timer<T, R>& obj);

protected:
	double total_;
	R start_;
	double max_;
	unsigned long num_epochs_;
	bool on_;
};

template<typename T, typename R>
std::ostream& operator<<(std::ostream& os, const Timer<T, R>& obj) {
	return static_cast<const T&>(obj).stream_out(os);
}

/**
 * This class uses CRTP (curiously recurring template pattern).
 *
 * T is the type of the subclass, R is the type returned by the get time method.
 * For MPI_Wtime, R is double, for LinuxTimer that uses time(), the type is time_t.
 */
template<typename T, typename R>
class TimerList {
public:
	TimerList(size_t size) :
			size_(size), total_(size, 0), start_(size, R()), max_(size,0), num_epochs_(
					size, 0), on_(size, false) {
	}
	void start(size_t index) {
		start_.at(index) = get_time();
		check(!on_.at(index), "starting timer already on");
		on_.at(index) = true;
	}

	R get_time() {
		return static_cast<T*>(this)->get_time_impl();
	}

	double diff(R start, R end) {
		return static_cast<T*>(this)->diff_impl(start, end);
	}

	void pause(size_t index) {
		double elapsed = diff(start_.at(index), get_time());
		total_.at(index) += elapsed;
		++num_epochs_.at(index);
		if (elapsed > max_.at(index)) {
			max_.at(index) = elapsed;
		}
		check(on_.at(index), "pausing time that is not on");
		on_.at(index) = false;
	}
	void inc(size_t index, double elapsed) {
		total_.at(index) += elapsed;
		++num_epochs_.at(index);
		if (elapsed > max_.at(index)) {
			max_.at(index) = elapsed;
		}
		check(!on_.at(index), "incrementing timer that is on");
	}

	void gather() {
		static_cast<T*>(this)->gather_impl();
	}
	void reduce() {
		static_cast<T*>(this)->reduce_impl();
	}

	void print_op_table_stats(std::ostream& os,
			const SipTables& sip_tables) const {
		static_cast<const T*>(this)->print_op_table_stats_impl(os, sip_tables);
	}

//	const std::vector<double>& get_total() {
//		return total_;
//	}
//	const std::vector<double>& get_max() {
//		return max_;
//	}
//	const std::vector<unsigned long> get_num_epochs() {
//		return num_epochs_;
//	}
//	size_t get_size() {
//		return size_;
//	}

	friend std::ostream& operator<<<>(std::ostream& os,
			const TimerList<T, R>& obj);

protected:
	std::vector<double> total_;
	std::vector<R> start_;
	std::vector<double> max_;
	std::vector<unsigned long> num_epochs_;
	std::vector<bool> on_;
	size_t size_;
};

template<typename T, typename R>
std::ostream& operator<<(std::ostream& os, const TimerList<T, R>& obj) {
	return static_cast<const T&>(obj).stream_out(os);
}

class LinuxTimerList: public TimerList<LinuxTimerList, std::time_t> {
public:
	LinuxTimerList(size_t size) :
			TimerList<LinuxTimerList, std::time_t>(size) {
	}
	~LinuxTimerList() {
	}
	std::time_t get_time_impl() {
		std::time_t now;
		time(&now);
		return now;
	}

	double diff_impl(std::time_t start, std::time_t end) {
		return std::difftime(end, start);
	}

	void gather_impl() {
	}
	void reduce_impl() {
	}

	std::ostream& stream_out(std::ostream& os) const {
		std::setprecision(30);
		os << "total_,  max_,  num_epochs_" << std::endl;
		std::vector<double>::const_iterator total_iter = total_.begin();
		;
		std::vector<double>::const_iterator max_iter = max_.begin();
		std::vector<unsigned long>::const_iterator epoch_iter =
				num_epochs_.begin();
		for (int i = 0; i < size_; ++i) {
			os << i << ',' << *(total_iter++) << ',' << *(max_iter++) << ','
					<< *(epoch_iter++) << std::endl;
		}
		return os;
	}

	void print_op_table_stats_impl(std::ostream& os,
			const SipTables& sip_tables) const {
		os << std::setprecision(30);
		os << "pc, line number, opcode, mean,  max,  num_epochs" << std::endl;
		std::vector<double>::const_iterator total_iter = total_.begin();
		;
		std::vector<double>::const_iterator max_iter = max_.begin();
		std::vector<unsigned long>::const_iterator epoch_iter =
				num_epochs_.begin();
		for (int i = 0; i < size_; ++i) {
			os << i << ',' << sip_tables.line_number(i) << ','
					<< sip_tables.opcode_name(i) << ',';
			double mean = *total_iter / *epoch_iter;
			os << mean << ',' << *max_iter << ',' << *epoch_iter;
			++total_iter;
			++max_iter;
			++epoch_iter;
		}
	}

};

class MPITimerList: public TimerList<MPITimerList, double> {
public:
	MPITimerList(const MPI_Comm& comm, size_t size) :
			TimerList<MPITimerList, double>(size), comm_(comm), gathered_total_(
					0), gathered_max_(0), gathered_num_epoch_(0), reduced_mean_(
					0), reduced_max_(0), gather_done_(false), reduce_done_(
					false) {
	}
	~MPITimerList() {
	}
	template<typename T, typename R> friend class TimerList;
	template<typename T, typename R> friend std::ostream& operator<<(
			std::ostream& os, const TimerList<T, R>& obj);

protected:
	double get_time_impl() {
		return MPI_Wtime();
	}
	double diff_impl(double start, double end) {
		return end - start;
	}

	const MPI_Comm& comm_;
	std::vector<double> gathered_total_;
	std::vector<double> gathered_max_;
	std::vector<unsigned long> gathered_num_epoch_;
	std::vector<double> reduced_mean_;
	std::vector<double> reduced_max_;
	std::vector<unsigned long> reduced_num_epoch_;
	bool gather_done_;
	bool reduce_done_;

	void gather_impl() {
		int rank; //this rank is relative to comm_
		int comm_size;
		MPI_Comm_rank(comm_, &rank);
		MPI_Comm_size(comm_, &comm_size);

		if (rank == 0) {
			size_t buffsize = size_ * comm_size;
			gathered_total_.resize(buffsize, 0.0);
			gathered_max_.resize(buffsize, 0.0);
			gathered_num_epoch_.resize(buffsize, 0.0);
			gather_done_ = true;
		}

		if (comm_size > 1) {
			MPI_Gather(total_.data(), size_, MPI_DOUBLE, gathered_total_.data(),
					size_, MPI_DOUBLE, 0, comm_);
			MPI_Gather(max_.data(), size_, MPI_DOUBLE, gathered_max_.data(),
					size_, MPI_DOUBLE, 0, comm_);
			MPI_Gather(num_epochs_.data(), size_, MPI_UNSIGNED_LONG,
					gathered_num_epoch_.data(), size_, MPI_UNSIGNED_LONG, 0,
					comm_);

		} else {
			if (rank == 0){
			gathered_total_ = total_;
			gathered_max_ = max_;
			gathered_num_epoch_ = num_epochs_;
			}
		}
	}

	void reduce_impl() {
		int rank;
		int comm_size;
		MPI_Comm_rank(comm_, &rank);
		MPI_Comm_size(comm_, &comm_size);

		if (rank == 0) {
			reduced_mean_.resize(size_, 0.0);
			reduced_max_.resize(size_, 0.0);
			reduced_num_epoch_.resize(size_, 0);
			reduce_done_=true;
		}
		//reduce total and calculate mean
		MPI_Reduce(total_.data(), reduced_mean_.data(), size_, MPI_DOUBLE,
				MPI_SUM, 0, comm_); //this is actually the total. need to divide by #epochs
		MPI_Reduce(max_.data(), reduced_max_.data(), size_, MPI_DOUBLE, MPI_MAX,
				0, comm_);
		MPI_Reduce(num_epochs_.data(), reduced_num_epoch_.data(), size_,
				MPI_UNSIGNED_LONG, MPI_SUM, 0, comm_);
		if (rank == 0) {
			std::vector<double>::iterator miter = reduced_mean_.begin();
			std::vector<unsigned long>::iterator niter =
					reduced_num_epoch_.begin();
			for (int j = 0; j < size_; ++j) {
				if (*niter > 0) {
					*miter = *miter / *niter;
				}
				++miter;
				++niter;
			}
		}
	}

	std::ostream& stream_out(std::ostream& os) const {

		if (!gather_done_ && !reduce_done_) {
			std::setprecision(30);
			//print own data
			int rank;
			MPI_Comm_rank(comm_, &rank);
			os << "rank=" << rank << std::endl;
			os << "total_  max_  num_epochs_" << std::endl;
			std::vector<double>::const_iterator total_iter = total_.begin();
			;
			std::vector<double>::const_iterator max_iter = max_.begin();
			std::vector<unsigned long>::const_iterator epoch_iter =
					num_epochs_.begin();
			for (int i = 0; i < size_; ++i) {
				os << i << ',' << *(total_iter++) << ',' << *(max_iter++) << ','
						<< *(epoch_iter++) << std::endl;
			}
			return os;
		}

		int comm_size;
		MPI_Comm_size(comm_, &comm_size);
		int num_double_vals = 2;  //total_ and max_
		int num_unsigned_long_vals = 1; //num_epochs_
		std::vector<double>::const_iterator it;

		if (gather_done_) {
			//output the gathered values in a comma separated list

			os << "total" << std::endl;
			;
			int i = 0;
			while (i < size_) {
				os << i;
				it = gathered_total_.begin() + i;
				for (int j = 0; j < comm_size; ++j) {
					os << ',' << *it;
					it += size_;
				}
				os << std::endl;
				++i;
			}
			//output the values for max
			os << "max" << std::endl;
			i = 0;
			while (i < size_) {
				os << i;
				it = gathered_max_.begin() + i;
				for (int j = 0; j < comm_size; ++j) {
					os << ',' << *it;
					it += size_;
				}
				os << std::endl;
				++i;
			}
		}

		if (reduce_done_) {
			std::vector<double>::const_iterator mean_iter =
					reduced_mean_.begin();
			std::vector<double>::const_iterator max_iter = reduced_max_.begin();
			std::vector<unsigned long>::const_iterator num_epoch_iter =
					reduced_num_epoch_.begin();
			os << "reduced values:  mean, max, num_epoch" << std::endl;
			for (int i = 0; i != size_; ++i) {
				os << i << ',' << *mean_iter << ',' << *max_iter << ','
						<< *num_epoch_iter++ << std::endl;
				++mean_iter;
				++max_iter;
				++num_epoch_iter;
			}
		}
		return os;
	}

	void print_op_table_stats_impl(std::ostream& os,
			const SipTables& sip_tables, bool filter = true) const {

		int comm_size;
		MPI_Comm_size(comm_, &comm_size);
		check(reduce_done_, "must call reduce before print_optable_stats");
		os << "pc, line number, opcode, mean,  max,  mean num_epochs" << std::endl;
		std::vector<double>::const_iterator mean_iter = reduced_mean_.begin();
		std::vector<double>::const_iterator max_iter = reduced_max_.begin();
		std::vector<unsigned long>::const_iterator num_epoch_iter =
				reduced_num_epoch_.begin();
		for (int i = 0; i != size_; ++i) {
			if (!filter || *num_epoch_iter > 0) {
				std::setprecision(30);
				os << i << ',' << sip_tables.line_number(i) << ','
						<< sip_tables.opcode_name(i) << ',';
				os << *mean_iter << ',' << *max_iter << ',';
				std::setprecision(3);
				os<< (double)(*num_epoch_iter) / comm_size;
				os << std::endl;
			}
			++mean_iter;
			++max_iter;
			++num_epoch_iter;
		}
	}

};

template<typename R>
class NoopTimerList {
public:
	NoopTimerList(size_t size) {}
	NoopTimerList(const MPI_Comm& comm, size_t size){}
	void start(size_t index) {
	}
	R get_time() {
		return 0;
	}

	double diff(time_t start, time_t end) {
		return 0;
	}

	void pause(size_t index) {
	}
	void inc(size_t index, double elapsed) {
	}

	double get_mean(){return 0.0;}

	void gather() {
	}
	void reduce() {
	}

	void print_op_table_stats(std::ostream& os,
			const SipTables& sip_tables) const {
	}

	const std::vector<double>& get_total() {
		return 0;
	}
	const std::vector<double>& get_max() {
		return 0;
	}
	const std::vector<unsigned long> get_num_epochs() {
		return 0;
	}
	size_t get_size() {
		return 0;
	}

	friend std::ostream& operator<<<>(std::ostream& os, const NoopTimerList<R>& obj);
};

template<typename T>
std::ostream& operator<<(std::ostream& os, const NoopTimerList<T>& obj) {
	return os;
}

class LinuxTimer: public Timer<LinuxTimer, time_t> {
public:
	LinuxTimer() {
	}
	~LinuxTimer() {
	}
protected:
	time_t get_time_impl() {
		time_t now;
		time(&now);
		return now;
	}

	double diff_impl(time_t start, time_t end) {
		return difftime(end, start);
	}
	double get_mean_impl(){return total_;}

	void gather_impl() {
	}
	void reduce_impl() {
	}
	std::ostream& stream_out(std::ostream& os) const {
		os << ", total_," << total_;
		os << ", max_," << max_;
		os << ", num_epochs_," << num_epochs_;
		os << std::endl;
		return os;
	}
};

class MPITimer: public Timer<MPITimer, double> {
public:
	MPITimer(const MPI_Comm& comm) :
			comm_(comm), gathered_total_(0), gathered_max_(0), gathered_num_epochs_(
					0), reduced_mean_(0), reduced_max_(0), gather_done_(false), reduce_done_(
					false) {
	}
	~MPITimer() {
	}

	template<typename T, typename R> friend class Timer;
	template<typename T, typename R> friend std::ostream& operator<<(
			std::ostream& os, const Timer<T, R>& obj);

protected:

	double get_time_impl() {
		return MPI_Wtime();
	}

	double diff_impl(double start, double end) {
		return end - start;
	}
	double get_mean_impl(){
		check(reduce_done_, "must call reduce before retrieving mean");
		return reduced_mean_;
	}
	const MPI_Comm& comm_;
	std::vector<double> gathered_total_;
	std::vector<double> gathered_max_;
	std::vector<unsigned long> gathered_num_epochs_;
	double reduced_mean_;
	double reduced_max_;
	bool gather_done_;
	bool reduce_done_;

	void gather_impl() {
		int rank;
		int comm_size;
		MPI_Comm_rank(comm_, &rank);
		MPI_Comm_size(comm_, &comm_size);

		//resize destination buffers
		if (rank == 0) {
			gathered_total_.resize(comm_size, 0.0);
			gathered_max_.resize(comm_size, 0.0);
			gathered_num_epochs_.resize(comm_size, 0);
			gather_done_ = true;
		}
		if (comm_size > 1) {
			MPI_Gather(&total_, 1, MPI_DOUBLE, gathered_total_.data(), 1,
					MPI_DOUBLE, 0, comm_);
			MPI_Gather(&max_, 1, MPI_DOUBLE, gathered_max_.data(), 1,
					MPI_DOUBLE, 0, comm_);
			MPI_Gather(&num_epochs_, 1, MPI_UNSIGNED_LONG,
					gathered_num_epochs_.data(), 1, MPI_UNSIGNED_LONG, 0,
					comm_);
		} else {
			gathered_total_[0] = total_;
			gathered_max_[0] = max_;
			gathered_num_epochs_[0] = num_epochs_;
		}

	}
	void reduce_impl() {
		int rank;
		int comm_size;
		MPI_Comm_rank(comm_, &rank);
		MPI_Comm_size(comm_, &comm_size);

		if (comm_size > 1) {
			//reduce total
			MPI_Reduce(&total_, &reduced_mean_, 1, MPI_DOUBLE, MPI_SUM, 0,
					comm_);
			//reduce max
			MPI_Reduce(&max_, &reduced_max_, 1, MPI_DOUBLE, MPI_MAX, 0, comm_);
		} else {
			reduced_mean_ = total_;
			reduced_max_ = max_;
		}
		if (rank == 0) {
			reduce_done_ = true;
			reduced_mean_ /= comm_size; //calculate mean
		}
	}

	std::ostream& stream_out(std::ostream& os) const {
		std::setprecision(30);
		if (!gather_done_ && !reduce_done_) {
			//print own data
			int rank;
			MPI_Comm_rank(comm_, &rank);
			os << "rank=" << rank;
			os << ", total_=" << total_;
			os << ", max_=" << max_;
			os << ", num_epochs_=" << num_epochs_;
			os << std::endl;
			return os;
		}

		int comm_size;
		MPI_Comm_size(comm_, &comm_size);
		int num_double_vals = 2;  //total_ and max_
		int num_unsigned_long_vals = 1; //num_epochs_
		std::vector<double>::const_iterator it;

		if (gather_done_) {
			//output the gathered values in a comma separated list
			os << "total";
			it = gathered_total_.begin();
			for (int j = 0; j < comm_size; ++j) {
				os << ',' << *it;
				++it;
			}
			os << std::endl;

			//output the values for max
			os << "max";
			it = gathered_max_.begin();
			for (int j = 0; j < comm_size; ++j) {
				os << ',' << *it;
				++it;
			}
			os << std::endl;
		}
		if (reduce_done_) {
			os << "mean," << reduced_mean_ << std::endl;
			os << "max epoch," << reduced_max_ << std::endl;
		}
		return os;
	}
private:
	DISALLOW_COPY_AND_ASSIGN(MPITimer);
};

} /* namespace sip */

#endif /* COUNTER_H_ */
