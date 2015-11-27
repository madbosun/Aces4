/*
 * block_shape.cpp
 *
 *  Created on: Apr 6, 2014
 *      Author: basbas
 */

#include <block_shape.h>
#include <iostream>
#include <limits>

using namespace std::rel_ops;

namespace sip {



BlockShape::BlockShape() {
	std::fill(segment_sizes_ + 0, segment_sizes_ + MAX_RANK, 1);
}

//BlockShape::BlockShape(const segment_size_array_t& segment_sizes) {
//	std::copy(segment_sizes + 0, segment_sizes + MAX_RANK, segment_sizes_ + 0);
//}

BlockShape::BlockShape(const segment_size_array_t& segment_sizes, int rank){
	std::copy(segment_sizes + 0, segment_sizes + rank, segment_sizes_+0);
	std::fill(segment_sizes_ + rank, segment_sizes_ + MAX_RANK, 1);
}

BlockShape::~BlockShape() {
}

std::ostream& operator<<(std::ostream& os, const BlockShape & shape) {
	os << "[";
	for (int i = 0; i < MAX_RANK; ++i) {
		os << (i == 0 ? "" : ",") << shape.segment_sizes_[i];
	}
	os << ']';

	return os;
}

/**
 * This method returns an int, but computes the block shape internally as
 * size_t and checks that it is in range.
 * TODO  This test should be done with a tool offline.
 * @return
 */
int BlockShape::num_elems() const{
	size_t num_elems = 1;
	for (int i = 0; i < MAX_RANK; i++) {
		num_elems *= segment_sizes_[i];
	}
	check (num_elems < std::numeric_limits<int>::max(), "block size too large");
	return num_elems;
}


bool BlockShape::operator==(const BlockShape& rhs) const {
	bool is_equal = true;
	for (int i = 0; is_equal && i < MAX_RANK; ++i) {
		is_equal = (segment_sizes_[i] == rhs.segment_sizes_[i]);
	}
	return is_equal;
}
bool BlockShape::operator<(const BlockShape& rhs) const {
	bool is_eq = true;
	bool is_leq = true;
	for (int i = 0; is_leq && i < MAX_RANK; ++i) {
		is_leq = (segment_sizes_[i] <= rhs.segment_sizes_[i]);
		is_eq = is_eq && (segment_sizes_[i] == rhs.segment_sizes_[i]);
	}
	return (is_leq && !is_eq);
}



} /* namespace sip */
