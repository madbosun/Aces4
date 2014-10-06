/*
 * profile_timer_store.h
 *
 *  Created on: Oct 4, 2014
 *      Author: jindal
 */

#ifndef PROFILE_TIMER_STORE_H_
#define PROFILE_TIMER_STORE_H_

#include <sqlite3.h>
#include <string>
#include <utility>

#include "profile_timer.h"

namespace sip {

class ProfileTimerStore {
public:
	ProfileTimerStore(const std::string &db_name);
	~ProfileTimerStore();

	/**
	 * Saves to a persistent database, the profiled time for an op with
	 * a given size and shape of its operands.
	 *
	 * @param opcode_operands opcode and its operand blocks' shapes & sizes
	 * @param time_count_pair pair of total time & number of times opcode was called
	 */
	void save_to_store(const ProfileTimer::Key& opcode_operands, const std::pair<double, int>& time_count_pair);

	/**
	 * Retrieves from a persistent database, the profiled total time & count
	 * for a given opcode and its operands.
	 * @param opcode_operands
	 * @return
	 */
	std::pair<double, int> get_from_store(const ProfileTimer::Key& opcode_operands);

	const static int MAX_BLOCK_OPERANDS;

private:
	std::string db_location_;
	sqlite3* db_;

	/**
	 * Utility function to print the sqlite3 error and throw an exception.
	 * @param rc
	 */
	void sip_sqlite3_error(int rc);

	/**
	 * Utility method to create a table for a given number of blocks
	 * as operands. Called by the constructor
	 * Tables are created for zero block operands, 1 block operand ... upto 5 block operands.
	 * @param num_blocks
	 */
	void create_table(int num_blocks);




	// SQLITE3 database table, columns names & prefixes.

	// TABLE Names
	const static std::string block_table[];

	// CHARACTER Columns
	const static std::string opcode_column;
	const static int max_opcode_size;

	// INTEGER Columns
	const static std::string block_prefix;
	const static std::string indices_prefix;
	const static std::string segment_prefix;

	// DOUBLE PRECISION COLUMN
	const static std::string tottime_column;

	// BIGINT COLUMN
	const static std::string count_column;

};

} /* namespace sip */

#endif /* PROFILE_TIMER_STORE_H_ */
