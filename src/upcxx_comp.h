/*
 * upcxx_comp.h
 *
 *  Created on: Jul 23, 2021
 *      Author: lizhen
 */

#ifndef SOURCE_DIRECTORY__SRC_UPCXX_UPCXX_COMP_H_
#define SOURCE_DIRECTORY__SRC_UPCXX_UPCXX_COMP_H_

#include <upcxx/upcxx.hpp>

#if UPCXX_SPEC_VERSION == 20210300L
#define upcxx_op_add upcxx::op_fast_add
#define upcxx_fatal_error upcxx::detail::fatal_error
#define upcxx_reduce_all upcxx::reduce_all
#define upcxx_broadcast_nontrivial upcxx::experimental::broadcast_nontrivial
#else
#define upcxx_op_add upcxx::op_add
#define upcxx_fatal_error upcxx::fatal_error
#define upcxx_broadcast_nontrivial upcxx::broadcast_nontrivial
#define upcxx_reduce_all upcxx::reduce_all
#endif




#endif /* SOURCE_DIRECTORY__SRC_UPCXX_UPCXX_COMP_H_ */
