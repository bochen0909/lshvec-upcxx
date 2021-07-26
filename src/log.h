/*
 * log.h
 *
 *  Created on: Jul 13, 2020
 *      Author:
 */

#ifndef SOURCE_DIRECTORY__SRC_SPARC_LOG_H_
#define SOURCE_DIRECTORY__SRC_SPARC_LOG_H_


void set_spdlog_pattern(const char* hostname, int rank);

int myinfo(const char *fmt, ...);

int mydebug(const char *fmt, ...);

int myerror(const char *fmt, ...);

int mywarn(const char *fmt, ...);


#endif /* SOURCE_DIRECTORY__SRC_SPARC_LOG_H_ */
