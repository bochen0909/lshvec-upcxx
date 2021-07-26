/*
 * log.cpp
 *
 *  Created on: Jul 13, 2020
 *      Author:
 */

#include <sstream>
#include "spdlog/spdlog.h"
#include "log.h"

using namespace std;

void set_spdlog_pattern(const char* hostname, int rank) {

	string tmp = "[%Y-%m-%d %H:%M:%S.%e] [thread %t][%^%l%$]";
	stringstream ss;
	ss << "[%Y-%m-%d %H:%M:%S.%e]";
	ss << " [" << hostname;
	ss << " r:" << rank << " thr:%t]";
	ss << " [%^%l%$] %v";

	spdlog::set_pattern(ss.str().c_str());
}

int mydebug(const char *fmt, ...) {
	char buffer[4096];
	va_list args;
	va_start(args, fmt);
	int rc = vsnprintf(buffer, sizeof(buffer), fmt, args);
	va_end(args);
	spdlog::debug(buffer);
	return rc;
}


int myinfo(const char *fmt, ...) {
	char buffer[4096];
	va_list args;
	va_start(args, fmt);
	int rc = vsnprintf(buffer, sizeof(buffer), fmt, args);
	va_end(args);
	spdlog::info(buffer);
	return rc;
}

int myerror(const char *fmt, ...) {
	char buffer[4096];
	va_list args;
	va_start(args, fmt);
	int rc = vsnprintf(buffer, sizeof(buffer), fmt, args);
	va_end(args);
	spdlog::error(buffer);
	return rc;
}

int mywarn(const char *fmt, ...) {
	char buffer[4096];
	va_list args;
	va_start(args, fmt);
	int rc = vsnprintf(buffer, sizeof(buffer), fmt, args);
	va_end(args);
	spdlog::warn(buffer);
	return rc;
}

