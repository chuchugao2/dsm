#ifndef UTILS_GLOBALS
#define UTILS_GLOBALS

#include <atomic>
#include "sys/types.h"

extern std::atomic<bool> reach_time_limit;
extern long long windowSize;
extern uint ts;
extern uint k;
extern float mw;

#endif //UTILS_GLOBALS
