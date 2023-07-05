//
// Created by ¸ß³þ³þ on 2023/7/5.
//

#include "Timer.h"
void Timer::StartTimer() {
    start_time=std::chrono::steady_clock::now();
}
void Timer::StopTimer() {
    std::chrono::steady_clock::time_point end_time=std::chrono::steady_clock::now();
    time_microseconds+=std::chrono::duration_cast<std::chrono::microseconds>(end_time-start_time);
}
long long Timer::GetTimer() {
    return time_microseconds.count();
}