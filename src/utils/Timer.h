//
// Created by ¸ß³þ³þ on 2023/7/5.
//

#ifndef BASELINE_TIMER_H
#define BASELINE_TIMER_H

#include <chrono>
#include <iostream>

class Timer {
private:
    std::chrono::steady_clock::time_point start_time;
    std::chrono::microseconds time_microseconds;
public:
    Timer() : time_microseconds(0) {};

    void StartTimer();

    void StopTimer();

    long long GetTimer();


};


#endif //BASELINE_TIMER_H
