#pragma once
#include <chrono>
#include <iostream>

class SimpleTimer
{
public:
	bool stop=false;
	SimpleTimer();
	float GetTime() {
		if (!stop)
			end = std::chrono::high_resolution_clock::now();
		duration = end - start;
		return duration.count();
	}
	void StopTime() {
		stop = true;
		end = std::chrono::high_resolution_clock::now();
	}

	~SimpleTimer() ;
private:
	std::chrono::time_point<std::chrono::steady_clock> start, end;
	std::chrono::duration<float> duration;
	float result;
};