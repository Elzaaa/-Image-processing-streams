#include "SimpleTimer.h"

SimpleTimer::SimpleTimer() {
	start = std::chrono::high_resolution_clock::now();

}


SimpleTimer::~SimpleTimer() {
}