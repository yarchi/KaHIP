#include <chrono>
#include <iostream>

#define TIME

#if defined(TIME)

#define CLOCK_START \
std::chrono::time_point<std::chrono::high_resolution_clock> __begin; \
do {__begin = std::chrono::high_resolution_clock::now();} while (false)

#define CLOCK_END(message) \
std::chrono::time_point<std::chrono::high_resolution_clock> __end \
= std::chrono::high_resolution_clock::now(); \
{std::cout << message << '\t' << std::chrono::duration_cast<std::chrono::microseconds>(__end - __begin).count() \
<< std::endl;} while (false)

#define CLOCK_START_N do {__begin = std::chrono::high_resolution_clock::now();} while (false)

#define CLOCK_END_N(message) \
__end = std::chrono::high_resolution_clock::now(); \
{std::cout << message << '\t' << std::chrono::duration_cast<std::chrono::microseconds>(__end - __begin).count() \
<< std::endl;} while (false)

#else
#define CLOCK_START
#define CLOCK_END(mes)
#define CLOCK_START_N
#define CLOCK_END_N(mes)
#endif