#include <chrono>
#include <iostream>

#define TIME

#if defined(TIME)

#define TIME_RANGE_TYPE

#define CLOCK \
std::chrono::high_resolution_clock::now()

#define CLOCK_START \
std::chrono::time_point<std::chrono::high_resolution_clock> __begin; \
do {__begin = std::chrono::high_resolution_clock::now();} while (false)

#define CLOCK_END(message) \
{std::cout << message << '\t' << std::chrono::duration<double>(std::chrono::high_resolution_clock::now() - __begin).count() \
<< std::endl;} while (false)

#define CLOCK_END_TIME \
std::chrono::duration<double>(CLOCK - __begin).count()

#define CLOCK_START_N do {__begin = std::chrono::high_resolution_clock::now();} while (false)

#define PRINT_TIME(message, time) \
{std::cout << message << '\t' << std::chrono::duration<double>(time).count() \
<< std::endl;} while (false)

#else
#define CLOCK
#define CLOCK_START
#define CLOCK_END(mes)
#define CLOCK_END_TIME
#define CLOCK_START_N
#define CLOCK_END_N(mes)
#endif