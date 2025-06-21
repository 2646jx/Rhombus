#ifndef RHOMBUS_DEFINE_H
#define RHOMBUS_DEFINE_H

#include <cstdint>

#define CREATE_THREAD_POOL(nths, len, func) \
do { \
    uint32_t thread_block = (len + nths - 1) / nths; \
    std::vector<std::thread> thread_pool(nths);\
    for (size_t i = 0; i < nths; ++i) {\
        size_t bgn = i * thread_block; \
        size_t end = std::min<size_t>(bgn + thread_block, len); \
        thread_pool[i] = std::thread(func, bgn, end); \
    } \
    std::for_each(thread_pool.begin(), thread_pool.end(), [](std::thread &t){ t.join(); }); \
} while(0)

#endif //RHOMBUS_DEFINE_H