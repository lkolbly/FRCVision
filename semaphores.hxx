#ifndef SEMAPHORES_HXX
#define SEMAPHORES_HXX

#include <pthread.h>

typedef struct threadData_t {
	pthread_mutex_t mutex;
	int var;
} threadData_t;

#endif
