#ifndef NETWORKING_HXX
#define NETWORKING_HXX 1

#include <pthread.h>

void *networkMain(void *arg);

typedef struct threadData_t {
	pthread_mutex_t mutex;
	int var;
} threadData_t;

#endif
