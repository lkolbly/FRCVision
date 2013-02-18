#ifndef KILLSERVER_HXX
#define KILLSERVER_HXX 1

#include <pthread.h>

typedef struct {
	unsigned char *secret;
	int secret_len;
	int needs_death;
	pthread_mutex_t mutex;
} killServerArg_t;

void *killServerMain(void *p);

#endif
