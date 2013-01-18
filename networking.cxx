#include <stdlib.h>
#include <stdio.h>
#include <curl/curl.h>
#include <pthread.h>
#include <windows.h>
#include "networking.hxx"

size_t write_func(void *buffer, size_t size, size_t nmemb, void *userp)
{
	//printf("Got %i bytes/%i memb.\n", size, nmemb);
	fwrite((unsigned char *)buffer, size, nmemb, *(FILE**)userp);
	return size*nmemb;
}

void *networkMain(void *arg)
{
	threadData_t *td = (threadData_t*)arg;
	printf("Starting networking thread.\n");
	curl_global_init(CURL_GLOBAL_DEFAULT);

	int last_downloaded = -1;
	
	CURL *c = curl_easy_init();
	curl_easy_setopt(c, CURLOPT_URL, "http://pillow.rscheme.org/robotics.jpg");
	curl_easy_setopt(c, CURLOPT_WRITEFUNCTION, write_func);

	while (1) {
		pthread_mutex_lock(&td->mutex);
		//printf("Networking: %i/%i\n", last_downloaded, td->var);
		if (last_downloaded < td->var) {
			printf("Downloading file.\n");
			FILE *f = fopen("out.jpg", "wb");
			if (!f) {
				fprintf(stderr, "We couldn't open out.jpg...\n");
				return NULL;
			}
			curl_easy_setopt(c, CURLOPT_WRITEDATA, &f);
			int success = curl_easy_perform(c);
			printf("Success was %s\n", curl_easy_strerror((CURLcode)success));
			fclose(f);
			last_downloaded = ++td->var;
		}
		pthread_mutex_unlock(&td->mutex);
		//Sleep(500);
	}

	curl_easy_cleanup(c);

	curl_global_cleanup();
	return 0;
}
