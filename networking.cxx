#include <stdlib.h>
#include <stdio.h>
#include <curl/curl.h>
#include "networking.hxx"

size_t write_func(void *buffer, size_t size, size_t nmemb, void *userp)
{
	printf("Got %i bytes/%i memb.\n", size, nmemb);
	fwrite((unsigned char *)buffer, size, nmemb, *(FILE**)userp);
	return size*nmemb;
}

void *networkMain(void *arg)
{
	printf("Starting networking thread.\n");
	curl_global_init(CURL_GLOBAL_DEFAULT);

	CURL *c = curl_easy_init();
	curl_easy_setopt(c, CURLOPT_URL, "http://pillow.rscheme.org/robotics.jpg");
	curl_easy_setopt(c, CURLOPT_WRITEFUNCTION, write_func);
	FILE *f = fopen("out.jpg", "wb");
	if (!f) {
		fprintf(stderr, "We couldn't open out.jpg...\n");
		return NULL;
	}
	curl_easy_setopt(c, CURLOPT_WRITEDATA, &f);
	int success = curl_easy_perform(c);
	printf("Success was %s\n", curl_easy_strerror((CURLcode)success));
	fclose(f);
	curl_easy_cleanup(c);

	curl_global_cleanup();
	return 0;
}
