#include <stdlib.h>
#include <stdio.h>
#include <curl/curl.h>
#include <pthread.h>
#include <windows.h>
#include "networking.hxx"
#include "semaphores.hxx"

size_t write_func(void *buffer, size_t size, size_t nmemb, void *userp)
{
	//printf("Got %i bytes/%i memb.\n", size, nmemb);
	fwrite((unsigned char *)buffer, size, nmemb, *(FILE**)userp);
	return size*nmemb;
}

void networkingDownloadImage(CURL *c)
{
	//printf("Downloading file.\n");
	FILE *f = fopen("network-tmp.jpg", "wb");
	if (!f) {
		fprintf(stderr, "We couldn't open tmp.jpg...\n");
		return;
	}
	curl_easy_setopt(c, CURLOPT_WRITEDATA, &f);
	int success = curl_easy_perform(c);
	//printf("Success was %s\n", curl_easy_strerror((CURLcode)success));
	fclose(f);
	//Sleep(1000);
}

void *networkMain(void *arg)
{
	threadData_t *td = (threadData_t*)arg;
	printf("Starting networking thread.\n");
	curl_global_init(CURL_GLOBAL_DEFAULT);

	int last_downloaded = -1;
	
	CURL *c = curl_easy_init();
	//curl_easy_setopt(c, CURLOPT_URL, "http://pillow.rscheme.org/robotics.jpg");
	//curl_easy_setopt(c, CURLOPT_URL, "http://10.4.18.12/axis-cgi/jpg/image.cgi");
	char camera_url[1024];
	snprintf(camera_url, 1024, "http://%s/axis-cgi/jpg/image.cgi", td->collection_cfg.camera_hostname);
	printf("Using camera URL '%s'\n", camera_url);
	curl_easy_setopt(c, CURLOPT_URL, camera_url);
	curl_easy_setopt(c, CURLOPT_USERPWD, "FRC:FRC");
	curl_easy_setopt(c, CURLOPT_WRITEFUNCTION, write_func);

	// Set our initial configuration
	pthread_mutex_lock(&td->image_file_lock);
	td->collection_cfg.image_filename = strdup("new.jpg");
	pthread_mutex_unlock(&td->image_file_lock);

	// TODO: Some sort of flow rate.
	// TODO: Multiple cameras.
	while (1) {
		networkingDownloadImage(c);

		// Move 'tmp.jpg' to the protected 'out.jpg'
		pthread_mutex_lock(&td->image_file_lock);
		//printf("DOWNLOADING IMAGE.\n");
		CopyFile("network-tmp.jpg", "storage-tmp.jpg", false);
		Sleep(50);
		td->has_processed_image = 0;
		td->collection_cfg.uid = rand();
		pthread_mutex_unlock(&td->image_file_lock);
	}

	curl_easy_cleanup(c);

	curl_global_cleanup();
	return 0;
}
