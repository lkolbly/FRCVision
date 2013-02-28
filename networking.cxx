#include <stdlib.h>
#include <stdio.h>
#include <curl/curl.h>
#include <pthread.h>
#include <windows.h>
#include "networking.hxx"
#include "semaphores.hxx"

FILE *network_Log_File;

// Get the seconds since an arbitrary time in the past
double netGetTime(void)
{
	return (double)GetTickCount() / 1000.0;
}

class RateLimiter
{
private:
	int m_BPS; // Bytes per Second. Everything's bytes.
	list<int> m_t_bytes;
	list<double> m_t_times;

	double m_average_bytes; // Average bytes per transaction
	int m_n_transactions; // Total number of transactions
	int m_recent_bytes;

public:
	RateLimiter(int BPS);
	void clean(void); // Removes bytes we don't care about
	int get_bytes(void); // Gets the bytes in the past 1 second
	void delay(void);
	void add_transaction(int bytes);
};

RateLimiter::RateLimiter(int BPS)
{
	m_BPS = BPS;
	m_average_bytes = 0.0;
	m_n_transactions = 0;
	m_recent_bytes = 0;
}

void RateLimiter::clean(void)
{
	double cur_time = netGetTime();
	double t = m_t_times.front();
	while (cur_time-1.0 > t) {
		//printf("%f %f\n", cur_time, m_t_times.front());
		if (m_t_bytes.size() == 0) return;
		m_recent_bytes -= m_t_bytes.front();
		m_t_times.pop_front();
		m_t_bytes.pop_front();
		t = m_t_times.front();
	}
}

int RateLimiter::get_bytes(void)
{
	return m_recent_bytes;
#if 0
	int bytes = 0;
	for (size_t i=0; i<m_t_bytes.size(); i++) {
		bytes += m_t_bytes[i];
	}
	return bytes;
#endif
}

void RateLimiter::delay(void)
{
	while (1) {
		clean();
		//printf("%i %i %f\n", m_BPS, get_bytes(), m_average_bytes);
		if (m_BPS - get_bytes() > m_average_bytes) {
			return;
		}
	}
}

void RateLimiter::add_transaction(int bytes)
{
	m_recent_bytes += bytes;
	m_t_bytes.push_back(bytes);
	m_t_times.push_back(netGetTime());
	if (m_n_transactions == 0) {
		m_average_bytes = bytes;
	} else {
		m_average_bytes = m_average_bytes + (double)bytes / m_n_transactions;
	}
}

struct netFile_t {
	unsigned char *data;
	unsigned int data_len;
	unsigned int alloced_len;
};

struct netFile_t network_Tmp_File = {NULL,0,0};

size_t write_func(void *buffer, size_t size, size_t nmemb, void *userp)
{
	if (network_Tmp_File.alloced_len == 0) {
		network_Tmp_File.data = new unsigned char[28000];
		network_Tmp_File.alloced_len = 28000;
		network_Tmp_File.data_len = 0;
	}
	if (network_Tmp_File.alloced_len-5 < network_Tmp_File.data_len+size*nmemb) {
		network_Tmp_File.alloced_len = (network_Tmp_File.data_len+size*nmemb)*2;
		unsigned char *tmp = new unsigned char[network_Tmp_File.alloced_len];
		memcpy(tmp, network_Tmp_File.data, network_Tmp_File.data_len);
		delete network_Tmp_File.data;
		network_Tmp_File.data = tmp;
	}

	//printf("%i %i %i\n", network_Tmp_File.alloced_len, nmemb, network_Tmp_File.data_len);
	memcpy(network_Tmp_File.data+network_Tmp_File.data_len, (unsigned char*)buffer, size*nmemb);
	network_Tmp_File.data_len += size*nmemb;
	//printf("Got %i bytes => %i.\n", size*nmemb, network_Tmp_File.data_len);

	//printf("Got %i bytes/%i memb.\n", size, nmemb);
	//fwrite((unsigned char *)buffer, size, nmemb, *(FILE**)userp);
	return size*nmemb;
}

int networkingDownloadImage(CURL *c)
{
	int success = curl_easy_perform(c);
	if (success) {
		fprintf(network_Log_File, "Success was %i: %s\n", success, curl_easy_strerror((CURLcode)success));
	}
	return success;
}

void networkCloseLogFile(void)
{
	fclose(network_Log_File);
}

void *networkMain(void *arg)
{
	network_Log_File = fopen("network.log", "w");

	threadData_t *td = (threadData_t*)arg;
	printf("Starting networking thread.\n");
	curl_global_init(CURL_GLOBAL_DEFAULT);

	//int last_downloaded = -1;
	RateLimiter limiter(1000000);

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
	
	network_Tmp_File.alloced_len = 0;

	// TODO: Some sort of flow rate.
	// TODO: Multiple cameras.
	while (1) {
		pthread_mutex_lock(&td->network_heartbeat_mutex);
		if (td->time_to_die) {
			pthread_mutex_unlock(&td->network_heartbeat_mutex);
			break;
		}
		pthread_mutex_unlock(&td->network_heartbeat_mutex);

		limiter.delay(); // Delay the transaction until we're safe
		if (networkingDownloadImage(c)) {
			fprintf(stderr, "An error has occurred while attempting to access the HTTP camera at '%s'\n", camera_url);
			break;
		}

		// Move 'tmp.jpg' to the protected 'out.jpg'
		pthread_mutex_lock(&td->image_file_lock);
		FILE *f = fopen("storage-tmp.jpg", "wb");
		int nbytes = fwrite(network_Tmp_File.data, 1, network_Tmp_File.data_len, f);
		limiter.add_transaction(nbytes);
		fclose(f);
		td->has_processed_image = 0;
		td->collection_cfg.uid = rand();
		pthread_mutex_unlock(&td->image_file_lock);

		delete network_Tmp_File.data;
		network_Tmp_File.data = NULL;
		network_Tmp_File.alloced_len = 0;
		network_Tmp_File.data_len = 0;
	}

	curl_easy_cleanup(c);

	curl_global_cleanup();
	pthread_mutex_lock(&td->network_heartbeat_mutex);
	td->networking_is_dead = 1;
	pthread_mutex_unlock(&td->network_heartbeat_mutex);
	return (void*)1;
}
