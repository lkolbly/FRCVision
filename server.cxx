#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <sys/types.h> 
#include <winsock.h>
#include <Winsock2.h>
#include <vector>
#include <opencv2/opencv.hpp>
#include "server.hxx"
#include "semaphores.hxx"
#include "processing.hxx"

#if 0
// compute sum of positive matrix elements
// (assuming that M is double-precision matrix)
double sum=0;
for(int i = 0; i < M.rows; i++)
{
    const double* Mi = M.ptr<double>(i);
    for(int j = 0; j < M.cols; j++)
        sum += std::max(Mi[j], 0.);
}
#endif

class Client
{
private:
	threadData_t *m_td;
	int m_fd;
	int m_conn_stage; // 0=Needs handshake, 1=has handshake
	int m_waiting_for_image;

	unsigned char *m_databuf;
	int m_databuf_len;
	
	int process_bytes(const unsigned char *buf, int nbytes);
	int send_raw_packet(unsigned int packet_type, const unsigned char *buf, int buflen);

public:
	Client();
	Client(int fd, threadData_t *td);

	int data_receivable(void);
	int update(void);
	int getFD(void);
};

Client::Client()
{
	m_fd = -1;
}

Client::Client(int fd, threadData_t *td)
{
	m_td = td;
	m_fd = fd;
	m_conn_stage = 0;
	m_databuf_len = 0;
	m_databuf = (unsigned char *)malloc(1);
}

int Client::send_raw_packet(unsigned int packet_type, const unsigned char *buf, int buflen)
{
	unsigned int l = buflen;
	send(m_fd, (const char *)&l, 4, 0);
	send(m_fd, (const char *)&packet_type, 4, 0);
	send(m_fd, (const char *)buf, buflen, 0);
	return 8+buflen;
}

// Returns the # of bytes processed
int Client::process_bytes(const unsigned char *buf, int nbytes)
{
	if (nbytes < 8) {
		return 0;
	}
	
	unsigned int len;
	memcpy(&len, buf, 4);

	unsigned int packet_type;
	memcpy(&packet_type, buf+4, 4);
	
	if (nbytes < len+8) {
		return 0;
	}

	int camera_ip;
	unsigned char outgoing[256];
	switch (packet_type) {
	case 0x00000001:
		// TODO: Actually process something here...
		outgoing[0] = 0x00;
		outgoing[1] = 0x00;
		send_raw_packet(0x80000001, outgoing, 2);
		return 12;
	case 0x01000003:
		// Get the camera IP address...
		memcpy(&camera_ip, buf+8, 4);

		// Set the outgoing packet
		memset(outgoing, 0, 6);
		send_raw_packet(0x81000003, outgoing, 6);
		return 12;
	case 0x02000001:
		// Oh boy...

#if 0 // This code will send stuff....
		double sum=0;
for(int i = 0; i < M.rows; i++)
{
    const double* Mi = M.ptr<double>(i);
    for(int j = 0; j < M.cols; j++)
        sum += std::max(Mi[j], 0.);
}
#endif
		m_waiting_for_image = 1;

		return 11;
	default:
		printf("Unknown packet!\n");
		return len+8; // Skip unknown packets...
	}

	return 0;
}

// Update is called when we're in fastloop mode (i.e. waiting for a mutex, or something else not selectable)
int Client::update(void)
{
	if (m_waiting_for_image) {
		if (pthread_mutex_trylock(&m_td->processed_data_lock) == 0) {
			
#if 0
			cv::Mat M = m_td->processing_result.img_data;
			unsigned char *outgoing = (unsigned char*)malloc(14+M.rows*M.cols);
			outgoing[0] = 0; // Camera ID
			long timestamp = 0;
			memcpy(outgoing+1, &timestamp, 8);
			short w=M.cols, h=M.rows;
			memcpy(outgoing+9, &w, 2);
			memcpy(outgoing+11, &h, 2);
			unsigned char type = 0x00;
			memcpy(outgoing+13, &type, 1);
			
			for(int i = 0; i < M.rows; i++)
			{
				const double* Mi = M.ptr<double>(i);
				for(int j = 0; j < M.cols; j++) {
					//sum += std::max(Mi[j], 0.);
					char pixel = Mi[j];
					outgoing[i*M.cols+j+14] = pixel;
				}
			}
#endif

			unsigned int datalen;
			unsigned char *contour_data = m_td->processing_result.render_contours(datalen);
			
			unsigned char *outgoing = (unsigned char*)malloc(10+datalen);
			short subframe_id = 0x0001;
			memcpy(outgoing, &subframe_id, 2);
			memcpy(outgoing+2, contour_data, datalen);
			
			send_raw_packet(0x82000001, outgoing, datalen+2);
			
			pthread_mutex_unlock(&m_td->processed_data_lock);
			m_waiting_for_image = 0;
		}
	}
	return 0;
}

// data_receivable is called when we know that data can be read.
int Client::data_receivable(void)
{
	if (m_fd == -1) {
		return 1;
	}
	printf("I have a m_fd! It's %i!\n", m_fd);
	char buf[1024];
	int nbytes = recv(m_fd, buf, 1024, 0);
	if (-1 == nbytes || 0 == nbytes) {
		printf("They've disconnected from us!\n");
		return 1;
	} else {
		printf("They sent us %i bytes!\n", nbytes);
		m_databuf = (unsigned char *)realloc((void*)m_databuf, m_databuf_len+nbytes);
		memcpy(m_databuf+m_databuf_len, buf, nbytes);
		m_databuf_len += nbytes;
		int processed = process_bytes(m_databuf, m_databuf_len);
		memmove(m_databuf, m_databuf+processed, m_databuf_len-processed);
		m_databuf_len -= processed;
	}
	if (m_waiting_for_image) {
		return 2;
	}
	return 0; // 0 is "OK". 2 is "We're waiting on an external event."
}

int Client::getFD(void)
{
	return m_fd;
}

void serverError(const char *msg)
{
    perror(msg);
    exit(1);
}

void *serverMain(void *arg)
{
	threadData_t *td = (threadData_t*)arg;
	std::vector<Client*> clients;

    int sockfd, newsockfd, portno;
    int clilen;
	// int socklen_
    char buffer[256];
    struct sockaddr_in serv_addr, cli_addr;
    int n;

    sockfd = socket(AF_INET, SOCK_STREAM, 0);
    // if (sockfd < 0) 
       // error("ERROR opening socket");
    // bzero((char *) &serv_addr, sizeof(serv_addr));
    portno = 4180;
    serv_addr.sin_family = AF_INET;
    serv_addr.sin_addr.s_addr = INADDR_ANY;
    serv_addr.sin_port = htons(portno);
    // if (bind(sockfd, (struct sockaddr *) &serv_addr,
             // sizeof(serv_addr)) < 0) 
             // error("ERROR on binding");
	bind(sockfd, (struct sockaddr *)&serv_addr, sizeof(serv_addr));
    listen(sockfd,5);
	
	int enable_fastloop = 0;
	while (1) {
		fd_set rfs;
		FD_ZERO(&rfs);
		FD_SET(sockfd, &rfs);

		for (std::vector<Client*>::iterator it=clients.begin(); it!=clients.end(); ++it) {
			FD_SET((*it)->getFD(),&rfs);
		}

		printf("Selecting...\n");
		int rval = select(1+clients.size(), &rfs, NULL, NULL, NULL);
		printf("%i\n", rval);
		if (rval == -1) {
			perror("select");
			printf("WSA: %i\n", WSAGetLastError());
			exit(0);
		}

		if (FD_ISSET(sockfd, &rfs)) {
			// There's someone new wanting a connection!
			printf("There's a new connection!\n");
			int new_fd = accept(sockfd, (struct sockaddr*)&cli_addr, &clilen);
			Client *c = new Client(new_fd, td);
			clients.push_back(c);
		}

		for (std::vector<Client*>::iterator it=clients.begin(); it!=clients.end(); ++it) {
			if (FD_ISSET((*it)->getFD(),&rfs)) {
				// Looks like we found someone!
				printf("Client can talk.\n");
				int retval = (*it)->data_receivable();
				if (retval == 1) {
					clients.erase(it);
					break;
				} else if (retval == 2) {
					enable_fastloop = 1;
				}
			}
			if (enable_fastloop) {
				int retval = (*it)->update();
				if (retval == 1) {
					clients.erase(it);
					break;
				} else if (retval == 2) {
					enable_fastloop = 1;
				}
			}
		}
	}

    return 0;
}
