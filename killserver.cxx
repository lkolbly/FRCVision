#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <sys/types.h> 
#include <winsock.h>
#include <Winsock2.h>
#include <vector>
#include <pthread.h>
#include "killserver.hxx"

void *killServerMain(void *p)
{
	//unsigned char *secret = (unsigned char *)arg;
	killServerArg_t *arg = (killServerArg_t*)p;

    int sockfd, portno;
	// int socklen_
    //char buffer[256];
    struct sockaddr_in serv_addr, cli_addr;
    int clilen = sizeof(cli_addr);
	
    sockfd = socket(AF_INET, SOCK_STREAM, 0);
    portno = 4181;
    serv_addr.sin_family = AF_INET;
    serv_addr.sin_addr.s_addr = INADDR_ANY;
    serv_addr.sin_port = htons(portno);
	bind(sockfd, (struct sockaddr *)&serv_addr, sizeof(serv_addr));
    listen(sockfd,5);

	while (1) {
		fd_set rfs;
		FD_ZERO(&rfs);
		FD_SET(sockfd, &rfs);
		
		struct timeval tv;
		tv.tv_sec = 0;
		tv.tv_usec = 0;

		//printf("Selecting...\n");
		int rval;
			rval = select(1, &rfs, NULL, NULL, NULL);
		//printf("%i\n", rval);
		if (rval == -1) {
			perror("select");
			printf("WSA: %i\n", WSAGetLastError());
			exit(0);
		}

		if (FD_ISSET(sockfd, &rfs)) {
			// There's someone new wanting a connection!
			int new_fd = accept(sockfd, (struct sockaddr*)&cli_addr, &clilen);
			printf("There's a new connection %i!\n", new_fd);
			pthread_mutex_lock(&arg->mutex);
			arg->needs_death = 1;
			pthread_mutex_unlock(&arg->mutex);
			return (void*)1;

			// Check to see if they have the right MD5sum
			char buf[arg->secret_len];
			int nread = recv(new_fd, buf, 1, arg->secret_len);
			if (nread < arg->secret_len) {
				shutdown(new_fd, 0);
				continue;
			}
			for (int i=0; i<arg->secret_len; i++) {
				if (buf[i] != arg->secret[i]) {
					shutdown(new_fd, 0);
					continue;
				}
			}
			pthread_mutex_lock(&arg->mutex);
			arg->needs_death = 1;
			pthread_mutex_unlock(&arg->mutex);
			return (void*)1;
		}
	}

    return 0;
}
