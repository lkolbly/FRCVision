#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <sys/types.h> 
#include <winsock.h>
#include <Winsock2.h>
#include <vector>
#include "server.hxx"

class Client
{
private:
	int m_fd;

public:
	Client();
	Client(int fd);
	
	int update(void);
	int getFD(void);
};

Client::Client()
{
	m_fd = -1;
}

Client::Client(int fd)
{
	m_fd = fd;
}

int Client::update(void)
{
	if (m_fd == -1) {
		return -1;
	}
	printf("I have a m_fd! It's %i!\n", m_fd);
	char buf[1024];
	int nbytes = recv(m_fd, buf, 1024, 0);
	if (-1 == nbytes || 0 == nbytes) {
		printf("They've disconnected from us!\n");
		return 1;
	} else {
		printf("They sent us %i bytes!\n", nbytes);
	}
	return 0;
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
			Client *c = new Client(new_fd);
			clients.push_back(c);
		}

		for (std::vector<Client*>::iterator it=clients.begin(); it!=clients.end(); ++it) {
			if (FD_ISSET((*it)->getFD(),&rfs)) {
				// Looks like we found someone!
				printf("Client can talk.\n");
				if ((*it)->update()) {
					clients.erase(it);
					break;
				}
			}
		}
	}

#if 0

	// Socket code courtesy of linuxhowtos.org, since I didn't feel like pasting together the code ;)
	/*
    int sockfd, newsockfd, portno;
    socklen_t clilen;
    char buffer[256];
    struct sockaddr_in serv_addr, cli_addr;
    int n;
	*/


    sockfd = socket(AF_INET, SOCK_STREAM, 0);
    if (sockfd < 0) 
       error("ERROR opening socket");
    bzero((char *) &serv_addr, sizeof(serv_addr));
    portno = 4180;
    serv_addr.sin_family = AF_INET;
    serv_addr.sin_addr.s_addr = INADDR_ANY;
    serv_addr.sin_port = htons(portno);
    if (bind(sockfd, (struct sockaddr *) &serv_addr,
             sizeof(serv_addr)) < 0) 
             error("ERROR on binding");
    listen(sockfd,5);


	while (1) {
		/*
		int e = epoll_create1(0);
		struct epoll_event ev, events[5];
		ev.events = EPOLLIN;
		ev.data.fd = sockfd;
		epoll_ctl(e, EPOLL_CTL_ADD, sockfd, &ev); // WARNING: Could equal -1! If so, return error message with perror.

		// Add the clients...
		for (std::set<int,Client>::iterator it=clients.begin(); it!=clients.end(); ++it) {
			ev.data.fd = (*it).getFD();
			epoll_ctl(e, EPOLL_CTL_ADD, (*it).getFD(), &ev);
		}

		// Now we wait, and respond...
		int nfds = epoll_wait(e, events, 5, -1);
		for (int i=0; i<nfds; i++) {
			if (events[i].data.fd == sockfd) {
				// Accept a new connection...
				int new_fd = accept(sockfd, (struct sockaddr*)&cli_addr, &clilen);
				Client c = new Client(new_fd);
				clients[new_fd] = c;
			} else {
				// It's a client!
				clients[events[i].data.fd].update();
			}
		}
		*/


		/*
		clilen = sizeof(cli_addr);
		newsockfd = accept(sockfd, (struct sockaddr *) &cli_addr, &clilen);
		if (newsockfd < 0) 
			error("ERROR on accept");
		bzero(buffer,256);
		n = read(newsockfd,buffer,255);
		if (n < 0) error("ERROR reading from socket");
		printf("Here is the message: %s\n",buffer);
		n = write(newsockfd,"I got your message",18);
		if (n < 0) error("ERROR writing to socket");
		*/
	}


    close(newsockfd);
    close(sockfd);
#endif
    return 0;
}
