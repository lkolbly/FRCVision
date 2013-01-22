#include <highgui.h>
#include "networking.hxx"
#include "processing.hxx"
#include "server.hxx"
#include "semaphores.hxx"
#include <pthread.h>
#include <stdio.h>

int main ( int argc, char **argv )
{
    WORD wVersionRequested;
    WSADATA wsaData;
    int err;

/* Use the MAKEWORD(lowbyte, highbyte) macro declared in Windef.h */
    wVersionRequested = MAKEWORD(2, 2);

    err = WSAStartup(wVersionRequested, &wsaData);
    if (err != 0) {
        /* Tell the user that we could not find a usable */
        /* Winsock DLL.                                  */
        printf("WSAStartup failed with error: %d\n", err);
        return 1;
    }


	pthread_mutex_t new_Image_Mutex;
	pthread_mutex_init(&new_Image_Mutex, NULL);
	threadData_t td;
	td.mutex = new_Image_Mutex;
	td.var = 0;
	pthread_mutex_init(&td.image_file_lock, NULL);
	pthread_mutex_init(&td.processed_data_lock, NULL);

	pthread_t networking_thread, processing_thread, server_thread;
	int rc;
	rc = pthread_create(&networking_thread, NULL, networkMain, &td);
	assert(0==rc);
	
	rc = pthread_create(&processing_thread, NULL, processingMain, &td);
	assert(0==rc);
	
	rc = pthread_create(&server_thread, NULL, serverMain, &td);
	assert(0==rc);
	
	rc = pthread_join(networking_thread, NULL);
	rc = pthread_join(processing_thread, NULL);
	rc = pthread_join(server_thread, NULL);
	return 0;

	networkMain(NULL);
	processingMain(NULL);
	return 0;

  cvNamedWindow( "My Window", 1 );
  IplImage *img = cvCreateImage( cvSize( 640, 480 ), IPL_DEPTH_8U, 1 );
  CvFont font;
  double hScale = 1.0;
  double vScale = 1.0;
  int lineWidth = 1;
  cvInitFont( &font, CV_FONT_HERSHEY_SIMPLEX | CV_FONT_ITALIC,
              hScale, vScale, 0, lineWidth );
  cvPutText( img, "Hello World!", cvPoint( 200, 400 ), &font,
             cvScalar( 255, 255, 0 ) );
  cvShowImage( "My Window", img );
  cvWaitKey();
  return 0;
}
