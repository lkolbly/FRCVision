#include <highgui.h>
#include "networking.hxx"
#include "processing.hxx"
#include "server.hxx"
#include "killserver.hxx"
#include "semaphores.hxx"
#include <pthread.h>
#include <stdio.h>
#include <expat.h>

typedef struct parsingState_t {
	int is_in_tracking;
	
	int is_in_object;
	trackingObject_t cur_obj;
	
	int is_parsing_text;
	char text[1024];
	
	threadData_t *td;
} parsingState_t;

const char *getAttribute(const char *key, const char **attrs)
{
	for (int i=0; attrs[i]!=NULL; i+=2) {
		if (strcmp(attrs[i], key) == 0) {
			return attrs[i+1];
		}
	}
	return NULL;
}

static void XMLCALL
cfgStartElem(void *data, const char *el, const char **attr)
{
	parsingState_t *o = (parsingState_t*)data;
	if (strcmp(el, "tracking") == 0 && !o->is_in_tracking) {
		o->is_in_tracking = 1;
	} else if (strcmp(el, "object") == 0) {
		//o->cur_obj.type = strdup(getAttribute("type", attr));
		const char *t = getAttribute("type", attr);
		if (strcmp(t, "rectangle") == 0) {
			o->cur_obj.type = OBJ_HOLLOW_RECT;
		}
		o->is_in_object = 1;
	} else if (strcmp(el, "width") == 0 && o->is_in_object) {
		o->is_parsing_text = 1;
	} else if (strcmp(el, "height") == 0 && o->is_in_object) {
		o->is_parsing_text = 1;
	} else if (strcmp(el, "camera") == 0) {
		o->td->collection_cfg.camera_hostname = strdup(getAttribute("ip", attr));
		printf("The camera hostname is set to '%s'\n", o->td->collection_cfg.camera_hostname);
	}
	return;
}

static void XMLCALL
cfgTextElem(void *data, const char *s, int len)
{
	parsingState_t *o = (parsingState_t*)data;
	if (o->is_parsing_text) {
		memcpy(o->text+strlen(o->text), s, len);
	}
}

static void XMLCALL
cfgEndElem(void *data, const char *el, const char **attr)
{
	parsingState_t *o = (parsingState_t*)data;
	if (strcmp(el, "tracking") == 0 && o->is_in_tracking) {
		o->is_in_tracking = 0;
	} else if (strcmp(el, "object") == 0) {
		o->is_in_object = 0;
		o->td->processing_config.targets.push_back(o->cur_obj);
	} else if (strcmp(el, "width") == 0) {
		if (o->is_in_object) {
			o->cur_obj.width = atof(o->text);
		}
	} else if (strcmp(el, "height") == 0) {
		if (o->is_in_object) {
			o->cur_obj.height = atof(o->text);
		}
	}
	
	if (o->is_parsing_text) {
		o->is_parsing_text = 0;
		memset(o->text, 0, 1024);
	}
	return;
}

// This is a bit of a hackish way to do it, since WE DON'T HAVE CONFIG FILES (Michael...)
int loadConfigFiles(threadData_t *td)
{
	FILE *f = fopen("camera.xml", "rb");
	if (!f) {
		fprintf(stderr, "Couldn't open 'camera.xml'\n");
		return 1;
	}

	XML_Parser p = XML_ParserCreate(NULL);
	XML_SetElementHandler(p, cfgStartElem, (XML_EndElementHandler)cfgEndElem);
	XML_SetCharacterDataHandler(p, cfgTextElem);
	parsingState_t state;
	state.td = td;
	XML_SetUserData(p, &state);
	for (;;) {
		char buf[1024];
		int done = 0;
		(int)fread(buf, 1, 1024, f);
		if (feof(f)) {
			done = 1;
		}
		
		XML_Parse(p, buf, 1024, done);
		if (done) {
			break;
		}
	}
	XML_ParserFree(p);
	fclose(f);
	return 0;
}

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
	td.processing_result = NULL;
	td.mutex = new_Image_Mutex;
	td.var = 0;
	if (loadConfigFiles(&td)) {
		fprintf(stderr, "There was an issue reading the config file.\n");
		return 2;
	}
	pthread_mutex_init(&td.image_file_lock, NULL);
	pthread_mutex_init(&td.processed_data_lock, NULL);
	pthread_mutex_init(&td.network_heartbeat_mutex, NULL);

	pthread_t networking_thread, processing_thread, server_thread, killserver_thread;
	int rc;
	rc = pthread_create(&networking_thread, NULL, networkMain, &td);
	assert(0==rc);
	
	rc = pthread_create(&processing_thread, NULL, processingMain, &td);
	assert(0==rc);
	
	rc = pthread_create(&server_thread, NULL, serverMain, &td);
	assert(0==rc);
	
	killServerArg_t ks_td;
	pthread_mutex_init(&ks_td.mutex, NULL);
	ks_td.needs_death = 0;
	rc = pthread_create(&killserver_thread, NULL, killServerMain, &ks_td);
	assert(0==rc);

	while (1) {
		// Check the networking thread
		pthread_mutex_lock(&td.network_heartbeat_mutex);
		if (td.networking_is_dead) {
			fprintf(stderr, "Networking had an issue and needs to close.\n");
			return 1;
		}
		pthread_mutex_unlock(&td.network_heartbeat_mutex);

		// Check the killserver thread
		pthread_mutex_lock(&ks_td.mutex);
		if (ks_td.needs_death) {
			fprintf(stderr, "Sombody connected to the death port.\n");
			return 2;
		}
		pthread_mutex_unlock(&ks_td.mutex);

		Sleep(100);
	}
	
	void *retval;
	rc = pthread_join(networking_thread, &retval);
	printf("%p\n", retval);
	if (retval != 0) {
		fprintf(stderr, "An error occurred in the networking thread. Exiting.\n");
		return 1; // 1 == error due to network issue
	}
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
