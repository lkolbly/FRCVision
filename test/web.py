import BaseHTTPServer, os, time

TEST_SET_DIR = "test-data\\set1"
TEST_SET_COUNT = 0

class MyRequestHandler(BaseHTTPServer.BaseHTTPRequestHandler):
    def do_GET(self):
        self.send_response(200)
        self.send_header("Content-type", "image/jpg")
        self.end_headers()

        files = os.listdir(TEST_SET_DIR)
        SWITCH_TIME = 1
        choice = int((time.time()%(SWITCH_TIME*len(files))) / SWITCH_TIME)
        #choice = TEST_SET_COUNT
        #TEST_SET_COUNT += 1
        #TEST_SET_COUNT = TEST_SET_COUNT % len(files)
        #data = open("%s\\%s"%(TEST_SET_DIR,files[choice]), "rb").read()
        data = open("%s\\%s"%(TEST_SET_DIR,"image.cgi.24"), "rb").read()
        print choice, files[choice], len(data)
        self.wfile.write(data)
        pass

def run(server_class=BaseHTTPServer.HTTPServer,
        handler_class=BaseHTTPServer.BaseHTTPRequestHandler):
    server_address = ('', 8006)
    httpd = server_class(server_address, handler_class)
    httpd.serve_forever()

run(handler_class=MyRequestHandler)
