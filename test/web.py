import BaseHTTPServer, os, time

class MyRequestHandler(BaseHTTPServer.BaseHTTPRequestHandler):
    def do_GET(self):
        self.send_response(200)
        self.send_header("Content-type", "image/jpg")
        self.end_headers()

        files = os.listdir("test-data")
        SWITCH_TIME = 10
        choice = int((time.time()%(SWITCH_TIME*len(files))) / SWITCH_TIME)
        self.wfile.write(open("test-data/%s"%files[choice]).read())
        pass

def run(server_class=BaseHTTPServer.HTTPServer,
        handler_class=BaseHTTPServer.BaseHTTPRequestHandler):
    server_address = ('', 8006)
    httpd = server_class(server_address, handler_class)
    httpd.serve_forever()

run(handler_class=MyRequestHandler)
