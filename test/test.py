# Tester for the FRC Vision.
# Connects to the server, and prints out the target data.
import socket, struct, time
import urllib2

def getData(s, nbytes):
        data = ""
        while len(data) < nbytes:
                data += sock.recv(nbytes-len(data))
        return data

# First, open the socket.
sock = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
sock.connect(("localhost", 4180))

sock.send(struct.pack("!III", 4, 1, 1)) # Handshake

# Read out the s->c handshake.
data = getData(sock, 4)
getData(sock, struct.unpack("!I", data)[0]+4)

cnt = 0
def getResponse():
        # Request a camera update frame.
        sock.send(struct.pack("!IIBH", 3, 0x02000001, 0x01, 1))

        # Wait for the response.
        data = getData(sock, 4)
        data = getData(sock, 4)
        if struct.unpack("!I", data)[0] != 0x82000002:
        	print "We were sent the wrong packet type"

        data = getData(sock, 2) # Eat the sub-type ID

        ntargets = struct.unpack("!H", getData(sock, 2))[0]
        print "Got %i targets."%ntargets
        for i in range(ntargets):
                #print i, ntargets
                values = struct.unpack("!hHHHHHH", getData(sock,14))
                print i, values
                f = open("unittest/%i"%cnt, "w");
                f.write("tgt %i, %s\n"%(i, values))
                f.close()

        f = urllib2.urlopen("http://localhost:8006/")
        open("unittest/test-case-%i.jpg"%cnt, "wb").write(f.read())
        f.close()

while 1:
        getResponse()
        time.sleep(2);
        cnt += 1

sock.shutdown(0)

# Real quick, let's go grab the test case we ran against

