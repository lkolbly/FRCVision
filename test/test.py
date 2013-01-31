# Tester for the FRC Vision.
# Connects to the server, and prints out the target data.
import socket, struct

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

# Request a camera update frame.
sock.send(struct.pack("!IIBH", 3, 0x02000001, 0x01, 1))

# Wait for the response.
data = getData(sock, 4)
data = getData(sock, 4)
print len(data)
if struct.unpack("!I", data)[0] != 0x82000002:
	print "We were sent the wrong packet type"

data = getData(sock, 2) # Eat the sub-type ID

ntargets = struct.unpack("!H", getData(sock, 2))[0]
for i in range(ntargets):
        #print i, ntargets
        values = struct.unpack("!hHHHHHH", getData(sock,14))
        print values

sock.shutdown(0)

# Real quick, let's go grab the test case we ran against
import urllib2

f = urllib2.urlopen("http://localhost:8006/")
open("test-case.jpg", "wb").write(f.read())
f.close()
