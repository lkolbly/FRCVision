# Tester for the FRC Vision.
# Connects to the server, and prints out the target data.
import socket, struct

# First, open the socket.
sock = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
sock.connect(("localhost", 4180))

sock.send(struct.pack("!III", 4, 1, 1)) # Handshake

# Read out the s->c handshake.
data = sock.recv(4)
sock.recv(struct.unpack("!I", data)[0])

# Request a camera update frame.
sock.send(struct.pack("!IIBH", 3, 0x02000001, 0x01, 1))

# Wait for the response.
data = sock.recv(4)
data = sock.recv(4)
if struct.unpack("!I", data)[0] != 0x82000002:
	print "We were sent the wrong packet type"

data = sock.recv(2) # Eat the sub-type ID

ntargets = struct.unpack("!I", sock.recv(2))[0]
for i in ntargets:
	values = struct.unpack("!hHHHHHH", sock.recv(14))
	print values
