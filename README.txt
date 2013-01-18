Build Instructions (known to work on 32-bit):
=============================================
Obviously, you'll need a working mingw setup.

You will also need OpenCV and libcurl. OpenCV needs to be compiled from source, at least I did, so download the source package from their website and run the following commands:
$ cmake -G "MinGW Makefiles" .
$ mingw32-make
This will create a directory "install" inside the opencv root folder. You need to copy all of the "dll" files within install\bin to the "bin" directory within the FRCVision source tree. Note that this bin directory does not exist when the code is checked out, so you must create it yourself.

Then you need to edit the CMakeLists.txt file within the FRCVision root. You need to change the "SET(OpenCV_DIR, '...')" to point to your OpenCV's install path. Note that you cannot delete OpenCV, or else you won't be able to compile against it.

The other library is libcurl, which is slightly easier. To do this, Simply download the Win32 zip file from their website, and don't choose the one labelled "binary".
