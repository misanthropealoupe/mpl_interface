SHELL = /bin/sh
PYTHON_INC = /usr/include/python2.7
NUMPY_CINC = /usr/local/lib/python2.7/dist-packages/numpy/core/include

.SUFFIXES:
.SUFFIXES: .cpp .hpp .o .so

all: mpl_interface.hpp
	g++ -std=c++11 -fPIC -fpermissive -I$(PYTHON_INC) -I$(NUMPY_CINC) -shared mpl_interface.cpp -o libmpl_interface.so -lpython2.7

test:	
	g++ -std=c++11 -fpermissive -I/usr/local/include -L/usr/local/lib -I$(PYTHON_INC) -I$(NUMPY_CINC) mpl_interface.cpp -o mpl_test.x -lmpl -lpython2.7

install:
	mv libmpl_interface.so /usr/lib64
	cp mpl_interface.hpp /usr/include
