include makefile.local

.SUFFIXES:
.SUFFIXES: .cpp .hpp .o .so

all: mpl_interface.hpp
	g++ -std=c++11 -fPIC -fpermissive -I$(PYTHON_INC) -I$(NUMPY_CINC) -shared mpl_interface.cpp -o libmpl_interface.so -lpython2.7

test:	
	g++ -std=c++11 -fpermissive -I$(INC) -L$(LIB) -I$(PYTHON_INC) -I$(NUMPY_CINC) mpl_interface.cpp -o mpl_test.x -lmpl -lpython2.7

install:
	mv libmpl_interface.so $(LIB)
	cp mpl_interface.hpp $(INC)