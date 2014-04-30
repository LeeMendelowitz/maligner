CXXFLAGS=-Wall -pedantic -std=c++0x -O3

PYTHON_INC = /fs/szattic-asmg7/lmendelo/packages/distribution/include/python2.7
PYTHON_LIB = /fs/szattic-asmg7/lmendelo/packages/distribution/lib/python2.7/config

BOOST_INC = /fs/szattic-asmg7/lmendelo/packages/distribution/include
BOOST_LIB = /fs/szattic-asmg7/lmendelo/packages/distribution/lib

maligner_lib: ../lib/maligner.a

all: ../lib/maligner.a ../lib/malignpy/malignpy.so

malignpy.so: malignpy.o ScoreCell.o ScoreMatrix.o globals.o align.o
	$(CXX) $(CXXFLAGS) -shared -Wl,--export-dynamic $^ -L$(BOOST_LIB) -lboost_python  -o $@
	
../lib/malignpy/malignpy.so: malignpy.so
	cp malignpy.so ../lib/malignpy/malignpy.so
	
malignpy.o: malignpy.cc align.h ScoreCell.h ScoreMatrix.h ScoreCell.h globals.h
	$(CXX) $(CXXFLAGS) -I$(PYTHON_INC) -I$(BOOST_INC) -c -fPIC $<

../lib/maligner.a: ScoreCell.o ScoreMatrix.o globals.o align.o
	ar rvs maligner.a ScoreCell.o ScoreMatrix.o globals.o align.o
	mv maligner.a ../lib/maligner.a

ScoreCell.o : ScoreCell.cc ScoreCell.h
	$(CXX) $(CXXFLAGS) -I$(BOOST_INC) -c -g -fPIC $<

ScoreMatrix.o : ScoreMatrix.cc ScoreMatrix.h
	$(CXX) $(CXXFLAGS) -I$(BOOST_INC) -c -g -fPIC $<

align.o : align.cc align.h ScoreCell.h ScoreMatrix.h
	$(CXX) $(CXXFLAGS) -I$(BOOST_INC) -c -g -fPIC $<

globals.o: globals.cc globals.h
	$(CXX) $(CXXFLAGS) -I$(BOOST_INC) -c -g -fPIC $<

clean:
	rm *.o *.so ../lib/maligner.a
