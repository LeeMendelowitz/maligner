CXXFLAGS=-Wall -pedantic -std=c++0x

PYTHON_INC = /fs/szattic-asmg7/lmendelo/packages/distribution/include/python2.7
PYTHON_LIB = /fs/szattic-asmg7/lmendelo/packages/distribution/lib/python2.7/config

BOOST_INC = /fs/szattic-asmg7/lmendelo/packages/distribution/include
BOOST_LIB = /fs/szattic-asmg7/lmendelo/packages/distribution/lib

maligner_lib: ../lib/maligner.a

all: boost_wrapper.so ../lib/maligner.a
	
boost_wrapper.so: boost_wrapper.o ScoreCell.o ScoreMatrix.o globals.o align.o
	$(CXX) $(CXXFLAGS) -shared -Wl,--export-dynamic $^ -L$(BOOST_LIB) -lboost_python  -o $@

boost_wrapper.o: boost_wrapper.cc
	$(CXX) $(CXXFLAGS) -I$(PYTHON_INC) -I$(BOOST_INC) -c -fPIC $<

malignpy.so: malignpy.o ScoreCell.o ScoreMatrix.o globals.o align.o
	$(CXX) $(CXXFLAGS) -shared -Wl,--export-dynamic $^ -L$(BOOST_LIB) -lboost_python  -o $@
	cp malignpy.so ../lib/malignpy.so
	
malignpy.o: malignpy.cc
	$(CXX) $(CXXFLAGS) -I$(PYTHON_INC) -I$(BOOST_INC) -c -fPIC $<

../lib/maligner.a: ScoreCell.o ScoreMatrix.o globals.o align.o
	ar rvs maligner.a ScoreCell.o ScoreMatrix.o globals.o align.o
	mv maligner.a ../lib/maligner.a

ScoreCell.o : ScoreCell.cc ScoreCell.h
	$(CXX) $(CXXFLAGS) -c -g -fPIC $<

ScoreMatrix.o : ScoreMatrix.cc ScoreMatrix.h
	$(CXX) $(CXXFLAGS) -c -g -fPIC $<

align.o : align.cc align.h ScoreCell.h ScoreMatrix.h
	$(CXX) $(CXXFLAGS) -c -g -fPIC $<

globals.o: globals.cc globals.h
	$(CXX) $(CXXFLAGS) -c -g -fPIC $<

clean:
	rm *.o *.so ../lib/maligner.a
