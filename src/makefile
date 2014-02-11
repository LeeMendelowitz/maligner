CXXFLAGS=-Wall -pedantic -std=c++0x

PYTHON_INC = /fs/szattic-asmg7/lmendelo/packages/distribution/include/python2.7
PYTHON_LIB = /fs/szattic-asmg7/lmendelo/packages/distribution/lib/python2.7/config

BOOST_INC = /fs/szattic-asmg7/lmendelo/packages/distribution/include
BOOST_LIB = /fs/szattic-asmg7/lmendelo/packages/distribution/lib

TARGET = boost_wrapper

boost_wrapper.so: boost_wrapper.o ScoreCell.o ScoreMatrix.o globals.o
	$(CXX) $(CXXFLAGS) -shared -Wl,--export-dynamic $^ -L$(BOOST_LIB) -lboost_python -o $@

boost_wrapper.o: boost_wrapper.cc
	$(CXX) $(CXXFLAGS) -I$(PYTHON_INC) -I$(BOOST_INC) -c -fPIC $<

ScoreCell.o : ScoreCell.cc ScoreCell.h
	$(CXX) $(CXXFLAGS) -I$(PYTHON_INC) -I$(BOOST_INC) -c -fPIC $<

ScoreMatrix.o : ScoreMatrix.cc ScoreMatrix.h
	$(CXX) $(CXXFLAGS) -I$(PYTHON_INC) -I$(BOOST_INC) -c -fPIC $<

globals.o: globals.cc globals.h
	$(CXX) $(CXXFLAGS) -I$(PYTHON_INC) -I$(BOOST_INC) -c -fPIC $<

clean:
	rm *.o *.so
