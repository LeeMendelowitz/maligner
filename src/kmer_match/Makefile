CXXFLAGS=-Wall -pedantic -std=c++11 -g
CXXFLAGS=-Wall -pedantic -std=c++11 -pg
CXXFLAGS=-Wall -pedantic -std=c++11 -O3
INC=-I/usr/local/include
LIB=-L/usr/local/lib


all: maligner

test: test_maligner_get_all test_maligner_get_best

map.o : map.cpp map.h
	$(CXX) $(CXXFLAGS) -c $<

map_frag.o : map_frag.cpp map.h
	$(CXX) $(CXXFLAGS) -c $<

error_model.o : error_model.cpp error_model.h 
	$(CXX) $(CXXFLAGS) -c $<

map_frag_db.o : map_frag_db.cpp map_frag_db.h map.h map_frag.h
	$(CXX) $(CXXFLAGS) -c $<

map_chunk.o: map_chunk.cpp map_chunk.h map.h 
	$(CXX) $(CXXFLAGS) -c $<

ref_alignment.o: ref_alignment.cpp ref_alignment.h map_chunk.h 
	$(CXX) $(CXXFLAGS) -c $<

map_chunk_db.o : map_chunk_db.cpp map.h map_chunk.h map_chunk_db.h ref_alignment.h
	$(CXX) $(CXXFLAGS) -c $<

kmer_match : main.cpp map.o map_frag.o map_frag_db.o error_model.o
	$(CXX) $(CXXFLAGS) $^ -o $@

kmer_match_count : main_count.cpp map.o map_frag.o map_frag_db.o error_model.o
	$(CXX) $(CXXFLAGS) $^ -o $@

kmer_match_no_boundary : main_no_boundary.cpp map.o map_frag.o map_frag_db.o error_model.o
	$(CXX) $(CXXFLAGS) $^ -o $@

maligner : maligner.cpp map_chunk.o map_chunk_db.o map.o map_frag.o map_frag_db.o error_model.o ref_alignment.o
	$(CXX) $(CXXFLAGS) $^ -o $@

test_maligner_get_all : test_maligner_get_all.cpp map_chunk.o map_chunk_db.o map.o map_frag.o map_frag_db.o error_model.o ref_alignment.o
	$(CXX) $(CXXFLAGS) $^ -o $@

test_maligner_get_best : test_maligner_get_best.cpp map_chunk.o map_chunk_db.o map.o map_frag.o map_frag_db.o error_model.o ref_alignment.o
	$(CXX) $(CXXFLAGS) $^ -o $@


clean:
	rm -rf *.o maligner kmer_match_no_boundary kmer_match kmer_match_count
