#include <boost/python.hpp>
#include <boost/python/class.hpp>
#include <boost/python/iterator.hpp>
#include <boost/python/stl_iterator.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>

#include "ScoreCell.h"
#include "ScoreMatrix.h"
#include "align.h"
#include "utils.h"

using namespace boost::python;
using self_ns::str;

template<typename T>
void vec_assign(std::vector<T>& v, object o) {
    // Turn a Python sequence into an STL input range
    stl_input_iterator<T> begin(o), end;
    v.assign(begin, end);
}

// Create a thin wrapper for the sum all function
template<typename T>
T sum_all(const std::vector<T>& v) {
  return sum(v);
}

// Construct a vector from a python iterable
template<typename T>
std::vector<T> make_vec(object o) {
    std::vector<T> vec;
    vec_assign(vec, o);
    return vec;
}


BOOST_PYTHON_MODULE(malignpy)
{

    class_<ScoreMatrix>("ScoreMatrix", init<size_t, size_t>())
        .def("resize", &ScoreMatrix::resize)
        .def("getSize", &ScoreMatrix::getSize)
        .def("getCapacity", &ScoreMatrix::getCapacity)
        .def("countFilledCells", &ScoreMatrix::countFilledCells)
        .def("getMaxScore", &ScoreMatrix::getMaxScore);

    class_<ScoreCell>("ScoreCell")
        .def_readonly("q", &ScoreCell::q_)
        .def_readonly("r", &ScoreCell::r_)
        .def_readonly("score", &ScoreCell::score_);

    class_<AlignOpts>("AlignOpts", init< double, double, int, int >())
        .def_readwrite("query_miss_penalty", &AlignOpts::query_miss_penalty)
        .def_readwrite("ref_miss_penalty", &AlignOpts::ref_miss_penalty)
        .def_readwrite("query_max_misses", &AlignOpts::query_max_misses)
        .def_readwrite("ref_max_misses", &AlignOpts::ref_max_misses);
        
    class_<AlignTask>("AlignTask", init<IntVec&, IntVec&, ScoreMatrix&, AlignOpts&>());

    class_< Chunk >("Chunk", init<int, int, int>())
                   .def_readwrite("start", &Chunk::start)
                   .def_readwrite("end", &Chunk::end)
                   .def_readwrite("size", &Chunk::size)
                   .def(self_ns::str(self_ns::self));

    class_< std::vector<Chunk> >("ChunkVec")
        .def(vector_indexing_suite< std::vector<Chunk> >());

    class_< ScoreCellPVec >("ScoreCellVec")
        .def(vector_indexing_suite< ScoreCellPVec >());

    class_< std::vector<int> >("IntVec")
        .def("assign", &vec_assign<int>)
        .def("sum", &sum_all<int>)
        .def(vector_indexing_suite< std::vector<int> >());
    
    // Functions
    def("fill_score_matrix", fill_score_matrix);
    def("build_trail", build_trail);
    def("build_chunk_trail", build_chunk_trail);
    def("get_best_alignment", get_best_alignment);



}
