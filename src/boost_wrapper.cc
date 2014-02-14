#include <boost/python.hpp>
#include <boost/python/iterator.hpp>
#include <boost/python/stl_iterator.hpp>
#include "ScoreCell.h"
#include "ScoreMatrix.h"
#include "align.h"
#include "utils.h"

using namespace boost::python;

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


BOOST_PYTHON_MODULE(boost_wrapper)
{

    class_<ScoreMatrix>("ScoreMatrix", init<size_t, size_t>())
        .def("resize", &ScoreMatrix::resize)
        .def("getSize", &ScoreMatrix::getSize)
        .def("getCapacity", &ScoreMatrix::getCapacity)
        .def("countFilledCells", &ScoreMatrix::countFilledCells)
        .def("getMaxScore", &ScoreMatrix::getMaxScore);

    class_<AlignOpts>("AlignOpts", init< double, double, int, int >())
        .def_readwrite("query_miss_penalty", &AlignOpts::query_miss_penalty)
        .def_readwrite("ref_miss_penalty", &AlignOpts::ref_miss_penalty)
        .def_readwrite("query_max_misses", &AlignOpts::query_max_misses)
        .def_readwrite("ref_max_misses", &AlignOpts::ref_max_misses);
        
    class_<AlignTask>("AlignTask", init<IntVec&, IntVec&, ScoreMatrix&, AlignOpts&>());

    class_< std::vector<int> >("IntVec")
        .def("assign", &vec_assign<int>)
        .def("sum", &sum_all<int>);
}
