#include <boost/python.hpp>
#include <boost/python/class.hpp>
#include <boost/python/iterator.hpp>
#include <boost/python/stl_iterator.hpp>
#include <boost/python/return_internal_reference.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>

#include "ScoreCell.h"
#include "ScoreMatrix.h"
#include "align.h"
#include "utils.h"
#include "types.h"

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

// Make a new Score Matrix, and return as a pointer.
// This has the important advantage that we avoid copies when passing to Python.
ScoreMatrix* make_score_matrix(size_t m, size_t n) {
    return new ScoreMatrix(m, n);
}

IntVec* make_int_vec() {
    return new IntVec();
}

void ptask(AlignTask& task) {
    print_align_task(std::cout, task) << std::endl;
}


BOOST_PYTHON_MODULE(malignpy)
{

    class_<ScoreMatrix>("ScoreMatrix", no_init)
        .def("resize", &ScoreMatrix::resize)
        .def("getSize", &ScoreMatrix::getSize)
        .def("getCapacity", &ScoreMatrix::getCapacity)
        .def("countFilledCells", &ScoreMatrix::countFilledCells)
        .def("percentFilled", &ScoreMatrix::percentFilled)
        .def("getMaxScore", &ScoreMatrix::getMaxScore);

    // C++ promises to never delete the ScoreMatrix. This is done by Python garbage
    // collector at shutdown or when matrix is explicitly deleted with del.
    def("make_score_matrix", make_score_matrix, return_value_policy<manage_new_object>());

    class_<ScoreCell>("ScoreCell")
        .def_readonly("q", &ScoreCell::q_)
        .def_readonly("r", &ScoreCell::r_)
        .def_readonly("score", &ScoreCell::score_);

    class_<AlignOpts>("AlignOpts", init< double, double, int, int, double >())
        .def_readwrite("query_miss_penalty", &AlignOpts::query_miss_penalty)
        .def_readwrite("ref_miss_penalty", &AlignOpts::ref_miss_penalty)
        .def_readwrite("query_max_misses", &AlignOpts::query_max_misses)
        .def_readwrite("ref_max_misses", &AlignOpts::ref_max_misses);
        
    class_<AlignTask>("AlignTask", init<IntVec&, IntVec&, ScoreMatrix *, AlignOpts&>());

    class_< Chunk >("Chunk", no_init)
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

    def("make_int_vec", make_int_vec, return_value_policy<manage_new_object>() );

    class_< Score >("Score")
        .def("total", &Score::total)
        .def_readonly("query_miss_score", &Score::query_miss_score)
        .def_readonly("ref_miss_score", &Score::ref_miss_score)
        .def_readonly("sizing_score", &Score::sizing_score);

    class_< Alignment >("Alignment")
        .def_readonly("matched_chunks", &Alignment::matched_chunks)
        .def_readonly("score", &Alignment::score)
        .def_readonly("num_matched_sites", &Alignment::num_matched_sites)
        .def_readonly("query_misses", &Alignment::query_misses)
        .def_readonly("ref_misses", &Alignment::ref_misses)
        .def_readonly("query_miss_rate", &Alignment::query_miss_rate)
        .def_readonly("ref_miss_rate", &Alignment::ref_miss_rate)
        .def_readonly("total_miss_rate", &Alignment::total_miss_rate)
        .def_readonly("query_interior_size", &Alignment::query_interior_size)
        .def_readonly("ref_interior_size", &Alignment::ref_interior_size)
        .def_readonly("interior_size_ratio", &Alignment::interior_size_ratio)                
        .def(self_ns::str(self_ns::self));
    
    // Functions
    def("fill_score_matrix", fill_score_matrix);
    def("build_trail", build_trail);
    def("build_chunk_trail", build_chunk_trail);
    def("get_best_alignment", get_best_alignment);
    def("alignment_from_trail", alignment_from_trail, return_value_policy<manage_new_object>());
    def("print_align_task", ptask);
    def("make_best_alignment", make_best_alignment, return_value_policy<manage_new_object>());

}
