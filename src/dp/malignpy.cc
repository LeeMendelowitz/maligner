#include <boost/python.hpp>
#include <boost/python/class.hpp>
#include <boost/python/iterator.hpp>
#include <boost/python/stl_iterator.hpp>
#include <boost/python/return_internal_reference.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>

#include <memory>
#include <string>

#include "ScoreCell.h"
#include "ScoreMatrix.h"
#include "align.h"
#include "utils.h"
#include "types.h"

using namespace boost::python;
using std::shared_ptr;
using self_ns::str;
using namespace maligner_dp;

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
ScoreMatrixPtr make_score_matrix(size_t m, size_t n) {
    return ScoreMatrixPtr(new ScoreMatrix(m, n));
}

IntVecPtr make_int_vec() {
    IntVecPtr p = IntVecPtr(new IntVec());
    return p;
}

BoolVecPtr make_bool_vec() {
    BoolVecPtr p = BoolVecPtr(new BoolVec());
    return p;
}

AlignmentPVecPtr make_alignment_vec() {
    AlignmentPVecPtr p = AlignmentPVecPtr(new AlignmentPVec());
    return p;
}

void print_align_task2(AlignTask& task) {
    print_align_task(std::cout, task) << std::endl;
}


BOOST_PYTHON_MODULE(malignpy)
{

    class_<ScoreMatrix, ScoreMatrixPtr >("ScoreMatrix", no_init)
        .def("resize", &ScoreMatrix::resize)
        .def("getSize", &ScoreMatrix::getSize)
        .def("getCapacity", &ScoreMatrix::getCapacity)
        .def("countFilledCells", &ScoreMatrix::countFilledCells)
        .def("countFilledByRow", &ScoreMatrix::countFilledByRow)
        .def("percentFilled", &ScoreMatrix::percentFilled)
        .def("getMaxScore", &ScoreMatrix::getMaxScore)
        .add_property("nrows", &ScoreMatrix::getNumRows)
        .add_property("m", &ScoreMatrix::getNumRows)
        .add_property("ncols", &ScoreMatrix::getNumCols)
        .add_property("n", &ScoreMatrix::getNumCols);

    // C++ promises to never delete the ScoreMatrix. This is done by Python garbage
    // collector at shutdown or when matrix is explicitly deleted with del.
    def("make_score_matrix", make_score_matrix);

    class_<ScoreCell>("ScoreCell")
        .def_readonly("q", &ScoreCell::q_)
        .def_readonly("r", &ScoreCell::r_)
        .def_readonly("score", &ScoreCell::score_);

    class_<AlignOpts>("AlignOpts", init< double, double, int, int, 
         double, double, double, int, int, int, bool>())
        .def_readwrite("query_miss_penalty", &AlignOpts::query_miss_penalty)
        .def_readwrite("ref_miss_penalty", &AlignOpts::ref_miss_penalty)
        .def_readwrite("query_max_misses", &AlignOpts::query_max_misses)
        .def_readwrite("ref_max_misses", &AlignOpts::ref_max_misses)
        .def_readwrite("sd_rate", &AlignOpts::sd_rate)
        .def_readwrite("max_chunk_sizing_error", &AlignOpts::max_chunk_sizing_error)
        .def_readwrite("min_sd", &AlignOpts::min_sd)
        .def_readwrite("alignments_per_reference", &AlignOpts::alignments_per_reference)
        .def_readwrite("min_alignment_spacing", &AlignOpts::min_alignment_spacing)
        .def_readwrite("neighbor_delta", &AlignOpts::neighbor_delta)
        .def_readwrite("query_is_bounded", &AlignOpts::query_is_bounded);
       
    class_< MapData, MapDataPtr >("MapData", 
            init< const std::string&, size_t, bool >() )
     .def_readwrite("map_name", &MapData::map_name_)
     .def_readwrite("num_frags", &MapData::num_frags_)
     .def_readwrite("is_bounded", &MapData::is_bounded_)
     .def("print_debug_stats", &MapData::print_debug_stats)
     .staticmethod("print_debug_stats");

    class_<AlignTask>("AlignTask",
                      init<MapDataPtr, MapDataPtr, IntVecPtr, IntVecPtr, PartialSumsPtr,
                           PartialSumsPtr , ScoreMatrixPtr,    
                           AlignmentPVecPtr,
                           AlignOpts&>()
                      )
        .def(init<MapDataPtr, MapDataPtr, IntVecPtr, IntVecPtr,
                           PartialSumsPtr, PartialSumsPtr,
                           int, ScoreMatrixPtr,    
                           AlignmentPVecPtr,
                           AlignOpts&>());

    class_< Chunk >("Chunk", no_init)
                   .def_readwrite("start", &Chunk::start)
                   .def_readwrite("end", &Chunk::end)
                   .def_readwrite("size", &Chunk::size)
                   .def_readonly("is_boundary", &Chunk::is_boundary)
                   .def(self_ns::str(self_ns::self));

    class_< std::vector<Chunk> >("ChunkVec")
        .def(vector_indexing_suite< std::vector<Chunk> >());

    class_< ScoreCellPVec >("ScoreCellVec")
        .def(vector_indexing_suite< ScoreCellPVec >());

    class_< IntVec, IntVecPtr >("IntVec")
        .def("assign", &vec_assign<int>)
        .def("sum", &sum_all<int>)
        .def(vector_indexing_suite< std::vector<int> >());

    // To correctly work with a vector of shared_ptr, the
    // second argument to VectorIndexingSuite must be true.
    // Black magic. See:
    //   http://stackoverflow.com/questions/5919170/boost-python-and-vectors-of-shared-ptr
    class_< AlignmentPVec, AlignmentPVecPtr >("AlignmentPVec")
        .def(vector_indexing_suite< AlignmentPVec, true >());

    class_< PartialSums, PartialSumsPtr >("PartialSums")
        .def(vector_indexing_suite< PartialSums >());

    class_< MatchedChunkVec >("MatchedChunkVec")
        .def(vector_indexing_suite< MatchedChunkVec >());


    class_< BoolVec, BoolVecPtr >("BoolVec")
        .def(vector_indexing_suite< BoolVec >());

    def("make_int_vec", make_int_vec);

    class_< Score >("Score")
        .add_property("total", &Score::total)
        .def_readonly("query_miss_score", &Score::query_miss_score)
        .def_readonly("ref_miss_score", &Score::ref_miss_score)
        .def_readonly("sizing_score", &Score::sizing_score);

    class_< Alignment, AlignmentPtr, boost::noncopyable >( "Alignment", no_init )
        .def_readonly("matched_chunks", &Alignment::matched_chunks)
        .def_readonly("rescaled_matched_chunks", &Alignment::rescaled_matched_chunks)
        .def_readonly("score", &Alignment::score)
        .def_readonly("rescaled_score", &Alignment::rescaled_score)
        .def_readonly("num_matched_sites", &Alignment::num_matched_sites)
        .def_readonly("query_misses", &Alignment::query_misses)
        .def_readonly("ref_misses", &Alignment::ref_misses)
        .def_readonly("query_miss_rate", &Alignment::query_miss_rate)
        .def_readonly("ref_miss_rate", &Alignment::ref_miss_rate)
        .def_readonly("total_miss_rate", &Alignment::total_miss_rate)
        .def_readonly("query_interior_size", &Alignment::query_interior_size)
        .def_readonly("ref_interior_size", &Alignment::ref_interior_size)
        .def_readonly("interior_size_ratio", &Alignment::interior_size_ratio)
        .def_readwrite("query_scaling_factor", &Alignment::query_scaling_factor)
        .def("rescale_matched_chunks", &Alignment::rescale_matched_chunks)                
        .def(self_ns::str(self_ns::self));

    class_< MatchedChunk >("MatchedChunk")
        .def_readonly("query_chunk", &MatchedChunk::query_chunk)
        .def_readonly("ref_chunk", &MatchedChunk::ref_chunk)
        .def_readonly("score", &MatchedChunk::score)
        .add_property("is_boundary", &MatchedChunk::is_boundary);
    
    // Functions
    def("fill_score_matrix", fill_score_matrix);
    def("build_trail", build_trail);
    def("build_chunk_trail", build_chunk_trail);
    def("get_best_alignment", get_best_alignment);
    def("alignment_from_trail", alignment_from_trail);
    def("print_align_task", print_align_task2);
    def("make_best_alignment", make_best_alignment);
    def("make_partial_sums", make_partial_sums_new);
    def("make_bool_vec", make_bool_vec);
    def("make_alignment_vec", make_alignment_vec);
    def("fill_score_matrix_using_partials", fill_score_matrix_using_partials);
    def("make_best_alignment_using_partials", make_best_alignment_using_partials);
    def("make_best_alignments_using_partials", make_best_alignments_using_partials);

}
