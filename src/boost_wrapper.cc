#include <boost/python.hpp>
#include "ScoreCell.h"
#include "ScoreMatrix.h"

using namespace boost::python;


BOOST_PYTHON_MODULE(boost_wrapper)
{
    class_<ScoreMatrix>("ScoreMatrix", init<size_t, size_t>())
        .def("resize", &ScoreMatrix::resize)
        .def("getSize", &ScoreMatrix::getSize)
        .def("getCapacity", &ScoreMatrix::getCapacity)
        .def("countFilledCells", &ScoreMatrix::countFilledCells)
        .def("getMaxScore", &ScoreMatrix::getMaxScore);
}
