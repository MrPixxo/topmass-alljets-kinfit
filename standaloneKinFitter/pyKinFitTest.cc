#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
namespace py = pybind11;

#include "KinFitTest.cc"

PYBIND11_MODULE(pyKinFitTest, m) {
    m.doc() = "pybind11 example plugin"; // optional module docstring

    m.def("setBestCombi", &setBestCombi, "A function that sets the best combination of jets based on the chi2 value of the kinematic fit", py::arg("inputdata"), py::arg("decayprodmats"));
    py::class_<Selection>(m, "Selection")
        .def(py::init<>())
        .def_readwrite("chi2", &Selection::chi2)
        .def_readwrite("bestPermutation", &Selection::bestPermutation)
        .def_readwrite("fitJets", &Selection::fitJets)
        .def_readwrite("combinationType", &Selection::combinationType);
}
