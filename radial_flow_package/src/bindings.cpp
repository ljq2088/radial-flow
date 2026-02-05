/**
 * @file bindings.cpp
 * @brief Python 绑定（使用 pybind11）
 */

#include <pybind11/pybind11.h>
#include <pybind11/complex.h>
#include <pybind11/stl.h>
#include "radial_flow.hpp"

namespace py = pybind11;
using namespace radial_flow;

PYBIND11_MODULE(_radial_flow_cpp, m) {
    m.doc() = "径向流计算 C++ 核心模块";

    // SchwarzschildParams
    py::class_<SchwarzschildParams>(m, "SchwarzschildParams")
        .def(py::init<>())
        .def_readwrite("r0", &SchwarzschildParams::r0)
        .def("f", &SchwarzschildParams::f)
        .def("f_prime", &SchwarzschildParams::f_prime);

    // WaveParams
    py::class_<WaveParams>(m, "WaveParams")
        .def(py::init([](Complex omega, int l, int m = 0, bool use_subleading_term = true) {
            WaveParams w;
            w.omega = omega;
            w.l = l;
            w.m = m;
            w.use_subleading_term = use_subleading_term;
            return w;
        }),
        py::arg("omega"),
        py::arg("l"),
        py::arg("m") = 0,
        py::arg("use_subleading_term") = true,
        "Initialize WaveParams with omega, l, optional m and use_subleading_term")
        .def_readwrite("omega", &WaveParams::omega)
        .def_readwrite("l", &WaveParams::l)
        .def_readwrite("m", &WaveParams::m)
        .def_readwrite("use_subleading_term", &WaveParams::use_subleading_term);

    // SolveResult
    py::class_<SolveResult>(m, "SolveResult")
        .def_readonly("r", &SolveResult::r)
        .def_readonly("sigma", &SolveResult::sigma)
        .def_readonly("reflection_coeff", &SolveResult::reflection_coeff)
        .def_readonly("success", &SolveResult::success);

    // RadialFlowSolver
    py::class_<RadialFlowSolver>(m, "RadialFlowSolver")
        .def(py::init<const SchwarzschildParams&, const WaveParams&>())
        .def("potential_V", &RadialFlowSolver::potential_V)
        .def("sigma_rhs", &RadialFlowSolver::sigma_rhs)
        .def("sigma_initial", &RadialFlowSolver::sigma_initial)
        .def("solve", &RadialFlowSolver::solve,
             py::arg("r_max") = 200.0,
             py::arg("epsilon") = 1e-6,
             py::arg("n_points") = 10000,
             py::arg("convergence_threshold") = 0.01,
             py::arg("n_periods") = 2,
             py::arg("max_r_factor") = 5.0)
        .def("compute_reflection_coefficient",
             &RadialFlowSolver::compute_reflection_coefficient,
             py::arg("r_max") = 200.0);
}