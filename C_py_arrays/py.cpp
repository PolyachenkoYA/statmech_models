#include "test.cpp"

py::array_t<double> compute_evolution(int N, double Temp, int Nt, int my_seed, int verbose)
{
    printf("Nt = %d\n", Nt);

    //    double* E = static_cast<double *>(info.ptr);
    auto E = py::array_t<double>(Nt);
    py::buffer_info info = E.request();
    double *E_ptr = static_cast<double *>(info.ptr);

    print_E(E_ptr, Nt, 'p');

    compute_evolution_C(N, Temp, Nt, my_seed, E_ptr);

    print_E(E_ptr, Nt, 'P');

    return E;
}

py::array_t<double> add_arrays(py::array_t<double> input1, py::array_t<double> input2)
{
    py::buffer_info buf1 = input1.request();
    py::buffer_info buf2 = input2.request();

    if (buf1.ndim != 1 || buf2.ndim != 1)
        throw std::runtime_error("Number of dimensions must be one");

    if (buf1.size != buf2.size)
        throw std::runtime_error("Input shapes must match");

    /* No pointer is passed, so NumPy will allocate the buffer */
    auto result = py::array_t<double>(buf1.size);

    py::buffer_info buf3 = result.request();

    double *ptr1 = static_cast<double *>(buf1.ptr);
    double *ptr2 = static_cast<double *>(buf2.ptr);
    double *ptr3 = static_cast<double *>(buf3.ptr);

	print_E(ptr1, buf1.size, 'a');
	print_E(ptr2, buf1.size, 'b');
	print_E(ptr3, buf1.size, 'c');

    for (size_t idx = 0; idx < buf1.shape[0]; idx++)
        ptr3[idx] = ptr1[idx] + ptr2[idx];

	print_E(ptr1, buf1.size, 'A');
	print_E(ptr2, buf1.size, 'B');
	print_E(ptr3, buf1.size, 'C');

    return result;
}

PYBIND11_MODULE(test, m) {
    m.def("add_arrays", &add_arrays, "Add two NumPy arrays");

    m.def("compute_evolution", &compute_evolution,
        "Fill preallocated array E with 'time'-evolution with Temp for Nt steps; init with my_seed (default=time(NULL))",
        py::arg("grid_size"),
        py::arg("Temp"),
        py::arg("N_steps"),
        py::arg("my_seed")=time(nullptr),
        py::arg("verbose")=0
    );

}
