#include <Python.h>
#include <numpy/arrayobject.h>
#include <stdio.h>
#include <gsl_randist.h>
#include <gsl_rng.h>
#include <sys/time.h>
#include <complex.h>
#include "switch.c"

static PyMethodDef SSAMethods[] = {
    {"fret_sim", fret_sim, METH_VARARGS, "SSA in C"},
    {NULL, NULL, 0, NULL},
};


static struct PyModuleDef _SSA = {
    PyModuleDef_HEAD_INIT,
    "_SSA",
    "Python interface for C",
    -1,
    SSAMethods
};

PyMODINIT_FUNC PyInit__SSA(void) {
    import_array();
    return PyModule_Create(&_SSA);
}
