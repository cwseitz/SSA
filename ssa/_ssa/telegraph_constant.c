#include <Python.h>
#include <numpy/arrayobject.h>
#include <stdio.h>
#include <gsl_randist.h>
#include <gsl_rng.h>
#include <sys/time.h>
#include "complex.h"

#define N         4		// number of reaction
#define M         3		// number of chemical species
	
void init(int x[]){

	// population of chemical species
	x[0] = 1;
	x[1] = 0;
	x[2] = 0;

}


void update_p(double t, double p[], double k_on, double k_off, double ksyn, double kdeg, int x[]){

	p[0] = k_on*x[0];
	p[1] = k_off*x[1];
	p[2] = ksyn*x[1];
	p[3] = kdeg;
	//printf("%f, %f, %f, %f",p[0],p[1],p[2],p[3]);
}


int select_reaction(double p[], int pn, double sum_propencity, double r){
	int reaction = -1;
	double sp = 0.0;
	int i;
	r = r * sum_propencity;
	for(i=0; i<pn; i++){
		sp += p[i];
		if(r < sp){
			reaction = i;
			break;
		}
	}
	return reaction;
}

void update_x(int x[], int reaction){
	if (reaction == 0){
	   x[0] = 0; 
	   x[1] = 1;
	} 
	else if (reaction == 1){
	   x[0] = 1;
	   x[1] = 0;
	   }
	else if (reaction == 2){
	    x[2] = x[2] + 1;
	}
	else if (reaction == 3){
	    x[2] = x[2] - 1;
	    if (x[2] < 0){
	       x[2] = 0;
	       }
	    }

}

double sum(double a[], int n){
	int i;
	double s=0.0;
	for(i=0; i<n; i++) 
		s += a[i];
	return(s);
}

int telegraph_constant_ssa(int* x1, int* x2, int* x3, double* times, double end_time, double k_on, double k_off, double ksyn, double kdeg){

	int x[M];			    // population of chemical species
	double p[N];			// propencities of reaction

	// initialization
	double sum_propencity = 0.0;	// sum of propencities
	double tau=0.0;			// step of time
	double t=0.0;			// time
	double r;			// random number
	int reaction;			// reaction number selected

	init(x);
    int i = 0;
    time_t current_time = time(NULL);
    double seed = (double)current_time + (double)clock() / CLOCKS_PER_SEC;
    srand((unsigned int)(seed * 1000000.0));
	x1[0] = 1;
	x2[0] = 0;
	x3[0] = 0;
		
	// main loop
	//printf("Time: 0 hours, x1=%d, x2=%d \n", x1[0], x2[0]);
	while(t < end_time){
		
		// update propencity
		update_p(t, p, k_on, k_off, ksyn, kdeg, x);
		sum_propencity = sum(p, N);
	    i += 1;

		// sample tau
		if(sum_propencity > 0){
			tau = -log((double)rand()/INT_MAX) / sum_propencity;
		}else{
			break;
		}
	
		// select reaction
		r = (double)rand()/INT_MAX;
		reaction = select_reaction(p, N, sum_propencity, r);
		// update chemical species
		update_x(x, reaction);
		x1[i] = x[0];
		x2[i] = x[1];
		x3[i] = x[2];
        
		// time
		t += tau;
		times[i] = t;
		//printf("Time: %f, Reaction: %d x1=%d, x2=%d, x3=%d\n", t, reaction, x1[i], x2[i], x3[i]);	
	}
	
	return(0);
}


static PyObject* telegraph_constant(PyObject* Py_UNUSED(self), PyObject* args) {

  PyObject* list;

  if (!PyArg_ParseTuple(args, "O!", &PyList_Type, &list))
    return NULL;
    
    double end_time = PyFloat_AsDouble(PyList_GetItem(list, 0));
    double k_on = PyFloat_AsDouble(PyList_GetItem(list, 1));
    double k_off = PyFloat_AsDouble(PyList_GetItem(list, 2));
    double ksyn = PyFloat_AsDouble(PyList_GetItem(list, 3));
    double kdeg = PyFloat_AsDouble(PyList_GetItem(list, 4));
    int Nt = PyLong_AsLong(PyList_GetItem(list, 5));
 
	int* x1 = calloc(Nt, sizeof(int));
	int* x2 = calloc(Nt, sizeof(int));
	int* x3 = calloc(Nt, sizeof(int));
	double* times = calloc(Nt, sizeof(double));
	telegraph_constant_ssa(x1,x2,x3,times,end_time,k_on,k_off,ksyn,kdeg);

	npy_intp dims[1] = {Nt};
	PyObject *x1_out = PyArray_SimpleNew(1, dims, NPY_INT);
	PyObject *x2_out = PyArray_SimpleNew(1, dims, NPY_INT);
	PyObject *x3_out = PyArray_SimpleNew(1, dims, NPY_INT);
	PyObject *times_out = PyArray_SimpleNew(1, dims, NPY_DOUBLE);

    memcpy(PyArray_DATA(x1_out), x1, Nt * sizeof(int));
	memcpy(PyArray_DATA(x2_out), x2, Nt * sizeof(int));
	memcpy(PyArray_DATA(x3_out), x3, Nt * sizeof(int));
	memcpy(PyArray_DATA(times_out), times, Nt * sizeof(double));

	return Py_BuildValue("(OOOO)", x1_out, x2_out, x3_out, times_out);
	free(x1);
	free(x2);
	free(x3);
	free(times);
}
