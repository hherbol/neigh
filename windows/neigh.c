/* neigh.c
A code to generate a neighbourhood list for python codes.

Given:
    - pos - List of 3D lists
    - R - Cutoff radius
    - PBC - List of x, y, z lengths if PBC is desired
Return:
    - neighbours - A list of lists, holding indices to the neighbours
                 - where indices are the respective positions in the
                 - given pos list.
*/
#include <Python.h>
#include <malloc.h>
#include <stdlib.h>

static PyObject* py_gen_neigh(PyObject* self, PyObject* args, PyObject* keywds) {
    // Initialize variables
    int R, R2, C, rows, cols, MAX_SIZE, TOO_CLOSE, DELTA, n, counter;
    int N_count, *N_neigh,*neighs;
    int *resized;
    int SKIN_A, passed;
    double *pos, pbc_pos[3], cutoff, rsq, dist, r, x, y, z, *PBC, *origin;
    static char *kwlist[] = {"list_obj", "cutoff","PBC","origin", NULL};
    PyObject * list_obj;
    PyObject * line;
    PyObject * neighbours;
    PyObject * neigh_line;
    PyObject * pbc = NULL;
    PyObject * origin_list = NULL;

    // Starting size for the number of neighbours we want, and the delta shift when approaching our limit
    resized = NULL;
    MAX_SIZE = 1000;
    TOO_CLOSE = 100;
    DELTA = 1000;
    N_count = 0;
    PBC = NULL;
    origin = NULL;

    // Ensure it is a list being passed
    // NOTE! The format string "O!" says that I am reading in an object and requiring it to be a certain type
    // Therefore, I give this function two variables.  The first is the type I want to read in and the second
    // Is the memory location for the list itself.
    // For more info, go here: https://docs.python.org/2/c-api/arg.html
    if ( !PyArg_ParseTupleAndKeywords( args, keywds, "O!d|O!O!", kwlist, &PyList_Type, &list_obj, &cutoff, &PyList_Type, &pbc, &PyList_Type, &origin_list) ) return NULL;
    if (cutoff < 0) return NULL;
    rsq = cutoff*cutoff;

    // Here we dont want a list of 0 rows / columns.
    // Get rows
    if ( (rows = PyList_Size(list_obj)) < 0) return NULL;
    // Get columns
    if ( (cols = PyList_Size(PyList_GetItem(list_obj, 0))) < 0) return NULL;

    // Generate local pos array
    // This is a 2D array compressed into 1D.  Grab indices by pos[R*cols + C] where R,C is the
    // 2D Row, Column index and cols is the number of values in each row.
    pos = (double *)malloc(rows*cols*sizeof(double));
    N_neigh = (int *)malloc(rows*sizeof(int));
    neighs = (int *)malloc(MAX_SIZE*sizeof(int));

    if (pos == NULL ||
        N_neigh == NULL ||
        neighs == NULL) {
        printf("ERROR - Unable to allocate memory...");
        exit(1);
    }

    if (pbc != NULL) {
        PBC = (double *)malloc(cols*sizeof(double));
        if (PBC == NULL) {
            printf("ERROR - Unable to allocate memory...");
            exit(1);
        }
        for (C = 0; C < cols; C++) PBC[C] = PyFloat_AsDouble(PyList_GetItem(pbc, C));
    } else PBC = NULL;

    if (origin_list != NULL) {
        origin = (double *)malloc(cols*sizeof(double));
        if (origin == NULL) {
            printf("ERROR - Unable to allocate memory...");
            exit(1);
        }
        for (C = 0; C < cols; C++) origin[C] = PyFloat_AsDouble(PyList_GetItem(origin_list, C));
    } else {
        origin = (double *)malloc(cols*sizeof(double));
        for (C = 0; C < cols; C++) origin[C] = 0.0;
    }

    for (R = 0; R < rows; R++) {
        N_neigh[R] = 0;
        line = PyList_GetItem(list_obj, R);
        for (C = 0; C < cols; C++) {
            pos[R*cols+C] = PyFloat_AsDouble(PyList_GetItem(line, C));
        }
    }

    // Now I have an array "pos" with "rows" rows and "cols" columns
    // Also reads in PBC and origin
    // Here is where I generate the neighbour list...
    // Loop through the different atomic positions
    for (R = 0; R < rows; R++) {

        if (PBC != NULL) {
            SKIN_A = 0;
            for (C = 0; C < cols; C++) {
                if ( (pos[R*cols+C]-cutoff < origin[C]) || (pos[R*cols+C]-cutoff > origin[C]) ) {
                    SKIN_A = 1;
                    break;
                }
            }
        }

        for (R2 = 0; R2 < rows; R2++) {
            if (R == R2) continue;
            dist = 0;
            for (C = 0; C < cols; C++) {
                r = (pos[R*cols+C]-pos[R2*cols+C]);
                dist += r*r;
            }
            // If distance is greater than rsq, see if we're on the 'skin' of our PBC
            // If so, then check the offset distance
            passed = 0;
            if (PBC != NULL && SKIN_A == 1 && dist > rsq) {
                x = pos[R2*cols];
                if (cols > 1) y = pos[R2*cols+1];
                if (cols > 2) z = pos[R2*cols+2];
                if (PBC != NULL && cols > 3) {
                    printf("Error - Code only works for PBC in dim <= 3.");
                    exit (1);
                }

                // DONT JUDGE ME!
                // Checking all combinations for periodic boundary
                // To gain a little extra speed, as soon as one is found to work, don't check the remainder.
                // There's probably a better way around this... but I don't care right now.
                if (passed == 0) { pbc_pos[0] = x+PBC[0]; if (cols > 1) pbc_pos[1] = y; if (cols > 2) pbc_pos[2] = z; }
                if (passed == 0) { dist = 0; for (C = 0; C < cols; C++) { r = (pos[R*cols+C]-pbc_pos[C]); dist += r*r;} if (dist < rsq) {passed=1;} }
                if (passed == 0) { pbc_pos[0] = x-PBC[0]; if (cols > 1) pbc_pos[1] = y; if (cols > 2) pbc_pos[2] = z; }
                if (passed == 0) { dist = 0; for (C = 0; C < cols; C++) { r = (pos[R*cols+C]-pbc_pos[C]); dist += r*r;} if (dist < rsq) {passed=1;} }
                if (passed == 0) { pbc_pos[0] = x; if (cols > 1) pbc_pos[1] = y+PBC[1]; if (cols > 2) pbc_pos[2] = z; }
                if (passed == 0) { dist = 0; for (C = 0; C < cols; C++) { r = (pos[R*cols+C]-pbc_pos[C]); dist += r*r;} if (dist < rsq) {passed=1;} }
                if (passed == 0) { pbc_pos[0] = x; if (cols > 1) pbc_pos[1] = y-PBC[1]; if (cols > 2) pbc_pos[2] = z; }
                if (passed == 0) { dist = 0; for (C = 0; C < cols; C++) { r = (pos[R*cols+C]-pbc_pos[C]); dist += r*r;} if (dist < rsq) {passed=1;} }
                if (passed == 0) { pbc_pos[0] = x; if (cols > 1) pbc_pos[1] = y; if (cols > 2) pbc_pos[2] = z+PBC[2]; }
                if (passed == 0) { dist = 0; for (C = 0; C < cols; C++) { r = (pos[R*cols+C]-pbc_pos[C]); dist += r*r;} if (dist < rsq) {passed=1;} }
                if (passed == 0) { pbc_pos[0] = x; if (cols > 1) pbc_pos[1] = y; if (cols > 2) pbc_pos[2] = z-PBC[2]; }
                if (passed == 0) { dist = 0; for (C = 0; C < cols; C++) { r = (pos[R*cols+C]-pbc_pos[C]); dist += r*r;} if (dist < rsq) {passed=1;} }

                if (passed == 0) { pbc_pos[0] = x+PBC[0]; if (cols > 1) pbc_pos[1] = y+PBC[1]; if (cols > 2) pbc_pos[2] = z; }
                if (passed == 0) { dist = 0; for (C = 0; C < cols; C++) { r = (pos[R*cols+C]-pbc_pos[C]); dist += r*r;} if (dist < rsq) {passed=1;} }
                if (passed == 0) { pbc_pos[0] = x+PBC[0]; if (cols > 1) pbc_pos[1] = y-PBC[1]; if (cols > 2) pbc_pos[2] = z; }
                if (passed == 0) { dist = 0; for (C = 0; C < cols; C++) { r = (pos[R*cols+C]-pbc_pos[C]); dist += r*r;} if (dist < rsq) {passed=1;} }
                if (passed == 0) { pbc_pos[0] = x-PBC[0]; if (cols > 1) pbc_pos[1] = y+PBC[1]; if (cols > 2) pbc_pos[2] = z; }
                if (passed == 0) { dist = 0; for (C = 0; C < cols; C++) { r = (pos[R*cols+C]-pbc_pos[C]); dist += r*r;} if (dist < rsq) {passed=1;} }
                if (passed == 0) { pbc_pos[0] = x-PBC[0]; if (cols > 1) pbc_pos[1] = y-PBC[1]; if (cols > 2) pbc_pos[2] = z; }
                if (passed == 0) { dist = 0; for (C = 0; C < cols; C++) { r = (pos[R*cols+C]-pbc_pos[C]); dist += r*r;} if (dist < rsq) {passed=1;} }

                if (passed == 0) { pbc_pos[0] = x; if (cols > 1) pbc_pos[1] = y+PBC[1]; if (cols > 2) pbc_pos[2] = z+PBC[2]; }
                if (passed == 0) { dist = 0; for (C = 0; C < cols; C++) { r = (pos[R*cols+C]-pbc_pos[C]); dist += r*r;} if (dist < rsq) {passed=1;} }
                if (passed == 0) { pbc_pos[0] = x; if (cols > 1) pbc_pos[1] = y-PBC[1]; if (cols > 2) pbc_pos[2] = z+PBC[2]; }
                if (passed == 0) { dist = 0; for (C = 0; C < cols; C++) { r = (pos[R*cols+C]-pbc_pos[C]); dist += r*r;} if (dist < rsq) {passed=1;} }
                if (passed == 0) { pbc_pos[0] = x; if (cols > 1) pbc_pos[1] = y+PBC[1]; if (cols > 2) pbc_pos[2] = z-PBC[2]; }
                if (passed == 0) { dist = 0; for (C = 0; C < cols; C++) { r = (pos[R*cols+C]-pbc_pos[C]); dist += r*r;} if (dist < rsq) {passed=1;} }
                if (passed == 0) { pbc_pos[0] = x; if (cols > 1) pbc_pos[1] = y-PBC[1]; if (cols > 2) pbc_pos[2] = z-PBC[2]; }
                if (passed == 0) { dist = 0; for (C = 0; C < cols; C++) { r = (pos[R*cols+C]-pbc_pos[C]); dist += r*r;} if (dist < rsq) {passed=1;} }

                if (passed == 0) { pbc_pos[0] = x+PBC[0]; if (cols > 1) pbc_pos[1] = y; if (cols > 2) pbc_pos[2] = z+PBC[2]; }
                if (passed == 0) { dist = 0; for (C = 0; C < cols; C++) { r = (pos[R*cols+C]-pbc_pos[C]); dist += r*r;} if (dist < rsq) {passed=1;} }
                if (passed == 0) { pbc_pos[0] = x-PBC[0]; if (cols > 1) pbc_pos[1] = y; if (cols > 2) pbc_pos[2] = z+PBC[2]; }
                if (passed == 0) { dist = 0; for (C = 0; C < cols; C++) { r = (pos[R*cols+C]-pbc_pos[C]); dist += r*r;} if (dist < rsq) {passed=1;} }
                if (passed == 0) { pbc_pos[0] = x+PBC[0]; if (cols > 1) pbc_pos[1] = y; if (cols > 2) pbc_pos[2] = z-PBC[2]; }
                if (passed == 0) { dist = 0; for (C = 0; C < cols; C++) { r = (pos[R*cols+C]-pbc_pos[C]); dist += r*r;} if (dist < rsq) {passed=1;} }
                if (passed == 0) { pbc_pos[0] = x-PBC[0]; if (cols > 1) pbc_pos[1] = y; if (cols > 2) pbc_pos[2] = z-PBC[2]; }
                if (passed == 0) { dist = 0; for (C = 0; C < cols; C++) { r = (pos[R*cols+C]-pbc_pos[C]); dist += r*r;} if (dist < rsq) {passed=1;} }

                if (passed == 0) { pbc_pos[0] = x+PBC[0]; if (cols > 1) pbc_pos[1] = y+PBC[1]; if (cols > 2) pbc_pos[2] = z+PBC[2]; }
                if (passed == 0) { dist = 0; for (C = 0; C < cols; C++) { r = (pos[R*cols+C]-pbc_pos[C]); dist += r*r;} if (dist < rsq) {passed=1;} }
                if (passed == 0) { pbc_pos[0] = x+PBC[0]; if (cols > 1) pbc_pos[1] = y+PBC[1]; if (cols > 2) pbc_pos[2] = z-PBC[2]; }
                if (passed == 0) { dist = 0; for (C = 0; C < cols; C++) { r = (pos[R*cols+C]-pbc_pos[C]); dist += r*r;} if (dist < rsq) {passed=1;} }
                if (passed == 0) { pbc_pos[0] = x+PBC[0]; if (cols > 1) pbc_pos[1] = y-PBC[1]; if (cols > 2) pbc_pos[2] = z+PBC[2]; }
                if (passed == 0) { dist = 0; for (C = 0; C < cols; C++) { r = (pos[R*cols+C]-pbc_pos[C]); dist += r*r;} if (dist < rsq) {passed=1;} }
                if (passed == 0) { pbc_pos[0] = x+PBC[0]; if (cols > 1) pbc_pos[1] = y-PBC[1]; if (cols > 2) pbc_pos[2] = z-PBC[2]; }
                if (passed == 0) { dist = 0; for (C = 0; C < cols; C++) { r = (pos[R*cols+C]-pbc_pos[C]); dist += r*r;} if (dist < rsq) {passed=1;} }
                if (passed == 0) { pbc_pos[0] = x-PBC[0]; if (cols > 1) pbc_pos[1] = y+PBC[1]; if (cols > 2) pbc_pos[2] = z+PBC[2]; }
                if (passed == 0) { dist = 0; for (C = 0; C < cols; C++) { r = (pos[R*cols+C]-pbc_pos[C]); dist += r*r;} if (dist < rsq) {passed=1;} }
                if (passed == 0) { pbc_pos[0] = x-PBC[0]; if (cols > 1) pbc_pos[1] = y+PBC[1]; if (cols > 2) pbc_pos[2] = z-PBC[2]; }
                if (passed == 0) { dist = 0; for (C = 0; C < cols; C++) { r = (pos[R*cols+C]-pbc_pos[C]); dist += r*r;} if (dist < rsq) {passed=1;} }
                if (passed == 0) { pbc_pos[0] = x-PBC[0]; if (cols > 1) pbc_pos[1] = y-PBC[1]; if (cols > 2) pbc_pos[2] = z+PBC[2]; }
                if (passed == 0) { dist = 0; for (C = 0; C < cols; C++) { r = (pos[R*cols+C]-pbc_pos[C]); dist += r*r;} if (dist < rsq) {passed=1;} }
                if (passed == 0) { pbc_pos[0] = x-PBC[0]; if (cols > 1) pbc_pos[1] = y-PBC[1]; if (cols > 2) pbc_pos[2] = z-PBC[2]; }
                if (passed == 0) { dist = 0; for (C = 0; C < cols; C++) { r = (pos[R*cols+C]-pbc_pos[C]); dist += r*r;} if (dist < rsq) {passed=1;} }
            }

            if (passed==1 || dist < rsq) {
                // Ensure we have memory...
                if ( (MAX_SIZE-N_count) < TOO_CLOSE) {
                    // Resize neighs
                    MAX_SIZE += DELTA;
                    resized = (int *)realloc(neighs, MAX_SIZE*sizeof(int));
                    if ( resized == NULL ) {
                        free(pos); free(N_neigh); free(neighs); free(origin);
                        if (PBC != NULL) free(PBC);
                        printf("\nERROR - Unable to resize arrays.");
                        exit(1);
                    }
                    neighs = resized;
                    resized = NULL;
                }
                N_neigh[R]++;
                neighs[N_count++] = R2;
            }
        }
    }

    // Regenerate a python list of lists for the various neighbours
    neighbours = PyList_New(rows);

    counter = 0;
    for (R = 0; R < rows; R++) {
        neigh_line = PyList_New(N_neigh[R]);
        for (n = 0; n < N_neigh[R]; n++) {
            PyList_SetItem(neigh_line,n,PyInt_FromSsize_t(neighs[counter++]));
        }
        PyList_SetItem(neighbours, R, neigh_line);
    }

    free(pos); free(N_neigh); free(neighs); free(origin);
    if (PBC != NULL) free(PBC);
    return neighbours;
    //return Py_BuildValue("d", 0);
}

// This binds our c and python function names
static PyMethodDef neigh_methods[] = {
  {"gen_neigh", py_gen_neigh, METH_VARARGS|METH_KEYWORDS},
  {NULL, NULL, NULL} // Sentinel, essentially look at what we're passing
};

// Python calls this to initialize module (essentially, this is our main function)
PyMODINIT_FUNC initneigh(void)
{
  (void) Py_InitModule("neigh", neigh_methods);
}
