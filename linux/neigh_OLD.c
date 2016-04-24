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

static PyObject* py_gen_neigh(PyObject* self, PyObject* args) {
    // Initialize variables
    int R, R2, C, rows, cols, MAX_SIZE, TOO_CLOSE, DELTA, n;
    int held, N_hold_low, N_count, N_hold_count, *N_neigh, *N_neigh_hold, *neighs, *neigh_hold;
    int *resized1, *resized2;
    double *pos, cutoff, rsq, dist, x;
    PyObject * list_obj;
    PyObject * line;
    PyObject * neighbours;
    PyObject * neigh_line;

    // Starting size for the number of neighbours we want, and the delta shift when approaching our limit
    resized1 = resized2 = NULL;
    MAX_SIZE = 1000;
    TOO_CLOSE = 100;
    DELTA = 1000;
    N_count = 0;
    N_hold_count = 0;

    // Ensure it is a list being passed
    // NOTE! The format string "O!" says that I am reading in an object and requiring it to be a certain type
    // Therefore, I give this function two variables.  The first is the type I want to read in and the second
    // Is the memory location for the list itself.
    // For more info, go here: https://docs.python.org/2/c-api/arg.html
    if ( !PyArg_ParseTuple( args, "O!d", &PyList_Type, &list_obj, &cutoff ) ) return NULL;
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
    N_neigh_hold = (int *)malloc(rows*sizeof(int));
    neighs = (int *)malloc(MAX_SIZE*sizeof(int));
    neigh_hold = (int *)malloc(MAX_SIZE*sizeof(int));

    if (pos == NULL ||
        N_neigh == NULL ||
        N_neigh_hold == NULL ||
        neighs == NULL ||
        neigh_hold == NULL) {
        printf("ERROR - Unable to allocate memory...");
        exit(1);
    }

    for (R = 0; R < rows; R++) {
        N_neigh[R] = 0;
        N_neigh_hold[R] = 0;
        line = PyList_GetItem(list_obj, R);
        for (C = 0; C < cols; C++) {
            pos[R*cols+C] = PyFloat_AsDouble(PyList_GetItem(line, C));
        }
    }

    // Now I have an array "pos" with "rows" rows and "cols" columns
    // Here is where I generate the neighbour list...
    // Loop through the different atomic positions
    for (R = 0; R < rows; R++) {
        // If R already in nearest neighbour of R2, we stored it in N_neigh_hold.
        // So, save it here
        if (N_neigh_hold[R] > 0) {
            N_hold_low = 0;
            for (held = 0; held < R; held++) N_hold_low += N_neigh_hold[held];
        }
        for (held = 0; held < N_neigh_hold[R]; held++) {
            if ( ( (MAX_SIZE-N_count) < TOO_CLOSE) || ( (MAX_SIZE-N_hold_count) < TOO_CLOSE) ) {
                // Resize neighs and neigh_hold
                MAX_SIZE += DELTA;
                resized1 = (int *)realloc(neighs, MAX_SIZE*sizeof(int));
                resized2 = (int *)realloc(neigh_hold, MAX_SIZE*sizeof(int));
                if ( (resized1 == NULL) || (resized2 == NULL) ) {
                    free(pos); free(N_neigh); free(neighs); free(N_neigh_hold); free(neigh_hold);
                    printf("\nERROR - Unable to resize arrays.");
                    exit(1);
                }
                neighs = resized1;
                neigh_hold = resized2;
                resized1 = resized2 = NULL;
            }

            N_neigh[R]++;
            neighs[N_count++] = neigh_hold[N_hold_low + held];
        }
        for (R2 = R+1; R2 < rows; R2++) {
            dist = 0;
            for (C = 0; C < cols; C++) {
                x = (pos[R*cols+C]-pos[R2*cols+C]);
                dist += x*x;
            }
            if (dist < rsq) {
                // Ensure we have memory...
                if ( ( (MAX_SIZE-N_count) < TOO_CLOSE) || ( (MAX_SIZE-N_hold_count) < TOO_CLOSE) ) {
                    // Resize neighs and neigh_hold
                    MAX_SIZE += DELTA;
                    resized1 = (int *)realloc(neighs, MAX_SIZE*sizeof(int));
                    resized2 = (int *)realloc(neigh_hold, MAX_SIZE*sizeof(int));
                    if ( (resized1 == NULL) || (resized2 == NULL) ) {
                        free(pos); free(N_neigh); free(neighs); free(N_neigh_hold); free(neigh_hold);
                        printf("\nERROR - Unable to resize arrays.");
                        exit(1);
                    }
                    neighs = resized1;
                    neigh_hold = resized2;
                    resized1 = resized2 = NULL;
                }
                N_neigh[R]++;
                neighs[N_count++] = R2;
                N_neigh_hold[R2]++;
                neigh_hold[N_hold_count++] = R;
            }
        }
        // If we're nearing our buffer limits, resize them
        if ( ( (MAX_SIZE-N_count) < TOO_CLOSE) || ( (MAX_SIZE-N_hold_count) < TOO_CLOSE) ) {
            // Resize neighs and neigh_hold
            MAX_SIZE += DELTA;
            resized1 = (int *)realloc(neighs, MAX_SIZE*sizeof(int));
            resized2 = (int *)realloc(neigh_hold, MAX_SIZE*sizeof(int));
            if ( (resized1 == NULL) || (resized2 == NULL) ) {
                free(pos); free(N_neigh); free(neighs); free(N_neigh_hold); free(neigh_hold);
                printf("\nERROR - Unable to resize arrays.");
                exit(1);
            }
            neighs = resized1;
            neigh_hold = resized2;
            resized1 = resized2 = NULL;
        }
    }

    // Regenerate a python list of lists for the various neighbours
    neighbours = PyList_New(rows);
    // As N_neigh and neighs holds everything, we can use N_hold_count freely
    N_hold_count = 0;
    for (R = 0; R < rows; R++) {
        neigh_line = PyList_New(N_neigh[R]);
        for (n = 0; n < N_neigh[R]; n++) {
            PyList_SetItem(neigh_line,n,PyInt_FromSsize_t(neighs[N_hold_count++]));
        }
        PyList_SetItem(neighbours, R, neigh_line);
    }

    free(pos); free(N_neigh); free(neighs); free(N_neigh_hold); free(neigh_hold);
    return neighbours;
    //return Py_BuildValue("d", 0);
}

// This binds our c and python function names
static PyMethodDef neigh_methods[] = {
  {"gen_neigh", py_gen_neigh, METH_VARARGS},
  {NULL} // Sentinel, essentially look at what we're passing
};

// Python calls this to initialize module (essentially, this is our main function)
void initneigh()
{
  (void) Py_InitModule("neigh", neigh_methods);
}
