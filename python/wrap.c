#include <Python.h>
#include <structmember.h>
#include <numpy/arrayobject.h>

#include <stdlib.h>
#include <inttypes.h>
#include <float.h>

#include "eb.h"

static PyObject *wrap_model (PyObject *self, PyObject *args, PyObject *keywds) {
  static char *kwlist[] = { "parm", "t", "typ",
                            "flags", "out", "ol1", "ol2",
                            NULL };
  PyObject *parmarg = NULL, *targ = NULL, *typarg = NULL;
  PyObject *parmarr = NULL, *tarr = NULL, *typarr = NULL;

  PyObject *outarg = NULL, *ol1arg = NULL, *ol2arg = NULL;
  PyObject *outarr = NULL, *ol1arr = NULL, *ol2arr = NULL;
  
  double *ol1 = NULL, *ol2 = NULL;

  int flags = 0, npt;

  /* Get arguments */
  if(!PyArg_ParseTupleAndKeywords(args, keywds, "OOO|iOOO", kwlist,
                                  &parmarg, &targ, &typarg,
                                  &flags, &outarg, &ol1arg, &ol2arg))
    goto error;

  /* Get array arguments */
  parmarr = PyArray_FROM_OTF(parmarg, NPY_DOUBLE, NPY_IN_ARRAY | NPY_FORCECAST);
  if(!parmarr)
    goto error;

  tarr = PyArray_FROM_OTF(targ, NPY_DOUBLE, NPY_IN_ARRAY | NPY_FORCECAST);
  if(!tarr)
    goto error;

  typarr = PyArray_FROM_OTF(typarg, NPY_UINT8, NPY_IN_ARRAY | NPY_FORCECAST);
  if(!typarr)
    goto error;

  /* Check sizes */
  if(PyArray_Size(parmarr) < EB_NPAR) {
    PyErr_SetString(PyExc_IndexError,
                    "Parameter vector is too short");
    goto error;
  }

  npt = PyArray_Size(tarr);

  if(PyArray_Size(typarr) != npt) {
    PyErr_SetString(PyExc_IndexError,
                    "Lengths of 't' and 'typ' arrays do not match");
    goto error;
  }

  /* Optional arguments */
  if(ol1arg) {
    ol1arr = PyArray_FROM_OTF(ol1arg, NPY_DOUBLE, NPY_IN_ARRAY | NPY_FORCECAST);
    if(!ol1arr)
      goto error;

    if(PyArray_Size(ol1arr) < npt) {
      PyErr_SetString(PyExc_IndexError,
                      "'ol1' vector is too short");
      goto error;
    }

    ol1 = (double *) PyArray_DATA(ol1arr);
  }

  if(ol2arg) {
    ol2arr = PyArray_FROM_OTF(ol2arg, NPY_DOUBLE, NPY_IN_ARRAY | NPY_FORCECAST);
    if(!ol2arr)
      goto error;

    if(PyArray_Size(ol2arr) < npt) {
      PyErr_SetString(PyExc_IndexError,
                      "'ol2' vector is too short");
      goto error;
    }

    ol2 = (double *) PyArray_DATA(ol2arr);
  }

  /* Output array */
  if(outarg) {
    outarr = PyArray_FROM_OTF(outarg, NPY_DOUBLE, NPY_OUT_ARRAY);
    if(!outarr)
      goto error;

    if(PyArray_Size(outarr) < npt) {
      PyErr_SetString(PyExc_IndexError,
                      "'out' vector is too short");
      goto error;
    }
  }
  else {
    outarr = PyArray_SimpleNew(PyArray_NDIM(tarr),
                               PyArray_DIMS(tarr),
                               NPY_DOUBLE);
    if(!outarr)
      goto error;
  }

  /* Compute result */
  eb_model_dbl((double *) PyArray_DATA(parmarr),
               (double *) PyArray_DATA(tarr),
               ol1, ol2,
               (uint8_t *) PyArray_DATA(typarr),
               (double *) PyArray_DATA(outarr),
               NULL,
               flags, npt);

  /* Done */
  Py_DECREF(parmarr);
  Py_DECREF(tarr);
  Py_DECREF(typarr);

  if(ol1arg) {
    Py_DECREF(ol1arr);
  }
  if(ol2arg) {
    Py_DECREF(ol2arr);
  }

  if(outarg) {
    Py_DECREF(outarr);
    return(Py_None);
  }
  else {
    return(PyArray_Return((PyArrayObject *) outarr));
  }

 error:
  Py_XDECREF(parmarr);
  Py_XDECREF(tarr);
  Py_XDECREF(typarr);
  Py_XDECREF(ol1arr);
  Py_XDECREF(ol2arr);
  PyArray_XDECREF_ERR((PyArrayObject *) outarr);

  return(NULL);
}

static PyObject *wrap_phisec (PyObject *self, PyObject *args) {
  int narg;

  PyObject *parmarg = NULL;
  PyObject *parmarr = NULL;
  double *parm;

  double esinw, ecosw;

  double phi = 0;

  /* How many non-keyword arguments? */
  narg = PyTuple_Size(args);

  if(narg == 1) {
    /* Parameter vector argument */
    if(!PyArg_ParseTuple(args, "O", &parmarg))
      return(NULL);
    
    /* Get array argument */
    parmarr = PyArray_FROM_OTF(parmarg, NPY_DOUBLE, NPY_IN_ARRAY | NPY_FORCECAST);
    if(!parmarr)
      goto error;
    
    /* Check sizes */
    if(PyArray_Size(parmarr) < EB_NPAR) {
      PyErr_SetString(PyExc_IndexError,
                      "Parameter vector is too short");
      goto error;
    }
    
    parm = (double *) PyArray_DATA(parmarr);
    
    phi = eb_phisec(parm[EB_PAR_ESINW], parm[EB_PAR_ECOSW]);
    
    Py_DECREF(parmarr);
  }
  else if(narg == 2) {
    /* esinw, ecosw given directly */
    if(!PyArg_ParseTuple(args, "dd", &esinw, &ecosw))
      return(NULL);

    phi = eb_phisec(esinw, ecosw);
  }
  else {
    /* Invalid */
    PyErr_SetString(PyExc_TypeError,
                    "Usage: phisec(esinw, ecosw) or phisec(parm)");
    goto error;
  }

  return(Py_BuildValue("d", phi));

 error:
  Py_XDECREF(parmarr);

  return(NULL);
}

static PyObject *wrap_phicont (PyObject *self, PyObject *args) {
  int narg;

  PyObject *parmarg = NULL;
  PyObject *parmarr = NULL;
  double *parm;

  double esinw, ecosw, cosi, rasum;

  double phi[4];

  /* How many non-keyword arguments? */
  narg = PyTuple_Size(args);

  if(narg == 1) {
    /* Parameter vector argument */
    if(!PyArg_ParseTuple(args, "O", &parmarg))
      return(NULL);
    
    /* Get array argument */
    parmarr = PyArray_FROM_OTF(parmarg, NPY_DOUBLE, NPY_IN_ARRAY | NPY_FORCECAST);
    if(!parmarr)
      goto error;
    
    /* Check sizes */
    if(PyArray_Size(parmarr) < EB_NPAR) {
      PyErr_SetString(PyExc_IndexError,
                      "Parameter vector is too short");
      goto error;
    }
    
    parm = (double *) PyArray_DATA(parmarr);
    
    eb_phicont(parm[EB_PAR_ESINW], parm[EB_PAR_ECOSW],
               parm[EB_PAR_COSI], parm[EB_PAR_RASUM],
               phi);
    
    Py_DECREF(parmarr);
  }
  else if(narg == 4) {
    /* Parameters given directly */
    if(!PyArg_ParseTuple(args, "dddd", &esinw, &ecosw, &cosi, &rasum))
      return(NULL);

    eb_phicont(esinw, ecosw, cosi, rasum, phi);
  }
  else {
    /* Invalid */
    PyErr_SetString(PyExc_TypeError,
                    "Usage: phisec(esinw, ecosw, cosi, rsum/a) or phisec(parm)");
    goto error;
  }

  return(Py_BuildValue("dddd", phi[0], phi[1], phi[2], phi[3]));

 error:
  Py_XDECREF(parmarr);

  return(NULL);
}

static PyObject *wrap_getvder (PyObject *self, PyObject *args, PyObject *keywds) {
  static char *kwlist[] = { "parm", "gamma", "ktot", NULL };
  PyObject *parmarg = NULL;
  PyObject *parmarr = NULL, *outarr = NULL;
  npy_intp outdims[1] = { EB_NDER };
  double gamma = 0.0, ktot = 0.0;

  /* Get arguments */
  if(!PyArg_ParseTupleAndKeywords(args, keywds, "O|dd", kwlist,
                                  &parmarg, &gamma, &ktot))
    goto error;

  /* Get array argument */
  parmarr = PyArray_FROM_OTF(parmarg, NPY_DOUBLE, NPY_IN_ARRAY | NPY_FORCECAST);
  if(!parmarr)
    goto error;

  /* Check sizes */
  if(PyArray_Size(parmarr) < EB_NPAR) {
    PyErr_SetString(PyExc_IndexError,
                    "Parameter vector is too short");
    goto error;
  }

  /* Create output array */
  outarr = PyArray_SimpleNew(1, outdims, NPY_DOUBLE);
  if(!outarr)
    goto error;

  /* Compute result */
  eb_getvder((double *) PyArray_DATA(parmarr),
             gamma, ktot,
             (double *) PyArray_DATA(outarr));

  /* Done */
  Py_DECREF(parmarr);

  return(PyArray_Return((PyArrayObject *) outarr));

 error:
  Py_XDECREF(parmarr);
  PyArray_XDECREF_ERR((PyArrayObject *) outarr);

  return(NULL);
}

static PyMethodDef eb_methods[] = {
  { "model", (PyCFunction) wrap_model, METH_VARARGS | METH_KEYWORDS,
    "y = model(parm, t, typ, flags=0)\n\n"
    "Main subroutine to compute the model.\n\n"
    "Required arguments:\n"
    "parm  -- double parameter vector, length NPAR (see PAR_*).\n"
    "t     -- double array of times at which to compute the model.\n"
    "typ   -- uint8  quantity (OBS_*) to compute for each time.\n\n"
    "Optional arguments:\n"
    "flags -- integer bitmask of FLAG_* values.\n"
    "out   -- store model in user-supplied double array.\n"
    "ol1   -- fractional adjustment to light of star 1.\n"
    "ol2   -- fractional adjustment to light of star 2.\n\n"
    "Arguments t and typ must be the same shape, and for best\n"
    "performance these arrays should use contiguous storage.\n"
    "Returned array y will be the same shape as t and contains\n"
    "the corresponding model values.\n\n"
    "The only requirement on units for the array 't' is they\n"
    "must be the same as parm[PAR_P] and parm[PAR_T0].  However,\n"
    "if you plan to use the routine getvder() on the same\n"
    "parameter vector, please note that the units of parm[PAR_P]\n"
    "must be Julian days.\n\n"
  },
  { "phisec", (PyCFunction) wrap_phisec, METH_VARARGS,
    "phi = phisec(esinw, ecosw)\n"
    "phi = phisec(parm)\n\n"
    "Compute [0,1] phase offset of superior conjunction."
  },
  { "phicont", (PyCFunction) wrap_phicont, METH_VARARGS,
    "(ps, pe, ss, se) = phicont(esinw, ecosw, cosi, (R_1+R_2)/a)\n"
    "(ps, pe, ss, se) = phicont(parm)\n\n"
    "Compute [0,1] phases of first and last contact points for\n"
    "both eclipses.  If there is no eclipse, the phases of the\n"
    "conjunctions are returned instead, so this condition can be\n"
    "detected by checking if the duration is zero.\n\n"
    "NOTE: light travel time is not taken into account.\n"
  },
  { "getvder", (PyCFunction) wrap_getvder, METH_VARARGS | METH_KEYWORDS,
    "vder = getvder(parm, gamma=0, ktot=0)\n\n"
    "Compute vector of derived parameters.\n\n"
    "Required arguments:\n"
    "parm  -- double parameter vector, length NPAR (see PAR_*).\n"
    "         parm[PAR_P] must be in Julian days.\n\n"
    "Optional arguments:\n"
    "gamma -- Barycentric radial velocity of center of mass (km/s).\n"
    "ktot  -- Sum of radial velocity semiamplitudes (km/s).\n\n"
    "Returns a vector of length NDER containing the derived\n"
    "parameters.  If the optional radial velocity parameters are\n"
    "not given, the parameters depending on mass and system scale\n"
    "(PAR_A - PAR_TCIRC) in the returned vector will be invalid.\n\n"
    "NOTE: for precise work, it is important to ensure consistent\n"
    "values are used for the physical constants throughout the\n"
    "analysis.  The values used by this routine are exposed in\n"
    "the module as constants GMSUN, AU, LIGHT and RSUN in MKS\n"
    "units and are the IAU 2009 TDB-compatible values for GMSUN\n"
    "and AU, and RSUN from Brown & Christensen-Dalsgaard 1998 (as\n"
    "adopted by Cox 2000 in Allen's Astrophysical Quantities,\n"
    "4th Ed.).\n"
  },
  { NULL, NULL, 0, NULL }
};

PyMODINIT_FUNC initeb (void) {
  PyObject *m, *o, *t;

#define MAKECONST(m) { #m, EB_##m }

  struct {
    char *name;
    int value;
  } iconst[] = {
    MAKECONST(PAR_J),
    MAKECONST(PAR_RASUM),
    MAKECONST(PAR_RR),
    MAKECONST(PAR_LDLIN1),
    MAKECONST(PAR_LDLIN2),
    MAKECONST(PAR_COSI),
    MAKECONST(PAR_ECOSW),
    MAKECONST(PAR_ESINW),
    MAKECONST(PAR_GD1),
    MAKECONST(PAR_GD2),
    MAKECONST(PAR_REFL1),
    MAKECONST(PAR_REFL2),
    MAKECONST(PAR_Q),
    MAKECONST(PAR_TIDANG),
    MAKECONST(PAR_L3),
    MAKECONST(PAR_PHI0),
    MAKECONST(PAR_M0),
    MAKECONST(PAR_INTEG),
    MAKECONST(PAR_P),
    MAKECONST(PAR_T0),
    MAKECONST(PAR_LDNON1),
    MAKECONST(PAR_LDNON2),
    MAKECONST(PAR_ROT1),
    MAKECONST(PAR_ROT2),
    MAKECONST(PAR_FSPOT1),
    MAKECONST(PAR_FSPOT2),
    MAKECONST(PAR_OOE1O),
    MAKECONST(PAR_OOE11A),
    MAKECONST(PAR_OOE11B),
    MAKECONST(PAR_OOE12A),
    MAKECONST(PAR_OOE12B),
    MAKECONST(PAR_OOE2O),
    MAKECONST(PAR_OOE21A),
    MAKECONST(PAR_OOE21B),
    MAKECONST(PAR_OOE22A),
    MAKECONST(PAR_OOE22B),
    MAKECONST(PAR_CLTT),
    MAKECONST(NPAR),
    MAKECONST(FLAG_REFL),
    MAKECONST(FLAG_PHI),
    MAKECONST(OBS_MAG),
    MAKECONST(OBS_LIGHT),
    MAKECONST(OBS_LRAT),
    MAKECONST(OBS_AVLR),
    MAKECONST(OBS_VRAD1),
    MAKECONST(OBS_VRAD2),
    MAKECONST(OBS_PSS),
    MAKECONST(OBS_A),
    MAKECONST(PAR_I),
    MAKECONST(PAR_R1A),
    MAKECONST(PAR_R2A),
    MAKECONST(PAR_E),
    MAKECONST(PAR_OMEGA),
    MAKECONST(PAR_A),
    MAKECONST(PAR_MTOT),
    MAKECONST(PAR_M1),
    MAKECONST(PAR_M2),
    MAKECONST(PAR_RTOT),
    MAKECONST(PAR_R1),
    MAKECONST(PAR_R2),
    MAKECONST(PAR_LOGG1),
    MAKECONST(PAR_LOGG2),
    MAKECONST(PAR_VSYNC1),
    MAKECONST(PAR_VSYNC2),
    MAKECONST(PAR_TSYNC),
    MAKECONST(PAR_TCIRC),
    MAKECONST(PAR_TSEC),
    MAKECONST(PAR_DURPRI),
    MAKECONST(PAR_DURSEC),
    MAKECONST(NDER),
  };

  struct {
    char *name;
    double value;
  } dconst[] = {
    MAKECONST(GMSUN),
    MAKECONST(AU),
    MAKECONST(LIGHT),
    MAKECONST(RSUN),
  };

  int c, nc;

  /* Init module */
  m = Py_InitModule("eb", eb_methods);
  if(!m)
    return;

  /* Import numpy */
  import_array();

  /* Create constants */
  nc = sizeof(iconst) / sizeof(iconst[0]);

  for(c = 0; c < nc; c++) {
    o = PyInt_FromLong(iconst[c].value);
    if(o)
      PyModule_AddObject(m, iconst[c].name, o);
  }

  nc = sizeof(dconst) / sizeof(dconst[0]);

  for(c = 0; c < nc; c++) {
    o = PyFloat_FromDouble(dconst[c].value);
    if(o)
      PyModule_AddObject(m, dconst[c].name, o);
  }

  /* Create tuples for name and unit lists */
  o = PyTuple_New(EB_NPAR);
  if(o) {
    for(c = 0; c < EB_NPAR; c++) {
      t = PyString_FromString(eb_parnames[c]);
      if(t)
        PyTuple_SetItem(o, c, t);
    }
    PyModule_AddObject(m, "parnames", o);
  }

  o = PyTuple_New(EB_NPAR);
  if(o) {
    for(c = 0; c < EB_NPAR; c++) {
      t = PyString_FromString(eb_parunits[c]);
      if(t)
        PyTuple_SetItem(o, c, t);
    }
    PyModule_AddObject(m, "parunits", o);
  }

  o = PyTuple_New(EB_NDER);
  if(o) {
    for(c = 0; c < EB_NDER; c++) {
      t = PyString_FromString(eb_dernames[c]);
      if(t)
        PyTuple_SetItem(o, c, t);
    }
    PyModule_AddObject(m, "dernames", o);
  }

  o = PyTuple_New(EB_NDER);
  if(o) {
    for(c = 0; c < EB_NDER; c++) {
      t = PyString_FromString(eb_derunits[c]);
      if(t)
        PyTuple_SetItem(o, c, t);
    }
    PyModule_AddObject(m, "derunits", o);
  }
}
