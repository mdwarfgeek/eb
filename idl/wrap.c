#include <stdio.h>
#include <stdlib.h>

#include <idl_export.h>

#include "eb.h"

struct wrap_eb_model_kw_result {
  IDL_KW_RESULT_FIRST_FIELD;

  IDL_LONG flags;

  IDL_VARIABLE ol1;
  int have_ol1;

  IDL_VARIABLE ol2;
  int have_ol2;

  IDL_VPTR out;
  int have_out;
};

static IDL_KW_PAR wrap_eb_model_kw_pars[] = {
  IDL_KW_FAST_SCAN,  /* must be in lexical order */
  { "FLAGS", IDL_TYP_LONG, 1, IDL_KW_ZERO,
    NULL,
    IDL_KW_OFFSETOF2(struct wrap_eb_model_kw_result, flags) },
  { "OL1", IDL_TYP_UNDEF, 1, IDL_KW_VIN,
    IDL_KW_OFFSETOF2(struct wrap_eb_model_kw_result, have_ol1),
    IDL_KW_OFFSETOF2(struct wrap_eb_model_kw_result, ol1) },
  { "OL2", IDL_TYP_UNDEF, 1, IDL_KW_VIN,
    IDL_KW_OFFSETOF2(struct wrap_eb_model_kw_result, have_ol2),
    IDL_KW_OFFSETOF2(struct wrap_eb_model_kw_result, ol2) },
  { "OUT", IDL_TYP_UNDEF, 1, IDL_KW_OUT,
    IDL_KW_OFFSETOF2(struct wrap_eb_model_kw_result, have_out),
    IDL_KW_OFFSETOF2(struct wrap_eb_model_kw_result, out) },
  { NULL }
};

IDL_VPTR wrap_eb_model (int argc, IDL_VPTR *argv, char *argk) {
  struct wrap_eb_model_kw_result kw;
  int nplain;

  IDL_VPTR parm, t, typ;
  double *parm_p, *t_p;
  unsigned char *typ_p;

  IDL_VPTR ol1, ol2;
  double *ol1_p = (double *) NULL;
  double *ol2_p = (double *) NULL;

  int npt;

  IDL_VPTR out;
  double *out_p;

  /* Process keywords */
  nplain = IDL_KWProcessByOffset(argc, argv, argk, wrap_eb_model_kw_pars,
                                 (IDL_VPTR *) NULL, 1, &kw);

  /* Check all required arguments are present */
  if(nplain != 3) {
    IDL_Message(IDL_M_NAMED_GENERIC, IDL_MSG_LONGJMP,
                "Usage: y = eb_model(parm, t, typ, [flags=, ol1=, ol2=])\n");
    return(IDL_GettmpInt(-1));
  }

  /* Unpack array arguments */
  parm = argv[0];
  IDL_ENSURE_SIMPLE(parm);
  IDL_ENSURE_ARRAY(parm);

  if(parm->type != IDL_TYP_DOUBLE)
    parm = IDL_CvtDbl(1, &parm);

  parm_p = (double *) parm->value.arr->data;

  t = argv[1];
  IDL_ENSURE_SIMPLE(t);
  IDL_ENSURE_ARRAY(t);

  if(t->type != IDL_TYP_DOUBLE)
    t = IDL_CvtDbl(1, &t);

  t_p = (double *) t->value.arr->data;

  typ = argv[2];
  IDL_ENSURE_SIMPLE(typ);
  IDL_ENSURE_ARRAY(typ);

  if(typ->type != IDL_TYP_BYTE)
    typ = IDL_CvtByte(1, &t);

  typ_p = (unsigned char *) typ->value.arr->data;

  /* Check sizes */
  if(parm->value.arr->n_elts < EB_NPAR) {
    IDL_Message(IDL_M_NAMED_GENERIC, IDL_MSG_LONGJMP,
                "Parameter vector is too short\n");
    goto error;
  }

  npt = t->value.arr->n_elts;

  if(typ->value.arr->n_elts != npt) {
    IDL_Message(IDL_M_NAMED_GENERIC, IDL_MSG_LONGJMP,
                "Lengths of 't' and 'typ' arrays do not match\n");
    goto error;
  }

  /* Optional arguments */
  if(kw.have_ol1) {
    ol1 = &(kw.ol1);

    IDL_ENSURE_SIMPLE(ol1);
    IDL_ENSURE_ARRAY(ol1);

    if(ol1->value.arr->n_elts < npt) {
      IDL_Message(IDL_M_NAMED_GENERIC, IDL_MSG_LONGJMP,
                  "Lengths of 't' and 'ol1' arrays do not match\n");
      goto error;
    }

    if(ol1->type != IDL_TYP_DOUBLE)
      ol1 = IDL_CvtDbl(1, &ol1);

    ol1_p = (double *) ol1->value.arr->data;
  }

  if(kw.have_ol2) {
    ol2 = &(kw.ol2);

    IDL_ENSURE_SIMPLE(ol2);
    IDL_ENSURE_ARRAY(ol2);

    if(ol2->value.arr->n_elts < npt) {
      IDL_Message(IDL_M_NAMED_GENERIC, IDL_MSG_LONGJMP,
                  "Lengths of 't' and 'ol2' arrays do not match\n");
      goto error;
    }

    if(ol2->type != IDL_TYP_DOUBLE)
      ol2 = IDL_CvtDbl(1, &ol2);

    ol2_p = (double *) ol2->value.arr->data;
  }

  if(kw.have_out) {
    out = kw.out;

    IDL_ENSURE_SIMPLE(out);
    IDL_ENSURE_ARRAY(out);

    if(out->value.arr->n_elts < npt) {
      IDL_Message(IDL_M_NAMED_GENERIC, IDL_MSG_LONGJMP,
                  "Lengths of 't' and 'out' arrays do not match\n");
      goto error;
    }

    if(out->type != IDL_TYP_DOUBLE) {
      IDL_Message(IDL_M_NAMED_GENERIC, IDL_MSG_LONGJMP,
                  "Type of 'out' array must be double\n");
      goto error;
    }

    out_p = (double *) out->value.arr->data;
  }
  else {
    /* Make output array (same shape as t) */
    out_p = (double *) IDL_MakeTempArray(IDL_TYP_DOUBLE,
                                         t->value.arr->n_dim,
                                         t->value.arr->dim,
                                         IDL_ARR_INI_NOP,
                                         &out);
  }

  /* Compute result */
  eb_model_dbl(parm_p,
               t_p,
               ol1_p, ol2_p,
               typ_p,
               out_p,
               NULL,
               kw.flags, npt);

  /* Free any temporary copies of arguments */
  if(parm != argv[0])
    IDL_DELTMP(parm);

  if(t != argv[1])
    IDL_DELTMP(t);

  if(typ != argv[2])
    IDL_DELTMP(typ);

  if(kw.have_ol1 && ol1 != &(kw.ol1))
    IDL_DELTMP(ol1);

  if(kw.have_ol2 && ol2 != &(kw.ol2))
    IDL_DELTMP(ol2);

  if(kw.have_out) {
    IDL_KW_FREE;
    return(IDL_GettmpInt(npt));
  }
  else {
    IDL_KW_FREE;
    return(out);
  }

 error:
  IDL_KW_FREE;

  return(IDL_GettmpInt(-1));
}

IDL_VPTR wrap_eb_phisec (int argc, IDL_VPTR *argv) {
  IDL_VPTR esinw, ecosw;
  double phi;

  /* Check all required arguments are present */
  if(argc != 2) {
    IDL_Message(IDL_M_NAMED_GENERIC, IDL_MSG_LONGJMP,
                "Usage: y = eb_phisec(esinw, ecosw)\n");
    return(IDL_GettmpInt(-1));
  }

  esinw = argv[0];
  
  IDL_ENSURE_SIMPLE(esinw);
  IDL_ENSURE_SCALAR(esinw);

  if(esinw->type != IDL_TYP_DOUBLE)
    esinw = IDL_CvtDbl(1, &esinw);

  ecosw = argv[1];
  
  IDL_ENSURE_SIMPLE(ecosw);
  IDL_ENSURE_SCALAR(ecosw);

  if(ecosw->type != IDL_TYP_DOUBLE)
    ecosw = IDL_CvtDbl(1, &ecosw);

  phi = eb_phisec(esinw->value.d, ecosw->value.d);

  if(esinw != argv[0])
    IDL_DELTMP(esinw);

  if(ecosw != argv[1])
    IDL_DELTMP(ecosw);

  return(IDL_GettmpDouble(phi));
}

IDL_VPTR wrap_eb_phicont (int argc, IDL_VPTR *argv) {
  IDL_VPTR esinw, ecosw, cosi, rasum, phi;
  double *phi_p;
  IDL_MEMINT dim[1] = { 4 };

  /* Check all required arguments are present */
  if(argc != 4) {
    IDL_Message(IDL_M_NAMED_GENERIC, IDL_MSG_LONGJMP,
                "Usage: y = eb_phicont(esinw, ecosw, cosi, rsum/a)\n");
    return(IDL_GettmpInt(-1));
  }

  esinw = argv[0];
  
  IDL_ENSURE_SIMPLE(esinw);
  IDL_ENSURE_SCALAR(esinw);

  if(esinw->type != IDL_TYP_DOUBLE)
    esinw = IDL_CvtDbl(1, &esinw);

  ecosw = argv[1];
  
  IDL_ENSURE_SIMPLE(ecosw);
  IDL_ENSURE_SCALAR(ecosw);

  if(ecosw->type != IDL_TYP_DOUBLE)
    ecosw = IDL_CvtDbl(1, &ecosw);

  cosi = argv[2];
  
  IDL_ENSURE_SIMPLE(cosi);
  IDL_ENSURE_SCALAR(cosi);

  if(cosi->type != IDL_TYP_DOUBLE)
    cosi = IDL_CvtDbl(1, &cosi);

  rasum = argv[3];
  
  IDL_ENSURE_SIMPLE(rasum);
  IDL_ENSURE_SCALAR(rasum);

  if(rasum->type != IDL_TYP_DOUBLE)
    rasum = IDL_CvtDbl(1, &rasum);

  phi_p = (double *) IDL_MakeTempArray(IDL_TYP_DOUBLE,
                                       1,
                                       dim,
                                       IDL_ARR_INI_NOP,
                                       &phi);

  eb_phicont(esinw->value.d, ecosw->value.d,
             cosi->value.d, rasum->value.d,
             phi_p);

  if(esinw != argv[0])
    IDL_DELTMP(esinw);

  if(ecosw != argv[1])
    IDL_DELTMP(ecosw);

  if(cosi != argv[2])
    IDL_DELTMP(cosi);

  if(rasum != argv[3])
    IDL_DELTMP(rasum);

  return(phi);
}

IDL_VPTR wrap_eb_getvder (int argc, IDL_VPTR *argv) {
  IDL_VPTR parm, gamma, ktot;
  double *parm_p;

  IDL_VPTR vder;
  double *vder_p;

  IDL_MEMINT dim[1] = { EB_NDER };

  /* Check all required arguments are present */
  if(argc != 3) {
    IDL_Message(IDL_M_NAMED_GENERIC, IDL_MSG_LONGJMP,
                "Usage: y = eb_phicont(parm, gamma, ktot)\n");
    return(IDL_GettmpInt(-1));
  }

  /* Unpack parameter vector */
  parm = argv[0];
  IDL_ENSURE_SIMPLE(parm);
  IDL_ENSURE_ARRAY(parm);

  if(parm->value.arr->n_elts < EB_NPAR) {
    IDL_Message(IDL_M_NAMED_GENERIC, IDL_MSG_LONGJMP,
                "Parameter vector is too short\n");
    goto error;
  }

  if(parm->type != IDL_TYP_DOUBLE)
    parm = IDL_CvtDbl(1, &parm);

  parm_p = (double *) parm->value.arr->data;

  /* Unpack scalar arguments */
  gamma = argv[1];
  
  IDL_ENSURE_SIMPLE(gamma);
  IDL_ENSURE_SCALAR(gamma);

  if(gamma->type != IDL_TYP_DOUBLE)
    gamma = IDL_CvtDbl(1, &gamma);

  ktot = argv[2];
  
  IDL_ENSURE_SIMPLE(ktot);
  IDL_ENSURE_SCALAR(ktot);

  if(ktot->type != IDL_TYP_DOUBLE)
    ktot = IDL_CvtDbl(1, &ktot);

  /* Output array */
  vder_p = (double *) IDL_MakeTempArray(IDL_TYP_DOUBLE,
                                        1,
                                        dim,
                                        IDL_ARR_INI_NOP,
                                        &vder);

  eb_getvder(parm_p,
             gamma->value.d, ktot->value.d,
             vder_p);

  if(parm != argv[0])
    IDL_DELTMP(parm);

  if(gamma != argv[1])
    IDL_DELTMP(gamma);

  if(ktot != argv[2])
    IDL_DELTMP(ktot);

  return(vder);

 error:
  return(IDL_GettmpInt(-1));
}

static IDL_SYSFUN_DEF2 eb_func_def[] = {
  { (IDL_SYSRTN_GENERIC) wrap_eb_model,
    "EB_MODEL",
    0,
    IDL_MAXPARAMS,
    IDL_SYSFUN_DEF_F_KEYWORDS,
    0 },
  { (IDL_SYSRTN_GENERIC) wrap_eb_phisec,
    "EB_PHISEC",
    0,
    IDL_MAXPARAMS,
    0,
    0 },
  { (IDL_SYSRTN_GENERIC) wrap_eb_phicont,
    "EB_PHICONT",
    0,
    IDL_MAXPARAMS,
    0,
    0 },
  { (IDL_SYSRTN_GENERIC) wrap_eb_getvder,
    "EB_GETVDER",
    0,
    IDL_MAXPARAMS,
    0,
    0 }
};

int IDL_Load (void) {
  return(IDL_SysRtnAdd(eb_func_def,
                       IDL_TRUE,  /* function, not procedure */
                       sizeof(eb_func_def)/sizeof(eb_func_def[0])));
}

