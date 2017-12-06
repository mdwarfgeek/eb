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

IDL_VPTR wrap_eb_phiperi (int argc, IDL_VPTR *argv) {
  IDL_VPTR esinw, ecosw;
  double phi;

  /* Check all required arguments are present */
  if(argc != 2) {
    IDL_Message(IDL_M_NAMED_GENERIC, IDL_MSG_LONGJMP,
                "Usage: y = eb_phiperi(esinw, ecosw)\n");
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

  phi = eb_phiperi(esinw->value.d, ecosw->value.d);

  if(esinw != argv[0])
    IDL_DELTMP(esinw);

  if(ecosw != argv[1])
    IDL_DELTMP(ecosw);

  return(IDL_GettmpDouble(phi));
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

static IDL_STRUCT_TAG_DEF eb_const_tags[] = {
  { "PAR_J",      NULL, (void *) IDL_TYP_LONG,   0 },
  { "PAR_RASUM",  NULL, (void *) IDL_TYP_LONG,   0 },
  { "PAR_RR",     NULL, (void *) IDL_TYP_LONG,   0 },
  { "PAR_LDLIN1", NULL, (void *) IDL_TYP_LONG,   0 },
  { "PAR_LDLIN2", NULL, (void *) IDL_TYP_LONG,   0 },
  { "PAR_COSI",   NULL, (void *) IDL_TYP_LONG,   0 },
  { "PAR_ECOSW",  NULL, (void *) IDL_TYP_LONG,   0 },
  { "PAR_ESINW",  NULL, (void *) IDL_TYP_LONG,   0 },
  { "PAR_GD1",    NULL, (void *) IDL_TYP_LONG,   0 },
  { "PAR_GD2",    NULL, (void *) IDL_TYP_LONG,   0 },
  { "PAR_REFL1",  NULL, (void *) IDL_TYP_LONG,   0 },
  { "PAR_REFL2",  NULL, (void *) IDL_TYP_LONG,   0 },
  { "PAR_Q",      NULL, (void *) IDL_TYP_LONG,   0 },
  { "PAR_TIDANG", NULL, (void *) IDL_TYP_LONG,   0 },
  { "PAR_L3",     NULL, (void *) IDL_TYP_LONG,   0 },
  { "PAR_PHI0",   NULL, (void *) IDL_TYP_LONG,   0 },
  { "PAR_M0",     NULL, (void *) IDL_TYP_LONG,   0 },
  { "PAR_INTEG",  NULL, (void *) IDL_TYP_LONG,   0 },
  { "PAR_P",      NULL, (void *) IDL_TYP_LONG,   0 },
  { "PAR_T0",     NULL, (void *) IDL_TYP_LONG,   0 },
  { "PAR_LDNON1", NULL, (void *) IDL_TYP_LONG,   0 },
  { "PAR_LDNON2", NULL, (void *) IDL_TYP_LONG,   0 },
  { "PAR_CLTT",   NULL, (void *) IDL_TYP_LONG,   0 },
  { "PAR_ROT1",   NULL, (void *) IDL_TYP_LONG,   0 },
  { "PAR_ROT2",   NULL, (void *) IDL_TYP_LONG,   0 },
  { "PAR_FSPOT1", NULL, (void *) IDL_TYP_LONG,   0 },
  { "PAR_FSPOT2", NULL, (void *) IDL_TYP_LONG,   0 },
  { "PAR_OOE1O",  NULL, (void *) IDL_TYP_LONG,   0 },
  { "PAR_OOE11A", NULL, (void *) IDL_TYP_LONG,   0 },
  { "PAR_OOE11B", NULL, (void *) IDL_TYP_LONG,   0 },
  { "PAR_OOE12A", NULL, (void *) IDL_TYP_LONG,   0 },
  { "PAR_OOE12B", NULL, (void *) IDL_TYP_LONG,   0 },
  { "PAR_OOE2O",  NULL, (void *) IDL_TYP_LONG,   0 },
  { "PAR_OOE21A", NULL, (void *) IDL_TYP_LONG,   0 },
  { "PAR_OOE21B", NULL, (void *) IDL_TYP_LONG,   0 },
  { "PAR_OOE22A", NULL, (void *) IDL_TYP_LONG,   0 },
  { "PAR_OOE22B", NULL, (void *) IDL_TYP_LONG,   0 },
  { "PAR_DWDT",   NULL, (void *) IDL_TYP_LONG,   0 },
  { "NPAR",       NULL, (void *) IDL_TYP_LONG,   0 },

  { "FLAG_REFL",  NULL, (void *) IDL_TYP_LONG,   0 },
  { "FLAG_PHI",   NULL, (void *) IDL_TYP_LONG,   0 },

  { "OBS_MAG",    NULL, (void *) IDL_TYP_BYTE,   0 },
  { "OBS_LIGHT",  NULL, (void *) IDL_TYP_BYTE,   0 },
  { "OBS_LRAT",   NULL, (void *) IDL_TYP_BYTE,   0 },
  { "OBS_AVLR",   NULL, (void *) IDL_TYP_BYTE,   0 },
  { "OBS_VRAD1",  NULL, (void *) IDL_TYP_BYTE,   0 },
  { "OBS_VRAD2",  NULL, (void *) IDL_TYP_BYTE,   0 },
  { "OBS_PSS",    NULL, (void *) IDL_TYP_BYTE,   0 },
  { "OBS_A",      NULL, (void *) IDL_TYP_BYTE,   0 },
  { "OBS_LSS",    NULL, (void *) IDL_TYP_BYTE,   0 },

  { "PAR_I",      NULL, (void *) IDL_TYP_LONG,   0 },
  { "PAR_R1A",    NULL, (void *) IDL_TYP_LONG,   0 },
  { "PAR_R2A",    NULL, (void *) IDL_TYP_LONG,   0 },
  { "PAR_E",      NULL, (void *) IDL_TYP_LONG,   0 },
  { "PAR_OMEGA",  NULL, (void *) IDL_TYP_LONG,   0 },
  { "PAR_A",      NULL, (void *) IDL_TYP_LONG,   0 },
  { "PAR_MTOT",   NULL, (void *) IDL_TYP_LONG,   0 },
  { "PAR_M1",     NULL, (void *) IDL_TYP_LONG,   0 },
  { "PAR_M2",     NULL, (void *) IDL_TYP_LONG,   0 },
  { "PAR_RTOT",   NULL, (void *) IDL_TYP_LONG,   0 },
  { "PAR_R1",     NULL, (void *) IDL_TYP_LONG,   0 },
  { "PAR_R2",     NULL, (void *) IDL_TYP_LONG,   0 },
  { "PAR_LOGG1",  NULL, (void *) IDL_TYP_LONG,   0 },
  { "PAR_LOGG2",  NULL, (void *) IDL_TYP_LONG,   0 },
  { "PAR_VSYNC1", NULL, (void *) IDL_TYP_LONG,   0 },
  { "PAR_VSYNC2", NULL, (void *) IDL_TYP_LONG,   0 },
  { "PAR_TSYNC",  NULL, (void *) IDL_TYP_LONG,   0 },
  { "PAR_TCIRC",  NULL, (void *) IDL_TYP_LONG,   0 },
  { "PAR_TSEC",   NULL, (void *) IDL_TYP_LONG,   0 },
  { "PAR_DURPRI", NULL, (void *) IDL_TYP_LONG,   0 },
  { "PAR_DURSEC", NULL, (void *) IDL_TYP_LONG,   0 },
  { "NDER",       NULL, (void *) IDL_TYP_LONG,   0 },

  { "GMSUN",      NULL, (void *) IDL_TYP_DOUBLE, 0 },
  { "AU",         NULL, (void *) IDL_TYP_DOUBLE, 0 },
  { "LIGHT",      NULL, (void *) IDL_TYP_DOUBLE, 0 },
  { "RSUN",       NULL, (void *) IDL_TYP_DOUBLE, 0 },

  { NULL, NULL, NULL, 0 }
};

static struct {
  IDL_LONG par_j;
  IDL_LONG par_rasum;
  IDL_LONG par_rr;
  IDL_LONG par_ldlin1;
  IDL_LONG par_ldlin2;
  IDL_LONG par_cosi;
  IDL_LONG par_ecosw;
  IDL_LONG par_esinw;
  IDL_LONG par_gd1;
  IDL_LONG par_gd2;
  IDL_LONG par_refl1;
  IDL_LONG par_refl2;
  IDL_LONG par_q;
  IDL_LONG par_tidang;
  IDL_LONG par_l3;
  IDL_LONG par_phi0;
  IDL_LONG par_m0;
  IDL_LONG par_integ;
  IDL_LONG par_p;
  IDL_LONG par_t0;
  IDL_LONG par_ldnon1;
  IDL_LONG par_ldnon2;
  IDL_LONG par_cltt;
  IDL_LONG par_rot1;
  IDL_LONG par_rot2;
  IDL_LONG par_fspot1;
  IDL_LONG par_fspot2;
  IDL_LONG par_ooe1o;
  IDL_LONG par_ooe11a;
  IDL_LONG par_ooe11b;
  IDL_LONG par_ooe12a;
  IDL_LONG par_ooe12b;
  IDL_LONG par_ooe2o;
  IDL_LONG par_ooe21a;
  IDL_LONG par_ooe21b;
  IDL_LONG par_ooe22a;
  IDL_LONG par_ooe22b;
  IDL_LONG par_dwdt;
  IDL_LONG npar;

  IDL_LONG flag_refl;
  IDL_LONG flag_phi;

  unsigned char obs_mag;
  unsigned char obs_light;
  unsigned char obs_lrat;
  unsigned char obs_avlr;
  unsigned char obs_vrad1;
  unsigned char obs_vrad2;
  unsigned char obs_pss;
  unsigned char obs_a;
  unsigned char obs_lss;

  IDL_LONG par_i;
  IDL_LONG par_r1a;
  IDL_LONG par_r2a;
  IDL_LONG par_e;
  IDL_LONG par_omega;
  IDL_LONG par_a;
  IDL_LONG par_mtot;
  IDL_LONG par_m1;
  IDL_LONG par_m2;
  IDL_LONG par_rtot;
  IDL_LONG par_r1;
  IDL_LONG par_r2;
  IDL_LONG par_logg1;
  IDL_LONG par_logg2;
  IDL_LONG par_vsync1;
  IDL_LONG par_vsync2;
  IDL_LONG par_tsync;
  IDL_LONG par_tcirc;
  IDL_LONG par_tsec;
  IDL_LONG par_durpri;
  IDL_LONG par_dursec;
  IDL_LONG nder;

  double gmsun;
  double au;
  double light;
  double rsun;

} eb_const_data = {
  EB_PAR_J,
  EB_PAR_RASUM,
  EB_PAR_RR,
  EB_PAR_LDLIN1,
  EB_PAR_LDLIN2,
  EB_PAR_COSI,
  EB_PAR_ECOSW,
  EB_PAR_ESINW,
  EB_PAR_GD1,
  EB_PAR_GD2,
  EB_PAR_REFL1,
  EB_PAR_REFL2,
  EB_PAR_Q,
  EB_PAR_TIDANG,
  EB_PAR_L3,
  EB_PAR_PHI0,
  EB_PAR_M0,
  EB_PAR_INTEG,
  EB_PAR_P,
  EB_PAR_T0,
  EB_PAR_LDNON1,
  EB_PAR_LDNON2,
  EB_PAR_CLTT,
  EB_PAR_ROT1,
  EB_PAR_ROT2,
  EB_PAR_FSPOT1,
  EB_PAR_FSPOT2,
  EB_PAR_OOE1O,
  EB_PAR_OOE11A,
  EB_PAR_OOE11B,
  EB_PAR_OOE12A,
  EB_PAR_OOE12B,
  EB_PAR_OOE2O,
  EB_PAR_OOE21A,
  EB_PAR_OOE21B,
  EB_PAR_OOE22A,
  EB_PAR_OOE22B,
  EB_PAR_DWDT,
  EB_NPAR,

  EB_FLAG_REFL,
  EB_FLAG_PHI,

  EB_OBS_MAG,
  EB_OBS_LIGHT,
  EB_OBS_LRAT,
  EB_OBS_AVLR,
  EB_OBS_VRAD1,
  EB_OBS_VRAD2,
  EB_OBS_PSS,
  EB_OBS_A,
  EB_OBS_LSS,

  EB_PAR_I,
  EB_PAR_R1A,
  EB_PAR_R2A,
  EB_PAR_E,
  EB_PAR_OMEGA,
  EB_PAR_A,
  EB_PAR_MTOT,
  EB_PAR_M1,
  EB_PAR_M2,
  EB_PAR_RTOT,
  EB_PAR_R1,
  EB_PAR_R2,
  EB_PAR_LOGG1,
  EB_PAR_LOGG2,
  EB_PAR_VSYNC1,
  EB_PAR_VSYNC2,
  EB_PAR_TSYNC,
  EB_PAR_TCIRC,
  EB_PAR_TSEC,
  EB_PAR_DURPRI,
  EB_PAR_DURSEC,
  EB_NDER,

  EB_GMSUN,
  EB_AU,
  EB_LIGHT,
  EB_RSUN
};

IDL_VPTR wrap_eb_const (int argc, IDL_VPTR *argv) {
  void *s;
  IDL_MEMINT dim[1] = { 1 };
  IDL_VPTR out;

  s = IDL_MakeStruct(NULL, eb_const_tags);

  out = IDL_ImportArray(1, dim, IDL_TYP_STRUCT,
                        (UCHAR *) &eb_const_data, NULL, s);

  return(out);
}

IDL_VPTR wrap_eb_parnames (int argc, IDL_VPTR *argv) {
  IDL_MEMINT dim[1] = { EB_NPAR };
  IDL_STRING *ss;
  IDL_VPTR out;
  int ipar;

  /* Make output array (same shape as t) */
  ss = (IDL_STRING *) IDL_MakeTempArray(IDL_TYP_STRING,
                                        1,
                                        dim,
                                        IDL_ARR_INI_NOP,
                                        &out);

  for(ipar = 0; ipar < EB_NPAR; ipar++) {
    ss[ipar].slen  = strlen(eb_parnames[ipar]);
    ss[ipar].stype = 0;
    ss[ipar].s     = eb_parnames[ipar];
  }

  return(out);
}

IDL_VPTR wrap_eb_parunits (int argc, IDL_VPTR *argv) {
  IDL_MEMINT dim[1] = { EB_NPAR };
  IDL_STRING *ss;
  IDL_VPTR out;
  int ipar;

  /* Make output array (same shape as t) */
  ss = (IDL_STRING *) IDL_MakeTempArray(IDL_TYP_STRING,
                                        1,
                                        dim,
                                        IDL_ARR_INI_NOP,
                                        &out);

  for(ipar = 0; ipar < EB_NPAR; ipar++) {
    ss[ipar].slen  = strlen(eb_parunits[ipar]);
    ss[ipar].stype = 0;
    ss[ipar].s     = eb_parunits[ipar];
  }

  return(out);
}

IDL_VPTR wrap_eb_dernames (int argc, IDL_VPTR *argv) {
  IDL_MEMINT dim[1] = { EB_NDER };
  IDL_STRING *ss;
  IDL_VPTR out;
  int ider;

  /* Make output array (same shape as t) */
  ss = (IDL_STRING *) IDL_MakeTempArray(IDL_TYP_STRING,
                                        1,
                                        dim,
                                        IDL_ARR_INI_NOP,
                                        &out);

  for(ider = 0; ider < EB_NDER; ider++) {
    ss[ider].slen  = strlen(eb_dernames[ider]);
    ss[ider].stype = 0;
    ss[ider].s     = eb_dernames[ider];
  }

  return(out);
}

IDL_VPTR wrap_eb_derunits (int argc, IDL_VPTR *argv) {
  IDL_MEMINT dim[1] = { EB_NDER };
  IDL_STRING *ss;
  IDL_VPTR out;
  int ider;

  /* Make output array (same shape as t) */
  ss = (IDL_STRING *) IDL_MakeTempArray(IDL_TYP_STRING,
                                        1,
                                        dim,
                                        IDL_ARR_INI_NOP,
                                        &out);

  for(ider = 0; ider < EB_NDER; ider++) {
    ss[ider].slen  = strlen(eb_derunits[ider]);
    ss[ider].stype = 0;
    ss[ider].s     = eb_derunits[ider];
  }

  return(out);
}

static IDL_SYSFUN_DEF2 eb_func_def[] = {
  { (IDL_SYSRTN_GENERIC) wrap_eb_model,
    "EB_MODEL",
    0,
    IDL_MAXPARAMS,
    IDL_SYSFUN_DEF_F_KEYWORDS,
    0 },
  { (IDL_SYSRTN_GENERIC) wrap_eb_phiperi,
    "EB_PHIPERI",
    0,
    IDL_MAXPARAMS,
    0,
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
    0 },
  { (IDL_SYSRTN_GENERIC) wrap_eb_const,
    "EB_CONST",
    0,
    IDL_MAXPARAMS,
    0,
    0 },
  { (IDL_SYSRTN_GENERIC) wrap_eb_parnames,
    "EB_PARNAMES",
    0,
    IDL_MAXPARAMS,
    0,
    0 },
  { (IDL_SYSRTN_GENERIC) wrap_eb_parunits,
    "EB_PARUNITS",
    0,
    IDL_MAXPARAMS,
    0,
    0 },
  { (IDL_SYSRTN_GENERIC) wrap_eb_dernames,
    "EB_DERNAMES",
    0,
    IDL_MAXPARAMS,
    0,
    0 },
  { (IDL_SYSRTN_GENERIC) wrap_eb_derunits,
    "EB_DERUNITS",
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

