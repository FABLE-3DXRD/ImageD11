

static char moduledocs [] =\
"/* These were Copyrighted by Enthought Inc */\n"\
" It is an interface to bispev from fitpack (Dierckx/netlib) \n"\
" see bisplev.py \n"\
"/* They  seemed  to  give permission  for  me to copy out the \n"\
"   parts I wanted. I promise that if it goes wrong it was not\n"\
"   their fault.  They are not responsible for\n"\
"   errors I might have introduced                  - JPW 2005  */"\
" Adpated to numpy (change of include file), JPW 2007 "\
" Change to use PyArray_SimpleNew etc depreciation warn, JPW 2009 ";

#include <Python.h>
/* upgrade to numpy
   #include "Numeric/arrayobject.h"  */
#include "numpy/arrayobject.h"  

#define BISPEV bispev_
#define PARDER parder_
#define SURFIT surfit_

void BISPEV(double*,int*,double*,int*,double*,int*,int*,double*,int*,
                      double*,int*,double*,double*,int*,int*,int*,int*);
void PARDER(double*,int*,double*,int*,double*,int*,int*,int*,int*,double*,
                  int*,double*,int*,double*,double*,int*,int*,int*,int*);
void SURFIT(int*,int*,double*,double*,double*,double*,double*,double*,double*,double*,
            int*,int*,double*,int*,int*,int*,double*,int*,double*,int*,double*,double*,
            double*,double*,int*,double*,int*,int*,int*,int*);

static char doc_bispev[] = " [z,ier] = _bispev(tx,ty,c,kx,ky,x,y,nux,nuy)";

static PyObject *_bispev(PyObject *dummy, PyObject *args) {
  int nx,ny,kx,ky,mxy,mx,my,lwrk,*iwrk,kwrk,ier,lwa,nux,nuy;
  npy_intp dmxy;
  double *tx,*ty,*c,*x,*y,*z,*wrk,*wa = NULL;
  PyArrayObject *ap_x = NULL,*ap_y = NULL,*ap_z = NULL,*ap_tx = NULL,\
    *ap_ty = NULL,*ap_c = NULL;
  PyObject *x_py = NULL,*y_py = NULL,*c_py = NULL,*tx_py = NULL,*ty_py = NULL;
  if (!PyArg_ParseTuple(args, "OOOiiOOii",&tx_py,&ty_py,&c_py,&kx,&ky,
                        &x_py,&y_py,&nux,&nuy))
    return NULL;
  ap_x = (PyArrayObject *)PyArray_ContiguousFromObject(x_py, PyArray_DOUBLE, 0, 1);
  ap_y = (PyArrayObject *)PyArray_ContiguousFromObject(y_py, PyArray_DOUBLE, 0, 1);
  ap_c = (PyArrayObject *)PyArray_ContiguousFromObject(c_py, PyArray_DOUBLE, 0, 1);
  ap_tx = (PyArrayObject *)PyArray_ContiguousFromObject(tx_py, PyArray_DOUBLE, 0, 1);
  ap_ty = (PyArrayObject *)PyArray_ContiguousFromObject(ty_py, PyArray_DOUBLE, 0, 1);
  if (ap_x == NULL || ap_y == NULL || ap_c == NULL || ap_tx == NULL \
      || ap_ty == NULL) goto fail;  
  x = (double *) ap_x->data;
  y = (double *) ap_y->data;
  c = (double *) ap_c->data;
  tx = (double *) ap_tx->data;
  ty = (double *) ap_ty->data;
  nx = ap_tx->dimensions[0];
  ny = ap_ty->dimensions[0];
  mx = ap_x->dimensions[0];
  my = ap_y->dimensions[0];
  mxy = mx*my;
  dmxy = mxy*1;
  ap_z = (PyArrayObject *)PyArray_SimpleNew(1,&dmxy,PyArray_DOUBLE);
  z = (double *) ap_z->data;
  if (nux || nuy) 
    lwrk = mx*(kx+1-nux)+my*(ky+1-nuy)+(nx-kx-1)*(ny-ky-1);
  else
    lwrk = mx*(kx+1)+my*(ky+1);
  kwrk = mx+my;
  lwa = lwrk+kwrk;
  if ((wa = (double *)malloc(lwa*sizeof(double)))==NULL) {
    PyErr_NoMemory();
    goto fail;
  }
  wrk = wa;
  iwrk = (int *)(wrk+lwrk);
  if (nux || nuy)
    PARDER(tx,&nx,ty,&ny,c,&kx,&ky,&nux,&nuy,x,&mx,y,&my,z,wrk,&lwrk,iwrk,&kwrk,&ier);
  else
    BISPEV(tx,&nx,ty,&ny,c,&kx,&ky,x,&mx,y,&my,z,wrk,&lwrk,iwrk,&kwrk,&ier);
  if (wa) free(wa);
  Py_DECREF(ap_x);
  Py_DECREF(ap_y);
  Py_DECREF(ap_c);
  Py_DECREF(ap_tx);
  Py_DECREF(ap_ty);
  return Py_BuildValue("Ni",PyArray_Return(ap_z),ier);
  fail:
  if (wa) free(wa);
  Py_XDECREF(ap_x);
  Py_XDECREF(ap_y);
  Py_XDECREF(ap_z);
  Py_XDECREF(ap_c);
  Py_XDECREF(ap_tx);
  Py_XDECREF(ap_ty);
  return NULL;
}

static char doc_surfit[] =
 " [tx,ty,c,o] = "
 "_surfit(x,y,z,w,xb,xe,yb,ye,kx,ky,iopt,s,eps,tx,ty,nxest,nyest,wrk,lwrk1,lwrk2)";

static PyObject *_surfit(PyObject *dummy, PyObject *args) {
  int iopt,m,kx,ky,nxest,nyest,lwrk1,lwrk2,*iwrk,kwrk,ier,lwa,nxo,nyo,\
    i,lcest,nmax,nx,ny,lc,n;
  npy_intp dnx,dny,dlc,dn;
  double *x,*y,*z,*w,xb,xe,yb,ye,s,*tx,*ty,*c,fp,*wrk1,*wrk2,*wa = NULL,eps;
  PyArrayObject *ap_x = NULL,*ap_y = NULL,*ap_z,*ap_w = NULL,\
    *ap_tx = NULL,*ap_ty = NULL,*ap_c = NULL;
  PyArrayObject *ap_wrk = NULL;
  PyObject *x_py = NULL,*y_py = NULL,*z_py = NULL,*w_py = NULL,\
    *tx_py = NULL,*ty_py = NULL;
  PyObject *wrk_py=NULL;
  nx=ny=ier=nxo=nyo=0;
  if (!PyArg_ParseTuple(args, "OOOOddddiiiddOOiiOii",\
                        &x_py,&y_py,&z_py,&w_py,&xb,&xe,\
                        &yb,&ye,&kx,&ky,&iopt,&s,&eps,&tx_py,&ty_py,&nxest,&nyest,\
                        &wrk_py,&lwrk1,&lwrk2)) return NULL;
  ap_x = (PyArrayObject *)PyArray_ContiguousFromObject(x_py, PyArray_DOUBLE, 0, 1);
  ap_y = (PyArrayObject *)PyArray_ContiguousFromObject(y_py, PyArray_DOUBLE, 0, 1);
  ap_z = (PyArrayObject *)PyArray_ContiguousFromObject(z_py, PyArray_DOUBLE, 0, 1);
  ap_w = (PyArrayObject *)PyArray_ContiguousFromObject(w_py, PyArray_DOUBLE, 0, 1);
  ap_wrk=(PyArrayObject *)PyArray_ContiguousFromObject(wrk_py, PyArray_DOUBLE, 0, 1);
  /*ap_iwrk=(PyArrayObject *)PyArray_ContiguousFromObject(iwrk_py, PyArray_INT, 0, 1);*/
  if (ap_x == NULL || ap_y == NULL || ap_z == NULL || ap_w == NULL \
      || ap_wrk == NULL) goto fail;
  x = (double *) ap_x->data;
  y = (double *) ap_y->data;
  z = (double *) ap_z->data;
  w = (double *) ap_w->data;
  m = ap_x->dimensions[0];
  nmax=nxest;
  if (nmax<nyest) nmax=nyest;
  lcest=(nxest-kx-1)*(nyest-ky-1);
  kwrk=m+(nxest-2*kx-1)*(nyest-2*ky-1);
  lwa = 2*nmax+lcest+lwrk1+lwrk2+kwrk;
  if ((wa = (double *)malloc(lwa*sizeof(double)))==NULL) {
    PyErr_NoMemory();
    goto fail;
  }
  tx = wa;
  ty = tx + nmax;
  c = ty + nmax;
  wrk1 = c + lcest;
  iwrk = (int *)(wrk1 + lwrk1);
  wrk2 = (double *)(iwrk+kwrk);
  if (iopt) {
    ap_tx=(PyArrayObject *)PyArray_ContiguousFromObject(tx_py, PyArray_DOUBLE, 0, 1);
    ap_ty=(PyArrayObject *)PyArray_ContiguousFromObject(ty_py, PyArray_DOUBLE, 0, 1);
    if (ap_tx == NULL || ap_ty == NULL) goto fail;
    nx = nxo = ap_tx->dimensions[0];
    ny = nyo = ap_ty->dimensions[0];
    memcpy(tx,ap_tx->data,nx*sizeof(double));
    memcpy(ty,ap_ty->data,ny*sizeof(double));
  }
  if (iopt==1) {
    lc = (nx-kx-1)*(ny-ky-1);
    memcpy(wrk1,ap_wrk->data,lc*sizeof(double));
    /*memcpy(iwrk,ap_iwrk->data,n*sizeof(int));*/
  }
  SURFIT(&iopt,&m,x,y,z,w,&xb,&xe,&yb,&ye,&kx,&ky,&s,&nxest,&nyest,&nmax,
        &eps,&nx,tx,&ny,ty,c,&fp,wrk1,&lwrk1,wrk2,&lwrk2,iwrk,&kwrk,&ier);
  i=0;
  while ((ier>10) && (i++<5)) {
    lwrk2=ier;
    if ((wrk2 = (double *)malloc(lwrk2*sizeof(double)))==NULL) {
      PyErr_NoMemory();
      goto fail;
    }
    SURFIT(&iopt,&m,x,y,z,w,&xb,&xe,&yb,&ye,&kx,&ky,&s,&nxest,&nyest,&nmax,
          &eps,&nx,tx,&ny,ty,c,&fp,wrk1,&lwrk1,wrk2,&lwrk2,iwrk,&kwrk,&ier);
    if (wrk2) free(wrk2);
  }
  if (ier==10) goto fail;
  lc = (nx-kx-1)*(ny-ky-1);
  dnx = 1*nx;
  dny = 1*ny;
  dlc = 1*lc;
  ap_tx = (PyArrayObject *)PyArray_SimpleNew(1,&dnx,PyArray_DOUBLE);
  ap_ty = (PyArrayObject *)PyArray_SimpleNew(1,&dny,PyArray_DOUBLE);
  ap_c = (PyArrayObject *)PyArray_SimpleNew(1,&dlc,PyArray_DOUBLE);
  if (ap_tx == NULL || ap_ty == NULL || ap_c == NULL) goto fail;
  if ((iopt==0)||(nx>nxo)||(ny>nyo)) {
    ap_wrk = (PyArrayObject *)PyArray_SimpleNew(1,&dlc,PyArray_DOUBLE);
    if (ap_wrk == NULL) goto fail;
    dn = 1*n;
    /*ap_iwrk = (PyArrayObject *)PyArray_SimpleNew(1,&dn,PyArray_INT);*/
  }
  memcpy(ap_tx->data,tx,nx*sizeof(double));
  memcpy(ap_ty->data,ty,ny*sizeof(double));
  memcpy(ap_c->data,c,lc*sizeof(double));
  memcpy(ap_wrk->data,wrk1,lc*sizeof(double));
  /*memcpy(ap_iwrk->data,iwrk,n*sizeof(int));*/
  if (wa) free(wa);
  Py_DECREF(ap_x);
  Py_DECREF(ap_y);
  Py_DECREF(ap_z);
  Py_DECREF(ap_w);
  return Py_BuildValue("NNN{s:N,s:i,s:d}",PyArray_Return(ap_tx),\
             PyArray_Return(ap_ty),PyArray_Return(ap_c),\
             "wrk",PyArray_Return(ap_wrk),\
             "ier",ier,"fp",fp);
  fail:
  if (wa) free(wa);
  Py_XDECREF(ap_x);
  Py_XDECREF(ap_y);
  Py_XDECREF(ap_z);
  Py_XDECREF(ap_w);
  Py_XDECREF(ap_tx);
  Py_XDECREF(ap_ty);
  Py_XDECREF(ap_wrk);
  /*Py_XDECREF(ap_iwrk);*/
  return NULL;
}



static struct PyMethodDef splines_methods[] = {
   {"_bispev", _bispev, METH_VARARGS, doc_bispev},
   {"_surfit", _surfit, METH_VARARGS, doc_surfit},

   {NULL,          NULL, 0, NULL}
   };

void init_splines(void) {
  PyObject *m, *d, *s, *ms;
  m = Py_InitModule("_splines", splines_methods);
  import_array();
  d = PyModule_GetDict(m);
  s = PyString_FromString(" 0.1 ");
  PyDict_SetItemString(d, "__version__", s);
  ms=PyString_FromString(moduledocs);
  PyDict_SetItemString(d,"__doc__",ms);

  Py_DECREF(s);
  Py_DECREF(ms);
  if (PyErr_Occurred())
    Py_FatalError("can't initialize module splines");
}
