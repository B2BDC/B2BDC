<<<<<<< Updated upstream
#include "mex.h"
#include <math.h>

void LOCALengineHR(int n, int L,
  double *M[], int *Aidx[], int Adim[],
  double *alpha, double *beta,
  double *x, double *v,
  double *xproj, double *vproj, double *xwork);

void LOCALengineGD(int n, int L,
  double *M[], int *Aidx[], int Adim[],
  double *alpha, double *beta,
  double *x, double *v,
  double *xproj, double *vproj, double *xwork);

void LOCALquadG2quadS(double *pP, double *ABC, double *x0,
    double *v, int n);

void box2initint(double *alpha, double *beta,
   double *x, double *v, int n,
   double *pI1, double *pI2, double BN,
   double *tBox);

int LOCALquadS2int(double *ABC, double *I, double *tPos);

void LOCALquadG2quadS(double *pP, double *ABC, double *x0,
    double *v, int n);

void  box2initint(double *alpha, double *beta,
   double *x, double *v, int n,
   double *pI1, double *pI2, double BN,
   double *tBox);

void mexFunction(int nLHS,  mxArray *pLHS[],
int nRHS,  const mxArray *pRHS[])
{
   mxArray *Mcell, *ActiveIndex, *H, *X, *V, *Xout;
   int L, **Aidx, *Adim, i,j;
   mwSize ndims, *dims, n;
   mxArray *Mi, *Ai;
   double **Mdata, *pAi, *pX, *pV;
   double *LowerBoundBox, *UpperBoundBox;
   double *Xwork, *Xproj, *Vproj;
   int *atmp;
   int usemethod;
   
   usemethod = 0;  /* DEFAULT: use H&R */
   
   /* Rename the inputs */
   Mcell = pRHS[0];
   ActiveIndex = pRHS[1];
   H = pRHS[2];
   X = pRHS[3];
   V = pRHS[4];
   if (nRHS==6) {
      usemethod = 1;  /* GD */
   }

   /* Extract parameter dimension from H (n-by-2) */
   n = mxGetM(H);
   
   /* Get pointers to X and V */
   pX = mxGetPr(X);
   pV = mxGetPr(V);
   
   /* Allocate space for Xwork, Xproj and Vproj */
   Xproj = (double *) mxCalloc(n,sizeof(double));
   Vproj = (double *) mxCalloc(n,sizeof(double));
   Xout = mxCreateDoubleMatrix(n,1,mxREAL);
   pLHS[0] = Xout;
   Xwork = mxGetPr(Xout);

   /* Get Lower/Upper bound Box constraints */
   LowerBoundBox = mxGetPr(H);
   UpperBoundBox = LowerBoundBox + n;
   
   /* Mcell and ActiveIndex are CELLS of the same dimension */
   ndims = mxGetNumberOfDimensions(Mcell);
   dims = mxCalloc(ndims,sizeof(mwSize));
   dims = mxGetDimensions(Mcell);
   L = 1;
   for (i=0;i<ndims;i++) {
      /* mexPrintf("DIMS %d\n",dims[i]); */
      L *= dims[i];
   }
   /* mexPrintf("BASIC %d\n",L); */
   Adim = (int *) mxCalloc(L,sizeof(int));
   Aidx = (int **) mxCalloc(L,sizeof(int *));
   Mdata = (double **) mxCalloc(L,sizeof(double *));
   for (i=0;i<L;i++) {
      Mi = mxGetCell(Mcell,i);
      Ai = mxGetCell(ActiveIndex,i);  /* row vector ???*/
      Adim[i] = mxGetM(Ai)*mxGetN(Ai); /* 1-by?  or ?-by-1 */
      Mdata[i] = mxGetPr(Mi);
      Aidx[i] = (int *) mxCalloc(Adim[i],sizeof(int));
      atmp = Aidx[i];
      pAi = mxGetPr(Ai);
      for (j=0;j<Adim[i];j++) {
         *atmp++ = (int) *pAi++;
      }
   }
   
   if (usemethod==0) {
      /* allocate Xwork bigger, to hold many samples;
      /* put this is FOR-loop
       *Each time, reset pX (starting point) to the column of Xwork that
       *just got created
       *Create a new random direction, filling pV */ 
      LOCALengineHR(n,L,Mdata,Aidx,Adim,
        LowerBoundBox,UpperBoundBox,
        pX,pV,Xproj,Vproj,Xwork);
   } else {
      LOCALengineGD(n,L,Mdata,Aidx,Adim,
        LowerBoundBox,UpperBoundBox,
        pX,pV,Xproj,Vproj,Xwork);
   }
   
   mxFree(Xproj);
   mxFree(Vproj);
}
     

void LOCALengineHR(int n, int L,
  double *M[], int *Aidx[], int Adim[],
  double *alpha, double *beta,
  double *x, double *v,
  double *xproj, double *vproj, double *xwork)
{
   int i, nactive, k, nsingle, *pidx;
   double scalarquad[3], Isingle[4];
   double I1, I2;
   double *Icurrent, *Inew;
   int ncurrent, nnew;
   int Ichoose;
   double tfac, BN, gamma, *tBox;
   double dummy;
      
   BN = mxGetInf();
   /* Handle the box contraints to get initial interval */
   /* tBox is unused in Hit&Run, but is used in GasDynamics */
   tBox = (double *) calloc(n, sizeof(double));
   box2initint(alpha,beta,x,v,n,&I1,&I2,BN,tBox);
   free(tBox);

   Icurrent = (double *) calloc(2*(L+1), sizeof(double));
   Icurrent[0] = I1;  Icurrent[1] = I2;
   ncurrent = 1;
   Inew = (double *) calloc(2*(L+1), sizeof(double));

   for (i=0;i<L;i++) {
      /* i'th Model */
      nactive = Adim[i];
      pidx = Aidx[i];
      /* Project down to active variables */
      for (k=0;k<nactive;k++) {
         xproj[k] = x[*pidx-1];
         vproj[k] = v[*pidx-1];
         pidx++;
      }
      LOCALquadG2quadS(M[i], scalarquad, xproj, vproj, nactive);
      nsingle = LOCALquadS2int(scalarquad, Isingle, &dummy);
      nnew = LOCALIntersect(Icurrent, ncurrent, Isingle, nsingle, Inew);
      /* Replace Current with New */
      ncurrent = nnew;
      for (k=0;k<2*ncurrent;k++) {
         Icurrent[k] = Inew[k];
      }
   }

   /* Choose an interval */
   Ichoose = (int) floor(0.9999*((double) rand())/((double) RAND_MAX)*ncurrent);
   gamma = 0.0001 + 0.9998*((double) rand())/((double) RAND_MAX);
   tfac = gamma*Icurrent[2*Ichoose] + (1-gamma)*Icurrent[2*Ichoose+1];

   /* Update Xwork */
   for (i=0;i<n;i++) {
      xwork[i] = x[i] + tfac*v[i];
   }
   /* Free memory */
   free(Icurrent);
   free(Inew);
}


void LOCALengineGD(int n, int L,
  double *M[], int *Aidx[], int Adim[],
  double *alpha, double *beta,
  double *x, double *v,
  double *xproj, double *vproj, double *xwork)
{
   int i, nactive, k, kk, nsingle, *pidx;
   double scalarquad[3], Isingle[4], MBox[4];
   double *vlocal, *xlocal, *ONproj, *xcontact;
   double I1, I2, *Muse;
   double BN;
   double *tBox, *tGeneral, *tAll;
   int *tIndex, *tIndexSave;
   double iprod, normout,TRemain;
   double BACKOFF, tPos;
   int ncols, minidx, idx, icnt;
   
   BACKOFF = 0.999;
   ncols = 0;
   TRemain = 0.5;
   
   xlocal = (double *) calloc(n, sizeof(double));
   xcontact = (double *) calloc(n, sizeof(double));
   vlocal = (double *) calloc(n, sizeof(double));
   ONproj = (double *) calloc(n, sizeof(double));
   tAll = (double *) calloc(n+L, sizeof(double));
   tIndex = (int *) calloc(n+L, sizeof(int));
   tIndexSave = tIndex;
   
   tBox = tAll;
   tGeneral = tAll + n;
   for (i=0;i<n+L;i++) {
      tIndex[i] = i;  /* fill with 0,...,n+L-1 */
   }
   for (i=0;i<n;i++) {
      vlocal[i] = v[i];
      xlocal[i] = x[i];
   }
   BN = mxGetInf();
   icnt = 0;
   while (TRemain>0 & icnt<6) {
      box2initint(alpha,beta,xlocal,vlocal,n,&I1,&I2,BN,tBox);
      for (i=0;i<L;i++) {
         nactive = Adim[i];
         pidx = Aidx[i];
         for (k=0;k<nactive;k++) {
            xproj[k] = xlocal[*pidx-1];
            vproj[k] = vlocal[*pidx-1];
            pidx++;
         }
         LOCALquadG2quadS(M[i], scalarquad, xproj, vproj, nactive);
         nsingle = LOCALquadS2int(scalarquad, Isingle, tGeneral);
         tGeneral++;
      }
      tGeneral -= L;
      /* LOCALmergesort(tAll, tIndex); */
      tPos = tAll[0];
      minidx = 0;
      for (k=1;k<n+L;k++) {
         if (tAll[k]<tPos) {
            tPos = tAll[k];
            minidx = k;
         }
      }
      if (tPos>TRemain) {
         /* Stop before hitting surface */
         for (i=0;i<n;i++) {
            xlocal[i] += TRemain*vlocal[i];
         }
         TRemain = 0.0;
      } else if (tPos<TRemain) {
         /* Collide with wall, continue on */
         /* update wall collision stats */ 
         ncols++;
         
         /* Xcontact, and Xcontact{BACKOFF} */
         for (i=0;i<n;i++) {
            xcontact[i] = xlocal[i] + tPos*vlocal[i];
            xlocal[i] += BACKOFF*tPos*vlocal[i];
         }
         
         /* Outward Normal at xcontact */
         /* NOTE: all nonactive variables will have 0 gradient */
         if (minidx<n) {
            nactive = 1;
            /* x<=BETA <=> -BETA+x<0 */
            MBox[0] = -beta[minidx];
            MBox[1] = 0.5;
            MBox[2] = 0.5;
            MBox[3] = 0.0;
            Muse = MBox;
            minidx += 1;  /* to agree with MATLAB indexing */
            pidx = &minidx;
         } else {
            minidx = minidx - n;  /* offset */
            nactive = Adim[minidx];
            pidx = Aidx[minidx];
            Muse = M[minidx];
         }
         /* Project down to active variables */
         for (k=0;k<nactive;k++) {
            xproj[k] = xcontact[*pidx-1];
            vproj[k] = vlocal[*pidx-1];
            ONproj[k] = -Muse[k+1];  /* initialize */
            pidx++;
         }
         pidx -= nactive;
         normout = 0;
         for (k=0;k<nactive;k++) {
            for (kk=0;kk<nactive;kk++) {
               idx = kk+1 + (k+1)*(nactive+1);
               ONproj[k] -= xproj[kk]*Muse[idx];
            }
            normout += ONproj[k]*ONproj[k];
         }
         normout = sqrt(normout);
         for (k=0;k<nactive;k++) {
            ONproj[k] = ONproj[k]/normout;
         }
         iprod = 0.0;
         for (k=0;k<nactive;k++) {
            iprod += ONproj[k]*vproj[k];
         }
         /* Update active components of V for reflection */
         for (k=0;k<nactive;k++) {
            vlocal[*pidx-1] -= 2.0*iprod*ONproj[k];
            pidx++;
         }
         pidx -= nactive;
         /* At this point, XLOCAL and VLOCAL are ready to use again */
         /* Update the remaining time, and go back to beginning */
         TRemain -= tPos;
      } else {
         /* Apparently will hit surface at exact time */
         /* Backoff by 10%, and call it done */
         for (i=0;i<n;i++) {
            xlocal[i] += 0.9*TRemain*vlocal[i];
         }
         TRemain = 0.0;
      }
      icnt++;
   }
   /* Update x */
   for (i=0;i<n;i++) {
      xwork[i] = xlocal[i];
   }
   /* Free memory */
   free(xlocal);
   free(xcontact);
   free(vlocal);
   free(ONproj);
   free(tAll);
   free(tIndexSave);
}


int LOCALquadS2int(double *ABC, double *I, double *tPos)
{
   double A, B, C, r, r1, r2, radi, sradi;
   double AINF;
   int nintervals;
   
   AINF = 10000000000000;

   A = *ABC++;
   B = *ABC++;
   C = *ABC++;
   if (A==0.0) {
      if (B==0.0) {
         I[0] = -AINF;
         I[1] = AINF;
         nintervals = 1;
     } else {
        r = -C/B;
        if (r<0) {
           I[0] = r;
           I[1] = AINF;
           nintervals = 1;
        } else if (r>0) {
           I[0] = -AINF;
           I[1] = r;
           nintervals = 1;
        } 
     }
      *tPos = I[1];
  } else {
     radi = B*B - 4*A*C;
     if (radi<0) {
        I[0] = -AINF;
        I[1] = AINF;
        nintervals = 1;
        *tPos = I[1];
     } else if (radi==0.0) {
        r = -B/(2*A);
        if (r<0) {
           I[0] = r;
           I[1] = AINF;
           nintervals = 1;
        } else if (r>0) {
           I[0] = -AINF;
           I[1] = r;
           nintervals = 1;
        }
        *tPos = I[1];
     } else {
        sradi = sqrt(radi);
        /* mexPrintf("RAD quals %g, SRAD is %g\n",radi,sradi); */
        if (A>0) {
            r1 = (-B - sradi)/(2*A);
            r2 = (-B + sradi)/(2*A);
        } else {
            r2 = (-B - sradi)/(2*A);
            r1 = (-B + sradi)/(2*A);
        }
        if ((r1<0) & (r2>0)) {
         I[0] = r1;
         I[1] = r2;
         nintervals = 1;
         *tPos = I[1];
        } else {
         I[0] = -AINF;
         I[1] = r1;
         I[2] = r2;
         I[3] = AINF;
         nintervals = 2;
         if (r1>0) {
            *tPos = r1;
         } else if (r2>0) {
            *tPos = r2;
         } else {
            *tPos = AINF;
         }
         }
      }
   }
   return nintervals;
}



void LOCALquadG2quadS(double *pP, double *ABC, double *x0,
    double *v, int n)
{
   int i, j, nsq;
   double *g22x, *g22v, g12v;
   double A, B, C;
   
   g22x = (double *) calloc(n,sizeof(double));
   g22v = (double *) calloc(n,sizeof(double));
   nsq = (1+n)*(1+n);
   g12v = 0.0;
   for (i=0;i<n;i++) {
      g12v += pP[i+1]*v[i];
      *g22x = 0.0;
      *g22v = 0.0;
      for (j=0;j<n;j++) {
         *g22v += pP[i+1+(j+1)*(n+1)]*v[j]; 
         *g22x += pP[i+1+(j+1)*(n+1)]*x0[j];
      }
      g22v++;
      g22x++;
   }
   g22v -= n;
   g22x -= n;
   A = 0.0;
   B = g12v;
   C = *pP;
   for (i=0;i<n;i++) {
      A += v[i]*g22v[i];
      B += v[i]*g22x[i];
      C += x0[i]*(2.0*pP[i+1] + g22x[i]);
   }
   B = 2.0*B; 
   *ABC++ = A;
   *ABC++ = B;
   *ABC++ = C; 
   free(g22v);
   free(g22x);
}

int LOCALIntersect(double *Icurrent, int ncurrent, 
   double *Isingle, int nsingle, double *Inew)
{
   double Lpt, Rpt;
   double L1i, R1i, L1j, R1j;
   int i, j, jstart, found, go;
   int nnew;

   if (ncurrent==1 & nsingle==1) {
      Lpt = Icurrent[0];
      if (Isingle[0]>Icurrent[0]) {
         Lpt = Isingle[0];
      }
      Rpt = Icurrent[1];
      if (Isingle[1]<Icurrent[1]) {
         Rpt = Isingle[1];
      }
      if (Lpt<Rpt) {
         nnew = 1;
         Inew[0] = Lpt;
         Inew[1] = Rpt;
      } else {
         nnew = 0;
      }
   } else {
      nnew = 0;
      jstart = 0;
      i = 0;
      j = 0;
      while (i<ncurrent) {
         L1i = Icurrent[2*i];
         R1i = Icurrent[2*i+1];
         go = 1;
         found = 0;
         j = jstart;
         while (j<nsingle & go==1) {
            L1j = Isingle[2*j];
            R1j = Isingle[2*j+1];
            Lpt = L1i;
            if (L1j>L1i) {
               Lpt = L1j;
            }
            Rpt = R1i;
            if (R1j<R1i) {
               Rpt = R1j;
            }
            /* Lpt = max(L1i,L1j);
            Rpt = min(R1i,R1j); */
            if (Lpt<Rpt) {
               Inew[2*nnew] = Lpt;
               Inew[2*nnew+1] = Rpt;
               nnew = nnew + 1;
               j++;
               found = 1;
            } else {
               if (found==1) {
                  go = 0;
                  jstart = j - 1;
               } else {
                  j++;
               }
            }
         }
         i++;
      }
   }
   return nnew;
}

void  box2initint(double *alpha, double *beta,
   double *x, double *v, int n,
   double *pI1, double *pI2, double BN,
   double *tBox)
{
   int i;
   double Right, Left;
   for (i=0;i<n;i++) {
      if (v[i]>0) {
         Right = (beta[i]-x[i])/v[i];
         Left = (alpha[i]-x[i])/v[i];
      } else if (v[i]<0) {
         Right = (alpha[i]-x[i])/v[i];
         Left = (beta[i]-x[i])/v[i];
      } else {
         Right = BN;
         Left = -BN;
      }
      *tBox++ = Right;
      if (i==0) {
         *pI1 = Left;
         *pI2 = Right;
      } else {
         if (Right<*pI2) {
            *pI2 = Right;
         }
         if (Left>*pI1) {
            *pI1 = Left;
         }
      }
   }
}

=======
#include "mex.h"
#include <math.h>

void LOCALengineHR(int n, int L,
  double *M[], int *Aidx[], int Adim[],
  double *alpha, double *beta,
  double *x, double *v,
  double *xproj, double *vproj, double *xwork);

void LOCALengineGD(int n, int L,
  double *M[], int *Aidx[], int Adim[],
  double *alpha, double *beta,
  double *x, double *v,
  double *xproj, double *vproj, double *xwork);

void LOCALquadG2quadS(double *pP, double *ABC, double *x0,
    double *v, int n);

void box2initint(double *alpha, double *beta,
   double *x, double *v, int n,
   double *pI1, double *pI2, double BN,
   double *tBox);

int LOCALquadS2int(double *ABC, double *I, double *tPos);

void LOCALquadG2quadS(double *pP, double *ABC, double *x0,
    double *v, int n);

void  box2initint(double *alpha, double *beta,
   double *x, double *v, int n,
   double *pI1, double *pI2, double BN,
   double *tBox);

void mexFunction(int nLHS,  mxArray *pLHS[],
int nRHS,  const mxArray *pRHS[])
{
   mxArray *Mcell, *ActiveIndex, *H, *X, *V, *Xout;
   int L, **Aidx, *Adim, i,j;
   mwSize ndims, *dims, n;
   mxArray *Mi, *Ai;
   double **Mdata, *pAi, *pX, *pV;
   double *LowerBoundBox, *UpperBoundBox;
   double *Xwork, *Xproj, *Vproj;
   int *atmp;
   int usemethod;
   
   usemethod = 0;  /* DEFAULT: use H&R */
   
   /* Rename the inputs */
   Mcell = pRHS[0];
   ActiveIndex = pRHS[1];
   H = pRHS[2];
   X = pRHS[3];
   V = pRHS[4];
   if (nRHS==6) {
      usemethod = 1;  /* GD */
   }

   /* Extract parameter dimension from H (n-by-2) */
   n = mxGetM(H);
   
   /* Get pointers to X and V */
   pX = mxGetPr(X);
   pV = mxGetPr(V);
   
   /* Allocate space for Xwork, Xproj and Vproj */
   Xproj = (double *) mxCalloc(n,sizeof(double));
   Vproj = (double *) mxCalloc(n,sizeof(double));
   Xout = mxCreateDoubleMatrix(n,1,mxREAL);
   pLHS[0] = Xout;
   Xwork = mxGetPr(Xout);

   /* Get Lower/Upper bound Box constraints */
   LowerBoundBox = mxGetPr(H);
   UpperBoundBox = LowerBoundBox + n;
   
   /* Mcell and ActiveIndex are CELLS of the same dimension */
   ndims = mxGetNumberOfDimensions(Mcell);
   dims = mxCalloc(ndims,sizeof(mwSize));
   dims = mxGetDimensions(Mcell);
   L = 1;
   for (i=0;i<ndims;i++) {
      /* mexPrintf("DIMS %d\n",dims[i]); */
      L *= dims[i];
   }
   /* mexPrintf("BASIC %d\n",L); */
   Adim = (int *) mxCalloc(L,sizeof(int));
   Aidx = (int **) mxCalloc(L,sizeof(int *));
   Mdata = (double **) mxCalloc(L,sizeof(double *));
   for (i=0;i<L;i++) {
      Mi = mxGetCell(Mcell,i);
      Ai = mxGetCell(ActiveIndex,i);  /* row vector ???*/
      Adim[i] = mxGetM(Ai)*mxGetN(Ai); /* 1-by?  or ?-by-1 */
      Mdata[i] = mxGetPr(Mi);
      Aidx[i] = (int *) mxCalloc(Adim[i],sizeof(int));
      atmp = Aidx[i];
      pAi = mxGetPr(Ai);
      for (j=0;j<Adim[i];j++) {
         *atmp++ = (int) *pAi++;
      }
   }
   
   if (usemethod==0) {
      /* allocate Xwork bigger, to hold many samples;
      /* put this is FOR-loop
       *Each time, reset pX (starting point) to the column of Xwork that
       *just got created
       *Create a new random direction, filling pV */ 
      LOCALengineHR(n,L,Mdata,Aidx,Adim,
        LowerBoundBox,UpperBoundBox,
        pX,pV,Xproj,Vproj,Xwork);
   } else {
      LOCALengineGD(n,L,Mdata,Aidx,Adim,
        LowerBoundBox,UpperBoundBox,
        pX,pV,Xproj,Vproj,Xwork);
   }
   
   mxFree(Xproj);
   mxFree(Vproj);
}
     

void LOCALengineHR(int n, int L,
  double *M[], int *Aidx[], int Adim[],
  double *alpha, double *beta,
  double *x, double *v,
  double *xproj, double *vproj, double *xwork)
{
   int i, nactive, k, nsingle, *pidx;
   double scalarquad[3], Isingle[4];
   double I1, I2;
   double *Icurrent, *Inew;
   int ncurrent, nnew;
   int Ichoose;
   double tfac, BN, gamma, *tBox;
   double dummy;
      
   BN = mxGetInf();
   /* Handle the box contraints to get initial interval */
   /* tBox is unused in Hit&Run, but is used in GasDynamics */
   tBox = (double *) calloc(n, sizeof(double));
   box2initint(alpha,beta,x,v,n,&I1,&I2,BN,tBox);
   free(tBox);

   Icurrent = (double *) calloc(2*(L+1), sizeof(double));
   Icurrent[0] = I1;  Icurrent[1] = I2;
   ncurrent = 1;
   Inew = (double *) calloc(2*(L+1), sizeof(double));

   for (i=0;i<L;i++) {
      /* i'th Model */
      nactive = Adim[i];
      pidx = Aidx[i];
      /* Project down to active variables */
      for (k=0;k<nactive;k++) {
         xproj[k] = x[*pidx-1];
         vproj[k] = v[*pidx-1];
         pidx++;
      }
      LOCALquadG2quadS(M[i], scalarquad, xproj, vproj, nactive);
      nsingle = LOCALquadS2int(scalarquad, Isingle, &dummy);
      nnew = LOCALIntersect(Icurrent, ncurrent, Isingle, nsingle, Inew);
      /* Replace Current with New */
      ncurrent = nnew;
      for (k=0;k<2*ncurrent;k++) {
         Icurrent[k] = Inew[k];
      }
   }

   /* Choose an interval */
   Ichoose = (int) floor(0.9999*((double) rand())/((double) RAND_MAX)*ncurrent);
   gamma = 0.0001 + 0.9998*((double) rand())/((double) RAND_MAX);
   tfac = gamma*Icurrent[2*Ichoose] + (1-gamma)*Icurrent[2*Ichoose+1];

   /* Update Xwork */
   for (i=0;i<n;i++) {
      xwork[i] = x[i] + tfac*v[i];
   }
   /* Free memory */
   free(Icurrent);
   free(Inew);
}


void LOCALengineGD(int n, int L,
  double *M[], int *Aidx[], int Adim[],
  double *alpha, double *beta,
  double *x, double *v,
  double *xproj, double *vproj, double *xwork)
{
   int i, nactive, k, kk, nsingle, *pidx;
   double scalarquad[3], Isingle[4], MBox[4];
   double *vlocal, *xlocal, *ONproj, *xcontact;
   double I1, I2, *Muse;
   double BN;
   double *tBox, *tGeneral, *tAll;
   int *tIndex, *tIndexSave;
   double iprod, normout,TRemain;
   double BACKOFF, tPos;
   int ncols, minidx, idx, icnt;
   
   BACKOFF = 0.999;
   ncols = 0;
   TRemain = 0.5;
   
   xlocal = (double *) calloc(n, sizeof(double));
   xcontact = (double *) calloc(n, sizeof(double));
   vlocal = (double *) calloc(n, sizeof(double));
   ONproj = (double *) calloc(n, sizeof(double));
   tAll = (double *) calloc(n+L, sizeof(double));
   tIndex = (int *) calloc(n+L, sizeof(int));
   tIndexSave = tIndex;
   
   tBox = tAll;
   tGeneral = tAll + n;
   for (i=0;i<n+L;i++) {
      tIndex[i] = i;  /* fill with 0,...,n+L-1 */
   }
   for (i=0;i<n;i++) {
      vlocal[i] = v[i];
      xlocal[i] = x[i];
   }
   BN = mxGetInf();
   icnt = 0;
   while (TRemain>0 & icnt<6) {
      box2initint(alpha,beta,xlocal,vlocal,n,&I1,&I2,BN,tBox);
      for (i=0;i<L;i++) {
         nactive = Adim[i];
         pidx = Aidx[i];
         for (k=0;k<nactive;k++) {
            xproj[k] = xlocal[*pidx-1];
            vproj[k] = vlocal[*pidx-1];
            pidx++;
         }
         LOCALquadG2quadS(M[i], scalarquad, xproj, vproj, nactive);
         nsingle = LOCALquadS2int(scalarquad, Isingle, tGeneral);
         tGeneral++;
      }
      tGeneral -= L;
      /* LOCALmergesort(tAll, tIndex); */
      tPos = tAll[0];
      minidx = 0;
      for (k=1;k<n+L;k++) {
         if (tAll[k]<tPos) {
            tPos = tAll[k];
            minidx = k;
         }
      }
      if (tPos>TRemain) {
         /* Stop before hitting surface */
         for (i=0;i<n;i++) {
            xlocal[i] += TRemain*vlocal[i];
         }
         TRemain = 0.0;
      } else if (tPos<TRemain) {
         /* Collide with wall, continue on */
         /* update wall collision stats */ 
         ncols++;
         
         /* Xcontact, and Xcontact{BACKOFF} */
         for (i=0;i<n;i++) {
            xcontact[i] = xlocal[i] + tPos*vlocal[i];
            xlocal[i] += BACKOFF*tPos*vlocal[i];
         }
         
         /* Outward Normal at xcontact */
         /* NOTE: all nonactive variables will have 0 gradient */
         if (minidx<n) {
            nactive = 1;
            /* x<=BETA <=> -BETA+x<0 */
            MBox[0] = -beta[minidx];
            MBox[1] = 0.5;
            MBox[2] = 0.5;
            MBox[3] = 0.0;
            Muse = MBox;
            minidx += 1;  /* to agree with MATLAB indexing */
            pidx = &minidx;
         } else {
            minidx = minidx - n;  /* offset */
            nactive = Adim[minidx];
            pidx = Aidx[minidx];
            Muse = M[minidx];
         }
         /* Project down to active variables */
         for (k=0;k<nactive;k++) {
            xproj[k] = xcontact[*pidx-1];
            vproj[k] = vlocal[*pidx-1];
            ONproj[k] = -Muse[k+1];  /* initialize */
            pidx++;
         }
         pidx -= nactive;
         normout = 0;
         for (k=0;k<nactive;k++) {
            for (kk=0;kk<nactive;kk++) {
               idx = kk+1 + (k+1)*(nactive+1);
               ONproj[k] -= xproj[kk]*Muse[idx];
            }
            normout += ONproj[k]*ONproj[k];
         }
         normout = sqrt(normout);
         for (k=0;k<nactive;k++) {
            ONproj[k] = ONproj[k]/normout;
         }
         iprod = 0.0;
         for (k=0;k<nactive;k++) {
            iprod += ONproj[k]*vproj[k];
         }
         /* Update active components of V for reflection */
         for (k=0;k<nactive;k++) {
            vlocal[*pidx-1] -= 2.0*iprod*ONproj[k];
            pidx++;
         }
         pidx -= nactive;
         /* At this point, XLOCAL and VLOCAL are ready to use again */
         /* Update the remaining time, and go back to beginning */
         TRemain -= tPos;
      } else {
         /* Apparently will hit surface at exact time */
         /* Backoff by 10%, and call it done */
         for (i=0;i<n;i++) {
            xlocal[i] += 0.9*TRemain*vlocal[i];
         }
         TRemain = 0.0;
      }
      icnt++;
   }
   /* Update x */
   for (i=0;i<n;i++) {
      xwork[i] = xlocal[i];
   }
   /* Free memory */
   free(xlocal);
   free(xcontact);
   free(vlocal);
   free(ONproj);
   free(tAll);
   free(tIndexSave);
}


int LOCALquadS2int(double *ABC, double *I, double *tPos)
{
   double A, B, C, r, r1, r2, radi, sradi;
   double AINF;
   int nintervals;
   
   AINF = 10000000000000;

   A = *ABC++;
   B = *ABC++;
   C = *ABC++;
   if (A==0.0) {
      if (B==0.0) {
         I[0] = -AINF;
         I[1] = AINF;
         nintervals = 1;
     } else {
        r = -C/B;
        if (r<0) {
           I[0] = r;
           I[1] = AINF;
           nintervals = 1;
        } else if (r>0) {
           I[0] = -AINF;
           I[1] = r;
           nintervals = 1;
        } 
     }
      *tPos = I[1];
  } else {
     radi = B*B - 4*A*C;
     if (radi<0) {
        I[0] = -AINF;
        I[1] = AINF;
        nintervals = 1;
        *tPos = I[1];
     } else if (radi==0.0) {
        r = -B/(2*A);
        if (r<0) {
           I[0] = r;
           I[1] = AINF;
           nintervals = 1;
        } else if (r>0) {
           I[0] = -AINF;
           I[1] = r;
           nintervals = 1;
        }
        *tPos = I[1];
     } else {
        sradi = sqrt(radi);
        /* mexPrintf("RAD quals %g, SRAD is %g\n",radi,sradi); */
        if (A>0) {
            r1 = (-B - sradi)/(2*A);
            r2 = (-B + sradi)/(2*A);
        } else {
            r2 = (-B - sradi)/(2*A);
            r1 = (-B + sradi)/(2*A);
        }
        if ((r1<0) & (r2>0)) {
         I[0] = r1;
         I[1] = r2;
         nintervals = 1;
         *tPos = I[1];
        } else {
         I[0] = -AINF;
         I[1] = r1;
         I[2] = r2;
         I[3] = AINF;
         nintervals = 2;
         if (r1>0) {
            *tPos = r1;
         } else if (r2>0) {
            *tPos = r2;
         } else {
            *tPos = AINF;
         }
         }
      }
   }
   return nintervals;
}



void LOCALquadG2quadS(double *pP, double *ABC, double *x0,
    double *v, int n)
{
   int i, j, nsq;
   double *g22x, *g22v, g12v;
   double A, B, C;
   
   g22x = (double *) calloc(n,sizeof(double));
   g22v = (double *) calloc(n,sizeof(double));
   nsq = (1+n)*(1+n);
   g12v = 0.0;
   for (i=0;i<n;i++) {
      g12v += pP[i+1]*v[i];
      *g22x = 0.0;
      *g22v = 0.0;
      for (j=0;j<n;j++) {
         *g22v += pP[i+1+(j+1)*(n+1)]*v[j]; 
         *g22x += pP[i+1+(j+1)*(n+1)]*x0[j];
      }
      g22v++;
      g22x++;
   }
   g22v -= n;
   g22x -= n;
   A = 0.0;
   B = g12v;
   C = *pP;
   for (i=0;i<n;i++) {
      A += v[i]*g22v[i];
      B += v[i]*g22x[i];
      C += x0[i]*(2.0*pP[i+1] + g22x[i]);
   }
   B = 2.0*B; 
   *ABC++ = A;
   *ABC++ = B;
   *ABC++ = C; 
   free(g22v);
   free(g22x);
}

int LOCALIntersect(double *Icurrent, int ncurrent, 
   double *Isingle, int nsingle, double *Inew)
{
   double Lpt, Rpt;
   double L1i, R1i, L1j, R1j;
   int i, j, jstart, found, go;
   int nnew;

   if (ncurrent==1 & nsingle==1) {
      Lpt = Icurrent[0];
      if (Isingle[0]>Icurrent[0]) {
         Lpt = Isingle[0];
      }
      Rpt = Icurrent[1];
      if (Isingle[1]<Icurrent[1]) {
         Rpt = Isingle[1];
      }
      if (Lpt<Rpt) {
         nnew = 1;
         Inew[0] = Lpt;
         Inew[1] = Rpt;
      } else {
         nnew = 0;
      }
   } else {
      nnew = 0;
      jstart = 0;
      i = 0;
      j = 0;
      while (i<ncurrent) {
         L1i = Icurrent[2*i];
         R1i = Icurrent[2*i+1];
         go = 1;
         found = 0;
         j = jstart;
         while (j<nsingle & go==1) {
            L1j = Isingle[2*j];
            R1j = Isingle[2*j+1];
            Lpt = L1i;
            if (L1j>L1i) {
               Lpt = L1j;
            }
            Rpt = R1i;
            if (R1j<R1i) {
               Rpt = R1j;
            }
            /* Lpt = max(L1i,L1j);
            Rpt = min(R1i,R1j); */
            if (Lpt<Rpt) {
               Inew[2*nnew] = Lpt;
               Inew[2*nnew+1] = Rpt;
               nnew = nnew + 1;
               j++;
               found = 1;
            } else {
               if (found==1) {
                  go = 0;
                  jstart = j - 1;
               } else {
                  j++;
               }
            }
         }
         i++;
      }
   }
   return nnew;
}

void  box2initint(double *alpha, double *beta,
   double *x, double *v, int n,
   double *pI1, double *pI2, double BN,
   double *tBox)
{
   int i;
   double Right, Left;
   for (i=0;i<n;i++) {
      if (v[i]>0) {
         Right = (beta[i]-x[i])/v[i];
         Left = (alpha[i]-x[i])/v[i];
      } else if (v[i]<0) {
         Right = (alpha[i]-x[i])/v[i];
         Left = (beta[i]-x[i])/v[i];
      } else {
         Right = BN;
         Left = -BN;
      }
      *tBox++ = Right;
      if (i==0) {
         *pI1 = Left;
         *pI2 = Right;
      } else {
         if (Right<*pI2) {
            *pI2 = Right;
         }
         if (Left>*pI1) {
            *pI1 = Left;
         }
      }
   }
}

>>>>>>> Stashed changes
