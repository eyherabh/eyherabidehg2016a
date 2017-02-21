/*
 This software is provided as supplementary material for the following publication:

 Eyherabide HG, Neural Stochastic Codes, Encoding and Decoding, arXiv (2017),
 https://arxiv.org/abs/1608.05501v2.

 Should you use this code, I kindly request you to cite the aforementioened publication.

 DESCRIPTION:

 Computes the communication information losses depicted in Figure 4 of the
 aforementioned publication, under the conditions stated in Section 3.6.

 This code is part of the Matlab class Fig4codeC.m. It requires to download 
 some additional libraries and to compile it within Matlab.

 Specifically, the code requires the following libraries

 - GSL (https://www.gnu.org/software/gsl/)
 - Cubature (http://ab-initio.mit.edu/wiki/index.php/Cubature)
 
 They should be installed wherever #include looks for headers, or
 else, the folders in the #include statements within the c-files
 (mex-files) should be modified.

 The code can be compiled as follows
 
   mex -v GCC='/usr/bin/gcc-4.7' -lm infoGauss.c

 where you should replace /usr/bin/gcc-4.7 for the appropriate folder
 and C compiler compatible with your Matlab installation.

 VERSION CONTROL

 V1.000 Hugo Gabriel Eyherabide (10 Feb 2017)

 Should you find bugs, please contact Hugo Gabriel Eyherabide (neuralinfo@eyherabidehg.com)

 LICENSE

 Copyright (c) 2017, Hugo Gabriel Eyherabide
 All rights reserved.

 Redistribution and use in source and binary forms, with or without modification,
 are permitted provided that the following conditions are met:

 1.  Redistributions of source code must retain the above copyright notice,
     this list of conditions and the following disclaimer.

 2.  Redistributions in binary form must reproduce the above copyright notice,
     this list of conditions and the following disclaimer in the documentation
     and/or other materials provided with the distribution.

 3.  Neither the name of the copyright holder nor the names of its contributors
     may be used to endorse or promote products derived from this software
     without specific prior written permission.

 THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.
 IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT,
 INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT
 NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
 PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
 WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY
 OF SUCH DAMAGE.
*/


#include<mex.h>
#include<math.h>
#include<ctype.h>
#include<cubature/cubature.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_min.h>

int NDIntegrand(unsigned xdim, unsigned numx, const double *x, void *par, unsigned didim, double *dival)
{
    unsigned    indx;
    unsigned    indd;
    unsigned    indmu;
    double      *params = (double*) par;
    
    double      psx;
    double      pisx;
    double      px;
    double      pix;
    double      dinow;
    double      xc;
    double      xc2sum;
    double      xcprod;
    double      mu[2] = {1,-1};
    
    for(indx=0; indx<numx; indx++)
    {
        dinow = 0;
        px  = 0;
        pix = 0;
        
        for(indmu=0;indmu<2;indmu++)
        {
            xc2sum = 0;
            xcprod = 1;
                
            for(indd=0;indd<xdim;indd++)
            {
                xc      = (*x) + mu[indmu];
                xcprod *= xc;
                xc2sum += xc*xc;
                x++;
            }
            x -= xdim;
            xc2sum *= -0.5;

            params += indmu;
            psx  = params[2]*exp(params[4]*xc2sum+params[6]*xcprod);
            pisx = params[0]*exp(params[8]*xc2sum);
            params -= indmu; 
            
            px  += psx;
            pix += pisx;   
        
            dinow += psx * log(psx/pisx);
            if(!isfinite(dinow)) break;
        }
        
        if(isfinite(dinow)) 
        {
            dinow -= px * log(px/pix);
            dival[indx] = isfinite(dinow) ? dinow : 0;
        }
        else
            dival[indx] = 0;
            
        x+=xdim;    
    }
    x -= numx*xdim;
    return 0;
}
        
double diTheta(double th, double * const par)
{
    double  dival2D;
    double  dival1D;
    double  errorval;
    double  params[10];
    double  xmin[2]={-5,-5};     
    double  xmax[2]={5,5};     
    
    params[0] = par[0];
    params[1] = 1-par[0];
    params[4] = 1.0/(1.0-par[1]*par[1]); 
    params[5] = 1.0/(1.0-par[2]*par[2]);
    params[2] = params[0]*sqrt(params[4])/(2.0*M_PI);
    params[3] = params[1]*sqrt(params[5])/(2.0*M_PI);
    params[6] = par[1]*params[4];
    params[7] = par[2]*params[5];
    params[8] = th; 
    params[9] = th; 
    
    hcubature_v(1, NDIntegrand, params, 2, xmin, xmax, 1000, 1E-6,1E-3,ERROR_INDIVIDUAL,&dival2D,&errorval);
    
    params[0] = params[1];
    params[1] = par[0];
    params[4] = 1; 
    params[5] = 1;
    params[2] = params[0]/sqrt(2.0*M_PI);
    params[3] = params[1]/sqrt(2.0*M_PI);
    params[6] = 0;
    params[7] = 0;
    
    hcubature_v(1, NDIntegrand, params, 1, xmin, xmax, 1000, 1E-6,1E-3,ERROR_INDIVIDUAL,&dival1D,&errorval);

    return dival1D+dival2D;
}
      
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    double  *par = (double*) mxGetPr(prhs[0]);
    plhs[0] = mxCreateDoubleScalar(mxREAL);
    double  *di = (double*) mxGetPr(plhs[0]);
    
    double thl = -0.5;
    double thm = 0.5;
    double thr = 1.5;
    double dil = diTheta(thl,par);
    double dim = diTheta(thm,par);
    double dir = diTheta(thr,par);
    
    /* Looking for lower limit of minimization interval*/
    while(dil<dim)
    {
        dir=dim; dim=dil;
        thr=thm; thm=thl; 
        dil = diTheta(thl*=2,par);
    }
    
    /* Looking for upper limit of minimization interval*/
    while(dir<dim)
    {
        dil=dim; dim=dir;
        thl=thm; thm=thr;
        dir=diTheta(thr*=2,par);
    }
    
    /* Minimizing the communication information loss*/
    gsl_function dith;
    dith.function = &diTheta;
    dith.params = par;
    int     iter = 0;
    int     max_iter = 1000;
    gsl_min_fminimizer *s = gsl_min_fminimizer_alloc (gsl_min_fminimizer_brent);
    
    gsl_min_fminimizer_set (s, &dith, thm, thl, thr);
    
    do
    {
        gsl_min_fminimizer_iterate (s);
        thl = gsl_min_fminimizer_x_lower(s);
        thr = gsl_min_fminimizer_x_upper(s);
    }
    while (gsl_min_test_interval(thl, thr, 1E-6, 1E-3) == GSL_CONTINUE && iter < max_iter);
    
    *di = gsl_min_fminimizer_f_minimum (s);
        
    gsl_min_fminimizer_free (s);
}

#include<cubature/hcubature.c>

