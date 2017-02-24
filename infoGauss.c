/*
 This software is provided as supplementary material for the following publication:

 Eyherabide HG, Disambiguating the role of noise correlations when decoding neural
 populations together, arXiv (2017), https://arxiv.org/abs/1608.05501v2.

 Should you use this code, I kindly request you to cite the aforementioened publication.

 DESCRIPTION:

 Computes the information for the examples depicted in Figure 4 of the
 aforementioned publication, under the conditions stated in Section 3.6.

 This code is part of the Matlab class Fig4codeC. It requires to download
 some additional libraries and compile it within Matlab.

 Specifically, the code requires the following library

 - Cubature (http://ab-initio.mit.edu/wiki/index.php/Cubature)
 
 It should be installed wherever #include looks for headers, or
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
 
 
int NDIntegrand(unsigned xdim, unsigned numx, const double *x, void *par, unsigned infodim, double *infoval)
{
    unsigned    indx;
    unsigned    indd;
    unsigned    indmu;
    double      *params = (double*) par;
    
    double      psx;
    double      px;
    double      infonow;
    double      xc;
    double      xc2sum;
    double      xcprod;
    double      mu[2] = {1,-1};
    
    for(indx=0; indx<numx; indx++)
    {
        infonow = 0;
        px  = 0;
        
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
            params -= indmu; 
            
            px  += psx;
        
            infonow += psx * log(psx);
            if(!isfinite(infonow)) break;
        }
        
        if(isfinite(infonow)) 
        {
            infonow -= px * log(px);
            infoval[indx] = isfinite(infonow) ? infonow : 0;
        }
        else
            infoval[indx] = 0;
            
        x+=xdim;    
    }
    x -= numx*xdim;
    return 0;
}
        
      
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    double  *par = (double*) mxGetPr(prhs[0]);
    int     *parnum = (int*) mxGetDimensions(prhs[0]);
    if(parnum[0]<parnum[1]) parnum[0]=parnum[1];
    plhs[0] = mxCreateDoubleScalar(mxREAL);
    double  *infoval = (double*) mxGetPr(plhs[0]);

    double  errorval;
    double  params[8];
    double  xmin[2]={-5,-5};     
    double  xmax[2]={5,5};     
    
    params[0] = par[0];
    params[1] = 1-par[0];
    if(parnum[0]==1)
    {
        params[4] = 1.0; 
        params[5] = 1.0;
        params[2] = params[0]/sqrt(2.0*M_PI);
        params[3] = params[1]/sqrt(2.0*M_PI);
        params[6] = 0;
        params[7] = 0;
    }
    else
    {
        params[4] = 1.0/(1.0-par[1]*par[1]); 
        params[5] = 1.0/(1.0-par[2]*par[2]);
        params[2] = params[0]*sqrt(params[4])/(2.0*M_PI);
        params[3] = params[1]*sqrt(params[5])/(2.0*M_PI);
        params[6] = par[1]*params[4];
        params[7] = par[2]*params[5];
        parnum[0] = 2;
    }
    
    
    hcubature_v(1, NDIntegrand, params, parnum[0], xmin, xmax, 1000, 1E-6,1E-3,ERROR_INDIVIDUAL,infoval,&errorval);
    (*infoval)+=-params[0]*log(params[0])-params[1]*log(params[1]);
}

#include<cubature/hcubature.c>

