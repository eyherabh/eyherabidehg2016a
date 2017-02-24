classdef Fig4codeC < handle
    % This software is provided as supplementary material for the following publication:
    %
    % Eyherabide HG, Disambiguating the role of noise correlations when decoding neural
    % populations together, arXiv (2017), https://arxiv.org/abs/1608.05501v2.
    %
    % Should you use this code, I kindly request you to cite the aforementioened publication.
    %
    % DESCRIPTION:
    %
    % Computes the communication information losses depicted in Figure 4 of the
    % aforementioned publication, under the conditions stated in Section 3.6.
    %
    % The computations are performed for different values of q and rho, where q denotes the
    % probability of boxes, and rho denotes the correlation coefficient.
    %
    % EXAMPLE:
    %
    % The following line initializes the object that will perform the computations.
    %
    % 	fig4 = Fig4codeC;
    %
    % To effectively compute the losses, you must provide values for q and rho, for example,
    % as follows
    %
    %   fig4.q = .1;
    %   fig4.rho = .3;
    %
    % The information measures can be obtained by invoking their correspondin properties.
    % For example, the communication information losses computed using all neurons in all
    % populations can be obtained as follows
    %
    % 	fig4.di12
    %
    % Repeating the above command does not make the object to recompute the value. Instead
    % the object returns the value computed before, which is stored internally by the object.
    %
    % REMARKS
    %
    % This code is analogous to Fig4code, but with the functions computing
    % information and information losses written in C. As a result, the code
    % is much faster, but it requires to download some additional libraries
    % and compile it within Matlab.
    %
    % Specifically, the code requires the following libraries
    %
    % - GSL (https://www.gnu.org/software/gsl/)
    % - Cubature (http://ab-initio.mit.edu/wiki/index.php/Cubature)
    % 
    % They should be installed wherever #include looks for headers, or
    % else, the folders in the #include statements within the c-files
    % (mex-files) should be modified.
    %
    % The functions that must be compiled are called
    %
    % - infoGauss.c
    % - dinidlGaussTheta.c
    %
    % These files can be compiled as follows
    %
    %   mex -v GCC='/usr/bin/gcc-4.7' -lm infoGauss.c
    %   mex -v GCC='/usr/bin/gcc-4.7' -lgsl -lgslcblas -lm dinidlGaussTheta.c
    %
    % where you should replace /usr/bin/gcc-4.7 for the appropriate folder
    % and C compiler compatible with your Matlab installation.
    %
    % VERSION CONTROL
    %
    % V1.000 Hugo Gabriel Eyherabide (10 Feb 2017)
    %
    % Should you find bugs, please contact Hugo Gabriel Eyherabide (neuralinfo@eyherabidehg.com)
    %
    % LICENSE
    %
    % Copyright (c) 2017, Hugo Gabriel Eyherabide
    % All rights reserved.
    %
    % Redistribution and use in source and binary forms, with or without modification,
    % are permitted provided that the following conditions are met:
    %
    % 1.  Redistributions of source code must retain the above copyright notice,
    %     this list of conditions and the following disclaimer.
    %
    % 2.  Redistributions in binary form must reproduce the above copyright notice,
    %     this list of conditions and the following disclaimer in the documentation
    %     and/or other materials provided with the distribution.
    %
    % 3.  Neither the name of the copyright holder nor the names of its contributors
    %     may be used to endorse or promote products derived from this software
    %     without specific prior written permission.
    %
    % THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
    % AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
    % WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.
    % IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT,
    % INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT
    % NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
    % PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
    % WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
    % ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY
    % OF SUCH DAMAGE.
    
    properties
        % Values of q
        q
        % Values of rho
        rho
    end
    
    properties(Dependent)
        % Communication information loss when decoding only frames
        di1
        % Communication information loss when decoding only letters
        di2
        % Communication information loss caused by parallel NI decoders
        di1p2
        % Communication information loss caused by joint NI decoders
        di12
        % Total transmmitted information
        info
        % Communication destructive interference
        dint12
    end
    
    properties(Access = private)
        di1_p
        di2_p
        di1p2_p
        di12_p
        info_p
        dint12_p
        isoctave
    end
    
    methods
        function this = Fig4codeC()
            this.di1_p   = [];
            this.di2_p   = [];
            this.di1p2_p = [];
            this.di12_p  = [];
            this.info_p  = [];
        end
        
        function set.q(this,val)
            if val<1 && val>0
                if isempty(this.q) || this.q~=val
                    this.q = val;
                    this.reset;
                end
            else
                error('The value must be greater than zero and less than unity');
            end
        end
        
        function reset(this)
            this.di1_p   = [];
            this.di2_p   = [];
            this.di12_p  = [];
            this.di1p2_p = [];
            this.info_p  = [];
        end
        
        function set.rho(this,val)
            if val<1 && val>-1
                if isempty(this.rho) || this.rho~=val
                    this.rho = val;
                    this.reset;
                end
            else
                error('The value must be greater than minus one and less than one');
            end
        end
        
        function val = get.di12(this)
            if isempty(this.di12_p)
                if isempty(this.q) || isempty(this.rho), error('Please specify "q" and "rho"'); end
                this.di12_p = dinidlGaussTheta([this.q,this.rho,this.rho]);
            end
            val = this.di12_p;
        end
        
        function val = get.info(this)
            if isempty(this.info_p)
                if isempty(this.q) || isempty(this.rho), error('Please specify "q" and "rho"'); end
                this.info_p = infoGauss([this.q,this.rho,this.rho])+infoGauss(1-this.q);
            end
            val = this.info_p;
        end
        
        function val = get.di1(this)
            if isempty(this.di1_p)
                if isempty(this.q) || isempty(this.rho), error('Please specify "q" and "rho"'); end
                this.di1_p = 0;
            end
            val = this.di1_p;
        end
        
        function val = get.di2(this)
            if isempty(this.di2_p)
                if isempty(this.q) || isempty(this.rho), error('Please specify "q" and "rho"'); end
                this.di2_p = 0;
            end
            val = this.di2_p;
        end
        
        function val = get.di1p2(this)
            if isempty(this.di1p2_p)
                if isempty(this.q) || isempty(this.rho), error('Please specify "q" and "rho"'); end
                this.di1p2_p = 0;
            end
            val = this.di1p2_p;
        end
        
    end
end
