classdef Fig3code < handle
    % This software is provided as supplementary material for the following publication:
    %
    % Eyherabide HG, Disambiguating the role of noise correlations when decoding neural
    % populations together, arXiv (2017), https://arxiv.org/abs/1608.05501v2.
    %
    % Should you use this code, I kindly request you to cite the aforementioened publication.
    %
    % DESCRIPTION:
    %
    % Computes the communication information losses depicted in Figure 3 of the
    % aforementioned publication, under the conditions stated in Section 3.5.
    %
    % The computations are based on analytical solutions of the equations defining the
    % information losses, as opposed to a direct numerical computation of these equations
    % using otpimization methods.
    %
    % These computations are performed for different values of q and alpha, where q denotes the 
    % probability of boxes, and alpha denotes the probability of the response [2,2] given a box.
    %
    % EXAMPLE:
    %
    % The following line initializes the object that will perform the computations.
    %	
    % 	fig3 = Fig3code;
    %
    % To effectively compute the losses, you must provide values for q and alpha, for example,
    % as follows
    %
    %   fig3.q = .1;
    %   fig3.alpha = .3;
    %
    % The information measures can be obtained by invoking their correspondin properties.
    % For example, the communication information losses computed using all neurons in all
    % populations can be obtained as follows
    %
    % 	fig3.di12
    %
    % Repeating the above command does not make the object to recompute the value. Instead
    % the object returns the value computed before, which is stored internally by the object.
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
        % Values of alpha
        alpha
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
    end
   
    methods
        function this = Fig3code()
            this.di1_p   = [];
            this.di2_p   = [];
            this.di12_p  = [];
            this.di1p2_p = [];
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
        
        function set.alpha(this,val)
            if val<1 && val>0
                if isempty(this.alpha) || this.alpha~=val
                    this.alpha = val;
                    this.reset;
                end
            else
                error('The value must be greater than zero and less than unity');
            end
        end
        
        function val = get.info(this)
            if isempty(this.info_p)
                if isempty(this.q) || isempty(this.alpha), error('Please specify "q" and "alpha"'); end
                qc          = 1-this.q;
                ac          = 1-this.alpha;
                qa          = this.q*this.alpha;
                qcac        = qc*ac;
                qasum       = qa + qcac;
                p1          = qa./qasum;
                p2          = qcac./qasum;
                this.info_p = 2*(-this.q.*log(this.q)-qc.*log(qc))-qasum.*(-p1.*log(p1)-p2.*log(p2));
            end
            val = this.info_p;
        end
        
        function val = get.di2(this)
            if isempty(this.di2_p)
                if isempty(this.q) || isempty(this.alpha), error('Please specify "q" and "alpha"'); end
                this.di2_p = 0;
            end
            val = this.di2_p;
        end
        
        function val = get.di1p2(this)
            if isempty(this.di1p2_p)
                this.di1p2_p = this.di1+this.di2;
            end
            val = this.di1p2_p;
        end
        
        function val = get.di1(this)
            if isempty(this.di1_p)
                if isempty(this.q) || isempty(this.alpha), error('Please specify "q" and "alpha"'); end
                if this.alpha==0.5
                    % Recall this di1_p is greater than zero only when alpha=0.5
                    this.di1_p = -this.q.*log(this.q)/2;
                else
                    this.di1_p = 0;
                end
            end
            val = this.di1_p;
        end
        
        function val = get.di12(this)
            if isempty(this.di12_p)
                if isempty(this.q) || isempty(this.alpha), error('Please specify "q" and "alpha"'); end
                if this.alpha==0.5
                    this.di12_p = -this.q.*log(this.q)/2;
                else
                    qc = 1-this.q;
                    ac = 1-this.alpha;
                    qa   = this.q*this.alpha;
                    qcac = qc*ac;
                    qasum = 2*qa + qcac;
                    p1 = qa./qasum;
                    p2 = 1-p1;
                    this.di12_p = qasum.*(-p1.*log(p1)-p2.*log(p2))-2*log(2)*qa;
                end
            end
            val = this.di12_p;
        end
        
        function val = get.dint12(this)
            if isempty(this.dint12_p)
                this.dint12_p = this.di12-this.di1p2;
            end
            val = this.dint12_p;
        end
    end
end
