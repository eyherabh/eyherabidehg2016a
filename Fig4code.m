classdef Fig4code < handle
    % This software is provided as supplementary material for the following publication:
    %
    % Eyherabide HG, Neural Stochastic Codes, Encoding and Decoding, arXiv (2017), 
    % https://arxiv.org/abs/1608.05501v2.
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
    % 	fig4 = Fig4code;
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
    end
    
    methods
        function this = Fig4code()
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
                this.di12_p = this.dinidlgauss2Dp1D(this.q,this.rho,this.rho,1);
            end
            val = this.di12_p;
        end
        
        function val = get.info(this)
            if isempty(this.info_p)
                if isempty(this.q) || isempty(this.rho), error('Please specify "q" and "rho"'); end
                this.info_p = this.infogauss2D(this.q,this.rho,this.rho) + this.infogauss1D(1-this.q);
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
    
    methods(Static)
        function amount = infogauss1D(ps1)
            % Computes the information carried by the response about two different stimuli
            % when the responses elicited by each stimulus have one-dimensional Gaussian
            % distributions
            %
            % ps1 is the probability of the first stimulus.
            
            ps2 = 1-ps1;
            pk1       = log(ps1/sqrt(2*pi));
            pk2       = log(ps2/sqrt(2*pi));
            
            opt = {-5,5,'RelTol',1E-3,'AbsTol',1E-6};
            fun2int = @(x)integrand_infogauss(x,pk1,pk2);
            amount = integral(fun2int,opt{:})-ps1*log(ps1)-ps2*log(ps2);
            
            function di = integrand_infogauss(x,k1,k2)
                aux1 = -(x+1).^2/2+k1; aux2 = -(x-1).^2/2+k2;
                ps1x = exp(aux1); ps2x = exp(aux2); ps12x  = ps1x+ps2x;
                di = ps1x.*aux1 + ps2x.*aux2 - ps12x.*log(ps12x);
                di(~isfinite(di)) = 0; di = sum(di,1);
            end
        end
        
        function amount = infogauss2D(ps1,rhos1,rhos2)
            % Computes the information carried by the response about two different stimuli
            % when the responses elicited by each stimulus have two-dimensional Gaussian
            % distributions
            %
            % ps1 is the probability of the first stimulus.
            
            ps2 = 1-ps1;
            params.c1invdet = 1./(1-rhos1.^2);
            params.c2invdet = 1./(1-rhos2.^2);
            params.rhos1mod = rhos1.*params.c1invdet;
            params.rhos2mod = rhos2.*params.c2invdet;
            params.k1       = ps1.*sqrt(params.c1invdet)/2/pi;
            params.k2       = ps2.*sqrt(params.c2invdet)/2/pi;
            
            opt = {-5,5,-5,5,'Method','iterated','RelTol',1E-3,'AbsTol',1E-6};
            fun2int = @(x,y)integrand_infogauss(x,y,params);
            amount = integral2(fun2int,opt{:})-ps1.'*log(ps1)-ps2.'*log(ps2);
            
            function di = integrand_infogauss(x,y,params)
                
                x1c = x+1; y1c = y+1; xy1sum = (x1c.^2 + y1c.^2)/2;
                x2c = x-1; y2c = y-1; xy2sum = (x2c.^2 + y2c.^2)/2;
                
                numx = numel(x);
                aux1 = -params.c1invdet*xy1sum+params.rhos1mod*(x1c.*y1c)...
                    + repmat(log(params.k1),1,numx);
                aux2 = -params.c2invdet*xy2sum+params.rhos2mod*(x2c.*y2c)...
                    + repmat(log(params.k2),1,numx);
                ps1x = exp(aux1); ps2x = exp(aux2); ps12x  = ps1x+ps2x;
                
                di = ps1x.*aux1 + ps2x.*aux2 - ps12x.*log(ps12x);
                di(~isfinite(di)) = 0; di = sum(di,1);
            end
        end
        
        function [di,thetanew] = dinidlgauss2Dp1D(ps,rhos1,rhos2,thetaini)
            % Computes the communication information loss for any two populations, one with two
            % neurons and one with one neuron, which transmit independent information as 
            % defined in the above publication.
            %
            % The responses of all neurons have a Gaussain distribution with unit variance.
            % The responses of the neurons in the first populations may be correlated.
            %
            % ps is the probability of boxes and the probability of B, as set in Section 3.6
            % rhos1 is the correlation coefficient of the responses of the neurons in the first
            % population elicited by boxes
            % rhos2 is the correlation coefficient of the responses of the neurons in the first
            % population elicited by circles (which in Section 3.6 were taken as equal).
            % thetaini is the initial value of theta for computing the communication information loss
            
            optionsmin = optimoptions('fminunc','Algorithm','quasi-newton','Display','off','GradObj','off');
            [thetanew,di] = fminunc(@(theta)dinidlgausstheta2D(ps,rhos1,rhos2,theta)+dinidlgausstheta1D(1-ps,theta),thetaini,optionsmin);
            
            function ditheta = dinidlgausstheta2D(ps,rhos1,rhos2,theta)
                
                params.c1invdet = 1/(1-rhos1.^2);
                params.c2invdet = 1/(1-rhos2.^2);
                params.rhos1mod = rhos1*params.c1invdet;
                params.rhos2mod = rhos2*params.c2invdet;
                params.ps       = ps;
                params.pcs      = 1-ps;
                params.k1       = params.ps*sqrt(params.c1invdet)/2/pi;
                params.k2       = params.pcs*sqrt(params.c2invdet)/2/pi;
                
                optionsint = {-5,5,-5,5,'Method','iterated','RelTol',1E-3,'AbsTol',1E-6};
                ditheta = integral2(@(x,y)integrand_dinidlthetagauss(x,y,params,theta),optionsint{:});
                
                function di = integrand_dinidlthetagauss(x,y,params,th)
                    x1c = x+1; y1c = y+1; xy1sum = (x1c.^2 + y1c.^2)/2;
                    x2c = x-1; y2c = y-1; xy2sum = (x2c.^2 + y2c.^2)/2;
                    numx = numel(x);
                    ps1x  = repmat(params.k1,1,numx).*exp(-params.c1invdet*xy1sum+params.rhos1mod*(x1c.*y1c));
                    ps2x  = repmat(params.k2,1,numx).*exp(-params.c2invdet*xy2sum+params.rhos2mod*(x2c.*y2c));
                    pis1x = params.ps*exp(-th*xy1sum);
                    pis2x = params.pcs*exp(-th*xy2sum);
                    ps12x = ps1x+ps2x;
                    di = ps1x.*log(ps1x./pis1x) + ps2x.*log(ps2x./pis2x) - ps12x.*log(ps12x./(pis1x+pis2x));
                    di(~isfinite(di)) = 0; di = sum(di,1);
                end
            end
            function ditheta = dinidlgausstheta1D(ps,theta)
                params.ps       = ps;
                params.pcs      = 1-ps;
                params.k1       = params.ps/sqrt(2*pi);
                params.k2       = params.pcs/sqrt(2*pi);
                
                optionsint = {-5,5,'RelTol',1E-3,'AbsTol',1E-6};
                ditheta = integral(@(x)integrand_dinidlthetagauss(x,params,theta),optionsint{:});
                
                function di = integrand_dinidlthetagauss(x,params,th)
                    x1c = x+1; x1sum = x1c.^2/2; x2c = x-1; x2sum = x2c.^2/2;
                    ps1x  = params.k1*exp(-x1sum); pis1x = params.ps*exp(-th*x1sum);
                    ps2x  = params.k2*exp(-x2sum); pis2x = params.pcs*exp(-th*x2sum);
                    ps12x = ps1x+ps2x;
                    di = ps1x.*log(ps1x./pis1x) + ps2x.*log(ps2x./pis2x) - ps12x.*log(ps12x./(pis1x+pis2x));
                    di(~isfinite(di)) = 0; di = sum(di,1);
                end
            end
        end
    end
end