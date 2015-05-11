%-------------------------------------------------------------------------%
% This program solves for the steady state in an S-period OLG model with
% i.i.d ability shocks
%
% Author: Richard W. Evans, Brigham Young University, 2010
%         For paper, "OLG Life Cycle Model Transition Paths: New Solution
%         Methods and Applications," joint with Kerk L. Phillips
%         June 2010
%-------------------------------------------------------------------------%

%-------------------------------------------------------------------------%
% Opening Commands
%-------------------------------------------------------------------------%
clear ;
%cd('/Users/rwe2/Documents/BYU Economics/Paper Ideas/OLG Soc Sec/MatLab/Abil') ;
cd('R:\EvansResearch\OLG Soc Sec\MatLab\Abil') ;
starttime = tic ; % start timer

%-------------------------------------------------------------------------%
% Set parameters
%-------------------------------------------------------------------------%
% S     = number of periods an individual lives
% beta  = discount factor (0.96 per year)
% sigma = coefficient of relative risk aversion
% alpha = capital share of income
% rho   = contraction parameter in steady state iteration process
%         representing the weight on the new distribution gamma_new
% A     = total factor productivity parameter in firms' production function
% n     = 1 x S vector of inelastic labor supply for each age s
% e     = S x J matrix of age dependent possible working abilities e_s
% f     = S x J matrix of age dependent discrete probability mass
%         function for e: f(e_s)
% J     = number of points in the support of e
% bmin  = minimum value of b (usually is bmin = 0)
% bmax  = maximum value of b
% bsize = scalar, number of discrete points in the support of b
% b     = 1 x bsize vector of possible values for initial wealth b and
%         savings b'
%-------------------------------------------------------------------------%
S     = 60 ;
beta  = 0.96^(60/S) ;
sigma = 3 ;
alpha = 0.35 ;
rho   = 0.20 ;
A     = 1 ;
n     = [(0.865:(1-0.865)/6:1) ones(1,33) ...
        (1+((0.465-1)/10):((0.465-1)/10):0.465) ...
        (0.465+((0.116-0.465)/5):((0.116-0.465)/5):0.116)...
        (0.116+((0.093-0.116)/5):((0.093-0.116)/5):0.093)] ;

% Set the support and distribution for ability
e = repmat([0.1 0.5 0.8 1 1.2 1.5 1.9],[S,1]) ;
f = repmat([0.04 0.09 0.2 0.34 0.2 0.09 0.04],[S,1]) ;
J = size(e,2) ;

% Set the discretized support of wealth
bmin  = 0 ;
bmax  = 15 ;
bsize = 350 ;
b = (bmin:(bmax-bmin)/(bsize-1):bmax) ;

%-------------------------------------------------------------------------%
% Find the steady state
%-------------------------------------------------------------------------%
% ssiter     = index value that tracks the iteration number
% ssmaxiter  = maximum number of iterations
% ssdist     = sup norm distance measure of the elements of the wealth
%              distribution gamma_init and gamma_new
% ssmindist  = minimum value of the sup norm distance measure below which
%              the process has converged to the steady state
% gamma_init = initial guess for steady state capital distribution:
%              (S-1) x J x bsize array
% gamma_new  = (S-1) x J x bsize array, new iteration distribution of
%              wealth computed from the old distribution gamma_init and the
%              policy rule bprime
% Kss        = steady state aggregate capital stock: scalar
% Nss        = steady state aggregate labor: scalar
% Yss        = steady state aggregate output: scalar
% wss        = steady state real wage: scalar
% rss        = steady state real rental rate: scalar
% phiind     = S x J x bsize policy function indicies for b' = phi(s,e,b).
%              The last row phiind(S,e,b) is ones and corresponds to
%              b_{S+1}. The first row phi(1,e,b) corresponds to b_2.
% phi        = S x J x bsize policy function values for b' = phi(s,e,b).
%              The last row phi(S,e,b) is zeros and corresponds to b_{S+1}.
%              The first row corresponds to b_2.
% Vinit      = bsize x J matrix of values of the state in the next period
%              V(b',e')
% sind       = year index variable in for-loop
% eind       = productivity shock index variable in for-loop
% bind       = wealth level index variable in for-loop
% c          = bsize x J x bsize matrix of values for consumption in the
%              current period: c(b,e,b')
% cposind    = bsize x J x bsize array of = 1 if c > 0 and = 0 if c <= 0
% cpos       = bsize x J x bsize c array with c <= 0 values replaces with
%              positive values very close to zero
% uc         = utility of consumption. The utility of c<0 is set to -10^8
% EVprime    = the expected value of the value function in the next period
% EVprimenew = the expected value of the value function in the next period
%              reshaped to be added to the utility of current consumption
% Vnewarray  = new value function in terms of b, e, and the unoptimized
%              possible values of b'
% Vnew       = the new value function when b' is chosen optimally
% bprimeind  = the index of the optimal
% gamma_ss   = (S-1) x J x bsize array of steady state distribution of
%              wealth
% phiind_ss  = S x J x bsize steady-state policy function indicies for
%              b' = phi(s,e,b)
% phi_ss     = S x J x bsize steady-state policy function values for
%              b' = phi(s,e,b)
%-------------------------------------------------------------------------%

ssiter    = 0 ;
ssmaxiter = 700 ;
ssdist    = 10 ;
ssmindist = 10^(-9) ;
gamma_init = repmat((f(1:S-1,:)./(S-1))./bsize,[1,1,bsize]) ;

while (ssiter < ssmaxiter) && (ssdist >= ssmindist)

   Kss = ((S-1)/S)*...
         sum(sum(sum(gamma_init.*repmat(reshape(b,[1,1,bsize]),...
         [S-1,J,1])))) ;
   Nss = (1/S)*sum(sum(f.*e.*repmat(n',[1,J]))) ;
   Yss = A*(Kss^alpha)*(Nss^(1-alpha)) ;
   wss = (1-alpha)*Yss/Nss ;
   rss = alpha*Yss/Kss ;

   phiind = zeros(S,J,bsize) ;
   phi    = zeros(S,J,bsize) ;
   Vinit  = zeros(bsize,J) ;
   for sind = 1:S
      if sind < S
         c = (1+rss)*repmat(b',[1,J,bsize]) + ...
             (wss*n(S-sind+1))*repmat(e(S-sind+1,:),[bsize,1,bsize]) ...
             - repmat(reshape(b,[1,1,bsize]),[bsize,J,1]) ;
      elseif sind == S
         c = (wss*n(S-sind+1))*repmat(e(S-sind+1,:),[bsize,1,bsize]) ...
             - repmat(reshape(b,[1,1,bsize]),[bsize,J,1]) ;
      end
      cposind = c > 0 ;
      cnonposind = ones(bsize,J,bsize) - cposind ;
      cpos    = c.*cposind + (10^(-8))*cnonposind ;
      uc      = ((cpos.^(1-sigma) - ones(bsize,J,bsize))./(1-sigma)).*...
                cposind - 10^(8)*cnonposind ;
      if sind==1
         EVprime = sum(Vinit,2) ;
      else
         EVprime   = sum(Vinit.*repmat(f(S-sind+2,:),[bsize,1]),2) ;
      end
      EVprimenew = repmat(reshape(EVprime,[1,1,bsize]),[bsize,J,1]) ;
      Vnewarray = uc + beta*(EVprimenew.*cposind) ;
      [Vnew,bprimeind] = max(Vnewarray,[],3) ;
      phiind(S-sind+1,:,:) = reshape(bprimeind',[1,J,bsize]) ;
      phi(S-sind+1,:,:) = reshape(b(bprimeind)',[1,J,bsize]) ;
      Vinit = Vnew ;
   end

   gamma_new = zeros(S-1,J,bsize) ;
   for sind = 1:S-1
      for eind = 1:J
         for bind = 1:bsize
            if sind == 1
               gamma_new(sind,eind,bind) = (1/(S-1))*f(sind+1,eind)*...
                                           sum((phiind(sind,:,1) == bind).*...
                                           f(sind,:)) ;
            else
               gamma_new(sind,eind,bind) = f(sind+1,eind)*...
                                           sum(sum((phiind(sind,:,:) ...
                                           == bind).*...
                                           gamma_init(sind-1,:,:))) ;
            end
         end
      end
   end
   ssiter = ssiter + 1
   ssdist = max(max(max(abs(gamma_new - gamma_init))))
   gamma_init = rho*gamma_new + (1-rho)*gamma_init ;
end

gamma_ss  = gamma_init ;
phiind_ss = phiind ;
phi_ss    = phi ;

runtime = toc(starttime) ;

%-------------------------------------------------------------------------%
% Generate graph of the steady-state distribution of wealth
%-------------------------------------------------------------------------%
% Kssvec     = 1 x S vector of the constant value of Kss
% bsavg      = 1 x S vector of the average wealth b of each age cohort
% b75        = 1 x S vector of the 75th percentile of wealth holdings of
%              each age cohort
% b25        = 1 x S vector of the 25th percentile of wealth holdings of
%              each age cohort
% wealthwgts = (S-1) x J x bsize array of distribution weights times the
%              value of wealth
% bpct       = (S-1) x 1 x bsize array of percent of population with each
%              wealth level for each age cohort
% bpctl      = (S-1) x 1 x bsize array of percentile of each wealth level
%              for each age cohort
% pctl_init  = (S-1) x 1 vector of zeros for initial percentile
%-------------------------------------------------------------------------%

Kssvec = repmat(Kss,[1,S]) ;

bsavg = zeros(1,S) ;
b90  = zeros(1,S) ;
b75  = zeros(1,S) ;
b50  = zeros(1,S) ;
b25  = zeros(1,S) ;
b10  = zeros(1,S) ;

wealthwgts = ((S-1)*gamma_ss).*repmat(reshape(b,[1,1,bsize]),[S-1,J,1]) ;
bpct = (S-1)*sum(gamma_ss,2) ;
bpctl = zeros(S-1,1,bsize) ;
pctl_init = zeros(S-1,1) ;
for bind = 1:bsize
   if bind == 1
      bpctl(:,1,bind) = bpct(:,1,bind) + pctl_init ;
   else
      bpctl(:,1,bind) = bpct(:,1,bind) + bpctl(:,1,bind-1) ;
   end
end
for sind = 2:S
   bsavg(1,sind) = sum(sum(wealthwgts(sind-1,:,:))) ;
   b90diffsq = (bpctl(sind-1,1,:) - 0.90*ones(1,1,bsize)).^2 ;
   b75diffsq = (bpctl(sind-1,1,:) - 0.75*ones(1,1,bsize)).^2 ;
   b50diffsq = (bpctl(sind-1,1,:) - 0.50*ones(1,1,bsize)).^2 ;
   b25diffsq = (bpctl(sind-1,1,:) - 0.25*ones(1,1,bsize)).^2 ;
   b10diffsq = (bpctl(sind-1,1,:) - 0.10*ones(1,1,bsize)).^2 ;
   b90minind = find(b90diffsq == min(b90diffsq),1,'first') ;
   b75minind = find(b75diffsq == min(b75diffsq),1,'first') ;
   b50minind = find(b50diffsq == min(b50diffsq),1,'first') ;
   b25minind = find(b25diffsq == min(b25diffsq),1,'first') ;
   b10minind = find(b10diffsq == min(b10diffsq),1,'first') ;
   b90(1,sind) = b(b90minind) ;
   b75(1,sind) = b(b75minind) ;
   b50(1,sind) = b(b50minind) ;
   b25(1,sind) = b(b25minind) ;
   b10(1,sind) = b(b10minind) ;
end

%-------------------------------------------------------------------------%
% Generate steady state multiplier values
%-------------------------------------------------------------------------%

ssbsavg  = (S-1)*sum(sum(gamma_ss.*...
           repmat(reshape(b,[1,1,bsize]),[S-1,J,1]),3),2) ;
esavg    = sum(e.*f,2) ;
nsavg    = n' ;
cssvec   = (1+rss)*[0;ssbsavg(1:S-2)] + wss*esavg(1:S-1).*nsavg(1:S-1) ...
           - ssbsavg ;
cp1ssvec = (1+rss)*ssbsavg + wss*esavg(2:S).*nsavg(2:S) - ...
           [ssbsavg(2:S-1);0] ;
gxbar    = (cssvec.^(-sigma))./((beta*(1+rss))*cp1ssvec.^(-sigma)) ;

b0array   = [zeros(1,J,bsize); repmat(reshape(b,[1,1,bsize]),[S-2,J,1])] ;
e0array   = repmat(e(1:S-1,:),[1,1,bsize]) ;
n0array   = repmat(n(1:S-1)',[1,J,bsize]) ;
b1array   = phi_ss(1:S-1,:,:) ;
c0array   = (1+rss)*b0array + wss*e0array.*n0array - b1array ;
mu0array  = c0array.^(-sigma) ;
Eb1array  = repmat(b1array,[1,1,1,J]) ;
Een1start = e(2:S,:).*repmat(n(2:S)',[1,J]) ;
Een1array = repmat(reshape(Een1start,[S-1,1,1,J]),[1,J,bsize,1]) ;
b2b1ind   = repmat(phiind_ss,[1,1,1,J]) ;
b2s1ind   = repmat((2:S)',[1,J,bsize,J]) ;
b2e1ind   = repmat(reshape((1:J),[1,1,1,J]),[S-1,J,bsize,1]) ;
Eb2array  = zeros(S-1,J,bsize,J) ;
for sind = 1:S-1
   for e1ind = 1:J
      for bind = 1:bsize
         for e2ind = 1:J
            Eb2array(sind,e1ind,bind,e2ind) = ...
               phi_ss(b2s1ind(sind,e1ind,bind,e2ind),...
               b2e1ind(sind,e1ind,bind,e2ind),...
               b2b1ind(sind,e1ind,bind,e2ind)) ;
         end
      end
   end
end
f1array     = repmat(reshape(f(2:S,:),[S-1,1,1,J]),[1,J,bsize,1]) ;
Ec1array    = (1+rss)*Eb1array + wss*Een1array - Eb2array ;
Emu1array   = sum((Ec1array.^(-sigma)).*f1array,4) ;
lamdif      = (mu0array - (beta*(1+rss))*Emu1array)./beta ;
lam1pos     = phiind_ss(1:S-1,:,:) == 1 ;
lamdifpos   = lamdif > 0 ;
lambda1     = lamdif.*lam1pos.*lamdifpos ;
lam2pos     = phiind(1:S-1,:,:) == bsize ;
lamdifneg   = lamdif < 0 ;
lambda2     = -lamdif.*lam2pos.*lamdifneg ;
lambda1sbar = sum(sum(((S-1)*lambda1.*gamma_ss),3),2) ;
lambda2sbar = sum(sum(((S-1)*lambda2.*gamma_ss),3),2) ;
lamgxbar    = (cssvec.^(-sigma) - beta*lambda1sbar + beta*lambda2sbar)./...
              ((beta*(1+rss))*cp1ssvec.^(-sigma)) ;


save OLGabilSS7.mat ;

figure(1)
plot(1:S,Kssvec,'k--',1:S,bsavg,'LineWidth',3)
axis([0 60 0 1.1*max([Kssvec bsavg])])
xlabel('s-age','FontSize',15)
ylabel('b wealth level','FontSize',15)
%title('Steady state wealth distribution: S = 60','FontSize',20)
legend('steady state K','average steady state b','Location','South')

figure(2)
plot(1:S-1,gxbar)
axis([0 60 0.9*min(gxbar) 1.1*max(gxbar)])
xlabel('s-age','FontSize',15)
ylabel('Euler error gxbar','FontSize',15)
%title('Euler errors: S = 60','FontSize',20)
