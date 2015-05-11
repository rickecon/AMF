%-------------------------------------------------------------------------%
% This program solves for transition path of the distribution of wealth and
% the aggregate capital stock using the time path iteration (TPI) method
%
% This m-file calls the following other file(s):
%             OLGabilSS7.mat
%
% Author: Richard W. Evans, Brigham Young University, 2010
%         For paper, "OLG Life Cycle Model Transition Paths: New Solution
%         Methods and Applications," joint with Kerk L. Phillips
%         April 2010 (updated April 2010)
%-------------------------------------------------------------------------%

%-------------------------------------------------------------------------%
% Opening Commands
%-------------------------------------------------------------------------%
clear ;
%cd('/Users/rwe2/Documents/BYU Economics/Paper Ideas/OLG Soc Sec/MatLab/Abil') ;
cd('R:\EvansResearch\OLG Soc Sec\MatLab\Abil') ;

%-------------------------------------------------------------------------%
% Import steady state distribution, parameters and other objects from
% steady state computation in OLGabilss.m
%-------------------------------------------------------------------------%
% S        = number of periods an individual lives
% beta     = discount factor (0.96 per year)
% sigma    = coefficient of relative risk aversion
% alpha    = capital share of income
% rho      = contraction parameter in steady state iteration process
%            representing the weight on the new distribution gamma_new
% A        = total factor productivity parameter in firms' production
%            function
% n        = 1 x S vector of inelastic labor supply for each age s
% e        = S x J matrix of age dependent possible working abilities e_s
% f        = S x J matrix of age dependent discrete probability mass
%            function for e: f(e_s)
% J        = number of points in the support of e
% bmin     = minimum value of b (usually is bmin = 0)
% bmax     = maximum value of b
% bsize    = scalar, number of discrete points in the support of b
% b        = 1 x bsize vector of possible values for initial wealth b and
%            savings b'
% gamma_ss = (S-1) x J x bsize array of steady state distribution of
%            wealth
% Kss      = steady state aggregate capital stock: scalar
% Nss      = steady state aggregate labor: scalar
% Yss      = steady state aggregate output: scalar
% wss      = steady state real wage: scalar
% rss      = steady state real rental rate: scalar
% phiind_ss = S x J x bsize steady-state policy function indicies for
%             b' = phi(s,e,b). The last row phiind(S,e,b) is ones and
%             corresponds to b_{S+1}. The first row phi(1,e,b) corresponds
%             to b_2.
% phi_ss    = S x J x bsize steady-state policy function values for
%             b' = phi(s,e,b). The last row phi(S,e,b) is zeros and
%             corresponds to b_{S+1}. The first row corresponds to b_2
%-------------------------------------------------------------------------%
load OLGabilSS7.mat S beta sigma alpha rho A n e f J bmin bmax bsize b ...
     gamma_ss Kss Nss Yss wss rss phiind_ss phi_ss ;

starttime = tic ; % start timer

%-------------------------------------------------------------------------%
% Set other parameters and objects
%-------------------------------------------------------------------------%
% T       = number of periods until the steady state
% gamma0  = (S-1) x J x bsize array, initial distribution of wealth
% K0      = initial aggregate capital stock as a function of the initial
%           distribution of wealth
% rho_TPI = contraction parameter in TPI process representing weight on new
%           time path of aggregate capital stock
%-------------------------------------------------------------------------%
T       = 60 ;
I0      = eye(bsize) ;
I0row   = I0(155,:) ;
gamma0  = repmat(reshape(I0row,[1,1,bsize]),[S-1,J,1]).*...
          repmat(f(1:S-1,:,:),[1,1,bsize])./(S-1) ;
K0      = ((S-1)/S)*...
          sum(sum(sum(gamma0.*repmat(reshape(b,[1,1,bsize]),...
          [S-1,J,1])))) ;
rho_TPI = 0.2 ;

%-------------------------------------------------------------------------%
% Solve for equilibrium transition path by TPI
%-------------------------------------------------------------------------%
% Kinit   = 1 x T vector, initial time path of aggregate capital stock
% Ninit   = 1 x T vector, initial time path of aggregate labor demand. This
%           is just equal to a 1 x T vector of Nss because labor is
%           supplied inelastically
% Yinit   = 1 x T vector, initial time path of aggregate output
% winit   = 1 x T vector, initial time path of real wage
% rinit   = 1 x T vector, initial time path of real interest rate
% gammat  = (S-1) x J x bsize x T array time path of the distribution of
%           capital
% phiindt = S x J x bsize x T array time path of policy function indicies
%           for b' = phi(s,e,b,t). The last row phiindt(S,e,b,t) is ones
%           and corresponds to b_{S+1}=0. The first row phi(1,e,b,t)
%           corresponds to b_2.
% phit    = S x J x bsize x T array time path of policy function values for
%           b' = phi(s,e,b,t). The last row phi(S,e,b,t) is zeros and
%           corresponds to b_{S+1}=0. The first row corresponds to b_2.
% p1aind  = index of period-1 age
% Vinit   = bsize x J matrix of values of the state in the next period:
%           V(b',e')
% sind    = index of age from period 1
% tind    = index of time period
%-------------------------------------------------------------------------%
% c          = bsize x J x bsize matrix of values for consumption in the
%              current period: c(b,e,b')
% cposind    = bsize x J x bsize array of = 1 if c > 0 and = 0 if c <= 0
% cpos       = bsize x J x bsize c array with c <= 0 values replaces with
%              positive values very close to zero
% bposind    = bsize x J x bsize array of = 1 if c >= 0 and = 0 if c < 0.
%              This matrix is important because it allows for b'=0 to be
%              the optimal choice when income equals zero
% uc         = utility of consumption. The utility of c<0 is set to -10^8
% EVprime    = the expected value of the value function in the next period
% EVprimenew = the expected value of the value function in the next period
%              reshaped to be added to the utility of current consumption
% Vnewarray  = new value function in terms of b, e, and the unoptimized
%              possible values of b'
% Vnew       = the new value function when b' is chosen optimally
% bprimeind  = the index of the optimal
%-------------------------------------------------------------------------%
Kinit = [(K0:(Kss-K0)/(T-1):Kss) Kss*ones(1,S-2)] ; % This is a linear initial time path
Ninit = Nss*ones(1,T+S-2) ;
Yinit = A*((Kinit.^alpha).*(Ninit.^(1-alpha))) ;
winit = (1-alpha)*(Yinit./Ninit) ;
rinit = alpha*(Yinit./Kinit) ;

TPIiter    = 0 ;
TPImaxiter = 500 ;
TPIdist    = 10 ;
TPImindist = 3.0*10^(-6) ; % finishes at iteration 55

while (TPIiter < TPImaxiter) && (TPIdist >= TPImindist)

   gammat  = zeros(S-1,J,bsize,T) ;
   gammat(:,:,:,1) = gamma0 ;
   phiindt = zeros(S,J,bsize,T+S-2) ;
   phit    = zeros(S,J,bsize,T+S-2) ;

   % Solve for the lifetime savings policy rules of each type of person at t=1 
   for p1aind = 1:S
      Vinit = zeros(bsize,J) ;
      for sind = 1:p1aind
         if sind < S
            c = (1+rinit(p1aind-sind+1))*repmat(b',[1,J,bsize]) + ...
                (winit(p1aind-sind+1)*n(S-sind+1))*repmat(e(S-sind+1,:),...
                [bsize,1,bsize]) - repmat(reshape(b,[1,1,bsize]),...
                [bsize,J,1]) ;
         elseif sind == S
            c = (winit(p1aind-sind+1)*n(S-sind+1))*repmat(e(S-sind+1,:),...
                [bsize,1,bsize]) - repmat(reshape(b,[1,1,bsize]),...
                [bsize,J,1]) ;
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
         phiindt(S-sind+1,:,:,p1aind-sind+1) = reshape(bprimeind',...
                                               [1,J,bsize,1]) ;
         phit(S-sind+1,:,:,p1aind-sind+1) = reshape(b(bprimeind)',...
                                            [1,J,bsize,1]) ;
         Vinit = Vnew ;
      end
   end

   % Solve for the lifetime savings policy rules of each age-1 person from t=2
   % to t=T-1
   for tind = 2:T-1
      Vinit = zeros(bsize,J) ;
      for sind = 1:S
         if sind < S
            c = (1+rinit(tind+S-sind))*repmat(b',[1,J,bsize]) + ...
                (winit(tind+S-sind)*n(S-sind+1))*repmat(e(S-sind+1,:),...
                [bsize,1,bsize]) - repmat(reshape(b,[1,1,bsize]),...
                [bsize,J,1]) ;
         elseif sind == S
            c = (winit(tind+S-sind)*n(S-sind+1))*repmat(e(S-sind+1,:),...
                [bsize,1,bsize]) - repmat(reshape(b,[1,1,bsize]),...
                [bsize,J,1]) ;
         end
         cposind = c > 0 ;
         cnonposind = ones(bsize,J,bsize) - cposind ;
         cpos    = c.*cposind + (10^(-8))*cnonposind ;
         uc      = ((cpos.^(1-sigma) - ones(bsize,J,bsize))./(1-sigma)).*...
                   cposind - 10^(8)*cnonposind ;
         if sind==1
            EVprime = sum(Vinit,2) ;
         else
            EVprime = sum(Vinit.*repmat(f(S-sind+2,:),[bsize,1]),2) ;
         end
         EVprimenew = repmat(reshape(EVprime,[1,1,bsize]),[bsize,J,1]) ;
         Vnewarray = uc + beta*EVprimenew.*cposind ;
         [Vnew,bprimeind] = max(Vnewarray,[],3) ;
         phiindt(S-sind+1,:,:,tind+S-sind) = reshape(bprimeind',...
                                             [1,J,bsize,1]) ;
         phit(S-sind+1,:,:,tind+S-sind) = reshape(b(bprimeind)',...
                                          [1,J,bsize,1]) ;
         Vinit = Vnew ;
      end
   end

   % Generate new time path for distribution of capital gammat and of the
   % aggregate capital stock K
   Knew = [K0 zeros(1,T-1) Kss*ones(1,S-2)] ;
   for tind = 2:T
      for sind = 1:S-1
         for eind = 1:J
            for bind = 1:bsize
               if sind == 1
                  gammat(sind,eind,bind,tind) = (1/(S-1))*f(sind+1,eind)*...
                                sum((phiindt(sind,:,1,tind-1) == bind).*...
                                f(sind,:)) ;
               else
                  gammat(sind,eind,bind,tind) = f(sind+1,eind)*...
                                      sum(sum((phiindt(sind,:,:,tind-1) ...
                                      == bind).*...
                                      gammat(sind-1,:,:,tind-1))) ;
               end
            end
         end
      end
      Knew(1,tind) = ((S-1)/S)*sum(sum(sum(gammat(:,:,:,tind).*...
                     repmat(reshape(b,[1,1,bsize]),[S-1,J,1])))) ;
   end
   TPIiter = TPIiter + 1
   TPIdist = max(abs(Knew - Kinit))
   Kinit = rho_TPI*Knew + (1-rho_TPI)*Kinit ;
   Yinit = A*((Kinit.^alpha).*(Ninit.^(1-alpha))) ;
   winit = (1-alpha)*(Yinit./Ninit) ;
   rinit = alpha*(Yinit./Kinit) ;
end

Kpath_TPI = Kinit ;
gammat_TPI = gammat ;

runtime = toc(starttime) ; % end timer

%-------------------------------------------------------------------------%
% Generate graph of the transition path of the aggregate capital stock
%-------------------------------------------------------------------------%
Kssvec = Kss*ones(1,T+10) ;
save OLGabilTPI9.mat ;

figure(1)
plot(1:T+10,Kpath_TPI(1:T+10),1:T+10,Kssvec,'k--','LineWidth',2)
axis([0 T+11 0.9*min(Kpath_TPI) 1.1*max(Kpath_TPI)])
xlabel('time t','FontSize',15)
ylabel('aggregate capital K','FontSize',15)
%title('Transition path of aggregate capital K: S = 60','FontSize',20)
legend('TPI time path K_t','steady-state K','Location','East')