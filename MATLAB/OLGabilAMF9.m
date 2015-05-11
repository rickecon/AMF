%-------------------------------------------------------------------------%
% This program solves for transition path of the distribution of wealth and
% the aggregate capital stock using the alternate model forecast (AMF)
% method
%
% This m-file calls the following other file(s):
%             OLGabilSS7.mat
%             OLGabilTPI9.mat
%
% Author: Richard W. Evans, Brigham Young University, 2010
%         For paper, "OLG Life Cycle Model Transition Paths: New Solution
%         Methods and Applications," joint with Kerk L. Phillips
%         April 2010 (updated June 2010)
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
% S         = number of periods an individual lives
% beta      = discount factor (0.96 per year)
% sigma     = coefficient of relative risk aversion
% alpha     = capital share of income
% rho       = contraction parameter in steady state iteration process
%             representing the weight on the new distribution gamma_new
% A         = total factor productivity parameter in firms' production
%             function
% n         = 1 x S vector of inelastic labor supply for each age s
% e         = S x J matrix of age dependent possible working abilities e_s
% f         = S x J matrix of age dependent discrete probability mass
%             function for e: f(e_s)
% J         = number of points in the support of e
% bmin      = minimum value of b (usually is bmin = 0)
% bmax      = maximum value of b
% bsize     = scalar, number of discrete points in the support of b
% b         = 1 x bsize vector of possible values for initial wealth b and
%             savings b'
% gamma_ss  = (S-1) x J x bsize array of steady state distribution of
%             wealth
% Kss       = steady state aggregate capital stock: scalar
% Nss       = steady state aggregate labor: scalar
% Yss       = steady state aggregate output: scalar
% wss       = steady state real wage: scalar
% rss       = steady state real rental rate: scalar
% phiind_ss = S x J x bsize steady-state policy function indicies for
%             b' = phi(s,e,b). The last row phiind(S,e,b) is ones and
%             corresponds to b_{S+1}. The first row phi(1,e,b) corresponds
%             to b_2.
% phi_ss    = S x J x bsize steady-state policy function values for
%             b' = phi(s,e,b). The last row phi(S,e,b) is zeros and
%             corresponds to b_{S+1}. The first row corresponds to b_2
% Kpath_TPI = computed time path of aggregate capital stock from TPI method
% T         = number of periods until the steady state
% gamma0    = (S-1) x J x bsize array, initial distribution of wealth
% K0        = initial aggregate capital stock as a function of the initial
%             distribution of wealth
%-------------------------------------------------------------------------%
load OLGabilSS7.mat S beta sigma alpha rho A n e f J bmin bmax bsize b ...
     gamma_ss Kss Nss Yss wss rss phiind_ss phi_ss ;
load OLGabilTPI9.mat Kpath_TPI T gamma0 K0 ;

starttime = tic ; % start timer

%-------------------------------------------------------------------------%
% Solve for equilibrium transition path by AMF method
%-------------------------------------------------------------------------%
% Kpath_AMF = transition path for aggregate capital by AMF method
% gammat    = (S-1) x J x bsize x T array time path of the distribution of
%             capital
% phiindt   = S x J x bsize x T+S-2 array forecasted time path of policy
%             function indicies
%             for b' = phi(s,e,b,t). The last row phiindt(S,e,b,t) is ones
%             and corresponds to b_{S+1}=0. The first row phi(1,e,b,t)
%             corresponds to b_2.
% phit      = S x J x bsize x T array time path of policy function values
%             for b' = phi(s,e,b,t). The last row phi(S,e,b,t) is zeros and
%             corresponds to b_{S+1}=0. The first row corresponds to b_2.

% Kpathf    = 1 x T+S-2 vector at initial forecast, forecast time path of
%             aggregate capital stock K
% Npathf    = 1 x T vector, initial time path of aggregate labor demand. This
%           is just equal to a 1 x T vector of Nss because labor is
%           supplied inelastically
% Yinit   = 1 x T vector, initial time path of aggregate output
% winit   = 1 x T vector, initial time path of real wage
% rinit   = 1 x T vector, initial time path of real interest rate


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
Kpath_AMF = [K0 zeros(1,T-1) Kss*ones(1,S-2)] ;
gammat  = zeros(S-1,J,bsize,T) ;
gammat(:,:,:,1) = gamma0 ;
phiindt = zeros(S,J,bsize,T+S-2) ;
phit    = zeros(S,J,bsize,T+S-2) ;

[1 K0]
for tind = 1:T-1
   % Alternative forecast is linear trend between K_t and Kss
   Kpathf = [(Kpath_AMF(tind):(Kss-Kpath_AMF(tind))/(T-tind):Kss) ...
            Kss*ones(1,S-2)] ;
   Npathf = Nss*ones(1,T+S-tind-1) ;
   Ypathf = A*((Kpathf.^alpha).*(Npathf.^(1-alpha))) ;
   wpathf = (1-alpha)*(Ypathf./Npathf) ;
   rpathf = alpha*(Ypathf./Kpathf) ;
   
   % Solve for the lifetime savings policy rules of each type of person at
   % time t
   for p1aind = 1:S
      Vinit = zeros(bsize,J) ;
      for sind = 1:p1aind
         if sind < S
            c = (1+rpathf(p1aind-sind+1))*repmat(b',[1,J,bsize]) + ...
                (wpathf(p1aind-sind+1)*n(S-sind+1))*repmat(e(S-sind+1,:),...
                [bsize,1,bsize]) - repmat(reshape(b,[1,1,bsize]),...
                [bsize,J,1]) ;
         elseif sind == S
            c = (wpathf(p1aind-sind+1)*n(S-sind+1))*repmat(e(S-sind+1,:),...
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
         Vnewarray = uc + beta*EVprimenew.*cposind ;
         [Vnew,bprimeind] = max(Vnewarray,[],3) ;
         phiindt(S-sind+1,:,:,tind+p1aind-sind) = reshape(bprimeind',...
                                               [1,J,bsize,1]) ;
         phit(S-sind+1,:,:,tind+p1aind-sind) = reshape(b(bprimeind)',...
                                               [1,J,bsize,1]) ;
         Vinit = Vnew ;
      end
   end
   
   % Generate next period's distribution of wealth from time-t policy rules
   for sind = 1:S-1
      for eind = 1:J
         for bind = 1:bsize
            if sind == 1
               gammat(sind,eind,bind,tind+1) = (1/(S-1))*f(sind+1,eind)*...
                             sum((phiindt(sind,:,1,tind) == bind).*...
                             f(sind,:)) ;
            else
               gammat(sind,eind,bind,tind+1) = f(sind+1,eind)*...
                                   sum(sum((phiindt(sind,:,:,tind) ...
                                   == bind).*...
                                   gammat(sind-1,:,:,tind))) ;
            end
         end
      end
   end
   Kpath_AMF(1,tind+1) = ((S-1)/S)*sum(sum(sum(gammat(:,:,:,tind+1).*...
                  repmat(reshape(b,[1,1,bsize]),[S-1,J,1])))) ;
   
   [tind+1 Kpath_AMF(tind+1)]
end

gammat_AMF = gammat ;

runtime = toc(starttime) ; % end AMF1 timer

starttime2 = tic ;

Kpath_AMF2 = [K0 zeros(1,T-1) Kss*ones(1,S-2)] ;
gammat  = zeros(S-1,J,bsize,T) ;
gammat(:,:,:,1) = gamma0 ;
phiindt = zeros(S,J,bsize,T+S-2) ;
phit    = zeros(S,J,bsize,T+S-2) ;

[1 K0]
for tind = 1:T-1
   % Alternative forecast is constant growth rate between K_t and Kss
   lnKpathf = [(log(Kpath_AMF2(tind)):(1/(T-tind))*...
              log(Kss/Kpath_AMF2(tind)):log(Kss)) log(Kss)*ones(1,S-2)] ;
   Kpathf = exp(lnKpathf) ;
   Npathf = Nss*ones(1,T+S-tind-1) ;
   Ypathf = A*((Kpathf.^alpha).*(Npathf.^(1-alpha))) ;
   wpathf = (1-alpha)*(Ypathf./Npathf) ;
   rpathf = alpha*(Ypathf./Kpathf) ;
   
   % Solve for the lifetime savings policy rules of each type of person at
   % time t
   for p1aind = 1:S
      Vinit = zeros(bsize,J) ;
      for sind = 1:p1aind
         if sind < S
            c = (1+rpathf(p1aind-sind+1))*repmat(b',[1,J,bsize]) + ...
                (wpathf(p1aind-sind+1)*n(S-sind+1))*repmat(e(S-sind+1,:),...
                [bsize,1,bsize]) - repmat(reshape(b,[1,1,bsize]),...
                [bsize,J,1]) ;
         elseif sind == S
            c = (wpathf(p1aind-sind+1)*n(S-sind+1))*repmat(e(S-sind+1,:),...
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
         Vnewarray = uc + beta*EVprimenew.*cposind ;
         [Vnew,bprimeind] = max(Vnewarray,[],3) ;
         phiindt(S-sind+1,:,:,tind+p1aind-sind) = reshape(bprimeind',...
                                               [1,J,bsize,1]) ;
         phit(S-sind+1,:,:,tind+p1aind-sind) = reshape(b(bprimeind)',...
                                               [1,J,bsize,1]) ;
         Vinit = Vnew ;
      end
   end
   
   % Generate next period's distribution of wealth from time-t policy rules
   for sind = 1:S-1
      for eind = 1:J
         for bind = 1:bsize
            if sind == 1
               gammat(sind,eind,bind,tind+1) = (1/(S-1))*f(sind+1,eind)*...
                             sum((phiindt(sind,:,1,tind) == bind).*...
                             f(sind,:)) ;
            else
               gammat(sind,eind,bind,tind+1) = f(sind+1,eind)*...
                                   sum(sum((phiindt(sind,:,:,tind) ...
                                   == bind).*...
                                   gammat(sind-1,:,:,tind))) ;
            end
         end
      end
   end
   Kpath_AMF2(1,tind+1) = ((S-1)/S)*sum(sum(sum(gammat(:,:,:,tind+1).*...
                  repmat(reshape(b,[1,1,bsize]),[S-1,J,1])))) ;
   
   [tind+1 Kpath_AMF2(tind+1)]
end

gammat_AMF2 = gammat ;

runtime2 = toc(starttime2) ; % end AMF2 timer

runtimetot = toc(starttime) ; % end total timer

MPD_AMF = (1/T)*100*...
          (((Kpath_AMF(1:T) - Kpath_TPI(1:T))./Kpath_TPI(1:T))*ones(T,1)) ;

%-------------------------------------------------------------------------%
% Generate graph of the transition path of the aggregate capital stock
%-------------------------------------------------------------------------%
Kssvec = Kss*ones(1,T+10) ;
save OLGabilAMF9.mat ;

figure(1) % Figure with just the linear forecast
plot(1:T+10,Kpath_TPI(1:T+10),1:T+10,Kpath_AMF(1:T+10),'-.',1:T+10,Kssvec,'k--','LineWidth',2)
axis([0 T+11 0.9*min([min(Kpath_TPI) min(Kpath_AMF)]) ...
     1.1*max([max(Kpath_TPI) max(Kpath_AMF)])])
xlabel('time t','FontSize',15)
ylabel('aggregate capital K','FontSize',15)
%title('AMF transition path of aggregate capital K: S = 60','FontSize',20)
legend('TPI time path K_t','AMF time path K_t','steady-state K','Location','Best')

% figure(2) % Figure with both the linear and log linear forecasts
% plot(1:T+10,Kpath_TPI(1:T+10),1:T+10,Kpath_AMF(1:T+10),'-.',1:T+10,Kpath_AMF2(1:T+10),':',1:T+10,Kssvec,'k--','LineWidth',2)
% axis([0 T+11 0.9*min([min(Kpath_TPI) min(Kpath_AMF) min(Kpath_AMF2)]) ...
%      1.1*max([max(Kpath_TPI) max(Kpath_AMF) max(Kpath_AMF2)])])
% xlabel('time t','FontSize',15)
% ylabel('aggregate capital K','FontSize',15)
% %title('AMF transition path of aggregate capital K: S = 60','FontSize',20)
% legend('TPI time path K_t','AMF time path K_t','AMF2 time path K_t','steady-state K','Location','East')