 function [vels, evels] = MatrixRV(V,IterN,C,SigThr,I0)  
 % ========================================================================== 
% FUNCTION: [vels, evels] = MatrixRV(V,IterN,C,SigThr,I0)  
% ========================================================================== 
% 
% This function recieves a matrix of velocities, derived from multi-order
% spectra, and calculated the velocity and its error.  See Engel et al. (2017) for a detailed description of the algorithm.
% Written by Sahar Shahaf, April 2018, as a part of UNICOR.
%
% Last update: *** By: ***
%
% --------------------------------------------------------------------------
%
% INPUT:
% V        - A matrix of velocities [ # observation , # order ]
% IterN    - An integer. The required number of iterations. Typically ~10
%            iterations are sufficient.
% C        - A strictly positive number. Minimal number of sigma to define
%            a systematic shift as "signigicant" (Default is 3).
%            ****
% SigThr   - The scater threshold, above which orders are rejected. 
%            Must be provided in the same units of V (Default is 1).
% I0       - Logical array ([1 , # orders]). Initial subset of orders to 
%            use from the matrix. Orders that are rejected by I0 will not
%            participate in the calculations. These should be orders that
%            are, for example, heavily polluted by tellurics or do not
%            contin significant spectral information of the target star.
% OUTPUT:
% vels     - A vector, containing the derived velocities.
% evels    - A vector, containing the derived errors.
%--------------------------------------------------------------------------
%
 

% 1) Initialize.
%  ------------- 

% Check input:
if nargin == 2
    C      = 5;
    SigThr = 1;  %Same units as V
    I0     = true(1,size(V,2)); 
elseif nargin == 4
    I0     = true(1,size(V,2));
elseif nargin < 5
    warning('Check Input.');
    vels   = NaN; 
    evels  = NaN; 
    return;
end

if C<=0 || SigThr<=0 
    warning('Check Input.');
    vels   = NaN; 
    evels  = NaN; 
    return;
end

% Set variables:
V          = V(:,I0);
Nobs       = size(V,1);
Nord       = size(V,2); 


% 2) Zeroth Iteration.  
% --------------------
% Let us denote by v(i,j) the RV obtained per order j for each exposure i. 
% We assume that 
%
%                     v(i,j) = v(i) + epsilon(i,j),
%
% where epsilon(i,j) is the error of the j-th order at the i-th exposure,
% and v(i) is the true stellar velocity.
%
% An initial estimate of the stellar velocity and its uncertainty, v0 & dv0, is
% obtained for each exposure by taking the median and median absolute deviation of this subset
% respectively. Hence,
v0        = trimmean(V,10,2,'weighted');
 
% The zero-iteration error matrix is
epsilon0  =  V - repmat(v0,1,Nord);
 
% With no systematic errors, we expect each column j  in the epsilon(i,j) matrix, which corresponds to the j-th
% order, to be scattered evenly around zero. An initial estimate of the systematic shift of an order
delta0    =  trimmean(epsilon0,10,1,'weighted');

% The scatter of order j, at the current iteration, is derived from the median absolute deviation of
% the corrected error matrix,
% sigma0    = 1.4826*mad(epsilon0 - repmat(delta0,Nobs,1),1,1);
sigma0    = 1.4826*mad(epsilon0,1,1);
d_delta0  = sigma0/sqrt(Nobs);

% We only correct significant systematics:
significant          = abs(delta0) > C * d_delta0;
delta0(~significant) = 0;


% 3) Iterate...
% --------------

% Initialize before loop:
delta_k = delta0;
I_k     = sigma0 <= SigThr;
 
for k = 1 :IterN
    % The process is repeated by using a new velocity matrix, V_k, from which
    % the previous estimate of the systematics was removed. 
    clear('V_k');
    V_k        = V(:,I_k) - repmat(delta_k(I_k),Nobs,1);
    Nord_k     = sum(real(I_k));
    % The K-thestimate of the stellar velocity, v_k & dv0,  
    v_k        = trimmean(V_k,10,2,'weighted');

    % The zero-iteration error matrix is
    epsilon_k  = V_k - repmat(v_k,1,Nord_k);
 
    % The k-th estimate of the systematic shift of an order
    delta_k    =  trimmean(epsilon_k,10,1,'weighted');

    % The scatter of order j, at the current iteration
    %sigma_k    = 1.4826*mad(epsilon_k - repmat(delta_k,Nobs,1),1,1);
    sigma_k    = 1.4826*mad(epsilon_k,1,1);
    d_delta_k  = sigma_k/sqrt(Nobs);

    % We only correct significant systematics:
    significant           = abs(delta_k) > C * d_delta_k;
    delta_k(~significant) = 0;
    I_k                   = sigma_k <= SigThr;
      
end

 
% 4) Mangage output:
% ------------------

vels  = v_k;
evels = 1.4826*mad(V_k,1,2);
 
% QA PLOT
% Plot initial scatter and shift
% ------------------------------
figure;
subplot(2,1,1)
errorbar(1:1:Nord,delta0,d_delta0,'.k','markerfacecolor','k'); grid on; hold on;
plot([1,Nord],[0 0],'--','color',[0.5 0.5 0.5],'linewidth',2); hold off;
xlabel('\#ORD','fontsize',12,'interpreter','latex');
ylabel('$\delta_0$' ,'fontsize',16,'interpreter','latex');
xlim([0,Nord+1]);
subplot(2,1,2)
plot(1:1:Nord,sigma0,'ok','markerfacecolor','k'); grid on; hold on;
plot([1,Nord],[SigThr,SigThr],'--r','linewidth',2); hold off;
xlabel('\#ORD','fontsize',12,'interpreter','latex');
ylabel('$\sigma_0$ ','fontsize',16,'interpreter','latex');
xlim([0,Nord+1]);
% plot the convergence.
% --------------------
% subplot(3,1,3)
% semilogy(1:1:IterN,vDiff(1:end),'-ok','markerfacecolor','k'); grid on;
% xlabel('\#ITERN','fontsize',12,'interpreter','latex');
% ylabel('$<|v_k-v_{k+1}|> $','fontsize',12,'interpreter','latex');
% xlim([1,IterN+1]);
 end  




