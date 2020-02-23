function [ BIS,eBIS ] = CalcBIS( CCF,  I0)
%CALCBIS Summary of this function goes here
%   Detailed explanation goes here
if nargin == 1
    Nord = min(size(CCF{1}.corr));
    I0   = 1:1:Nord;
else
    Nord = length(I0);
end
for i = 1:length(CCF)
    BISvec = zeros(1,Nord);
    
    for j = 1:Nord
    c            = CCF{i}.corr(:,I0(j));
    v            = CCF{i}.vels(:,I0(j)); 
    
    [cmax,cmaxI] = max(c);    
    c            = c/cmax;
    
    
    vR           = v(cmaxI+1:end);
    vL           = v(1:cmaxI-1);
    cR           = c(cmaxI+1:end);
    cL           = c(1:cmaxI-1);    
    
    cL(vL<vL(1)-40) = [];
    vL(vL<vL(1)-40) = [];
    cR(vR>vR(1)+40) = [];
    vR(vR>vR(1)+40) = [];
    
    BIS_low      = (mean(vR(cR > 0.25 & cR < 0.4)) + mean(vL(cL > 0.25 & cL < 0.4)))/2;
    BIS_high     = (mean(vR(cR > 0.75 & cR < 0.9)) + mean(vL(cL > 0.75 & cL < 0.9)))/2;
    
    BISvec(j)    = (BIS_low - BIS_high);
    end
  BIS(i)  = median(BISvec);
  eBIS(i) = mad(BISvec,1)*1.48/sqrt(Nord);
end

 

end

