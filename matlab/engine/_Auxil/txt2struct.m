function [spec] = txt2struct(filename)
%
% Convert .txt files that contain spectrum to the format recognized
% by TODCOR and return the information. 
% 
% The spectrum must appear as in the following example:
% Star: HD89269A            
% JD: 2454538.41948964
% HCV: -13.70837964
% Velocity: -0031.1999
% number_of_orders: 39
% wl(1,1) sp(1,1) wl(2,1) sp(2,1) ...
% wl(1,2) sp(1,2) wl(2,2) sp(2,2) ...
% wl(1,3) sp(1,3) wl(2,3) sp(2,3) ...
%   :        :      :       :
% where the first index is the order and the second is the running index
% within the order.
% origin - 'string' that contains the name of the spectrograph
%
% output:
% wv,sp: cell arrays of wavelength and spectrum that correspond to the orders vector.
% name: the name to print in the file. 
% vel: initial guess for the velocity (if there is such)
% hcv: Helio-Centric Velocity
% jd: Julian Date
% orders: the orders of the spectrum
%
% input:
% filename: full string of the .txt file
% 
% Created: 18.11.09 LT
% Modified: 16.12.09 LT

fid=fopen(filename,'r');

% retrieving the name from the .txt file:
lin1=textscan(fid,'%s%s',1);
name1=lin1{2};
name=name1{1};

% retrieving the JD from the .txt file:
fclose(fid); % the first line may contain more information than the name
fid=fopen(filename,'r');
lin2=textscan(fid,'%s%f64',1,'Headerlines',1);
jd=lin2{2};

% retrieving the HCV from the .txt file:
lin3=textscan(fid,'%s%f64',1);
hcv=lin3{2};

% retrieving the Velocity (if there is such) from the .txt file:
lin4=textscan(fid,'%s%f64',1);
vel=lin4{2};

% retrieving the number of orders from the .txt file:
lin5=textscan(fid,'%s%d',1);
num_orders=lin5{2};
orders = 1:num_orders;

fclose(fid);
fid=fopen(filename,'r');
% Reading the spectrum from the text file:
M = dlmread(filename,'\t',5,0);
pixnum = size(M,1);
wv = zeros(pixnum,num_orders);
sp = wv;
for i=1:num_orders
   wv(:,i) = M(:,2*i-1); 
   sp(:,i) = M(:,2*i);
end
fclose(fid);

spec = struct('wv',wv,'sp',sp,'bjd',jd,'bcv',hcv,'vel',vel,'name',name,'orderN',num_orders);