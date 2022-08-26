function matrixOut = smooth2(matrixIn, Nr, Nc)

% SMOOTH2.M: Smooths matrix data.

%			MATRIXOUT=SMOOTH2(MATRIXIN,Nr,Nc) smooths the data in MATRIXIN 
%           using a running mean over N successive points. N=5 by default.
%           At the ends of the series skewed or one-sided means are used.  
%
%           Inputs: matrixIn - original matrix
%                   Nr - number of points used to smooth rows
%                   Nc - number of points to smooth columns
%           Outputs:matrixOut - smoothed version of original matrix
%
%           Remark: By default, if Nc is omitted, Nc = Nr.
%
%           Written by Yun Xu, Feb 2014
%


%Initial error statements and definitions
if nargin<2, N(1)=5; N(2)=5; end

N(1) = (Nr-1)/2; 
if nargin<3 
    N(2) = N(1); 
else
    N(2) = (Nc-1)/2;
end

if length(N(1))~=1, error('Nr must be a scalar!'), end
if length(N(2))~=1, error('Nc must be a scalar!'), end

[row,col]=size(matrixIn);
temp=nan(row+Nr-1,col+Nc-1);

temp(N(1)+1:N(1)+row,N(2)+1:N(2)+col)=matrixIn;

matrixOut=nan(row,col);

for i=1:row
 matrixOut(i,:)=nanmean(temp(i:i+Nr-1,N(2)+1:N(2)+col),1);
end
 aaa=find(isnan(matrixIn));matrixOut(aaa)=nan;
 temp(N(1)+1:N(1)+row,N(2)+1:N(2)+col)=matrixOut;
for i=1:col
 matrixOut(:,i)=nanmean(temp(N(1)+1:N(1)+row,i:i+Nc-1),2);
end
 aaa=find(isnan(matrixIn));matrixOut(aaa)=nan;

return


