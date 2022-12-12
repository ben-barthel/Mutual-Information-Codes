function [I] = mutualInfo_MC(vec1,vec2,nbins)
% Computes mutual information between 2 signals using kernel density
% estimation to estimate the probability density

% Written by Benedikt Barthel


vec1 = squeeze(vec1); vec2 = squeeze(vec2);
n = length(vec1);
Q = [vec1',vec2'];
% Joint Probability Density function
if nbins ==0 
    
r1 = iqr(vec1); r2 = iqr(vec2);
w1 = 2*r1/(n^0.333); w2 = 2*r2/(n^0.333); 
L1 = max(vec1)-min(vec1); L2 = max(vec2) - min(vec2);
nbins1 = round(L1/w1); nbins2 = round(L2/w2);
% disp(['Optimal Number of Bins:',num2str(nbins1),', ',num2str(nbins2)])
[f12,x12] = hist3(Q,'Nbins',[nbins1,nbins2]);
else
% disp(['User Specified Number of Bins:',num2str(nbins),', ',num2str(nbins+1)])
[f12,x12] = hist3(Q,'Nbins',[nbins,nbins+1]);
end

x1 = x12{1};
x2 = x12{2};
% Normalize
P12 = trapz(x1,f12,1);P12 = trapz(x2,P12,2);f12 = f12/P12;
% Marginal Probability
f1 = trapz(x2,f12,2)';
f2 = trapz(x1,f12,1)';


% figure
% surf(x1,x2,f12')

% Mutual Information Integration
% Joint Entropy Calculation
H1 = -f1.*log(f1); H1(isnan(H1)) = 0;H1(isinf(H1)) = 0;
H1 = trapz(x1,H1);

H2 = -f2.*log(f2); H2(isnan(H2)) = 0;H2(isinf(H2)) = 0;
H2 = trapz(x2,H2);

H12 = -f12.*log(f12);

H12(isnan(H12)) = 0;
H12(isinf(H12)) = 0;
H12(f12<10^-9) = 0;
H12 = trapz(x1,H12,1);H12 = trapz(x2,H12,2);

% compute Mutual information
I = H1 + H2 - H12;



end
