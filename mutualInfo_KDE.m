function [I] = mutualInfo_KDE(vec1,vec2,nbins)
% Computes mutual information between 2 signals using kernel density
% estimation to estimate the probability density

% Written by Benedikt Barthel


vec1 = squeeze(vec1); vec2 = squeeze(vec2);


n = length(vec1);
Q = [vec1',vec2'];

% Joint Probability Density function

%%5
if nbins == 0
    r1 = iqr(vec1); r2 = iqr(vec2);
    w1 = 2*r1/(n^0.333); w2 = 2*r2/(n^0.333);
    L1 = max(vec1)-min(vec1); L2 = max(vec2) - min(vec2);
    nbins1 = round(L1/w1); nbins2 = round(L2/w2);
    pts1 = linspace(min(vec1),max(vec1),nbins1);
    pts2 = linspace(min(vec2),max(vec2),nbins2);

    % probabilty density function of vec1
    [f1,x1,bw1] = ksdensity(vec1,pts1);
    np1 = length(x1);
    % probabilty density function of vec2
    [f2,x2,bw2] = ksdensity(vec2,pts2); %
    np2 = length(x2);
%     disp(['Automatic Number of Points (FD rule):',num2str(nbins1),', ',num2str(nbins2)])
else
    pts1 = linspace(min(vec1),max(vec1),nbins);
    pts2 = linspace(min(vec2),max(vec2),nbins+1);
    % probabilty density function of vec1
    [f1,x1,bw1] = ksdensity(vec1,pts1);
    np1 = length(x1);
    % probabilty density function of vec2
    [f2,x2,bw2] = ksdensity(vec2,pts2); %
    np2 = length(x2);
%     disp(['User Specified Number of Points: ', num2str(nbins)])
end


%%%
% Joint PDF Domain
[x,y] = meshgrid(x1,x2);
x = reshape(x,[],1);y = reshape(y,[],1);
pts = [x,y];
vec12 = [vec1;vec2]'; % 2D variable
[f12,x12,bw12] = ksdensity(vec12,pts);% joint probability density
f12 = reshape(f12,[np2,np1])'; % f_xy(x,y)

% % Marginal Probability
% f1 = trapz(x2,f12,2)';
% f2 = trapz(x1,f12,1)';


% Normalize
P12 = trapz(x1,f12,1);P12 = trapz(x2,P12,2);f12 = f12/P12;

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
