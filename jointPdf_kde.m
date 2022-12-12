function [prob] = jointPdf_kde(vec1,vec2,nbins)
% Computes mutual information between 2 signals using kernel density
% estimation to estimate the probability density

% Written by Benedikt Barthel

vec1 = squeeze(vec1); vec2 = squeeze(vec2);
pts1 = linspace(min(vec1),max(vec1),nbins);
pts2 = linspace(min(vec2),max(vec2),nbins+2);

if nbins ==0
    % probabilty density function of vec1
    [f1,x1] = ksdensity(vec1);
    np1 = length(x1);
    % probabilty density function of vec2
    [f2,x2] = ksdensity(vec2); %
    np2 = length(x2);
else
    % probabilty density function of vec1
    [f1,x1] = ksdensity(vec1,pts1);
    np1 = length(x1);
    % probabilty density function of vec2
    [f2,x2] = ksdensity(vec2,pts2); %
    np2 = length(x2);
end
% Joint Probability
[x,y] = meshgrid(x1,x2);
x = reshape(x,[],1);
y = reshape(y,[],1);
pts = [x,y];


vec12 = [vec1;vec2]'; % 2D variable
[f12,x12] = ksdensity(vec12,pts);% joint probability density
f12 = reshape(f12,[np2,np1])'; % f_xy(x,y)

% Normalize
P12 = trapz(x1,f12,1);P12 = trapz(x2,P12,2);
f12 = f12/P12;

% assign outputs
prob.f12 = f12;
prob.f1 = f1;
prob.f2 = f2;
prob.x1 = x1;
prob.x2 = x2;
end

