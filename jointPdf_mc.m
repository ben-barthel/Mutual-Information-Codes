function [prob] = jointPdf_mc(vec1,vec2,nbins)
Q = [vec1',vec2'];
[f12,x12] = hist3(Q,'Nbins',[nbins,nbins+1]);
x1 = x12{1};x2 = x12{2};
% Normalize
P12 = trapz(x1,f12,1);P12 = trapz(x2,P12,2);f12 = f12/P12;
% Marginal Probability
f1 = trapz(x2,f12,2)';f2 = trapz(x1,f12,1)';

prob.f1=f1;
prob.f2=f2;
prob.f12=f12;
prob.x1=x1;
prob.x2=x2;
end