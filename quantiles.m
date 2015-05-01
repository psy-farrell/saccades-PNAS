function [cdf]=quantiles(data,p)
%computes a cumulative distribution function consisting of length(p)
%quantiles specified by p: 0~=p(1)<p(2)...p(m)~=1. Algorithm described in
%Heathcote, Brown, & Mewhort (2002).

N=length(data);
data=sort(data);
for i=1:length(p)
    if(i==1) %compute number of data points within quantile i
        cdf(i,2)=p(i)*N;
    else cdf(i,2)=(p(i)-p(i-1))*N;
    end
    Index=p(i)*N+0.5;
    Index_min=floor(Index);
    Index_plus=ceil(Index);
    cdf(i,1)=data(Index_min)+(data(Index_plus)-data(Index_min))*(Index-Index_min);
end
% cdf=[cdf;max(data) length(data)-sum(cdf(:,2))];