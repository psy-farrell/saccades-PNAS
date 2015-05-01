function n1p=LBA_pdf(t,A,b,v,sdv)

% Returns the predicted probability density of latency t under the Linear
% Ballistic Accumulator model of Brown and Heathcote. The input parameters
% A, b, v, sdv are, respectively, upper bounderary of uniform starting
% point variability, response threshold, mean drift rate(s), and the shared
% standard deviation of the drift rates.
% The pdf is the defective pdf; that is, it returns the probability density
% of accumulator 1 winning the race at time t.

N=length(v); % number of accumulators
if(length(A)<N)
    A=repmat(A,1,N);
end

if(N>2)
    tmp=zeros(length(t),N-1);
    for(i=2:N) %probability of accumulator i having reached threshold by time t
        tmp(:,i-1)=fptcdf(t,A(i),b(i),v(i),sdv);
    end
    G=prod(1-tmp'); %joint probability of all accumulators bar the first one NOT having reached threshold by time t
else
    G=1-fptcdf(t,A(2),b(2),v(2),sdv);
end
n1p=G'.*fptpdf(t,A(1),b(1),v(1),sdv); %equation 3 of Brown & Heathcote


function nicdf=fptcdf(t,A,b,v,sdv)
%computes equation 1 from Brown & Heathcote
ts=t*sdv;
tv=t*v;
bmintv=b-tv;
bminAmintv=bmintv-A;
x1=bmintv./ts;
x2=bminAmintv./ts;
tmp1=ts.*(normpdf(x2,0,1)-normpdf(x1,0,1));
tmp2=(bminAmintv.*normcdf(x2,0,1))-(bmintv.*normcdf(x1,0,1));
nicdf=1+(tmp1+tmp2)/A;


function nipdf=fptpdf(t,A,b,v,sdv)
%computes equation 2 of Brown & Heathcote
ts=t*sdv;
tv=t*v;
bmintv=b-tv;
bminAmintv=bmintv-A;
x1=bmintv./ts;
x2=bminAmintv./ts;
nipdf=(v*(normcdf(x1,0,1)-normcdf(x2,0,1))+sdv*(normpdf(x2,0,1)-normpdf(x1,0,1)))/A;