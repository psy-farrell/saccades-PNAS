function [Gsq,ret_cdf,new_cdf,err,ploty]=LBA_ISR_cdf(parms,retSRT,newSRT,errors,getploty)

global fixPs

%Fit a 7 parameter LBA model to defective 2nd saccade latency distributions. The
%parameters are, in order: 1) v_t: mean rate of new target accumulator; 2)
%v_n: mean rate of new non-target accumulator; 3) alpha: mean rate
%multiplier: may reduce the rate of the accumulator associated with the
%return location; 4) sd: shared standard deviation of the accumulation
%rates; 5) A: upper boundary of the uniform distribution of starting
%points; 6) beta: multiplier of the threshold to enable a different
%threshold for return locations; 7) constant non-decisional delay
%(intercept). The threshold for a non-return location is fixed to 1.
% 'errors' is a 1x2 vector with the number of errors for return and
% non-return conditions respectively.

%assign parameters
v_t=parms(1);
v_n=parms(2);
alpha=parms(3);
sd=parms(4);
A=parms(5);
beta=parms(6);
Tr=parms(7);
theta=1;

%determine quantile distributions
q=[.1 .3 .5 .7 .9]; %same quantiles as Ratcliff & Smith (2004)
ret_cdf=quantiles(retSRT,q);
new_cdf=quantiles(newSRT,q);
%add the final quantile
ret_cdf(end+1,:)=[inf (1-q(end))*length(retSRT)];
new_cdf(end+1,:)=[inf (1-q(end))*length(newSRT)];
ret_cdf(end+1,:)=errors(1,1:2);
new_cdf((end+1):(end+2),:)=[errors(2,1:2); errors(3,1:2)];

% get proportions
ret_cdf(:,3)=ret_cdf(:,2)/nansum(ret_cdf(:,2)); %observed proportions in each bin
new_cdf(:,3)=new_cdf(:,2)/nansum(new_cdf(:,2));
maxSRT=1000; %need this to compute asymptotic performance level

%For a given set of parameters, compute the model predicted CDFs
%Return saccades
model_ret_cdf=LBA_cdf(ret_cdf(1:length(q),1)-Tr,A,[beta*theta theta theta],[alpha*v_t v_n v_n],sd);
model_ret_cdf(end+1)=LBA_cdf(maxSRT-Tr,A,[beta*theta theta theta],[alpha*v_t v_n v_n],sd);
%Non-return saccades
model_new_cdf=LBA_cdf(new_cdf(1:length(q),1)-Tr,A,[theta theta beta*theta],[v_t v_n alpha*v_n],sd);
model_new_cdf(end+1)=LBA_cdf(maxSRT,A,[theta theta beta*theta],[v_t v_n alpha*v_n],sd);

% calculate model predicted probabilities
model_ret_p = diff([0; model_ret_cdf]);
model_ret_p(end+1) = 1-model_ret_cdf(end);

%this prevents numerical problems
model_ret_p(model_ret_p < realmin) = realmin;
model_ret_p = model_ret_p./sum(model_ret_p);

% note for the above, the alternative would be to get the error proportions
% explicitly. This is just the quicker way of doing it, as we can't
% distinguish between the two other types of errors
model_new_p = diff([0; model_new_cdf]);
if fixPs
    endP = new_cdf((end-1):end,2);
    endP = endP./sum(endP);
    model_new_p((end+1):(end+2)) = (1-model_new_cdf(end)).*endP;
else
    model_new_p(end+1) = LBA_cdf(maxSRT,A,[beta*theta theta theta],[alpha*v_n v_n v_t],sd);
    model_new_p(end+1) = LBA_cdf(maxSRT,A,[theta beta*theta theta],[v_n alpha*v_n v_t],sd);
end
model_new_p(model_new_p < realmin) = realmin;
model_new_p = model_new_p./sum(model_new_p);

%Everything is in place to compute the G^2 error statistic...
if (ret_cdf(end,2)>0) %only evaluate the "error" cells if there are any errors
    Gsq_ret=2*nansum(ret_cdf(:,2).*log(ret_cdf(:,3)./model_ret_p));
else
    rind = ret_cdf(:,2)>0;
    Gsq_ret=2*nansum(ret_cdf(rind,2).*log(ret_cdf(rind,3)./model_ret_p(rind)));
end

if (all(new_cdf((end-1):end,2)>0))
    Gsq_new=2*nansum(new_cdf(:,2).*log(new_cdf(:,3)./model_new_p));
else
    nind = new_cdf(:,2)>0;
    Gsq_new=2*nansum(new_cdf(nind,2).*log(new_cdf(nind,3)./model_new_p(nind)));
end
Gsq=Gsq_ret+Gsq_new;

if getploty
    if (ret_cdf(end,2)>0) %only evaluate the "error" cells if there are any errors
        2*(ret_cdf(:,2).*log(ret_cdf(:,3)./model_ret_p));
    else
        rind = ret_cdf(:,2)>0;
        2*(ret_cdf(rind,2).*log(ret_cdf(rind,3)./model_ret_p(rind)));
    end
    
    if (all(new_cdf((end-1):end,2)>0))
        2*(new_cdf(:,2).*log(new_cdf(:,3)./model_new_p));
    else
        nind = new_cdf(:,2)>0;
        2*(new_cdf(nind,2).*log(new_cdf(nind,3)./model_new_p(nind)));
    end
end

%%
%The following code may be used to generate a goodness-of-fit plot: it
%shows the predicted cdfs for correct saccades to return and non-return
%targets (errors are so rare that there is not much point in plotting
%them for individuals).
err = [model_ret_p(end); model_new_p((end-1):end)];

%Not a good idea to trigger this when fitting!
if ~getploty
    ploty=0;
else
    plott=1:1000; %to be used for g-o-f plots
    ploty=zeros(length(plott),2); %model densities for each of the 5 possible outcomes
    ploty(:,1)=LBA_pdf(plott'-Tr,A,[beta*theta theta theta],[alpha*v_t v_n v_n],sd);
    ploty(:,2)=LBA_pdf(plott'-Tr,A,[theta theta beta*theta],[v_t v_n alpha*v_n],sd);
    ploty(:,3)=LBA_pdf(plott'-Tr,A,[theta theta beta*theta],[v_n v_n alpha*v_t],sd);
    ploty(:,4)=LBA_pdf(plott'-Tr,A,[beta*theta theta theta],[alpha*v_n v_n v_t],sd);
    ploty(:,5)=LBA_pdf(plott'-Tr,A,[theta beta*theta theta],[v_n alpha*v_n v_t],sd);
    ploty=cumsum(ploty);
end

% uncomment this if you actually want the plot right now
% figure;hold on;
% plot(ret_cdf(1:length(q),1),cumsum(ret_cdf(1:length(q),3)),'ks','MarkerFaceColor','k','linewidth',2); %observed correct, return location
% plot(new_cdf(1:length(q),1),cumsum(new_cdf(1:length(q),3)),'kv','MarkerFaceColor','w','linewidth',2); %observed correct, non-return location
% plot(plott,ploty(:,1),'k-','linewidth',2); %predicted correct, return location
% plot(plott,ploty(:,2),'k--','linewidth',2); %predicted correct,
% non-return location
