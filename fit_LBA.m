function [fparms,minGsq,ret_cdf,new_cdf, err,ret_ploty, new_ploty]=fit_LBA(ret,nonret,errors)

%Fits five competing models to account for the differences between two
%latency distributions. The input are the two latency distributions, and a
%starting parameter vector. The latter contains, in order, mean rate 1,
%mean rate 2, sd rate 1, sd rate 2, threshold, shift1, shift2. IC contains the -2*log
%likelihood error term for the 4 different models (column 1), the AIC and
%BIC derived from these quantities (columns 2 and 3), and the associated
%weights (that add up to 1 - columns 4 and 5).

% guide start values
% [.005 .005 .001 .001 1 30 10]

if(min(ret)<1) %simple check for the timescale; if latency is specified in seconds, convert to ms
    ret=ret*1000;
    nonret=nonret*1000;
end
minT=min([min(ret) min(nonret)]);

% parameters are:

% target drift
% non-target drift
% ISR drift mult
% sd drift
% start range
% ISR boundary
% Tr

wval = 100000000;
for a=[.003 .005]
    for b=[1e-10 .001]
        for c=[.7 .9]
            for d=[.0005 .002]
                for e=[.2 .7]
                    for f=[.8 1.2];
                        for g=[50 150]
                            startp=[a b c d e f g] %#ok<NOPRT>
                            %[Gsq,ret_cdf,new_cdf,ploty] = LBA_ISR_cdf(parms,retSRT,newSRT,errors)
                            [tparms, fval, tflag] = fminsearchbnd(@LBA_ISR_cdf,startp,[0 0 0 eps 0 0 50],[1 1 Inf Inf Inf Inf minT],...
                                optimset('Display','none', 'TolX', .1),ret,nonret,errors, 0); %#ok<NASGU>
                            if fval<wval
                                wval=fval;
                                fparms=tparms %#ok<NOPRT>
                                minGsq=fval %#ok<NOPRT>
                            end
                        end
                    end
                end
            end
        end
    end
end

%startp=[.005 .0025 .8 .001 .2 .8 60] %#ok<NOPRT>
%[Gsq,ret_cdf,new_cdf,ploty] = LBA_ISR_cdf(parms,retSRT,newSRT,errors)
% [fparms, wval, tflag] = fminsearchbnd(@LBA_ISR_cdf,startp,[0 0 0 eps 0 0 50],[1 1 Inf Inf Inf Inf minT],...
%     optimset('Display','none', 'TolX',.1),ret,nonret,errors, 0); %#ok<NASGU>
% minGsq=wval

% Having fit the model, run it one last time with fitted parameters to get
% final predictions
[Gsq,ret_cdf,new_cdf,err, ploty] = LBA_ISR_cdf(fparms,ret,nonret,errors, 1);
ret_cdf = ret_cdf(:,1)';
new_cdf = new_cdf(:,1)';
ret_ploty = ploty(:,1);
new_ploty = ploty(:,2);