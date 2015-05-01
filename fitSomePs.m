global fixPs

fixPs = 0;

ESs = [1, 10, 13, 16, 20, 23, 26, 29, 32, 35, 4, 7];
NSs = [11, 14, 17, 2, 21, 24, 27, 30, 33, 36, 5, 8];
RSs = [12, 15, 18, 19, 22, 25, 28, 3, 31, 34, 6, 9];

E = 'R'; % set this to 'E', 'N', or 'R' to fit participants in 
        %'E'qual, 'N'ew, or 'R'eturn conditions 

tSs = eval([E 'Ss']);

for Ss=tSs(1)
    i=1;
    sind = find(tSs==Ss)
    
    load(['P7' E 'fil.mat']);
    tmat = eval(['p' num2str(Ss) 'secm']);
    
    % each matrix in the data file belongs to a single participant
    % the matrix columns are:
    % 1) session (1 to 3)
    % 2) event; a number from 1 to 5 indicating what happened:
    %       1: new target, correct response
    %       2: return target, correct response
    %       3: return target, error (saccade to new)
    %       4: new target, error (saccade to return)
    %       5: new targer, error (saccade to other new)
    % 3) saccadic reaction time, in milliseconds
    
    
    for cSess = 1:3
        sess = tmat(:,1);
        event = tmat(:,2);
        RTs = tmat(:,3);
        event = event(sess==cSess);
        RTs = RTs(sess==cSess);

        cSess
        
        if ~isempty(event)
            event=event(2:end);
            RTs=RTs(2:end);

            % this makes the complete errors structure required by the fitting routine
            errors = zeros(3,3);
            errors(:,1) = NaN;
            errors(1,2) = sum(event==3);
            errors(2,2) = sum(event==4);
            errors(3,2) = sum(event==5);
            errors(1,3) = sum(event==3)./sum(event==3| event==2);
            errors(2,3) = sum(event==4)./sum(event==1| event==4| event==5);
            errors(3,3) = sum(event==5)./sum(event==1| event==4| event==5);

            disp(errors);
            ret = RTs(event==2); % correct return latencies
            nonret = RTs(event==1); % correct new latencies

            % cut out outliers
            ret = ret(ret>80 & ret<1000);
            nonret = nonret(nonret>80 & nonret<1000);

            tic
            [fparms,minGsq,ret_cdf,new_cdf,err, ret_ploty, new_ploty]=fit_LBA(ret,nonret,errors);
            toc
            err
            P(sind).parms{i}=fparms;
            P(sind).Gsq(i) = minGsq;
            P(sind).ret_cdf(i,:)=ret_cdf(1:5);
            P(sind).new_cdf(i,:)=new_cdf(1:5);
            P(sind).moderr(i,:) = err;
            P(sind).daterr(i,:) = errors(:,3);
            P(sind).ret_ploty(i,:) = ret_ploty;
            P(sind).new_ploty(i,:) = new_ploty;
        else
            P(sind).parms{i}=repmat(NaN,1,7);
            P(sind).Gsq(i) = NaN;
            P(sind).ret_cdf(i,:)=repmat(NaN,1,5);
            P(sind).new_cdf(i,:)=repmat(NaN,1,5);
            P(sind).err(i,:) = [NaN NaN NaN];
            P(sind).ret_ploty(i,:) = repmat(NaN,1,1000);
            P(sind).new_ploty(i,:) = repmat(NaN,1,1000);
        end
        i=i+1; % i indexes session
    end
end

savname = [E 'P2015b.mat'];
save(savname, 'P');