function [M all_loglike]=MetroHast(X0, mu, cov,NumSamples, flag)
%Metropolis-Hasting sampler function for Gaussian distribution observations.
%Outputs: 
%   M - if flag=0, result modle matrix (by each column); flag=1, file name of output M matrix.
%   all_log - flag=0,log of likelihood vector; flag=1, filename of output
%   loglikelihood vector.
%
%Inputs:
%   X0 -  Starting sampling point/starting model;
%   mu - mean of the Gaussian distributed model;
%   cov - model covariance matrix;
%   NumSamples - number of samples.
%   flag - 0, return values to matrix M and all_loglike; 1, output result model matrix to
%       file M.txt & loglikelihood vector to LogLikelihood.txt.
%
%  Xiaotao Yang   Nov 10, 2012 based on Kaj Johnson's template.

X=X0;

%generate files to store results
%first write files (create file)
if flag == 1
    
    M_fname='M.txt'; all_loglike_fname='LogLikelihood.txt';
    
    fid = fopen(M_fname,'w'); fclose(fid);  %store model parameters
%    fid = fopen('DHAT.txt','w'); fclose(fid);   %store fit to data
    fid = fopen(all_loglike_fname,'w'); fclose(fid);    %store log-likelihood
%now open files to append
    fidM = fopen(M_fname,'a');
%fidD = fopen('DHAT.txt','a');
    fidL = fopen(all_loglike_fname,'a');
    
elseif flag==0
    
else
    error('Worng flag value, which should be either 1 or 0.')
end

% Initial value of loglike.
loglike=-0.5*(X - mu)'*inv(cov)*(X - mu);

%intialize radnom number generator
rand('state', sum(100*clock)); %reset random number generator

%define 'previous values'
Xprev=X;
loglikeprev=loglike;

if flag==0
    M=zeros(length(X),NumSamples);
    all_loglike=zeros(1,NumSamples);
end

for loop=1:NumSamples

    %take a random step in model space    
    X = Xprev + rand(size(Xprev))*2 - 1;

    %compute log-likelihood
    loglike=-0.5*(X - mu)'*inv(cov)*(X - mu);
    
    %use Metropolis rule to decide whether or not to accept model
    rat=exp(loglike-loglikeprev);
    if rat>1
       accept=1;
    else
       r=rand;
       if r<rat
          accept=1;
       else
          accept=0;
       end
    end

    
    if accept==1;
        loglikeprev = loglike;
        Xprev = X;
    else
        loglike = loglikeprev;
        X = Xprev;
    end
    
    if flag==1
        %save every kth sample
        k=10;
        if mod(loop,k)==0
            fprintf(fidM,'%6.8f\t',X'); fprintf(fidM,'\n',' ');
            fprintf(fidL,'%6.8f\t',loglike); fprintf(fidL,'\n',' ');
        end
    else  %flag ==0
        M(:,loop)=X;
        all_loglike(loop)=loglike;
    end
    
end %loop

end