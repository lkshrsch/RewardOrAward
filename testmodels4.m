% Causal inference analysis to answer the question: 
%
% Do citations "cause" funding (H1), or does funding "cause" citations (H2)? 
%
% Based on NIH data originally released here:
%
% Lauer et al., bioRxiv, May. 26, 2017; doi: http://dx.doi.org/10.1101/142554 
%
% get latest version of this code at http://parralab.org/reward-or-award/showdata.m
% get corresponding data file at     http://parralab.org/reward-or-award/RCR_table.mat
%
% (c) June 7, 2017, Lucas C. Parra, parra@ccny.cuny.edu

clear all; clf; 

R_grant_recipients_only = 1; % flag to load corresponding data table

if R_grant_recipients_only
    load('Lauer_2017_data_R.mat','product_R','funding_R'); 
    %rename data so it can be stored all the way through
    productData = product_R;
    fundingData = funding_R;
else
    load('Lauer_2017_data.mat','product','funding');
    %rename data so it can be stored all the way through
    productData = product;
    fundingData = funding;
end

% set some parameters for shuffle
Nrand=1000; % number of shuffles when testing hypothesis that conditionals are the same except for mean and std
showconditionals = 0; % to check if the conditionals looks like properly adjusted accross bins

summaryTable = (zeros(35,4));
l=1; % index for saving the results in summaryTable



for p = 1:7
    for f = 1:5

        p
        f
% pick a measure of productivity and funding (see converttable.m). 
product = productData(p);
funding = fundingData(f);

% filter for low values as those were hard-coded for zero publications/funidng
indx=find(funding.data<=funding.min); funding.data(indx)=[]; product.data(indx)=[]; 
indx=find(product.data<=product.min); funding.data(indx)=[]; product.data(indx)=[]; 

figure('position', [0, 0, 1800, 1000])
set(gcf,'Visible', 'off'); 
subplot(2,1,1);
loglog(funding.data,product.data,'.','markersize',4);
grid on; axis tight; 
ylabel(product.name)
xlabel(funding.name)
legend(['PI (N=' num2str(length(product.data)) ')'],'Location','southeast')
title ('Research Productivity versus Funding')

for i=1:2
    
    switch i % selcect H1 or H2
        case 1, X=(product.data); Y=(funding.data); xstr=product.acronym; ystr=funding.acronym; 
        case 2, X=(funding.data); Y=(product.data); xstr=funding.acronym; ystr=product.acronym; 
    end
    
    % select bins edges with equal number of points per bin
    edge=prctile(X,0:4:99);
    
    % compute conditional statistics in each bin of conditioning variable
    for k=length(edge)-1:-1:1
        indx = find(edge(k)<=X & X<edge(k+1));
        count(k,1) = length(indx);          % points per bin
        cond(k,1)  = (edge(k)+edge(k+1))/2; % center of conditioning bin
        mu(k,1)    = mean(Y(indx));         % conditional mean of dependent variable
        sigma(k,1) = std(Y(indx));          % conditional std of dependent variable
        data{k}    = Y(indx);               % keep for shuffling bins later
    end
    % remove bins that did not have enough data to computer reliable stats
    k=find(count>10);mu=mu(k);sigma=sigma(k);count=count(k);cond=cond(k);data=data(k); 
    
    % number of conditioning bins that had some data in them
    K = length(mu);
    
    % some displays
    subplot(4,4,5+i*4)
    scatter(cond,mu,[], cond,'filled'); 
    xlabel(xstr); ylabel(['H' num2str(i) ': ' xstr '->' ystr 10 'mean ' ystr]); axis tight
    
    subplot(4,4,6+i*4)
    scatter(cond,sigma,[],cond,'filled'); 
    xlabel(xstr); ylabel(['std ' ystr]); axis tight
    
    subplot(4,4,7+i*4)
    scatter(mu,sigma,[], cond,'filled'); 
    xlabel(['mean ' ystr]); ylabel(['std ' ystr]); axis tight

    % estimate the std from the mean assuming linear model
    sest=polyval(polyfit(mu,sigma,1),mu);    
    hold on; plot(mu,sest); hold off;
    rsquare(i) = 1-sum((sest-sigma).^2)/sum((sigma-mean(sigma)).^2);
    ax=axis; t=text(ax(1)+0.5*(ax(2)-ax(1)),ax(3)+0.15*(ax(4)-ax(3)),['R^2=' num2str(rsquare(i),2)]);
    set(t,'BackgroundColor', [0.95 0.95 0.95])
  
    drawnow

    % Prepare to sample by shuffling bins. The model assumes conditionals are the
    % same across bins except for mean and std determined from linear model
    
    % convert to lognormal and adjust mean and estimated std of each bin assuming lognormal
    for k=1:K, data{k} = log(data{k}); end
    sest = sqrt(log(sest.^2./mu.^2+1));     % Using sest from linear regression to calculate new sest (lognormal)
    mu = log(mu)-(sest.^2/2);
    
    % generate standartized version of the data
    for k=1:K, zdata{k} = (data{k}-mean(data{k}))/std(data{k}); end

    % storing all information in stats, meanstats and zscore. But iteration
    % over k bins changes k depending on which hypothesis we are testing 
    % i = 1 has different k bins than i = 2. Set stats to zero after an
    % hypothesis was tested.
    stats = [0,0];
    
    % do shuffle stats
    for n=Nrand:-1:1  % n random samples (something around 10000?)
        
        for k=1:K   % for each iteration, go through all bins. Compare selected bin k against randomly selected one "sample1". Compare two random bin with each other too "sample1" vs "sample2" (why?)
            
            % sample byst shuffling bins of conditioning variable and adjusting mean and estimated std
            krand = [1:k-1 k+1:K]; krand=krand(randperm(K-1,2)); % 2 random bins excluding kth bin
            sample1 = zdata{krand(1)}*sest(k)+mu(k); % de-standardize using std and mean from bin k for both samples
            sample2 = zdata{krand(2)}*sest(k)+mu(k);
            
            % compute stats for current bin of conditioning variable
            % (kstest2 test decision H0 = data in samples come from same
            % distribution. stats is the test statistic which is the max
            % absolute difference between the cdfs
            [~,~,stats(k,1)] = kstest2(sample1,sample2);
            [~,~,stats(k,2)] = kstest2(sample1,data{k});
            
            if showconditionals  % plot histograms for each bin
                [~,bin]=hist([sample1;sample2]);
                h1=hist(sample1,bin)/count(krand(1));
                h2=hist(data{k},bin)/count(krand(2));
                clf; bar(bin,[h1' h2']); title([i n k]); pause;
            end
                           
        end  % end of one iteration over all bins
        
        % agregate stats over conditioning bins. stats has k rows (bins)
        % and 2 columns (1 = in-between sample statistic, 2 = selected bin
        % vs sample.
        meanstats(n,:) = mean(stats);  % Append iteration to meanstats.
        % 2 values: mean stats for sample vs sample (1) and for data vs
        % sample (2) over all bins k
        % mean(stats) computes mean over columns
        
    end  % end of n resampling iterations 
   
    % evaluate shuffle stats
    % i can be 1 = H1 or 2 = H2 
    % zscore has two values. One for H1 and one for H2
    zscore(i) = diff(mean(meanstats))/mean(std(meanstats)); %substract the mean between both averages
    
    subplot(4,4,8+i*4)
    hist(meanstats,Nrand/10); hold off; axis tight
    ax=axis; t=text(ax(1)+0.5*(ax(2)-ax(1)),ax(3)+0.85*(ax(4)-ax(3)),['z=' num2str(zscore(i),2)]); 
    set(t,'BackgroundColor', [0.95 0.95 0.95])
    xlabel('KS stats'); 
    
end % end testing hypothesis, go to next.

colormap cool
%sublabel;

saveas(gca,strcat('KS_results/',product.name, funding.name,'.png'))

summaryTable(l,1) = p;
summaryTable(l,2) = f;
summaryTable(l,5:6) = zscore;
summaryTable(l,3:4) = rsquare;
%summaryTable(l,7) = pstrlong;
%summaryTable(l,8) = fstrlong;
summaryTable(l,:)
l = l+1;
%closeall
    end  
end
export(dataset(summaryTable),'File' ,'KS_results/summary_table_Rsq_z_scores.csv','Delimiter',',');
