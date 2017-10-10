% Causal inference analysis to answer the question:
%
% Do citations "cause" funding (H1), or does funding "cause" citations (H2)?
%
% Based on NIH data originally released here:
%
% Lauer et al., bioRxiv, May. 26, 2017; doi: http://dx.doi.org/10.1101/142554
%
% latest version of this code
%   http://parralab.org/reward-or-award/testmodel.m
% data files:
%   http://parralab.org/reward-or-award/Lauer_2017_data.mat
%   http://parralab.org/reward-or-award/Lauer_2017_data_R.mat
%
% The second file only inludes R grant recipients, excluding in particular
% P grants which have been questioned for the way that publications were
% attributed and GSI are counted.
%
% (c) June 15, 2017, Lucas C. Parra, parra@ccny.cuny.edu

%clear all; clf;

% product
%  1   'Mean Relative Citation Ratio (RCR)'
%  2   'Sum publications'
%  3   'Sum Relative Citation Ratio (RCR)'
%  4   'Annual Publications (Pubs)'
%  5   'Median Relative Citation Ratio (RCR)'
%  6   'Annual Relative Citation Ratio (RCR)'
%  7   'Max Relative Citation Ratio (RCR)'
% funding
%  1   'Years of funding (GSI)'
%  2   'Total Grant Support Index (GSI)'
%  3   'Annual Dollars ($$)'
%  4   'Total Dollar ($$)'
%  5   'Annual Grant Support Index (GSI)'

% set some parameters for shuffle
Nrand=100; % number of shuffles when testing hypothesis that conditionals are the same except for mean and std
showconditionals = 0; % to check if the conditionals looks like properly adjusted accross bins

% for pnr=[1:7], for fnr=[1:5]
figure(1)
for pnr=[1 4:7], for fnr=[3]

        load('Lauer_2017_data.mat','product','funding'); % all award recipients
       % load('Lauer_2017_data_R.mat','product','funding'); % only R grant recipients.

        % pick a measure of productivity and funding (see converttable.m).
        product = product(pnr);
        funding = funding(fnr);
        
        
        % filter for low values as those were hard-coded for zero publications/funidng
        indx=find(funding.data<=funding.min); funding.data(indx)=[]; product.data(indx)=[];
        indx=find(product.data<=product.min); funding.data(indx)=[]; product.data(indx)=[];
        
        subplot(2,1,1);
        loglog(funding.data,product.data,'.','markersize',4);
        grid on; axis tight;
        ylabel(product.name)
        xlabel(funding.name)
        legend(['PI (N=' num2str(length(product.data)) ')'],'Location','southeast')
        title ('Research Productivity versus Funding')
        
        clear meanstats
        for i=1:2
            
            switch i % selcect H1 or H2
                case 1, X=(product.data); Y=(funding.data); xstr=product.acronym; ystr=funding.acronym;
                case 2, X=(funding.data); Y=(product.data); xstr=funding.acronym; ystr=product.acronym;
            end
            
            % select bins edges with equal number of points per bin
            edge=prctile(X,0:4:99);
            
            % compute conditional statistics in each bin of conditioning variable
            for k=length(edge)-1:-1:1
                
                %binnoise = 0; randn(size(X))*(edge(k+1)-edge(k))*2; 
                %indx = find(edge(k)<=(X+binnoise) & (X+binnoise)<edge(k+1));
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
            scatter(cond,mu,'filled');
            xlabel(xstr); ylabel(['H' num2str(i) ': ' xstr '->' ystr 10 'mean ' ystr]); axis tight
            
            subplot(4,4,6+i*4)
            scatter(cond,sigma,'filled');
            xlabel(xstr); ylabel(['std ' ystr]); axis tight
            
            subplot(4,4,7+i*4)
            scatter(mu,sigma,'filled');
            xlabel(['mean ' ystr]); ylabel(['std ' ystr]); axis tight
            
            % estimate the std from the mean assuming linear model
            sest=polyval(polyfit(mu,sigma,2),mu);  % with offset
            
            hold on; [~,indx]=sort(mu); plot(mu(indx),sest(indx)); hold off;
            rsquare(pnr,fnr,i) = 1-sum((sest-sigma).^2)/sum((sigma-mean(sigma)).^2);
            ax=axis; t=text(ax(1)+0.5*(ax(2)-ax(1)),ax(3)+0.15*(ax(4)-ax(3)),['R^2=' num2str(rsquare(pnr,fnr,i),2)]);
            set(t,'BackgroundColor', [0.95 0.95 0.95])
            
            drawnow
            
            % generate standartized version of the data only using mu, and sest derived from mu
            for k=1:K, zdata{k} = (data{k}-mu(k))/sest(k); end
            
            % do shuffle stats
            for n=Nrand:-1:1
        
                clear stats
                for k=1:K
                    
                    % sample by shuffling bins of conditioning variable and adjusting mean and estimated std
                    krand = [1:k-1 k+1:K]; krand=krand(randperm(K-1,2)); % 2 random bins excluding kth bin
                    sample1 = zdata{krand(1)};
                    sample2 = zdata{krand(2)};
                    
                    % compute stats for current bin of conditioning variable
                    [~,~,stats(k)] = kstest2(sample1,sample2);
                    
                    if showconditionals
                        [~,bin]=hist([sample1;sample2]);
                        h1=hist(sample1,bin)/count(krand(1));
                        h2=hist(sample2,bin)/count(krand(2));
                        clf; bar(bin,[h1' h2']); title([i n k]);
                        [sigma(k) sest(k); mean(data{k}) mean(sample1); std(data{k}) std(sample1)]
                        pause;
                        
                    end
                    
                end
                
                % agregate stats over conditioning bins
                meanstats(n,i) = mean(stats);
                
            end
            
            
        end
        
        ksmean (pnr,fnr,:) = mean(meanstats);
        ksdiff(pnr,fnr) = diff(mean(meanstats))/mean(std(meanstats)); % positive if H1 more uniform across bins
        
        subplot(2,4,8)
        [h,bin]=hist(meanstats,10); 
        bar(bin,h/Nrand,1); axis tight
        legend('H1','H2')
        ax=axis; t=text(ax(1)+0.5*(ax(2)-ax(1)),ax(3)+0.15*(ax(4)-ax(3)),['d''=' num2str(ksdiff(pnr,fnr),2)]);
        set(t,'BackgroundColor', [0.95 0.95 0.95])
        xlabel('KS stats');
        
        
        colormap cool
        % sublabel; % get this function on the mathworks community repository.
        
        %saveas(gca,'funding_vs_productivity_model.png')
        
    end
end


figure(2)
subplot(2,2,1); plot(reshape(rsquare([1 4 5 6 7],[3 ],:),5,2)'); ylabel('R^2')
set(gca,'xtick',[1 2]); xlim([0.5 2.5]); set(gca,'xticklabel',{'H1','H2'})
subplot(2,2,2); plot(reshape(rsquare([2 3],[2 4],:),4,2)'); ylabel('R^2')
set(gca,'xtick',[1 2]); xlim([0.5 2.5]); set(gca,'xticklabel',{'H1','H2'})
subplot(2,2,3); plot(reshape(ksmean([1 4 5 6 7],[3 ],:),5,2)'); ylabel('KS')
set(gca,'xtick',[1 2]); xlim([0.5 2.5]); set(gca,'xticklabel',{'H1','H2'})
subplot(2,2,4); plot(reshape(ksmean([2 3],[2 4],:),4,2)'); ylabel('KS')
set(gca,'xtick',[1 2]); xlim([0.5 2.5]); set(gca,'xticklabel',{'H1','H2'})

