% Tests on retrieval of direction of causality as inferred by method of
% nonlinear causal inference (additive noise n??)
%
% bivariate models
%
% Simulations of model y = f(x), p(y|f(x)) 
%
% Given x -> y 
% Choose: 
% p(x) [uniform, normal, logNormal..] 
% y = f(x) [polynomial degree 1, 2 , ...]
% N [ number of samples ]
% g(n) [some noise function, additive noise, normal]

close all; clf

N = 70000;

% p(x)

% uniform 0 1
%x = rand(1,N);

%normal
mu = 3; var = 1.4;
x = (randn(N,1)+mu)*var;

%log normal
x = exp(x);

%measurement error P(x)
%x = x + randn(size(x)) ;

% additive noise
n = randn(size(x))*200000;

% deterministic linear function
b = 0.01;
a = 0.08;
% calculate y from x
y = a*x.^2 + b*x + n;

%y = 0.4*x + n;

subplot(1,2,1);
loglog(x,y,'.','markersize',8);
xlabel('X -> Y');
grid on; axis tight;
subplot(1,2,2);
loglog(y,x,'.','markersize',8);
xlabel('Y -> X');
grid on; axis tight;



%%

% Run causal inference on data

Nrand=100; % number of shuffles when testing hypothesis that conditionals are the same except for mean and std
showconditionals = 0; % to check if the conditionals looks like properly adjusted accross bins

    
        subplot(2,1,1);
        loglog(x,y,'.','markersize',6);
        grid on; axis tight;
        ylabel('y')
        xlabel('x')
        legend(['N=' num2str(length(x))],'Location','southeast')
        title ('X vs Y')
        
        xlim([min(x)-0.5,max(x)+0.5]);
        %set(gca,'XTick',[min(x)-0.5 : 0.15 : max(x)+0.5]);
        
        clear meanstats
        for i=1:2
            
            switch i % selcect H1 or H2
                case 1, X=x; Y=y; xstr='X'; ystr='Y';
                case 2, X=y; Y=x; xstr='Y'; ystr='X';
            end
            
            % select bins edges with equal number of points per bin
            edge=prctile(X,0:4:99);
            
            % compute conditional statistics in each bin of conditioning variable
            for k=length(edge)-1:-1:1
                %binnoise = 0; randn(size(X))*(edge(k+1)-edge(k))*2; 
                %indx = find(edge(k)<=(X+binnoise) & (X+binnoise)<edge(k+1));
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
            scatter(cond,mu,'filled');
            xlabel(xstr); ylabel(['H' num2str(i) ': ' xstr '->' ystr 10 'mean ' ystr]); axis tight
            
            subplot(4,4,6+i*4)
            scatter(cond,sigma,'filled');
            xlabel(xstr); ylabel(['std ' ystr]); axis tight
            
            subplot(4,4,7+i*4)
            scatter(mu,sigma,'filled');
            xlabel(['mean ' ystr]); ylabel(['std ' ystr]); axis tight
            
            % estimate the std from the mean assuming linear model
            sest=polyval(polyfit(mu,sigma,3),mu);  % with offset
            
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
        [h,bin]=hist(meanstats,50); 
        bar(bin,h/Nrand,1); axis tight
        legend('H1','H2')
        ax=axis; t=text(ax(1)+0.5*(ax(2)-ax(1)),ax(3)+0.15*(ax(4)-ax(3)),['d''=' num2str(ksdiff(pnr,fnr),2)]);
        set(t,'BackgroundColor', [0.95 0.95 0.95])
        xlabel('KS stats');
        
        
        colormap cool
        % sublabel; % get this function on the mathworks community repository.
        
        %saveas(gca,'funding_vs_productivity_model.png')
        
   

