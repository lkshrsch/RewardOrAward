x=csvread('KS_results/summary_table_Rsq_z_scores.csv',1,0);

close all; clc

GSI_Marker = 'O+*pX';
%RCR_Color = 'kymcrgb';
RCR_Color = {[0 0 0],[1 0 1],[0 204 204]/255, [1 0 0],[0 1 0],[0 0 1], [153 0 153]/255};

% 1=mean RCR; 2=sum pubs; 3=sum RCR; 4=annual pubs; 
% 5=median RCR; 6=annual RCR; 7=max RCR
variablesOfInterestRCR = [1,4,5,6,7];

% 1=years of funding; 2=total rpg $$; 3=annual $$;
% 4=total rpg rci; 5=annual rci
variablesOfInterestGSI = [1,3,5];

RCRNames = {'mean RCR','sum pubs','sum RCR','annual pubs','median RCR','annual RCR','max RCR'};
GSINames = {'years of funding','total rpg $$','annual $$','total rpg rci','annual rci'};

% R² plot
rsq = subplot(1,6,1);
hold on;
for i = 1:length(x)
    if any(x(i,1) == variablesOfInterestRCR) && any(x(i,2)==variablesOfInterestGSI); plot(x(i,3:4)',strcat(GSI_Marker(x(i,2)),'-'),'color',cell2mat(RCR_Color(x(i,1)))); end
grid on;
end

ylabel('R^2')
set(gca,'xtick',[1 2]);
set(gca,'xticklabel',{'H1','H2'})
hold off;

xlim([0.5 2.5])

% d' plot 
subplot(1,6,2)
hold on;
for i = 1:length(x)
if any(x(i,1) == variablesOfInterestRCR) && any(x(i,2)==variablesOfInterestGSI); plot(x(i,5:6)',strcat(GSI_Marker(x(i,2)),'-'),'color',cell2mat(RCR_Color(x(i,1)))); end
grid on;
end
ylabel('d''')
set(gca,'xtick',[1 2]);
set(gca,'xticklabel',{'H1','H2'})
xlim([0.5 2.5])

hold off;

text(-0.9, 72, 'Annual','Fontsize',16);
text(5.6, 72, 'Total','Fontsize',16);


% Plot Totals -----------------------------------

% 1=mean RCR; 2=sum pubs; 3=sum RCR; 4=annual pubs; 
% 5=median RCR; 6=annual RCR; 7=max RCR
variablesOfInterestRCR = [2,3];

% 1=years of funding; 2=total rpg $$; 3=annual $$;
% 4=total rpg rci; 5=annual rci
variablesOfInterestGSI = [2,4];

RCRNames = {'mean RCR','sum pubs','sum RCR','annual pubs','median RCR','annual RCR','max RCR'};
GSINames = {'years of funding','total rpg $$','annual $$','total rpg rci','annual rci'};

% R² plot
rsq = subplot(1,6,3);
hold on;
for i = 1:length(x)
    if any(x(i,1) == variablesOfInterestRCR) && any(x(i,2)==variablesOfInterestGSI); plot(x(i,3:4)',strcat(GSI_Marker(x(i,2)),'-'),'color',cell2mat(RCR_Color(x(i,1)))); end
grid on;
end

ylabel('R^2')
set(gca,'xtick',[1 2]);
set(gca,'xticklabel',{'H1','H2'})
hold off;

xlim([0.5 2.5])

% d' plot 
subplot(1,6,4)
hold on;
for i = 1:length(x)
if any(x(i,1) == variablesOfInterestRCR) && any(x(i,2)==variablesOfInterestGSI); plot(x(i,5:6)',strcat(GSI_Marker(x(i,2)),'-'),'color',cell2mat(RCR_Color(x(i,1)))); end
grid on;
end
ylabel('d''')
set(gca,'xtick',[1 2]);
set(gca,'xticklabel',{'H1','H2'})
xlim([0.5 2.5])
hold off;

% dummy plot for legend
hsub = subplot(1,6,5);
hold on;
set(hsub, 'Visible','off');
k = 1;
for j = [1,4,5,6,7,2,3]
hsub(k) = plot(nan,nan,'-','color',cell2mat(RCR_Color(j)));
k = k+1;
end
hold off;

columnlegend(1,hsub,'Location','North', legend(RCRNames(1,[1,4,5,6,7,2,3])));


hsub2 = subplot(1,6,6);
hold on;
set(hsub2, 'Visible','off');
k = 1;
for j = [1,3,5,2,4]
hsub2(k) = plot(nan,nan,strcat(GSI_Marker(j),'-'), 'color','k');
k = k+1;
end
hold off;

columnlegend(1, hsub2,'Location','North', legend(GSINames(1,[1,3,5,2,4])));

%columnlegend(2,hsub, 'Location','north',legend([RCRNames,GSINames]));

%set(gcf, 'PaperUnits', 'centimeters');
%set(gcf, 'PaperPosition', [0 0 26 10]); %x_width=10cm y_width=15cm
%saveas(gcf,'KS_results/summaryfigure.pdf')