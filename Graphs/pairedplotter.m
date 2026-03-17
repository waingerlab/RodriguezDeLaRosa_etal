function [graph] = pairedplotter(datacolumns,lim,points,pointlabels,pointsize)
%This function plots paired data
%
%Inputs are:
%datacolumns = cell array where columns are genotypes and rows are the
%pairs (should only be two rows)
%lim = the y axis lower and upper limits
%points = the points on the y axis to plot
%pointlabels = a cell array of how the y axis points should be labeled
%pointsize = how large individual points should be
%
%
%Example inputs are:
%datacolumns = {wtns,c9ns,tdpwtns,tdpmutns; wtnsnacl,c9nsnacl,tdpwtnsnacl,tdpmutnsnacl};
%lim = [0 250];
%points = [0 50 100 150 200 250];
%pointlabels = {'0','50','100','150','200','250'};
%pointsize = 50;


groupnum= numel(datacolumns(1,:));
graph = figure('Color',[1 1 1],'Position',[100, 100,(125*groupnum)+300,500]);

data=datacolumns;

for n=1:groupnum
    %statistics
    avg1= mean(data{1,n});
    counts1 = numel(data{1,n});
    SEM1 = std(data{1,n})/sqrt(counts1);
    
    avg2= mean(data{2,n});
    counts2 = numel(data{2,n});
    SEM2 = std(data{2,n})/sqrt(counts2);
    
    
    nscale1 = (1.2*n)-.75; %controls how closely genotypes are grouped
    nscale2 = (1.2*n)-.25; %controls how closely genotypes are grouped
    
    %format the data for plotting
    %scatter(nscale*ones(numel(data{n}),1),data{n}, 50, 'x','k') %plots x and y values, with x representing group number
    position1 = nscale1*ones(counts1,1);
    position2 = nscale2*ones(counts2,1);
    
    %the if statement allows a specific column's data point color to be
    %changed
    if n ==4
    scatter(position1,sort(data{1,n}),pointsize,[.2 .2 .2],'filled');
    else
    scatter(position1,sort(data{1,n}),pointsize,[.2 .2 .2],'filled');
    end
    hold on
    
    %this next section draws lines between datapoints
    for nn =1:counts1
        line([nscale1,nscale2],[data{1,n}(nn),data{2,n}(nn)],'color','k','linewidth',1.5);
    end
    
    
    %this next bit specifies bar colors and widths
    meanwidth = .2;
    SEMwidth = meanwidth/2; 
    
    %this next section plots the first of each pair
    if n == 4 %this first if statement can be used to change the color of a specific bar
    color = 'r';
    line([nscale1-meanwidth, nscale1+meanwidth],[avg1,avg1],'color',color,'linewidth',2.5); %mean
    line([nscale1-SEMwidth, nscale1+SEMwidth],[avg1-SEM1,avg1-SEM1],'color',color,'linewidth',2); %lower limit
    line([nscale1-SEMwidth, nscale1+SEMwidth],[avg1+SEM1,avg1+SEM1],'color',color,'linewidth',2); %upper limit
    line([nscale1, nscale1],[avg1+SEM1,avg1-SEM1],'color',color,'linewidth',2); %connects upper and lower limit
    else
    line([nscale1-meanwidth, nscale1+meanwidth],[avg1,avg1],'color','r','linewidth',2.5); %mean
    line([nscale1-SEMwidth, nscale1+SEMwidth],[avg1-SEM1,avg1-SEM1],'color','r','linewidth',1.5); %lower limit
    line([nscale1-SEMwidth, nscale1+SEMwidth],[avg1+SEM1,avg1+SEM1],'color','r','linewidth',1.5); %upper limit
    line([nscale1, nscale1],[avg1+SEM1,avg1-SEM1],'color','r','linewidth',1.5); %connects upper and lower limit
    end

    %this section plots the second of each pair
    if n ==4
    scatter(position2,sort(data{2,n}),pointsize,[.2 .2 .2],'filled');
    else
    scatter(position2,sort(data{2,n}),pointsize,[.2 .2 .2],'filled');
    end
    hold on
    
    
    if n == 4 %this first if statement can be used to change the color of a specific bar
    color = 'r';
    line([nscale2-meanwidth, nscale2+meanwidth],[avg2,avg2],'color',color,'linewidth',2.5); %mean
    line([nscale2-SEMwidth, nscale2+SEMwidth],[avg2-SEM2,avg2-SEM2],'color',color,'linewidth',2); %lower limit
    line([nscale2-SEMwidth, nscale2+SEMwidth],[avg2+SEM2,avg2+SEM2],'color',color,'linewidth',2); %upper limit
    line([nscale2, nscale2],[avg2+SEM2,avg2-SEM2],'color',color,'linewidth',2); %connects upper and lower limit
    else
    line([nscale2-meanwidth, nscale2+meanwidth],[avg2,avg2],'color','r','linewidth',2.5); %mean
    line([nscale2-SEMwidth, nscale2+SEMwidth],[avg2-SEM2,avg2-SEM2],'color','r','linewidth',1.5); %lower limit
    line([nscale2-SEMwidth, nscale2+SEMwidth],[avg2+SEM2,avg2+SEM2],'color','r','linewidth',1.5); %upper limit
    line([nscale2, nscale2],[avg2+SEM2,avg2-SEM2],'color','r','linewidth',1.5); %connects upper and lower limit
    end

     
end

ax.YAxis.LineWidth = 10000;
ax.XAxis.LineWidth = 10000;
xlim([0 nscale2+(meanwidth*1.5)]);
%these next few lines control the y axis
ylim(lim)
%ylim([0 mean(data{1}*2)]);
set(gca,'Ytick',points,'YTickLabel',pointlabels);
%this next line normalizes the y label to 100% of 1st column
%set(gca,'YTickLabel',{'0','50','100','150','200'},'YTick',[0,mean(data{1})/2,mean(data{1}),mean(data{1})*1.5,mean(data{1}*2)]);%this line is specifically for locomotion
set(gca,'XTickLabel',{},'XTick',[]);
set(gca,'FontSize',20);
set(gca,'FontName','Arial');
set(gca,'Linewidth',2);

hold off
end