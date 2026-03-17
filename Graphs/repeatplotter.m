%% Repeated measures plotter
%data should be cell array where each column is a genotype and each row is
%a timepoint

%for example:




function graph = repeatplotter(data,xpoints, xlimit, ypoints, ylimit, linecolors, shadecolors);

groupnum = size(data);

average = cellfun(@mean,data);
stdev = cellfun(@std,data);
count = cellfun(@numel,data);
SEM = stdev./sqrt(count);
high = average+SEM;
low = average-SEM;

hold on
for nn=1:groupnum(2)
    for n=1:(groupnum(1)-1)
        f=fill([xpoints(n) xpoints(n) xpoints(n+1) xpoints(n+1)],[high(n,nn),low(n,nn),  low(n+1,nn), high(n+1,nn)],shadecolors{nn},'linestyle','none');
        set(f,'facealpha',.33);
        line([xpoints(n);xpoints(n+1)],[average(n,nn); average(n+1,nn)],'color',linecolors{nn});
    end
end
xlim(xlimit);
ylim(ylimit);
set(gca,'Xtick',xpoints,'Ytick',ypoints);
hold off

        