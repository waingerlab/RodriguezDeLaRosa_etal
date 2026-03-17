%% Cumulative plotter takes a cell array where rows are replicates and columns are genotypes. meancolor is the averaged distribution and replicatecolor is the color of the replicates (these are matricies)

function cumulativeplotter(cellarray_in, meancolor, replicatecolor,replicateoption)

if exist('replicateoption','var') & strcmp(replicateoption,'showreplicates') %this allows replicates to be plotted
    options.color = true;
else options.color = false;
end

dimensions = size(cellarray_in);
cellarray_in = cellfun(@sort, cellarray_in,'uniformoutput',0);
objectnumber = cellfun(@numel, cellarray_in);
minimum = min(objectnumber,[],'all');

figure; hold on;
for n=1:dimensions(2) %number of genotypes
    minobservations = min(objectnumber(:,n));
    stored_distributions = zeros(minobservations,dimensions(1));
    for nn=1:dimensions(1) %number of replicates
        workingdata = cellarray_in{nn,n};
        interval1 = numel(workingdata)/minobservations; %used for storing data for averaging later
        stored_distributions(:,nn) = workingdata(floor(1:interval1:numel(workingdata)));
        
        interval2 = minimum/numel(workingdata); %used for plotting data
        if exist('replicateoption','var') & strcmp(replicateoption,'showreplicates') %this allows replicates to be plotted
            p=plot(workingdata,[interval2:interval2:minimum]','color',replicatecolor{n});
            p.Color(4)=0.33;
        end
    end
    averageddata = mean(stored_distributions,2);
    interval3 = minimum/numel(averageddata);
    plot(averageddata,[interval3:interval3:minimum]','color',meancolor{n},'LineWidth',2);
end

ylim([0 minimum]);
yticks([0 minimum/2 minimum]);
yticklabels({'0', '50', '100'});

hold off