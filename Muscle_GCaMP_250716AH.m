%% Batch script for muscle GCAMP analysis
%To use this script, navigate to the folder above all the timepoint
%folders. Then run the first block.

% This script takes the IMX output from Pram's calcium imaging and then:
% 1) compiles it into a video
% 2) masks off the YFP channel
% 3) measures the calcium activity in each object at each timepoint
% 4) determines if the calcium response is from a cell based on amplitude
% 5) writes out a mask, traces for each object, and a spreadsheet of data

addpath('F:/Scripts/Universal_functions')
options.color = false; %for writing out 3d tif

% delete all thumb files
folders = glob('*/Timepoint_*'); %this should point to the folder above where the timepoint folders are
for n=1:numel(folders)
    cd(folders{n})
    delete('*_thumb*.tif'); %deletes thumb files
    cd ../..
end

%identify all wells to analyze
filenames = glob('*/Timepoint_1/*.tif'); %need to add w1* so that only 1 channel is capatured for filenames
wells = cell(numel(filenames),1);
for n=1:numel(filenames)
    prelimwell = regexp(filenames{n},'-[0-9][0-9][0-9][0-9][0-9][0-9]_[A-Z][0-9][0-9]','match');
    wells(n) = {prelimwell{1}((end-2):end)};
end

%make list of frames to keep in order
framenumber = numel(folders);
[~,endindex] = regexp(folders{1},'TimePoint_');
basename = folders{1}(1:endindex);
folders = cell(framenumber,1);
for n=1:framenumber;
    folders{n}=strcat(basename,int2str(n));
end

%make stacks for each well and do background subtraction
mkdir('stacks')

for nn=1:numel(wells)
    stack = zeros(1024,1024,numel(folders),'uint16'); %3rd dimension must be 3 or 4, so last one is just black
    for n=1:numel(folders)
        %read in frame 1
        frame = glob(strcat(folders{n},'\*[0-9][0-9][0-9]_',wells{nn},'*'));
        frame = tifread(frame{1});
        %stack(:,:,n) = frame;
        background = imgaussfilt(frame,50);
        stack(:,:,n) = frame-background;

    end
    saveastiff(stack,strcat('stacks/',wells{nn}(1:3),'.tif'),options);
    
end

%% video analysis
mkdir('traces');
mkdir('masks');
videonames = glob('stacks/*.tif');
options.color = true;


%final output values - preallocation
frequencyout = NaN(numel(videonames),1);
amplitudeout = NaN(numel(videonames),1);
objectsout = NaN(numel(videonames),1);

for n=1:numel(videonames)
    n
    video = squeeze(tifread(videonames{n})); %load video
    framenumber = size(video);
    framenumber = framenumber(end);
    
    %normalize intensity across time

   %find average change in intensity over time
    frameintensity = mean(video,[1 2]);
    frameintensity = squeeze(frameintensity);
    
    difference = double(frameintensity./min(frameintensity)); %find minimal intensity frame and make normalization factor (usually last frame), make uint16 for multiplying out
    
    %normalize fluorescence of each frame
    normarray = ones(1024,1024,framenumber,'uint16'); %preallocate memory, add 10,000 to keep all values above 0 when normarlizing, may not be necessary with signed integers
    
    for nn=1:framenumber
        normarray(:,:,nn) = video(:,:,nn)./difference(nn);
    end
    
    clear video %eliminates video to free up RAM

    %finding areas that fluctuate
    int=5; %frame interval
    interval = floor(framenumber/int);
    diffstack = zeros(1024,1024,interval,'uint16'); %this is contrast stack
    for nn=1:interval
        ind1 = (nn-1)*int+1;
        ind2 = (nn-1)*int+int;
        sample = normarray(:,:,ind1:ind2);
        maximum = max(sample,[],3);
        minimum = min(sample,[],3);
        diff= maximum-minimum;
        diffstack(:,:,nn)=diff;
    end

    diffstack2 = diffstack - 700;
    diffsummary = sum(diffstack2,3);

    clear diffstack
    clear diffstack2

    diffbw = im2bw(diffsummary);
    diffareafilt = bwareafilt(diffbw,[40 20000]);
    difflabel = bwlabel(diffareafilt);


    %% now go through stacks and save information for each area
    objectnumber = max(difflabel,[],'all'); %will be used for object loop
    fluorescence = zeros(objectnumber,framenumber);
    
    
    
    %rewriting loop to be more efficient
    dataextraction = cell(objectnumber,1);
    
    for nn=1:objectnumber
        dataextraction{nn}=find(difflabel==nn);
    end
    
    
    %go through objects and pull values of interest for both masks; use 3d
    %matricies
    threedimindex = 0:numel(difflabel):(numel(difflabel)*(framenumber-1));
    
    for nn=1:objectnumber
        positions = dataextraction{nn}; %find positions in all frames
        threedimindex2 = repmat(threedimindex,[numel(positions),1]);
        threedimpositions = repmat(positions,[1,framenumber]);
        threedimpositions = threedimpositions+threedimindex2;
        values = normarray(threedimpositions);
        ROIintensity = sum(values,1)./numel(positions);
        
        fluorescence(nn,:) = ROIintensity;
    end

    %% This section calculates peaks using dF/F methodology
    tracepeaks = cell(objectnumber,1); %make cell array to store peak information in each trace
    dFarray = zeros(objectnumber,framenumber);
    amplitude = NaN(objectnumber,1);
    frequency = NaN(objectnumber,1);

    for nn = 1:objectnumber
        trace = fluorescence(nn,:);
        trace = trace - min(trace)+1; %normalization to make all values >=1
        med = median(trace);
        dF = trace/med; %normalization to dF/F
        dFarray(nn,:)=dF;
        
        nnn = 13;
        peaks = [];
        baseline = [];
        while nnn<(framenumber-7)
            start = nnn-4;
            [maxval, maxind] = max(dF(start:nnn));
            testpeakind = start + maxind-1;
            [minval, minind] = min(dF((testpeakind-5):(testpeakind)));
            minind = testpeakind-6+minind; %should be -5, but have to account for counting within
            meanmin = mean(dF(minind-2:minind)); %could also possibly determine local standard deviation and use to compute peak value below
            delta = maxval-meanmin;
            if delta>(.8)
                %refine the peak (make sure not just catching rising phase or decay phase in
                %window)
                [rmaxval, rmaxind] = max(dF(testpeakind-2:testpeakind+7));
                rpeakind = testpeakind + rmaxind -3;
                peaks = [peaks; rpeakind, dF(rpeakind)-dF(minind)];
                nnn=rpeakind+7;
            else
                baseline = [baseline; std(dF(start:nnn))];
            end
            nnn=nnn+4;
        end
        if 4*mean(baseline)<=.5 %if initial baseline < noise, write
            tracepeaks{nn} = peaks;  
        else %if it is not less than noise, use the refined baseline
            nnn = 13;
            peaks = [];
            baseline = 4*mean(baseline);
            while nnn<(framenumber-7)
                start = nnn-4;
                [maxval, maxind] = max(dF(start:nnn));
                testpeakind = start + maxind-1;
                [minval, minind] = min(dF((testpeakind-5):(testpeakind)));
                minind = testpeakind-6+minind; %should be -5, but have to account for counting within
                meanmin = mean(dF(minind-2:minind)); %could also possibly determine local standard deviation and use to compute peak value below
                delta = maxval-meanmin;
                if delta>(baseline)
                    %refine the peak (make sure not just catching rising phase or decay phase in
                    %window)
                    [rmaxval, rmaxind] = max(dF(testpeakind-2:testpeakind+7));
                    rpeakind = testpeakind + rmaxind -3;
                    peaks = [peaks; rpeakind, dF(rpeakind)-dF(minind)];
                    nnn=rpeakind+7;
                end
                nnn=nnn+4;
            end
            tracepeaks{nn} = peaks; 
        end
        if numel(peaks>0)
            amplitude(nn,1) = mean(peaks(:,2));
            frequency(nn,1) = numel(peaks(:,1));
        end
    end

    amplitude(isnan(amplitude))=[]; %delete non-responding objects
    frequency(isnan(frequency))=[];

    %% plot traces
    f=figure('Position',[0 0 1024 1024], 'visible','off'); %make figure; not visible
    colorwheel = {'r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m','r','g','b','m'};
    for nn = 1:objectnumber
        trace = dFarray(nn,:)+(-10*nn);%the 2*n factor is to offset the traces
        
        plot(trace,colorwheel{nn});
        hold on

        if numel(tracepeaks{nn})>0
            peaks = tracepeaks{nn}(:,1);
            scatter(peaks,trace(peaks),10,'k');
        end
        text(10,trace(1)+3,int2str(nn),'Fontsize',15);
    end
    hold off


    %% write out information
    set(f, 'CreateFcn', 'set(gcbo,''Visible'',''on'')');
    saveas(f,strcat('traces/',videonames{n}(8:10)),'fig'); 
    close(f);

    saveastiff(label2rgb(difflabel,'jet','k','shuffle'),strcat('masks/',videonames{n}(8:10),'.tif'),options);

    amplitudeout(n)=mean(amplitude);
    frequencyout(n)=mean(frequency);
    objectsout(n)=numel(amplitude);


end

%write out data
xlswrite('muscle_analysis.xlsx',videonames,'Summary','A2');
xlswrite('muscle_analysis.xlsx',{'files','objectnumber','amplitude (dF/F)','events/min'},'Summary','A1');
xlswrite('muscle_analysis.xlsx',[objectsout,amplitudeout,frequencyout],'Summary','B2');

    
%% plot data - template, could be modified to graph your data

addpath('F:/Scripts/Graphs')

file = 'muscle_analysis_250714.xlsx';
sheet = 'Combined';

pre_muscle_objects = xlsread(file,sheet,'b2:b11');
post_muscle_objects = xlsread(file,sheet,'e2:e11');

pre_MN_objects = xlsread(file,sheet,'b12:b24');
post_MN_objects = xlsread(file,sheet,'e12:e24');

datacolumns = {pre_MN_objects, pre_muscle_objects;post_MN_objects,post_muscle_objects};
lim = [0 300];
points = [0:100:300];
pointlabels = {'0','100','200','300'};
pointsize = 50;

pairedplotter(datacolumns,lim,points,pointlabels,pointsize)