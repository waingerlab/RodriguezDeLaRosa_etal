%% Function that transforms data from IEI into APs/min

function dataout = ieitransf(file,sheet,coordinates)

data = xlsread(file,sheet,coordinates);
data = data(:);
data(isnan(data))=[];
dataout = 60./data;
