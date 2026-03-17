%% This script loads data for statistical analysis and plotting

%load data into matlab
file = 'AH_190422_IEI_data_organized.xlsx';

%11a and mgh5b section
sheet = '190220_11a_mgh5b';

r111ad11 = xlsread(file,sheet,'b17:ee22'); %replicate1, genotype =11a, day 11
r111ad11 = r111ad11(:);
r111ad11(isnan(r111ad11))=[];
r111ad11=60./r111ad11; %transforms into events/s

r1mgh5bd11 = xlsread(file,sheet,'b24:ee29');
r1mgh5bd11 = r1mgh5bd11(:);
r1mgh5bd11(isnan(r1mgh5bd11))=[];
r1mgh5bd11=60./r1mgh5bd11;

r111ad15 = xlsread(file,sheet,'b33:ee44');
r111ad15 = r111ad15(:);
r111ad15(isnan(r111ad15))=[];
r111ad15=60./r111ad15;

r1mgh5bd15 = xlsread(file,sheet,'b46:ee57');
r1mgh5bd15 = r1mgh5bd15(:);
r1mgh5bd15(isnan(r1mgh5bd15))=[];
r1mgh5bd15=60./r1mgh5bd15;

r111ad20 = xlsread(file,sheet,'b61:ee67');
r111ad20 = r111ad20(:);
r111ad20(isnan(r111ad20))=[];
r111ad20=60./r111ad20;

r1mgh5bd20 = xlsread(file,sheet,'b69:ee75');
r1mgh5bd20 = r1mgh5bd20(:);
r1mgh5bd20(isnan(r1mgh5bd20))=[];
r1mgh5bd20=60./r1mgh5bd20;

r111ad27 = xlsread(file,sheet,'b79:ee84');
r111ad27 = r111ad27(:);
r111ad27(isnan(r111ad27))=[];
r111ad27=60./r111ad27;

r1mgh5bd27 = xlsread(file,sheet,'b86:ee91');
r1mgh5bd27 = r1mgh5bd27(:);
r1mgh5bd27(isnan(r1mgh5bd27))=[];
r1mgh5bd27=60./r1mgh5bd27;

r111ad33 = xlsread(file,sheet,'b95:ee103');
r111ad33 = r111ad33(:);
r111ad33(isnan(r111ad33))=[];
r111ad33=60./r111ad33;

r1mgh5bd33 = xlsread(file,sheet,'b105:ee113');
r1mgh5bd33 = r1mgh5bd33(:);
r1mgh5bd33(isnan(r1mgh5bd33))=[];
r1mgh5bd33=60./r1mgh5bd33;

r111ad39 = xlsread(file,sheet,'b117:ee122');
r111ad39 = r111ad39(:);
r111ad39(isnan(r111ad39))=[];
r111ad39=60./r111ad39;

r1mgh5bd39 = xlsread(file,sheet,'b124:ee129');
r1mgh5bd39 = r1mgh5bd39(:);
r1mgh5bd39(isnan(r1mgh5bd39))=[];
r1mgh5bd39=60./r1mgh5bd39;

r111ad49 = xlsread(file,sheet,'b133:b144');
r111ad49 = r111ad49(:);
r111ad49(isnan(r111ad49))=[];
r111ad49=60./r111ad49;

r1mgh5bd49 = xlsread(file,sheet,'b146:ee157');
r1mgh5bd49 = r1mgh5bd49(:);
r1mgh5bd49(isnan(r1mgh5bd49))=[];
r1mgh5bd49=60./r1mgh5bd49;

r111ad57 = xlsread(file,sheet,'b161:ee172');
r111ad57=r111ad57(:);
r111ad57(isnan(r111ad57))=[];
r111ad57=60./r111ad57;

r1mgh5bd57 = xlsread(file,sheet,'b174:ee185');
r1mgh5bd57 = r1mgh5bd57(:);
r1mgh5bd57(isnan(r1mgh5bd57))=[];
r1mgh5bd57=60./r1mgh5bd57;

r111ad64 = xlsread(file,sheet,'b189:ee197');
r111ad64 = r111ad64(:);
r111ad64(isnan(r111ad64))=[];
r111ad64=60./r111ad64;

r1mgh5bd64 = xlsread(file,sheet,'b200:ee208');
r1mgh5bd64 = r1mgh5bd64(:);
r1mgh5bd64(isnan(r1mgh5bd64))=[];
r1mgh5bd64=60./r1mgh5bd64;

%11a and 19f section (I accidentally named 19f as mgh5b in the code)
sheet = '190314_11a_19f';

r111anomatd18 = xlsread(file,sheet,'b30:ee35');
r111anomatd18 = r111anomatd18(:);
r111anomatd18(isnan(r111anomatd18))=[];
r111anomatd18=60./r111anomatd18;

r111amatd18 = xlsread(file,sheet,'b24:ee28');
r111amatd18 = r111amatd18(:);
r111amatd18(isnan(r111amatd18))=[];
r111amatd18=60./r111amatd18;

r1mgh5bnomatd18 = xlsread(file,sheet,'b37:ee42');
r1mgh5bnomatd18 = r1mgh5bnomatd18(:);
r1mgh5bnomatd18(isnan(r1mgh5bnomatd18))=[];
r1mgh5bnomatd18=60./r1mgh5bnomatd18;

r1mgh5bmatd18 = xlsread(file,sheet,'b44:ee49');
r1mgh5bmatd18 = r1mgh5bmatd18(:);
r1mgh5bmatd18(isnan(r1mgh5bmatd18))=[];
r1mgh5bmatd18=60./r1mgh5bmatd18;

r111amatd26 = xlsread(file,sheet,'b53:ee60');
r111amatd26 = r111amatd26(:);
r111amatd26(isnan(r111amatd26))=[];
r111amatd26=60./r111amatd26;

r111anomatd26 = xlsread(file,sheet,'b62:ee67');
r111anomatd26 = r111anomatd26(:);
r111anomatd26(isnan(r111anomatd26))=[];
r111anomatd26=60./r111anomatd26;

r1mgh5bnomatd26 = xlsread(file,sheet,'b69:ee74');
r1mgh5bnomatd26 = r1mgh5bnomatd26(:);
r1mgh5bnomatd26(isnan(r1mgh5bnomatd26))=[];
r1mgh5bnomatd26=60./r1mgh5bnomatd26;

r1mgh5bmatd26 = xlsread(file,sheet,'b76:ee81');
r1mgh5bmatd26 = r1mgh5bmatd26(:);
r1mgh5bmatd26(isnan(r1mgh5bmatd26))=[];
r1mgh5bmatd26=60./r1mgh5bmatd26;

r111amatd33 = xlsread(file,sheet,'b85:ee90');
r111amatd33 = r111amatd33(:);
r111amatd33(isnan(r111amatd33))=[];
r111amatd33=60./r111amatd33;

r111anomatd33 = xlsread(file,sheet,'b92:ee97');
r111anomatd33 = r111anomatd33(:);
r111anomatd33(isnan(r111anomatd33))=[];
r111anomatd33=60./r111anomatd33;

r1mgh5bnomatd33 = xlsread(file,sheet,'b99:ee104');
r1mgh5bnomatd33 = r1mgh5bnomatd33(:);
r1mgh5bnomatd33(isnan(r1mgh5bnomatd33))=[];
r1mgh5bnomatd33=60./r1mgh5bnomatd33;

r1mgh5bmatd33 = xlsread(file,sheet,'b106:ee111');
r1mgh5bmatd33 = r1mgh5bmatd33(:);
r1mgh5bmatd33(isnan(r1mgh5bmatd33))=[];
r1mgh5bmatd33=60./r1mgh5bmatd33;

r111amatd40 = xlsread(file,sheet,'b115:ee120');
r111amatd40 = r111amatd40(:);
r111amatd40(isnan(r111amatd40))=[];
r111amatd40=60./r111amatd40;

r111anomatd40 = xlsread(file,sheet,'b122:ee127');
r111anomatd40 = r111anomatd40(:);
r111anomatd40(isnan(r111anomatd40))=[];
r111anomatd40=60./r111anomatd40;

r1mgh5bnomatd40 = xlsread(file,sheet,'b129:ee134');
r1mgh5bnomatd40 = r1mgh5bnomatd40(:);
r1mgh5bnomatd40(isnan(r1mgh5bnomatd40))=[];
r1mgh5bnomatd40=60./r1mgh5bnomatd40;

r1mgh5bmatd40 = xlsread(file,sheet,'b136:ee141');
r1mgh5bmatd40 = r1mgh5bmatd40(:);
r1mgh5bmatd40(isnan(r1mgh5bmatd40))=[];
r1mgh5bmatd40=60./r1mgh5bmatd40;

%11a and mgh5b + kynurenic acid section
sheet = '190403_11a_mgh5b';

r111ad19 = xlsread(file,sheet,'b21:ee26');
r111ad19 = r111ad19(:);
r111ad19(isnan(r111ad19))=[];
r111ad19=60./r111ad19;

r1mgh5bd19 = xlsread(file,sheet,'b29:ee34');
r1mgh5bd19 = r1mgh5bd19(:);
r1mgh5bd19(isnan(r1mgh5bd19))=[];
r1mgh5bd19=60./r1mgh5bd19;

r111ad19kyur = xlsread(file,sheet,'b37:ee42');
r111ad19kyur = r111ad19kyur(:);
r111ad19kyur(isnan(r111ad19kyur))=[];
r111ad19kyur=60./r111ad19kyur;

r1mgh5bd19kyur = xlsread(file,sheet,'b44:ee49');
r1mgh5bd19kyur = r1mgh5bd19kyur(:);
r1mgh5bd19kyur(isnan(r1mgh5bd19kyur))=[];
r1mgh5bd19kyur=60./r1mgh5bd19kyur;

%a2 and e4 section
sheet = '190403_a2_e4';

r1a2d19 = xlsread(file,sheet,'b21:ee26');
r1a2d19 = r1a2d19(:);
r1a2d19(isnan(r1a2d19))=[];
r1a2d19=60./r1a2d19;

r1e4d19 = xlsread(file,sheet,'b29:ee34');
r1e4d19 = r1e4d19(:);
r1e4d19(isnan(r1e4d19))=[];
r1e4d19=60./r1e4d19;

r1a2d19kyur = xlsread(file,sheet,'b37:ee42');
r1a2d19kyur = r1a2d19kyur(:);
r1a2d19kyur(isnan(r1a2d19kyur))=[];
r1a2d19kyur=60./r1a2d19kyur;

r1e4d19kyur = xlsread(file,sheet,'b46:ee51');
r1e4d19kyur = r1e4d19kyur(:);
r1e4d19kyur(isnan(r1e4d19kyur))=[];
r1e4d19kyur = 60./r1e4d19kyur;


%TDP43 section
sheet = '190424_TDP43_glutamate';

tdpwtglutlong = xlsread(file,sheet,'b2:ee13');
tdpwtglutlong = tdpwtglutlong(:);
tdpwtglutlong(isnan(tdpwtglutlong))=[];
tdpwtglutlong =60./tdpwtglutlong;

tdpmutglutlong = ieitransf(file,sheet,'b16:ee27');


tdpwtglutshort = ieitransf(file,sheet,'b30:ee41');
tdpmutglutshort = ieitransf(file,sheet,'b44:ee55');

tdpwtsodiumars = ieitransf(file,sheet,'b58:ee69');
tdpmutsodiumars = ieitransf(file,sheet,'b72:ee83');

tdpwtctrl = ieitransf(file,sheet,'b86:ee97');
tdpmutctrl = ieitransf(file,sheet,'b100:ee111');

%% Section comparing different replicates
datacolumns = {r111ad20, r111ad19,r111anomatd18}
proximity = 0.1;
lim = 8;
points = [0 2 4 6 8];
pointlabels = {'0','2','4','6','8'};

beeswarmbar3(datacolumns,proximity,lim,points,pointlabels)

%% section looking at effect of matrigel on iPSC-MN activity over time
datacolumns = {r111anomatd18, r111anomatd26, r111anomatd33, r111anomatd40};
beeswarmbar3(datacolumns,proximity,lim,points,pointlabels)

datacolumns = {r1mgh5bnomatd18, r1mgh5bnomatd26, r1mgh5bnomatd33, r1mgh5bnomatd40};
beeswarmbar3(datacolumns,proximity,lim,points,pointlabels)

datacolumns = {r111amatd18, r111amatd26, r111amatd33, r111amatd40};
beeswarmbar3(datacolumns,proximity,lim,points,pointlabels)

datacolumns = {r1mgh5bmatd18, r1mgh5bmatd26, r1mgh5bmatd33, r1mgh5bmatd40};
beeswarmbar3(datacolumns,proximity,lim,points,pointlabels)

data = {r111anomatd18, r111amatd18; r111anomatd26, r111amatd26; r111anomatd33, r111amatd33; r111anomatd40, r111amatd40}; 

xpoints = [18 26 33 40];
linecolors = {[0 0 0],[1 0 0]};
shadecolors = {[.5 .5 .5],[.5 0 0]};
ypoints = [0 2 4 6 8];
xlimit = [0 40];
ylimit = [0 4];

repeatplotter(data,xpoints, xlimit, ypoints, ylimit, linecolors, shadecolors);

data = {r1mgh5bmatd18, r1mgh5bnomatd18; r1mgh5bmatd26, r1mgh5bnomatd26; r1mgh5bmatd33, r1mgh5bnomatd33; r1mgh5bmatd40, r1mgh5bnomatd40};

repeatplotter(data,xpoints, xlimit, ypoints, ylimit, linecolors, shadecolors);

% section looking at network vs intrinsic
datacolumns = {r111ad19, r111ad19kyur}

proximity = 0.08;
lim = 8;
points = [0 2 4 6 8];
pointlabels = {'0','2','4','6','8'};

beeswarmbar3(datacolumns,proximity,lim,points,pointlabels)

%% section looking at different genotypes
%11a vs 19f
data = {r111amatd18,r1mgh5bmatd18; r111amatd26, r1mgh5bmatd26; r111amatd33, r1mgh5bmatd33; r111amatd40, r1mgh5bmatd40};

xpoints = [18 26 33 40];
linecolors = {[0 0 0],[1 0 0]};
shadecolors = {[.5 .5 .5],[.5 0 0]};
ypoints = [0 2 4 6 8];
xlimit = [0 40];
ylimit = [0 4];

repeatplotter(data,xpoints, xlimit, ypoints, ylimit, linecolors, shadecolors);

%11a vs mgh5b
data = {r111ad11, r1mgh5bd11; r111ad15, r1mgh5bd15; r111ad27, r1mgh5bd27; r111ad33, r1mgh5bd33; r111ad39, r1mgh5bd39; r111ad49, r1mgh5bd49; r111ad57, r1mgh5bd57; r111ad64, r1mgh5bd64};

xlimit = [0 64];
xpoints = [11 15 27 33 39 49 57 64];
ylimit = [0 6];
ypoints = [0 2 4 6];
linecolors = {[0 0 0],[1 0 0]};
shadecolors = {[.5 .5 .5],[.5 0 0]};

repeatplotter(data,xpoints, xlimit, ypoints, ylimit, linecolors, shadecolors);

% A2 vs E4
datacolumns = {r1a2d19, r1e4d19}
proximity = 0.08;
lim = 8;
points = [0 2 4 6 8];
pointlabels = {'0','2','4','6','8'};

beeswarmbar3(datacolumns,proximity,lim,points,pointlabels)

%11a over time
datacolumns = {r111ad11,r111ad15,r111ad27, r111ad33, r111ad39, r111ad49,r111ad57};
proximity = 0.08;
lim = 8;
points = [0 2 4 6 8];
pointlabels = {'0','2','4','6','8'};
beeswarmbar3(datacolumns,proximity,lim,points,pointlabels)

%mgh5b over time
datacolumns = {r1mgh5bd11,r1mgh5bd15,r1mgh5bd27,r1mgh5bd33,r1mgh5bd39,r1mgh5bd49,r1mgh5bd57}
beeswarmbar3(datacolumns,proximity,lim,points,pointlabels)

%tdp43wt control, mut control, tdp43wt glutshort, tdp43mut glutshort,
%tdp43wt glutlong, tdp43mut glutlong, tdp43wt naars, tdp43mut naars
datacolumns = {tdpwtctrl, tdpmutctrl, tdpwtglutshort, tdpmutglutshort, tdpwtglutlong, tdpmutglutlong, tdpwtsodiumars, tdpmutsodiumars};
proximity = 0.08;
lim = 8;
points = [0 2 4 6 8];
pointlabels = {'0','2','4','6','8'};
beeswarmbar3(datacolumns,proximity,lim,points,pointlabels)






