function [] = AnimationGen_Experimental(tn,sn)
% ---------------------------- Begin Header ---------------------------- %
% FILE:     StrainAnalysis.m
% AUTHOR:   Darian Emmett
% DATE:   07/06/2020
% 
% PURPOSE: This function takes the pos data from each trial and generates
% an animation of the markers moving in 3d from muliple perspectives
%
%
% INPUTS: tn (trial name, written exactly as the name of the corresponding file),
% sn (subject number, should be formatted 'Sxx', where 'S01' corresponds to subject 1),
% a MATLAB data file (ext .mat) with the same name as tn
% 
% 
% OUTPUTS: External .avi file that animates the strain calculations at each moment of the trial
%
% NOTES: 
%
%
% VERSION HISTORY
% V1 - Original.
% ----------------------------- End Header ----------------------------- %

%% Data Organization %%
name={strcat(tn,'.mat'),tn};
input=load(name{1},name{2});

Fs=input.(name{2}).FrameRate; % Sample Frequency %
ns=input.(name{2}).Frames; % total number of data points
nm=length(input.(name{2}).Trajectories.Labeled.Labels); % total number of markers
labels=input.(name{2}).Trajectories.Labeled.Labels;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Marker Layout %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%                         T11L(34)              T11C(35)              T11R(36)                        %
%                         T12L(31)              T12C(32)              T12R(33)                        %
%                         L1L(28)               L1C(29)               L1R(30)                         %
%                         L2L(25)               L2C(26)               L2R(27)                         %
%                         L3L(22)               L3C(23)               L3R(24)                         %
%                         L4L(19)               L4C(20)               L4R(21)                         %
% L5L4(10)     L5L3(11)   L5L2(12)   L5L1(13)   L5C(14)    L5R1(15)   L5R2(16)   L5R3(17)   L5R4(18)  %
% PSISL3(1)    PSISL2(2)  PSISL1(3)  S2L(4)     S2C(5)     S2R(6)     PSISR1(7)  PSISR2(8)  PSISR3(9) %

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% If some markers were not reported, notify and quit
if length(labels)~=36
    return
end
% Create matrix to hold every marker's location organized according to above layout
fd=16; % Window for filter
BIG=nan(3,ns,nm); % initialize the matrix to hold all of the markers, ordered.
Tran=nan(3,ns,nm); % initialize the matrix used to hold the markers after rotational transformation
h = fir1(fd, 6/(Fs/2)); % N is filter length, Fs is sampling frequency (Hz)
% Name the markers
%           1        2        3      4     5     6      7       8        9       10     11     12     13    14     15     16     17     18    19    20    21    22    23    24    25    26    27    28    29    30     31     32     33     34     35     36
markers={'PSISL3','PSISL2','PSISL','S2L','S2C','S2R','PSISR','PSISR2','PSISR3','L5L4','L5L3','L5L2','L5L1','L5C','L5R1','L5R2','L5R3','L5R4','L4L','L4C','L4R','L3L','L3C','L3R','L2L','L2C','L2R','L1L','L1C','L1R','T12L','T12C','T12R','T11L','T11C','T11R'};
% Populate the ordered matrix of marker xyz values, using the names of each marker to index througth the names of the markers 
% in the data and reorder them in BIG
for i=1:nm
    BIG(:,:,i)=input.(name{2}).Trajectories.Labeled.Data(strcmp(labels,(markers{i})),1:3,:);
    BIG(1,:,i)=filter(h, 1, BIG(1,:,i)); % Filter x-coord
    BIG(2,:,i)=filter(h, 1, BIG(2,:,i)); % Filter y-coord
    BIG(3,:,i)=filter(h, 1, BIG(3,:,i)); % Filter z-coord
end
% Cut out area messed up by filter and transpose the BIG matrix
BIG=BIG(:,fd+1:end,:); % Cut off filter problem zone
ns=ns-fd; % Adjust number of samples to account for the cut-off area due to filter window
BIG=BIG-BIG(:,:,1);
clear input;

% % % % Generate Video/Animation % % % %
X=nan(nm,3);
X2=nan(nm,3);

set(0, 'DefaulttextInterpreter', 'none'); % Prevent subscript generation after _

myVideo = VideoWriter(['Anim2_' sn '_' tn]); % open/create video file
myVideo.FrameRate = Fs;  % set frame rate
open(myVideo); % open the video file for writing

for i=1:ns
    uy=(BIG(:,i,7)-BIG(:,i,3))/norm(BIG(:,i,7)-BIG(:,i,3)); % Get the unit vector from PSISL1 to PSISR1 for each step
    uz=(BIG(:,i,14)-BIG(:,i,5))/norm(BIG(:,i,14)-BIG(:,i,5)); % Get the unit vector from S2 to L5 for each step
    ux=cross(uy,uz);
    uy=cross(uz,ux);
    for im=1:nm
        Tran(1,i,im)=dot(BIG(:,i,im),ux);
        Tran(2,i,im)=dot(BIG(:,i,im),uy);
        Tran(3,i,im)=dot(BIG(:,i,im),uz);
    end
end

n=[1,0,0];
v=[0,0,1];
b=[0,1,0];

v2=(BIG(:,1,14)-BIG(:,1,5))'/norm((BIG(:,1,14)-BIG(:,1,5))'); % Get vertical vector pointing for S2-L5 for non-transformed data
v2=(v2+[.01,0,0])/norm((v2+[.01,0,0])); % Make c2 a unit vector
n2=cross((BIG(:,1,7)-BIG(:,1,3))',(BIG(:,1,14)-BIG(:,1,5))'); % Get normal vector for S2-L5 for non-transformed data
n2=n2/norm(n2); % Make n2 a unit vector
b2=cross(n2,v2); % get binormal vector for S2-L5 (normal to both vertical and normal vector) for non-transformed data
orth=(n2+v2+b2)/norm(n2+v2+b2); % Get orthogonal vector for non-transformed data

for is=1:ns
    for i=1:nm
        X(i,1)=Tran(1,is,i);
        X(i,2)=Tran(2,is,i);
        X(i,3)=Tran(3,is,i);
        X2(i,1)=BIG(1,is,i);
        X2(i,2)=BIG(2,is,i);
        X2(i,3)=BIG(3,is,i);
    end
    if is==1
        % Top Left plot: Top View
        subplot(2,2,1);
        TV=scatter3(X(:,1),X(:,2),X(:,3),10,'filled');
        title(strcat(sn,'_',tn," Top View"),'FontSize',8);
        axis equal off
        view(v)
        
        % Top Right plot: Ortho View
        subplot(2,2,2);
        OV=scatter3(X2(:,1),X2(:,2),X2(:,3),10,'filled');
        title(strcat(sn,'_',tn," Ortho View"),'FontSize',8);
        axis equal off
        view(orth)
        
        % Bottom Left plot: Back View
        subplot(2,2,3);
        BV=scatter3(X(:,1),X(:,2),X(:,3),10,'filled');
        title(strcat(sn,'_',tn," Back View"),'FontSize',8);
        axis equal off
        view(n)
        
        % Bottom Right plot: Side View (From Right)
        subplot(2,2,4);
        SV=scatter3(X(:,1),X(:,2),X(:,3),10,'filled');
        title(strcat(sn,'_',tn," Side View"),'FontSize',8);
        axis equal off
        view(b)
        
        pause(1/Fs)
    else
        TV.XData=X(:,1); TV.YData=X(:,2); TV.ZData=X(:,3);
        OV.XData=X2(:,1); OV.YData=X2(:,2); OV.ZData=X2(:,3);
        BV.XData=X(:,1); BV.YData=X(:,2); BV.ZData=X(:,3);
        SV.XData=X(:,1); SV.YData=X(:,2); SV.ZData=X(:,3);
        pause(1/Fs)
    end
    set(gcf, 'Units', 'centimeters', 'Position', [0 0 20 20]); % Set the anchor for the figure
    set(gcf, 'PaperUnits', 'centimeters', 'PaperSize', [20 20]); % Set the size of the figure
    frame = getframe(gcf); %get frame
    writeVideo(myVideo, frame);
end
close(myVideo);