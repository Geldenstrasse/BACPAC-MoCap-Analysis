function [] = AnimationGen(tn,sn)
% ---------------------------- Begin Header ---------------------------- %
% FILE:     StrainAnalysis.m
% AUTHOR:   Darian Emmett
% DATE:   07/06/2020
% 
% PURPOSE: This function takes the pos data from each trial and generates
% an animation of the markers moving in 3D from muliple perspectives
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

%% Data Import and Organization %%
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
nc=nm*3; % gives the number of columns in BIG (the total number of markers x3)
BIG=nan(nc,ns); % initialize the matrix to hold all of the markers, ordered.
h = fir1(fd, 6/(Fs/2)); % N is filter length, Fs is sampling frequency (Hz)

% Name the markers
%           1        2        3      4     5     6      7       8        9       10     11     12     13    14     15     16     17     18    19    20    21    22    23    24    25    26    27    28    29    30     31     32     33     34     35     36
markers={'PSISL3','PSISL2','PSISL','S2L','S2C','S2R','PSISR','PSISR2','PSISR3','L5L4','L5L3','L5L2','L5L1','L5C','L5R1','L5R2','L5R3','L5R4','L4L','L4C','L4R','L3L','L3C','L3R','L2L','L2C','L2R','L1L','L1C','L1R','T12L','T12C','T12R','T11L','T11C','T11R'};

% Populate the ordered matrix of marker xyz values, using the names of each marker to index through the names of the markers 
% in the data and reorder them in BIG
for i=1:nm
    BIG(3*i-2:3*i,:)=input.(name{2}).Trajectories.Labeled.Data(strcmp(labels,(markers{i})),1:3,:);
    BIG(3*i-2,:)=filter(h, 1, BIG(3*i-2,:)); % Filter x-coord
    BIG(3*i-1,:)=filter(h, 1, BIG(3*i-1,:)); % Filter y-coord
    BIG(3*i,:)=filter(h, 1, BIG(3*i,:)); % Filter z-coord
end

clear input; % Get rid of unused inputs from the .mat file

% Cut out area messed up by filter and transpose the BIG matrix
BIG=BIG(:,fd+1:end)'; % The data comes in with the samples/frames on the horizontal axis, this adjusts to have the markers specified by column, the sample specified by row
counter=1:ns;
BIG=BIG(min(counter(not(isnan(mean(BIG,2))))):end,:); % Cut off any NaN entries in the start
% ^ This works by finding the minimum time step (housed in 'counter') at
% which BIG does not have a NaN value for any marker positional coordinate
ns=size(BIG,1); % Adjust number of samples to account for the cut off frames

% % % % Generate Video/Animation % % % %
set(0, 'DefaulttextInterpreter', 'none'); % Prevent subscript generation after _

myVideo = VideoWriter(['Anim_' sn '_' tn]); % open/create video file
myVideo.FrameRate = Fs;  % set frame rate
open(myVideo); % open the video file for writing

X=nan(nm,3); % Initialize a matric for holding the current time-step of the marker positions

for is=1:ns % Do this for every time step
    for i=1:nm % Do this for every marker
        X(i,1)=BIG(is,3*i-2); % Fill marker 'x' coord with x data from BIG for current marker
        X(i,2)=BIG(is,3*i-1); % Fill marker 'y' coord with y data from BIG for current marker
        X(i,3)=BIG(is,3*i); % Fill marker 'z' coord with z data from BIG for current marker
    end
    n=cross((X(7,:)-X(3,:)),(X(14,:)-X(5,:))); % Get normal vector for S2-L5 segment
    n=n/norm(n); % Make n a unit vector
    v=(X(14,:)-X(5,:))/norm((X(14,:)-X(5,:))); % Get vertical vector for S2-L5
    v=v/norm(v); % Make c a unit vector
    b=cross(n,v); % Get binormal vector for S2-L5 (normal to both vertical and normal vector)
    orth=(n+v+b)/norm(n+v+b); % Get orthogonal view vector (average of the normal, vertical, and binormal)
    if is==1 % Only do the following for the first time-step
        subplot(2,2,1); % In the top-left subplot:
        TV=scatter3(X(:,1),X(:,2),X(:,3),10,'filled'); % Plot the marker positions in 3D space
        title(strcat(sn,'_',tn," Top View"),'FontSize',8); % Title the top view plot
        axis equal off % Remove visual axes and ensure equal axis tick spacing
        view(v) % Set view to be from above
        
        subplot(2,2,2); % In the top-right subplot:
        OV=scatter3(X(:,1),X(:,2),X(:,3),10,'filled'); % Plot the marker positions in 3D space
        title(strcat(sn,'_',tn," Ortho View"),'FontSize',8); % Title the ortho view plot
        axis equal off % Remove visual axes and ensure equal axis tick spacing
        view(orth) % Set view to be from ortho perspective
        
        subplot(2,2,3);
        BV=scatter3(X(:,1),X(:,2),X(:,3),10,'filled'); % Plot the marker positions in 3D space
        title(strcat(sn,'_',tn," Back View"),'FontSize',8); % Title the back view plot
        axis equal off % Remove visual axes and ensure equal axis tick spacing
        view(n) % Set view to be from the back
        
        subplot(2,2,4);
        SV=scatter3(X(:,1),X(:,2),X(:,3),10,'filled'); % Plot the marker positions in 3D space
        title(strcat(sn,'_',tn," Side View"),'FontSize',8); % Title the side view plot
        axis equal off % Remove visual axes and ensure equal axis tick spacing
        view(b) % Set view to be from the side
        
        pause(1/Fs) % Pause to give framerate to video
    else % For all frames after the first time-step
        TV.XData=X(:,1); TV.YData=X(:,2); TV.ZData=X(:,3); % Update the top-view plot
        OV.XData=X(:,1); OV.YData=X(:,2); OV.ZData=X(:,3); % Update the ortho-view plot
        BV.XData=X(:,1); BV.YData=X(:,2); BV.ZData=X(:,3); % Update the back-view plot
        SV.XData=X(:,1); SV.YData=X(:,2); SV.ZData=X(:,3); % Update the side-view plot
        pause(1/Fs) % Pause to give framerate to video
    end
    set(gcf, 'Units', 'centimeters', 'Position', [0 0 20 20]); % Set the position of the figure
    set(gcf, 'PaperUnits', 'centimeters', 'PaperSize', [20 20]); % Set the size of the figure
    frame = getframe(gcf); % Get frame
    writeVideo(myVideo, frame); % Write the frame to the video
end % End time-step loop

close(myVideo); % Close the video once finished saving

end % End the Function