function [romL,romT] = RoMK_F(tn,sn)
%%%% SUMMARY %%%%
% RoMK_F reads Qualisys mocap data for a set of 36 markers, laid out as
% depicted below, then uses the positional information of the markers to
% calculate angular position of segments of the back, relative to the 
% lowest segment, at each time step.

%%%% OUTPUTS %%%%
% Lum = Time-series of angular position (in degrees) for gross lumbar sement
% Tho = Time-series of angular position (in degrees) for gross thoracic segment
% romL = Max/min for each angle (in degrees) of the gross lumbar segment. Type: Double Array
% romT = Max/min for each angle (in degrees) of the gross thoracic segment. Type: Double Array
% Pos = nx3x6 matrix (where n= number of frames) holding the time-series for the
%       3-component angular position (in radians) of each spine segment.
%       1st index = frame, 2nd = angle (1=forward, 2=side, 3=twist),
%       3rd = segment (1=L5->L4, ... , 6=T12->T11)
% Vel = nx3x6 matrix holding the time-series for the 3-component angular velocity
%       (in rad/sec) of each spine segment. Formatted the same as Pos
% Acc = nx3x6 matrix holding the time-series for the 3-component angular
%       acceleration (in rad/s^2) of each spine segment. Formatted the same as Pos

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

%%%% Data Read-in %%%%
name={strcat(tn,'.mat'),tn}; % MoCap data Filename and main directory
input=load(name{1},name{2}); % Load MoCap data

Fs=input.(name{2}).FrameRate; % Pull Sample Frequency
ns=input.(name{2}).Frames; % Pull total number of data points
nm=length(input.(name{2}).Trajectories.Labeled.Labels); % Pull total number of markers
labels=input.(name{2}).Trajectories.Labeled.Labels; % Extract marker labels

% If some markers were not reported, notify and quit
if length(input.(name{2}).Trajectories.Labeled.Labels)~=36
    romL=0; % Report output variables as zeros (to notify program using function that it went wrong)
    romT=0;
    return
end

% Create matrix to hold every marker's location organized according to above layout
BIG=nan(nm,3,ns); % initialize the matrix to hold all of the markers' data, ordered.
% Name the markers
%           1        2        3      4     5     6      7       8        9       10     11     12     13    14     15     16     17     18    19    20    21    22    23    24    25    26    27    28    29    30     31     32     33     34     35     36
markers={'PSISL3','PSISL2','PSISL','S2L','S2C','S2R','PSISR','PSISR2','PSISR3','L5L4','L5L3','L5L2','L5L1','L5C','L5R1','L5R2','L5R3','L5R4','L4L','L4C','L4R','L3L','L3C','L3R','L2L','L2C','L2R','L1L','L1C','L1R','T12L','T12C','T12R','T11L','T11C','T11R'};
% Populate the ordered matrix of marker xyz values, using the names of each marker to index througth the names of the markers 
% in the data and reorder them in BIG
[b,a] = butter(6, 6/(Fs/2)); % Build a 6th order butterworth filter with a cutoff frequecy of 6 Hz
for i=1:nm % iterate through the markers (in the desired order)
    BIG(i,:,:)=input.(name{2}).Trajectories.Labeled.Data(strcmp(labels,(markers{i})),1:3,:); % Copy marker trajectory into BIG
    BIG(i,1,:) = filter(b, a, BIG(i,1,:)); % Filter x-coord
    BIG(i,2,:) = filter(b, a, BIG(i,2,:)); % Filter y-coord
    BIG(i,3,:) = filter(b, a, BIG(i,3,:)); % Filter z-coord
end
clear input; % Get rid of 'input' variables, since relevant data has been copied over already.

%%%% Initialize variables for RoM calculations %%%%
SC=nan(3,ns); SN=nan(3,ns); SB=nan(3,ns); % SC, SN, & SB are the coincident, normal, and binormal vectors for the sacral segment (S2-L5)
LC=nan(3,ns); LN=nan(3,ns); LB=nan(3,ns); % LC, LN, & LB are the coincident, normal, and binormal vectors for the lumbar segment (L5-L1)
L5C=nan(3,ns); L5N=nan(3,ns); L5B=nan(3,ns); % Coincident, normal, and binormal vectors for L5-L4
L4C=nan(3,ns); L4N=nan(3,ns); L4B=nan(3,ns); % Coincident, normal, and binormal vectors for L4-L3
L3C=nan(3,ns); L3N=nan(3,ns); L3B=nan(3,ns); % Coincident, normal, and binormal vectors for L3-L2
L2C=nan(3,ns); L2N=nan(3,ns); L2B=nan(3,ns); % Coincident, normal, and binormal vectors for L2-L1
L1C=nan(3,ns); L1N=nan(3,ns); L1B=nan(3,ns); % Coincident, normal, and binormal vectors for L1-T12
TC=nan(3,ns); TN=nan(3,ns); TB=nan(3,ns); % TC, TN, & TB are the coincident, normal, and binormal vectors for the thoracic segment (T12-T11)
Pos=nan(ns,3,6); % Holds pitch, yaw, and roll for all functional segments (bottom to top) at each frame
                 % ^ The dimensions are (number of samples, degrees of freedom, segment (1=L5, 2=L4, 3=L3, etc.) 
Lum=nan(ns,3); % Holds pitch, yaw, and roll for lumbar spine at each frame
Tho=nan(ns,3); % Holds pitch, yaw, and roll for thoracic spine at each frame

%%%% Calculate angles for RoM %%%%
for i=1:ns
    % Get descriptive vectors for S2
    SC(:,i)=BIG(14,:,i)-BIG(5,:,i); % Get the vector going from S2 to L5
    SC(:,i)=SC(:,i)/norm(SC(:,i)); % Make the above vector a unit vector
    SN(:,i)=cross(BIG(7,:,i)-BIG(3,:,i),SC(:,i)); % Get normal vector for S2
    SN(:,i)=SN(:,i)/norm(SN(:,i)); % Make sacral normal vector a unit vector
    SB(:,i)=cross(SC(:,i),SN(:,i)); % Get binormal vector for S2
    SB(:,i)=SB(:,i)/norm(SB(:,i)); % Make sacral binormal vector a unit vector
    % Get descriptive vectors for L5
    L5C(:,i)=BIG(20,:,i)-BIG(14,:,i); % Get vector going from L5 to L4
    L5C(:,i)=L5C(:,i)/norm(L5C(:,i)); % Normalize coincident vector
    L5N(:,i)=cross(BIG(16,:,i)-BIG(12,:,i),L5C(:,i)); % Get normal vector
    L5N(:,i)=L5N(:,i)/norm(L5N(:,i)); % Convert normal vector to unit normal
    L5B(:,i)=cross(L5C(:,i),L5N(:,i)); % Get binormal vector for L5
    L5B(:,i)=L5B(:,i)/norm(L5B(:,i)); % Make L5 binormal vector a unit vector
    % Get descriptive vectors for L4
    L4C(:,i)=BIG(23,:,i)-BIG(20,:,i); % Get vector going from L4 to L3
    L4C(:,i)=L4C(:,i)/norm(L4C(:,i)); % Normalize parallel vector
    L4N(:,i)=cross(BIG(21,:,i)-BIG(19,:,i),L4C(:,i)); % Get normal vector
    L4N(:,i)=L4N(:,i)/norm(L4N(:,i)); % Convert normal vector to unit normal
    L4B(:,i)=cross(L4C(:,i),L4N(:,i)); % Get binormal vector
    L4B(:,i)=L4B(:,i)/norm(L4B(:,i)); % Make binormal vector a unit vector
    % Get descriptive vectors for L3
    L3C(:,i)=BIG(26,:,i)-BIG(23,:,i); % Get vector going from L3 to L2
    L3C(:,i)=L3C(:,i)/norm(L3C(:,i)); % Normalize coincident vector
    L3N(:,i)=cross(BIG(24,:,i)-BIG(22,:,i),L3C(:,i)); % Get normal vector
    L3N(:,i)=L3N(:,i)/norm(L3N(:,i)); % Convert normal vector to unit normal
    L3B(:,i)=cross(L3C(:,i),L3N(:,i)); % Get binormal vector
    L3B(:,i)=L3B(:,i)/norm(L3B(:,i)); % Make binormal vector a unit vector
    % Get descriptive vectors for L2
    L2C(:,i)=BIG(29,:,i)-BIG(26,:,i); % Get vector going from L2 to L1
    L2C(:,i)=L2C(:,i)/norm(L2C(:,i)); % Normalize coincident vector
    L2N(:,i)=cross(BIG(27,:,i)-BIG(25,:,i),L2C(:,i)); % Get normal vector
    L2N(:,i)=L2N(:,i)/norm(L2N(:,i)); % Convert normal vector to unit normal
    L2B(:,i)=cross(L2C(:,i),L2N(:,i)); % Get binormal vector
    L2B(:,i)=L2B(:,i)/norm(L2B(:,i)); % Make binormal vector a unit vector
    % Get descriptive vectors for L1
    L1C(:,i)=BIG(32,:,i)-BIG(29,:,i); % Get vector going from L1 to T11
    L1C(:,i)=L1C(:,i)/norm(L1C(:,i)); % Normalize coincident vector
    L1N(:,i)=cross(BIG(30,:,i)-BIG(28,:,i),L1C(:,i)); % Get normal vector
    L1N(:,i)=L1N(:,i)/norm(L1N(:,i)); % Convert normal vector to unit normal
    L1B(:,i)=cross(L1C(:,i),L1N(:,i)); % Get binormal vector
    L1B(:,i)=L1B(:,i)/norm(L1B(:,i)); % Make binormal vector a unit vector
    % Get descriptive vectors for Thoracic spine
    TC(:,i)=BIG(35,:,i)-BIG(32,:,i); % Get vector going from T12 to T11
    TC(:,i)=TC(:,i)/norm(TC(:,i)); % Normalize coincident vector
    TN(:,i)=cross(BIG(36,:,i)-BIG(34,:,i),TC(:,i)); % Get normal vector
    TN(:,i)=TN(:,i)/norm(TN(:,i));
    TB(:,i)=cross(TC(:,i),TN(:,i));
    TB(:,i)=TB(:,i)/norm(TB(:,i));
    % get descriptive vectors for Lumbar spine as one segment
    LC(:,i)=BIG(29,:,i)-BIG(14,:,i); % Get vector going through lumbar vertebrae
    LC(:,i)=LC(:,i)/norm(LC(:,i)); % Normalize parallel vector
    LNB=cross(BIG(16,:,i)-BIG(12,:,i),LC(:,i)); % Get normal vector using bottom plane bias
    LNT=cross(BIG(30,:,i)-BIG(28,:,i),LC(:,i)); % Get normal vector using top plane bias
    LN(:,i)=(LNB+LNT); % Get approximate normal vector for lumbar back plane
    LN(:,i)=LN(:,i)/norm(LN(:,i)); % Convert normal vector to unit normal
    LB(:,i)=cross(LC(:,i),LN(:,i)); % Get binormal vector for lumbar segment
    LB(:,i)=LB(:,i)/norm(LB(:,i)); % Make lumbar binormal vector a unit vector
    % Get angles for lumbar and thoracic segments of the back, relative to the sacral spine
    Lum(i,1)=-atan(dot(LC(:,i),SN(:,i))/dot(LC(:,i),SC(:,i))); % project LC into plane [SN,SC], calculate angle relative to SC
    Lum(i,2)=-atan(dot(LB(:,i),SC(:,i))/dot(LB(:,i),SB(:,i))); % project LB into plane [SC,SB], calculate angle relative to SB
    Lum(i,3)=atan(dot(LB(:,i),SN(:,i))/dot(LB(:,i),SB(:,i))); % project LB into plane [SN,SB], calculate angle relative to SN
    Tho(i,1)=-atan(dot(TC(:,i),SN(:,i))/dot(TC(:,i),SC(:,i))); % project TC into plane [SN,SC], calculate angle relative to SC
    Tho(i,2)=-atan(dot(TB(:,i),SC(:,i))/dot(TB(:,i),SB(:,i))); % project TB into plane [SC,SB], calculate angle relative to SB
    Tho(i,3)=atan(dot(TB(:,i),SN(:,i))/dot(TB(:,i),SB(:,i))); % project TB into plane [SN,SC], calculate angle relative to SB
    % Get angles of each segment of the spine at each moment, relative to S2-L5 segment
    Pos(i,1,1)=-atan(dot(L5C(:,i),SN(:,i))/dot(L5C(:,i),SC(:,i))); % Flex/Ext angle for L5-L4 segment. '-' is to make Forward positive
    Pos(i,2,1)=-atan(dot(L5B(:,i),SC(:,i))/dot(L5B(:,i),SB(:,i))); % Lateral angle for L5-L4 segment. '-' is to make Right positive
    Pos(i,3,1)=atan(dot(L5B(:,i),SN(:,i))/dot(L5B(:,i),SB(:,i))); % Axial angle for L5-L4 segment. No '-' since Right is already positive
    Pos(i,1,2)=-atan(dot(L4C(:,i),SN(:,i))/dot(L4C(:,i),SC(:,i))); % Flex/Ext angle for L4-L3 segment
    Pos(i,2,2)=-atan(dot(L4B(:,i),SC(:,i))/dot(L4B(:,i),SB(:,i))); % Lateral angle for L4-L3 segment
    Pos(i,3,2)=atan(dot(L4B(:,i),SN(:,i))/dot(L4B(:,i),SB(:,i))); % Axial angle for L4-L3 segment
    Pos(i,1,3)=-atan(dot(L3C(:,i),SN(:,i))/dot(L3C(:,i),SC(:,i))); % Flex/Ext angle for L3-L2 segment
    Pos(i,2,3)=-atan(dot(L3B(:,i),SC(:,i))/dot(L3B(:,i),SB(:,i))); % Lateral angle for L3-L2 segment
    Pos(i,3,3)=atan(dot(L3B(:,i),SN(:,i))/dot(L3B(:,i),SB(:,i))); % Axial angle for L3-L2 segment
    Pos(i,1,4)=-atan(dot(L2C(:,i),SN(:,i))/dot(L2C(:,i),SC(:,i))); % Flex/Ext angle for L2-L1 segment
    Pos(i,2,4)=-atan(dot(L2B(:,i),SC(:,i))/dot(L2B(:,i),SB(:,i))); % Lateral angle for L2-L1 segment
    Pos(i,3,4)=atan(dot(L2B(:,i),SN(:,i))/dot(L2B(:,i),SB(:,i))); % Axial angle for L2-L1 segment
    Pos(i,1,5)=-atan(dot(L1C(:,i),SN(:,i))/dot(L1C(:,i),SC(:,i))); % Flex/Ext angle for L1-T12 segment
    Pos(i,2,5)=-atan(dot(L1B(:,i),SC(:,i))/dot(L1B(:,i),SB(:,i))); % Lateral angle for L1-T12 segment
    Pos(i,3,5)=atan(dot(L1B(:,i),SN(:,i))/dot(L1B(:,i),SB(:,i))); % Axial angle for L1-T12 segment
    Pos(i,:,6)=Tho(i,:); % Copy thoracic segment (T12-T11) angles over (already calculated for Tho)
end

% Normalize position vector based on relaxed position
Pos = (Pos-Pos(1,:,:));
% Get velocity and acceleration vectors
Vel=nan(ns,3,6);
Vel(:,:,1)= [gradient(Pos(:,1,1),1/Fs),gradient(Pos(:,2,1),1/Fs),gradient(Pos(:,3,1),1/Fs)];
Vel(:,:,2)= [gradient(Pos(:,1,2),1/Fs),gradient(Pos(:,2,2),1/Fs),gradient(Pos(:,3,2),1/Fs)];
Vel(:,:,3)= [gradient(Pos(:,1,3),1/Fs),gradient(Pos(:,2,3),1/Fs),gradient(Pos(:,3,3),1/Fs)];
Vel(:,:,4)= [gradient(Pos(:,1,4),1/Fs),gradient(Pos(:,2,4),1/Fs),gradient(Pos(:,3,4),1/Fs)];
Vel(:,:,5)= [gradient(Pos(:,1,5),1/Fs),gradient(Pos(:,2,5),1/Fs),gradient(Pos(:,3,5),1/Fs)];
Vel(:,:,6)= [gradient(Pos(:,1,6),1/Fs),gradient(Pos(:,2,6),1/Fs),gradient(Pos(:,3,6),1/Fs)];
Acc=nan(ns,3,6);
Acc(:,:,1)= [gradient(Vel(:,1,1),1/Fs),gradient(Vel(:,2,1),1/Fs),gradient(Vel(:,3,1),1/Fs)];
Acc(:,:,2)= [gradient(Vel(:,1,2),1/Fs),gradient(Vel(:,2,2),1/Fs),gradient(Vel(:,3,2),1/Fs)];
Acc(:,:,3)= [gradient(Vel(:,1,3),1/Fs),gradient(Vel(:,2,3),1/Fs),gradient(Vel(:,3,3),1/Fs)];
Acc(:,:,4)= [gradient(Vel(:,1,4),1/Fs),gradient(Vel(:,2,4),1/Fs),gradient(Vel(:,3,4),1/Fs)];
Acc(:,:,5)= [gradient(Vel(:,1,5),1/Fs),gradient(Vel(:,2,5),1/Fs),gradient(Vel(:,3,5),1/Fs)];
Acc(:,:,6)= [gradient(Vel(:,1,6),1/Fs),gradient(Vel(:,2,6),1/Fs),gradient(Vel(:,3,6),1/Fs)];

% Normalize Lum and Tho angles to be w/ respect to inital position
Lum=(Lum-Lum(1,:))*180/pi();
Tho=(Tho-Tho(1,:))*180/pi();

[~,sm]=max(abs(Pos),[],'all','linear'); % Get frame at which max deflection is reached
[sm,~,~] = ind2sub(size(Pos), sm); % [time step, angle, segment]

romL(1,:)=[min(Lum(:,1)),max(Lum(:,1))]; % Lumbar Pitch
romL(2,:)=[min(Lum(:,2)),max(Lum(:,2))]; % Lumbar Yaw
romL(3,:)=[min(Lum(:,3)),max(Lum(:,3))]; % Lumbar Roll
romT(1,:)=[min(Tho(:,1)),max(Tho(:,1))]; % Thoracic Pitch
romT(2,:)=[min(Tho(:,2)),max(Tho(:,2))]; % Thoracic Yaw
romT(3,:)=[min(Tho(:,3)),max(Tho(:,3))]; % Thoracic Roll

%%%% Export Outputs %%%%
datafile=strcat(sn,'_',tn,'_Data.mat'); % Name of data export file
VariableExplanation=["Lum = nx3 matrix (where n= number of frames) holding time-series of angular position (in degrees) for gross lumbar sement"
                     "Tho = nx3 matrix (where n= number of frames) holding time-series of angular position (in degrees) for gross thoracic segment"
                     "romL = 3x2 matrix holding max/min for each angle (in degrees) of the gross lumbar segment. 1st index specifies axis (1=flex, 2=lat, 3=axial), 2nd index specifies min(1)/max(2)"
                     "romT = 3x2 matrix holding max/min for each angle (in degrees) of the gross thoracic segment. 1st index specifies axis (1=flex, 2=lat, 3=axial), 2nd index specifies min(1)/max(2)"
                     "Pos = nx3x6 matrix (where n= number of frames) holding the time-series for the 3-component angular position (in radians) of each spine segment. 1st index = frame, 2nd = angle (1=forward, 2=side, 3=twist), 3rd = segment (1=L5->L4, ... , 6=T12->T11)"
                     "Vel = nx3x6 matrix holding the time-series for the 3-component angular velocity (in rad/sec) of each spine segment. Formatted the same as Pos"
                     "Acc = nx3x6 matrix holding the time-series for the 3-component angular acceleration (in rad/s^2) of each spine segment. Formatted the same as Pos"];
save(datafile,'Lum','Tho','romL','romT','Pos','Vel','Acc','sm','ns','VariableExplanation');

end
