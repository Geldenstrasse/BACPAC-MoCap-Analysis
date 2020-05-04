function [] = StrainAnalysis(tn,sn)
% ---------------------------- Begin Header ---------------------------- %
% FILE:     StrainAnalysis.m
% AUTHOR:   Darian Emmett
% DATE:   03/07/2020
% 
% PURPOSE: This function runs through the strain analysis code for each
% dataset name passed into it
%
%
% INPUTS: tn (trial name, written exactly as the name of the corresponding file), sn (subject number, should be formatted 'Sxx_', where 'S01_' corresponds to subject 1)
% 
% 
% OUTPUTS: generates an external .avi file that animates the strain
% calculations at each moment of the trial
%
%
% NOTES: 
%
%
% VERSION HISTORY
% V1 - Original
% ----------------------------- End Header ----------------------------- %

name={strcat(tn,'.mat'),tn};
input=load(name{1},name{2});

Fs=input.(name{2}).FrameRate; % Sample Frequency %
ns=input.(name{2}).Frames; % total number of data points
nm=length(input.(name{2}).Trajectories.Labeled.Labels); % total number of markers
labels=input.(name{2}).Trajectories.Labeled.Labels;
P=ns/Fs; % Overall data collection interval time Period (s) %
time=0:1/Fs:P-1/Fs;

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
if length(input.(name{2}).Trajectories.Labeled.Labels)~=36
    DataFile=strcat('MError',sn,tn,'_Data.mat');
    save(DataFile,'Fs','ns','P');
    return
end
% create matrix to hold every marker's location
nc=nm*3; % gives the number of columns in BIG (the total number of markers x3)
BIG=nan(nc,ns); % initialize the matrix to hold all of the markers, ordered.
%% Name the markers
%           1        2        3      4     5     6      7       8        9       10     11     12     13    14     15     16     17     18    19    20    21    22    23    24    25    26    27    28    29    30     31     32     33     34     35     36
markers={'PSISL3','PSISL2','PSISL','S2L','S2C','S2R','PSISR','PSISR2','PSISR3','L5L4','L5L3','L5L2','L5L1','L5C','L5R1','L5R2','L5R3','L5R4','L4L','L4C','L4R','L3L','L3C','L3R','L2L','L2C','L2R','L1L','L1C','L1R','T12L','T12C','T12R','T11L','T11C','T11R'};
% Populate the ordered matrix of marker xyz values, using the names of each marker to index througth the names of the markers 
% in the data and reorder them in BIG
for i=1:nm
    BIG(3*i-2:3*i,:)=input.(name{2}).Trajectories.Labeled.Data(strcmp(labels,(markers{i})),1:3,:);
end
BIG=BIG'; % The data comes in with the samples/frames on the horizontal axis, this adjusts to have the markers specified by column, the sample specified by row
clear input;

% % Generate parameters for 3rd order low-pass Butterworth filter
% [b,a]=butter(3,6/(Fs/2));
% % Filter BIG to clean the data a bit
% for i=1:nm*3
%     BIG(:,i)=filter(b,a,BIG(:,i));
% end

% Generate matrix for 2D coordinates
D2=zeros(ns,2,36); % (sample, xy, marker)

refp=1; % set reference frame/sample for the strain calculations

for s=1:20
    % Set coord for each point on the square (BL=bottom left, TR=top right, etc)
    if s<9
        BL=[3*s-2,3*s-1,3*s];
        BR=[3*s+1,3*s+2,3*s+3];
        TR=[3*(s+9)+1,3*(s+9)+2,3*(s+9)+3];
        TL=[3*(s+9)-2,3*(s+9)-1,3*(s+9)];
    elseif s==9
        BL=[3*12-2,3*12-1,3*12];
        BR=[3*14-2,3*14-1,3*14];
        TL=[3*19-2,3*19-1,3*19];
        TR=[3*20-2,3*20-1,3*20];
    elseif s==10
        BL=[3*14-2,3*14-1,3*14];
        BR=[3*16-2,3*16-1,3*16];
        TL=[3*20-2,3*20-1,3*20];
        TR=[3*21-2,3*21-1,3*21];
    elseif s>10
        temp=.5*mod(s,2);
        BL=[3*(s*3/2+2+temp)-2,3*(s*3/2+2+temp)-1,3*(s*3/2+2+temp)];
        BR=[3*(s*3/2+2+temp)+1,3*(s*3/2+2+temp)+2,3*(s*3/2+2+temp)+3];
        TR=[3*(s*3/2+5+temp)+1,3*(s*3/2+5+temp)+2,3*(s*3/2+5+temp)+3];
        TL=[3*(s*3/2+5+temp)-2,3*(s*3/2+5+temp)-1,3*(s*3/2+5+temp)];
    end
    
    m1=BL(3)/3; % marker 1 in following calc will be defined by the BL of the current square
    m2=BR(3)/3; % marker 2 in following calc will be defined by the BR of the current square
    m3=TL(3)/3; % marker 3 in following calc will be defined by the TL of the current square
    m4=TR(3)/3; % marker 4 in following calc will be defined by the TR of the current square

    % Generate matrix of 2D coordinates
    for is=1:ns % for this for each sample
        AB=BIG(is,BR)-BIG(is,BL); % % Define vector from BL to BR
        AC=BIG(is,TL)-BIG(is,BL); % % Define vector from BL to TL
        AD=BIG(is,TR)-BIG(is,BL); % % Define vector from BL to TR
        BD=BIG(is,TR)-BIG(is,BR); % % Define vector from BR to TR
        
        uX=AB/norm(AB);
        uY=cross(cross(AB,AC),AB);
        uY=uY/norm(uY);

        if m1==1 % only do this for the very first marker. Each other BL will be defined as the TL of the square beneath it
            D2(is,:,m1)=[0,0]; % set coord for first marker to (0,0)
        end
            D2(is,:,m2)=D2(is,:,m1)+[dot(AB,uX),dot(AB,uY)]; % set coord for BR marker in current square
        if or(m3==10,mod(m3-19,3)==0) % only do this for left-most markers of each row
            D2(is,:,m3)=D2(is,:,m1)+[dot(AC,uX),dot(AC,uY)]; % set coord for TL marker in current square
        end
        uX=AB/norm(AB);
        uY=cross(cross(AB,BD),AB);
        uY=uY/norm(uY);
        D2(is,:,m4)=D2(is,:,m1)+[dot(AD,uX),dot(AD,uY)]; % set coord for TR marker in current square
    end
end


% Square: Starts with bottom left and moves left-to-right, bottom-to-top
%                        __         __
% Triangle: 1 - |\ , 2 - \ | , 3 - | / , 4 - /|
%               |_\       \|       |/       /_|
% Line: Horiz always 1, vert always 2, diag always 3
ds0=[nan;nan;nan]; % [segment1, segment2, segment3] initial delta-s^2 of each segment of current triangle
ds=[nan;nan;nan]; % [segment1, segment2, segment3] current delta-s^2 of each segment of current triangle
dds=[nan;nan;nan]; % [segment1, segment2, segment3] difference in ds from ds0 for current triangle at current sampling
dX=nan(3,3); % (component(x,y,xy), segment)
E=nan(3,4,20,ns); % ((Eyy,Exx,Exy), triangle, square, sample) matrix for the E values of the strain map
M=[0,0,0,0,0,0,0,0]; % [BL,BR,TL,TR,BL,BR,TL,TR]
for s=1:20
    % Set coord for each point on the square (BL=bottom left, TR=top right, etc)
    if s<9
        M(1)=s;
        M(2)=s+1;
        M(3)=s+9;
        M(4)=s+10;
    elseif s==9
        M(1)=12;
        M(2)=14;
        M(3)=19;
        M(4)=20;
    elseif s==10
        M(1)=14;
        M(2)=16;
        M(3)=20;
        M(4)=21;
    elseif s>10
        temp=2+.5*mod(s,2);
        M(1)=s*3/2+temp;
        M(2)=s*3/2+temp+1;
        M(3)=s*3/2+temp+4;
        M(4)=s*3/2+temp+3;
    end
    M(5)=M(1);M(6)=M(2);M(7)=M(3);M(8)=M(4);
    for it=1:4
        % Initial distance components for segments of triangle 1
        dX(1:2,1)=(D2(refp,:,M(it))-D2(refp,:,M(it+1))); % dx, dy for segment 1
        dX(1:2,2)=(D2(refp,:,M(it))-D2(refp,:,M(it+2))); % dx, dy for segment 2
        dX(1:2,3)=(D2(refp,:,M(it+1))-D2(refp,:,M(it+2))); % dx, dy for segment 3
        dX(3,1)=dX(1,1)*dX(2,1); % dxdy for segment 1
        dX(3,2)=dX(1,2)*dX(2,2); % dxdy for segment 2
        dX(3,3)=dX(1,3)*dX(2,3); % dxdy for segment 3
        % matrix of dX^2 stuffs
        dXM=[2*dX(1,1)^2,2*dX(2,1)^2,4*dX(3,1);2*dX(1,2)^2,2*dX(2,2)^2,4*dX(3,2);2*dX(1,3)^2,2*dX(2,3)^2,4*dX(3,3)];
        % Initial squared distance magnitudes for segments
        ds0(1)=dX(1,1)^2+dX(2,1)^2;
        ds0(2)=dX(1,2)^2+dX(2,2)^2;
        ds0(3)=dX(1,3)^2+dX(2,3)^2;
        % Current squared distance magnitudes for segments
        for is=1:ns
            ds(1)=dot((D2(is,:,M(it))-D2(is,:,M(it+1))),(D2(is,:,M(it))-D2(is,:,M(it+1))));
            ds(2)=dot((D2(is,:,M(it))-D2(is,:,M(it+2))),(D2(is,:,M(it))-D2(is,:,M(it+2))));
            ds(3)=dot((D2(is,:,M(it+1))-D2(is,:,M(it+2))),(D2(is,:,M(it+1))-D2(is,:,M(it+2))));
            % Get difference from ds0 to ds:
            dds(1)=ds(1)-ds0(1);
            dds(2)=ds(2)-ds0(2);
            dds(3)=ds(3)-ds0(3);
            % Calculation of Exx, Eyy, Exy
            E(:,it,s,is)=mldivide(dXM,dds);
        end
    end
end

Ex=nan(4,20,ns); Ex(:,:,:)=E(1,:,:,:);
Ey=nan(4,20,ns); Ey(:,:,:)=E(2,:,:,:);
Exy=nan(4,20,ns); Exy(:,:,:)=E(3,:,:,:);

% clean any crazy outliers
Ex(abs(Ex)>2)=nan; % Clear big outliers for Exx
Ey(abs(Ey)>2)=nan; % Clear big outliers for Eyy
Exy(abs(Exy)>2)=nan; % Clear big outliers for Exy

Holes=max(isnan(Ex),[],'all')+max(isnan(Ey),[],'all')+max(isnan(Exy),[],'all');

%                        __         __
% Triangle: 1 - |\ , 2 - \ | , 3 - | / , 4 - /|
%               |_\       \|       |/       /_|

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Marker Layout %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%                         T11L(34)              T11C(35)              T11R(36)                        % 28
%                         T12L(31)              T12C(32)              T12R(33)                        % 24
%                         L1L(28)               L1C(29)               L1R(30)                         % 20
%                         L2L(25)               L2C(26)               L2R(27)                         % 16
%                         L3L(22)               L3C(23)               L3R(24)                         % 12
%                         L4L(19)               L4C(20)               L4R(21)                         % 8
% L5L4(10)     L5L3(11)   L5L2(12)   L5L1(13)   L5C(14)    L5R1(15)   L5R2(16)   L5R3(17)   L5R4(18)  % 4
% PSISL3(1)    PSISL2(2)  PSISL1(3)  S2L(4)     S2C(5)     S2R(6)     PSISR1(7)  PSISR2(8)  PSISR3(9) % 0
%     0     2      4    6    8     9   10    11   12    13   14    15    16    18   20    22   24
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% T1:           T2:           T3:           T4:
%      3           2          1________2     2
%     /\          /|           \      /      |\
%    /  \       3/ |            \    /       | \3
%   /    \       \ |             \  /        | /
% 1/______\2      \|1             \/3       1|/

% initialize xpatch coordinates matrix (each column gives the x-coord for the 3 points of each triangle)
xpatch=zeros(3,80);
% xpatch for triangle 1 of each square (this is different from the triangle 1 of strain calculating section; see diagram directly above)
xpatch(1,1:20)=[0, 4, 8,10,12,14,16,20,8, 12,8, 12,8, 12,8, 12,8, 12,8, 12];
xpatch(2,1:20)=[4, 8,10,12,14,16,20,24,12,16,12,16,12,16,12,16,12,16,12,16];
xpatch(3,1:20)=[2, 6, 9,11,13,15,18,22,10,14,10,14,10,14,10,14,10,14,10,14];
% xpatch for triangle 2 of each square
xpatch(1,21:40)=[4, 8,10,12,14,16,20,24,12,16,12,16,12,16,12,16,12,16,12,16];
xpatch(2,21:40)=[4, 8,10,12,14,16,20,24,12,16,12,16,12,16,12,16,12,16,12,16];
xpatch(3,21:40)=[2, 6, 9,11,13,15,18,22,10,14,10,14,10,14,10,14,10,14,10,14];
% xpatch for triangle 3 of each square
xpatch(1,41:60)=[0, 4, 8,10,12,14,16,20,8, 12,8, 12,8, 12,8, 12,8, 12,8, 12];
xpatch(2,41:60)=[4, 8,10,12,14,16,20,24,12,16,12,16,12,16,12,16,12,16,12,16];
xpatch(3,41:60)=[2, 6, 9,11,13,15,18,22,10,14,10,14,10,14,10,14,10,14,10,14];
% xpatch for triangle 4 of each square
xpatch(1,61:80)=[0, 4, 8,10,12,14,16,20, 8,12, 8,12, 8,12, 8,12, 8,12, 8,12];
xpatch(2,61:80)=[0, 4, 8,10,12,14,16,20, 8,12, 8,12, 8,12, 8,12, 8,12, 8,12];
xpatch(3,61:80)=[2, 6, 9,11,13,15,18,22,10,14,10,14,10,14,10,14,10,14,10,14];
% xpatch for calibration / color key
xpatch(1,81:82)=[28,28]; % x-val for the bottom left point in each triangle
xpatch(2,81:82)=[30,30]; % x-val for the bottom right point in each triangle
xpatch(3,81:82)=[29,29]; % x-val for the top middle point in each triangle

% initialize ypatch coordinates matrix (each column gives the y-coord for the 3 points of each triangle)
ypatch=zeros(3,80);
% ypatch for triangle 1 of each square (this is different from the triangle 1 of strain calculating section; see diagram directly above)
ypatch(1,1:20)=[0, 0, 0, 0, 0, 0, 0, 0, 4, 4, 8, 8,12,12,16,16,20,20,24,24];
ypatch(2,1:20)=[0, 0, 0, 0, 0, 0, 0, 0, 4, 4, 8, 8,12,12,16,16,20,20,24,24];
ypatch(3,1:20)=[2, 2, 2, 2, 2, 2, 2, 2, 6, 6,10,10,14,14,18,18,22,22,26,26];
% ypatch for triangle 2 of each square
ypatch(1,21:40)=[0, 0, 0, 0, 0, 0, 0, 0, 4, 4, 8, 8,12,12,16,16,20,20,24,24];
ypatch(2,21:40)=[4, 4, 4, 4, 4, 4, 4, 4, 8, 8,12,12,16,16,20,20,24,24,28,28];
ypatch(3,21:40)=[2, 2, 2, 2, 2, 2, 2, 2, 6, 6,10,10,14,14,18,18,22,22,26,26];
% ypatch for triangle 3 of each square
ypatch(1,41:60)=[4, 4, 4, 4, 4, 4, 4, 4, 8, 8,12,12,16,16,20,20,24,24,28,28];
ypatch(2,41:60)=[4, 4, 4, 4, 4, 4, 4, 4, 8, 8,12,12,16,16,20,20,24,24,28,28];
ypatch(3,41:60)=[2, 2, 2, 2, 2, 2, 2, 2, 6, 6,10,10,14,14,18,18,22,22,26,26];
% ypatch for triangle 4 of each square
ypatch(1,61:80)=[0, 0, 0, 0, 0, 0, 0, 0, 4, 4, 8, 8,12,12,16,16,20,20,24,24];
ypatch(2,61:80)=[4, 4, 4, 4, 4, 4, 4, 4, 8, 8,12,12,16,16,20,20,24,24,28,28];
ypatch(3,61:80)=[2, 2, 2, 2, 2, 2, 2, 2, 6, 6,10,10,14,14,18,18,22,22,26,26];
% ypatch for calibration / color key
ypatch(1,81:82)=[1,27]; % y-val for the bottom left point in each triangle
ypatch(2,81:82)=[1,27]; % y-val for the bottom right point in each triangle
ypatch(3,81:82)=[0,28]; % y-val for the top middle point in each triangle

% E(1,2,3,4) = ((Exx,Eyy,Exy), triangle, square, sample)
Cx=zeros(1,82); Cy=zeros(1,82); Cxy=zeros(1,82); % initialize the color coding variables for the patch plot
Cx(81)=-max(abs(Ex),[],'all'); Cx(82)=-Cx(81);
Cy(81)=-max(abs(Ey),[],'all'); Cy(82)=-Cy(81);
Cxy(81)=-max(abs(Exy),[],'all'); Cxy(82)=-Cxy(81);

X=nan(nc/3,3);

CMap=zeros(200,3);
CMap(1:50,1)=zeros(50,1); CMap(51:100,1)=linspace(0,1,50); CMap(101:200,1)=ones(100,1);
CMap(1:50,2)=linspace(0,1,50); CMap(51:100,2)=ones(50,1); CMap(101:200,2)=linspace(1,0,100);
CMap(1:50,3)=linspace(1,0,50); CMap(51:100,3)=zeros(50,1); CMap(101:200,3)=zeros(100,1);
colormap(CMap);

%% Initialize video

myVideo = VideoWriter([sn tn '_Anim']); % open/create video file
myVideo.FrameRate = Fs;  % set frame rate
open(myVideo) % open the video file for writing

set(0, 'DefaulttextInterpreter', 'none');

for is=1:ns
    for i=1:nm
        X(i,1)=BIG(is,3*i-2);
        X(i,2)=BIG(is,3*i-1);
        X(i,3)=BIG(is,3*i);
    end
    for s=1:20 % for each square, make the smaller triangles the average of the larger overlapping ones
        % first do the color coding for Exx
        Cx(s)=(Ex(1,s,is)+Ex(4,s,is))/2;
        Cx(s+20)=(Ex(2,s,is)+Ex(4,s,is))/2;
        Cx(s+40)=(Ex(2,s,is)+Ex(3,s,is))/2;
        Cx(s+60)=(Ex(1,s,is)+Ex(3,s,is))/2;
        % next color code for Eyy
        Cy(s)=(Ey(1,s,is)+Ey(4,s,is))/2;
        Cy(s+20)=(Ey(2,s,is)+Ey(4,s,is))/2;
        Cy(s+40)=(Ey(2,s,is)+Ey(3,s,is))/2;
        Cy(s+60)=(Ey(1,s,is)+Ey(3,s,is))/2;
        % next color code for Exy
        Cxy(s)=(Exy(1,s,is)+Exy(4,s,is))/2;
        Cxy(s+20)=(Exy(2,s,is)+Exy(4,s,is))/2;
        Cxy(s+40)=(Exy(2,s,is)+Exy(3,s,is))/2;
        Cxy(s+60)=(Exy(1,s,is)+Exy(3,s,is))/2;
    end
    if is==1
        tiledlayout(2,2);
        nexttile
        EyP=patch(xpatch,ypatch,Cy);
        colorbar
        ptitle=strcat("Eyy ",tn);
        title(ptitle,'FontSize',8);
        xlabel('X','FontSize',8); ylabel('Y','FontSize',8);
        nexttile
        sp=scatter3(X(:,1),X(:,2),X(:,3),10,'filled');
        title(tn,'FontSize',8);
        axis equal off
        view(45,30)
        nexttile
        ExyP=patch(xpatch,ypatch,Cxy);
        colorbar
        ptitle=strcat("Exy ",tn);
        title(ptitle,'FontSize',8);
        xlabel('X','FontSize',8); ylabel('Y','FontSize',8);
        nexttile
        ExP=patch(xpatch,ypatch,Cx);
        colorbar
        ptitle=strcat("Exx ",tn);
        title(ptitle,'FontSize',8);
        xlabel('X','FontSize',8); ylabel('Y','FontSize',8);
        pause(1/Fs)
    else
        sp.XData=X(:,1);
        sp.YData=X(:,2);
        sp.ZData=X(:,3);
        ExyP.CData=Cxy;
        ExP.CData=Cx;
        EyP.CData=Cy;
        pause(1/Fs)
    end
    frame = getframe(gcf); %get frame
    writeVideo(myVideo, frame);
end
close(myVideo)

close % close any figures that might still be open

% Find maximum strain moment and identify the corresponding indices
% (triangle square, sample)
[~,B]=max(Ey,[],'all','linear');
[it,s,is] = ind2sub(size(Ey), B);numsq=['Square' ' ' num2str(s)];

TimeMaxStrain=(is-1)/Fs;
mt=string(TimeMaxStrain);

ExMAX=zeros(1,82); EyMAX=zeros(1,82); ExyMAX=zeros(1,82); % initialize the color coding variables for the patch plot
ExMAX(81)=-max(abs(Ex),[],'all'); ExMAX(82)=-ExMAX(81);
EyMAX(81)=-max(abs(Ey),[],'all'); EyMAX(82)=-EyMAX(81);
ExyMAX(81)=-max(abs(Exy),[],'all'); ExyMAX(82)=-ExyMAX(81);

% Now recalculate C values for the max-strain time step
for s=1:20 % for each square, make the smaller triangles the average of the larger overlapping ones
    % first do the color coding for Exx
    ExMAX(s)=(Ex(1,s,is)+Ex(4,s,is))/2;
    ExMAX(s+20)=(Ex(2,s,is)+Ex(4,s,is))/2;
    ExMAX(s+40)=(Ex(2,s,is)+Ex(3,s,is))/2;
    ExMAX(s+60)=(Ex(1,s,is)+Ex(3,s,is))/2;
    % next color code for Eyy
    EyMAX(s)=(Ey(1,s,is)+Ey(4,s,is))/2;
    EyMAX(s+20)=(Ey(2,s,is)+Ey(4,s,is))/2;
    EyMAX(s+40)=(Ey(2,s,is)+Ey(3,s,is))/2;
    EyMAX(s+60)=(Ey(1,s,is)+Ey(3,s,is))/2;
    % next color code for Exy
    ExyMAX(s)=(Exy(1,s,is)+Exy(4,s,is))/2;
    ExyMAX(s+20)=(Exy(2,s,is)+Exy(4,s,is))/2;
    ExyMAX(s+40)=(Exy(2,s,is)+Exy(3,s,is))/2;
    ExyMAX(s+60)=(Exy(1,s,is)+Exy(3,s,is))/2;
end

switch it
    case 1
        numtr="Bottom Triangle";
    case 2
        numtr="Right Triangle";
    case 3
        numtr="Top Triangle";
    case 4
        numtr="Left Triangle";
end
% Generate vectors of the strain in the max strain triangle over the full
% sample period
StrX=nan(ns,1); StrX(:)=Ex(it,s,:);
StrY=nan(ns,1); StrY(:)=Ey(it,s,:);
StrXY=nan(ns,1); StrXY(:)=Exy(it,s,:);

colormap(CMap);

tiledlayout(2,2);
nexttile
patch(xpatch,ypatch,EyMAX);
colorbar
ptitle=strcat("Eyy for ",tn," at ",mt,"s");
title(ptitle,'FontSize',8);
xlabel('X','FontSize',8); ylabel('Y','FontSize',8);
nexttile
plot(time,StrX);
ptitle=strcat("Strain v Time for ",numtr," of ",numsq);
title(ptitle,'FontSize',8);
xlabel('Time (s)','FontSize',8)
ylabel('Strain','FontSize',8)
hold on
plot(time,StrY);
plot(time,StrXY);
legend('Normal X', 'Normal Y', 'Shear','FontSize',5);
nexttile
patch(xpatch,ypatch,ExyMAX);
colorbar
ptitle=strcat("Exy for ",tn," at ",mt,"s");
title(ptitle,'FontSize',8);
xlabel('X','FontSize',8); ylabel('Y','FontSize',8);
nexttile
patch(xpatch,ypatch,ExMAX);
colorbar
ptitle=strcat("Exx for ",tn," at ",mt,"s");
title(ptitle,'FontSize',8);
xlabel('X','FontSize',8); ylabel('Y','FontSize',8);

set(gcf,'PaperUnits','inches');
set(gcf,'PaperPosition',[0 0 8.5 6.5]);
set(gca,'FontSize',8);
set(gca,'FontName','Calibri');

print('-r300',strcat(sn,tn,'_EvT'),'-dpng');

DataFile=strcat(sn,tn,'_Data.mat');
save(DataFile,'Fs','ns','P','E','Ex','Ey','Exy','ExMAX','EyMAX','ExyMAX','TimeMaxStrain','StrX','StrY','StrXY','Holes');
end
