function [Marker,MarkerFreq,EMG,AnalogFreq,Measures]=load_c3d(filename)

% load_VICON_c3d function read data from .c3d file and extract trajectory (markers positions),
% forces and other analog data. Center of pression is also calculated
%
% Syntax:
% [Nframes,Nmarkers,Marker,MarkerName,MarkerFreq,Nchannels,Analog,AnalogName,AnalogFreq]=load_VICON_c3d(filename)
%
% Input:        filename          - Name of .c3d file
%               bodyModel         - Model used on markers to construct body
%                                       =1 if lower limb model
%                                       =2 if full body model
%
% Output:       Nframes           - Number of frames
%               FirstFrame        - Number of the first frame of c3d data
%               Nmarkers          - Number of markers
%               MarkerFreq        - Sampling frequency of trajectory measurements
%               Marker            - Structure where a field is the name of one marker and contains the x,y,z positions at each time
%                                   -> mkr.'MakerName'.x or .y or .z
%               Nchannels         - Number of analog data
%               Analog            - Structure where a field is the name of the analog data and contains the analog value at each time
%                                   -> analog.'AnalogName'
%               AnalogFreq        - Sampling frequency of analog data
%               Measures          - Reported metrics and weight of the subject
%

if strcmp(filename(1:3),'Hea')
    nEMGs = 16;
elseif strcmp(filename(1:3),'SCI')
    nEMGs = 28;
end

DEBUG = 0;

factor=4;

Side = struct('L',struct,'R',struct);
Marker=struct('Pos',Side,'Seg',Side,'Ang',Side);
Analog=struct('EMG',Side);
Measures=struct;

fid=fopen(filename);
if fid<0, msgbox('ERROR when opening the .c3d file!','Load_c3d','error'), return, end;


%% C3D HEADER
parameter_section=fread(fid,1,'int8'); % Points to the first block of the parameter section
key=fread(fid,1,'int8');               % Key indicating file format
if key~=80, msgbox({'ERROR when opening the .c3d file!';'The file format is not .c3d'},'Load_VICON_c3d','error'), return, end;

Nmarkers=fread(fid,1,'int16');   % Number of trajectories stored within the file
Nanalog=fread(fid,1,'int16');    % Number of analog measurements recorded in the file (Nchannels*SamplesPerChannels)
FirstFrame=fread(fid,1,'int16'); % Location of the start of the parameter records within the file (1 based, not 0!)
LastFrame=fread(fid,1,'int16');  % Location of the end of the parameter records within the file
if LastFrame>0, Nframes=LastFrame-FirstFrame+1;
else msgbox({'ERROR when opening the .c3d file!';'The file does not contain markers data'},'Load_VICON_c3d','error'), return, end
gap_size=fread(fid,1,'int16');    % Maximum interpolation gap in 3D frames
SCALE=fread(fid,1,'real*4');      % If the header SCALE value is negative then the 3D and analog data is always stored in floating-point format
if SCALE>0, msgbox({'ERROR when opening the .c3d file!';'The data must be rescaled'},'Load_VICON_c3d','error'), return, end;
data_start=fread(fid,1,'int16');     % Location of the first block of the trajectories and analog data section
analog_samples=fread(fid,1,'int16'); % Number of analog samples per 3D frame
MarkerFreq=fread(fid,1,'real*4');  % Trajectories frame rate in Hz
AnalogFreq=MarkerFreq*analog_samples; % Analog measurements frame rate in Hz
if Nanalog==0, Nchannels=1; else Nchannels=Nanalog/analog_samples; end

%% C3D DATA

fseek(fid,512*(data_start-1),-1);
DATAmatrix=fread(fid,[Nmarkers*factor+Nanalog,Nframes],'real*4'); % Matrix containing all data
MARKERS=DATAmatrix(1:Nmarkers*factor,:);     % Matrix containing markers data
ANALOG=reshape(DATAmatrix(Nmarkers*factor+1:end,:),Nchannels,[]); % Matrix containing analog data


%% C3D PARAMETERS

fseek(fid,512*(parameter_section-1),'bof');
nothing=fread(fid,1,'int16'); % First 2 bytes
nothing=fread(fid,1,'int8');  % Number of parameter blocks to follow
nothing=fread(fid,1,'int8');  % Processor type
Nchar=fread(fid,1,'int8');    % Number of characters in the Group name
IDgroup=fread(fid,1,'int8');  % Group ID number (always negative)

while Nchar > 0
    if IDgroup < 0
        IDgroup=abs(IDgroup);
        groupName=fread(fid,[1,Nchar],'char'); % Group name (ASCII characters – upper case A-Z, 0-9 and underscore _ only)
        offset=fread(fid,1,'int16');    % A signed integer offset in bytes pointing to the start of the next group/parameter.
        Nchar_desc=fread(fid,1,'int8'); % Number of characters in the Group description.
        desc=fread(fid,[1,Nchar_desc],'char'); % Group description (ASCII characters – mixed case).
        disp(sprintf('%s, %s',groupName,desc));
        
        fseek(fid,offset-3-Nchar_desc,'cof');
        group.(strcat('GROUP',(num2str(IDgroup))))=groupName;
        
    else
        groupName=group.(strcat('GROUP',(num2str(IDgroup))));
        paramName=fread(fid,[1,Nchar],'char'); % Parameter name (ASCII characters – normally upper case numeric or underscore only)
        disp(sprintf('- %s',paramName));
        
        offset=fread(fid,1,'int16'); % A signed integer offset in bytes pointing to the start of the next group/parameter.
        nextPos=ftell(fid)+offset(1)-2; % Position of next param bloc
        
        
        %% TRIAL BLOC
        if strcmp(sprintf('%c',groupName),'TRIAL') & strcmp(sprintf('%c',paramName),'ACTUAL_START_FIELD'),
            datatype = whichType(fid); % Data type (real, integer, ...)
            Ndim=fread(fid,1,'int8'); % Number of dimensions (0-7) of the parameter – 0 if the parameter is scalar
            dimension=fread(fid,[1,Ndim],'int8');  % Parameter dimensions
            if FirstFrame<0,
                FirstFrame=fread(fid,1,'int32');
            else
                paramData=fread(fid,dimension,datatype); % Parameter data
                FirstFrame=paramData(1);
            end
            
        elseif strcmp(sprintf('%c',groupName),'TRIAL') & strcmp(sprintf('%c',paramName),'ACTUAL_END_FIELD'),
            datatype = whichType(fid); % Data type (real, integer, ...)
            Ndim=fread(fid,1,'int8'); % Number of dimensions (0-7) of the parameter – 0 if the parameter is scalar
            dimension=fread(fid,[1,Ndim],'int8');  % Parameter dimensions
            if LastFrame<0,
                LastFrame=fread(fid,1,'int32');
            else
                paramData=fread(fid,dimension,datatype); % Parameter data
                LastFrame=paramData(1);
            end
            
        elseif strcmp(sprintf('%c',groupName),'POINT') & strcmp(sprintf('%c',paramName),'LABELS'),
            datatype = whichType(fid); % Data type (real, integer, ...)
            Ndim=fread(fid,1,'int8');  % Number of dimensions (0-7) of the parameter – 0 if the parameter is scalar
            dimension=fread(fid,[1,Ndim],'int8');    % Parameter dimensions
            
            if dimension(2) < 0, dimension(2) = 2^7 - dimension(2); end
            paramData=fread(fid,dimension,datatype); % Parameter data
            
            ih=1;
            for jj=1:dimension(2),
                labels=sprintf('%c',paramData(:,jj));
                disp(sprintf('-- %s (%d)',labels,jj))
                
                if ~isempty(strfind(labels,'Power')), break, end
                if isempty(strfind(labels,'*'))
                    [Marker,name]=create_marker(labels,MARKERS,jj,Marker);
                    ih=ih+1;
                end
            end
            
        elseif strcmp(sprintf('%c',groupName),'ANALOG') & strcmp(sprintf('%c',paramName),'GEN_SCALE'),
            datatype = whichType(fid); % Data type (real, integer, ...)
            Ndim=fread(fid,1,'int8');  % Number of dimensions (0-7) of the parameter – 0 if the parameter is scalar
            dimension=fread(fid,[1,Ndim],'int8');    % Parameter dimensions
            if Ndim == 0,
                paramData=fread(fid,1,datatype); % Parameter data
                param.GEN_SCALE=1;
            else
                paramData=fread(fid,dimension,datatype); % Parameter data
                param.GEN_SCALE=paramData;
            end
            
        elseif strcmp(sprintf('%c',groupName),'ANALOG') & strcmp(sprintf('%c',paramName),'SCALE'),
            datatype = whichType(fid); % Data type (real, integer, ...)
            Ndim=fread(fid,1,'int8');  % Number of dimensions (0-7) of the parameter – 0 if the parameter is scalar
            dimension=fread(fid,[1,Ndim],'int8');    % Parameter dimensions
            paramData=fread(fid,dimension,datatype); % Parameter data
            param.SCALE=paramData;
            
        elseif strcmp(sprintf('%c',groupName),'ANALOG') & strcmp(sprintf('%c',paramName),'OFFSET'),
            datatype = whichType(fid); % Data type (real, integer, ...)
            Ndim=fread(fid,1,'int8');  % Number of dimensions (0-7) of the parameter – 0 if the parameter is scalar
            dimension=fread(fid,[1,Ndim],'int8');    % Parameter dimensions
            paramData=fread(fid,dimension,datatype); % Parameter data
            param.OFFSET=paramData;
            
            
        elseif strcmp(sprintf('%c',groupName),'ANALOG') && strcmp(sprintf('%c',paramName),'LABELS'),
            datatype = whichType(fid); % Data type (real, integer, ...)
            Ndim=fread(fid,1,'int8');  % Number of dimensions (0-7) of the parameter – 0 if the parameter is scalar
            dimension=fread(fid,[1,Ndim],'int8');    % Parameter dimensions
            
            paramData=fread(fid,dimension,datatype); % Parameter data
            param.LABELS=paramData;
            analogTemp=ANALOG;
            
            ip=1;
            chan=1; jj = 1;
            while chan <= nEMGs
                labels=sprintf('%c',paramData(:,jj));
                disp(sprintf('-- %s (%d)', labels,chan));
                if isempty(strfind(labels,'Pin')) & isempty(strfind(labels,'MAIN'))
                    analogTemp(jj,:)=(ANALOG(jj,:)-abs(param.OFFSET(jj)))*param.SCALE(jj)*param.GEN_SCALE;
                    [Analog,name]=create_analog(labels,analogTemp,jj,Analog);
                    chan = chan+1;
                end
                jj = jj+1;
            end
            
            
        elseif strcmp(sprintf('%c',groupName),'PROCESSING')
            datatype = whichType(fid); % Data type (real, integer, ...)
            Ndim=fread(fid,1,'int8');  % Number of dimensions (0-7) of the parameter – 0 if the parameter is scalar
            dimension=fread(fid,[1,Ndim],'int8');    % Parameter dimensions
            paramData=fread(fid,dimension,datatype); % Parameter data
            
            switch strcat(paramName)
                case 'Bodymass', Measures.Mass = paramData;
                case 'Height', Measures.Height = paramData;
                case 'InterAsisDistance', Measures.InterAsis = paramData;
                case 'LLegLength', Measures.LLeg = paramData;
                case 'LAsisTrocanterDistance', Measures.LAsisTrocanter = paramData;
                case 'LKneeWidth', Measures.LKnee = paramData;
                case 'LAnkleWidth', Measures.LAnkle = paramData;
                case 'RLegLength', Measures.RLeg = paramData;
                case 'RAsisTrocanterDistance', Measures.RAsisTrocanter = paramData;
                case 'RKneeWidth', Measures.RKnee = paramData;
                case 'RAnkleWidth', Measures.RAnkle = paramData;
            end
        end
        
        fseek(fid,nextPos,'bof');
    end
    
    Nchar=fread(fid,1,'int8');    % Number of characters in the Group name
    IDgroup=fread(fid,1,'int8');  % Group ID number (always negative)
end


fclose(fid);

if ip==1, AnalogName=' '; end
Nchannels = length(AnalogName);

EMG = Analog.EMG;
if ~DEBUG, clc, end


% REVERSE is TRUE if subject is walking in the reverse direction
if Marker.Pos.L.ASI(2,1) > Marker.Pos.L.ASI(2,end)
     REVERSE = 1;
else REVERSE = 0; end

names = fieldnames(Marker.Pos.L);
for i = 1:length(names)
    Marker.Pos.L.(char(names(i))) = Marker.Pos.L.(char(names(i)))(1:3,:)./10; % convert mm to cm
    if REVERSE, Marker.Pos.L.(char(names(i)))(1:2,:) = -Marker.Pos.L.(char(names(i)))(1:2,:); end
end
names = fieldnames(Marker.Pos.R);
for i = 1:length(names)
    Marker.Pos.R.(char(names(i))) = Marker.Pos.R.(char(names(i)))(1:3,:)./10; % convert mm to cm
    if REVERSE, Marker.Pos.R.(char(names(i)))(1:2,:) = -Marker.Pos.R.(char(names(i)))(1:2,:); end
end
% SEGMENTS data
names = fieldnames(Marker.Seg.L);
for i = 1:length(names)
    Marker.Seg.L.(char(names(i))) = Marker.Seg.L.(char(names(i)))(1:3,:)./10; % convert mm to cm
    if REVERSE, Marker.Seg.L.(char(names(i)))(1:2,:) = -Marker.Seg.L.(char(names(i)))(1:2,:); end
end
names = fieldnames(Marker.Seg.R);
for i = 1:length(names)
    Marker.Seg.R.(char(names(i))) = Marker.Seg.R.(char(names(i)))(1:3,:)./10; % convert mm to cm
    if REVERSE, Marker.Seg.R.(char(names(i)))(1:2,:) = -Marker.Seg.R.(char(names(i)))(1:2,:); end
end
% ANGLES data
names = fieldnames(Marker.Ang.L);
for i = 1:length(names)
    Marker.Ang.L.(char(names(i))) = Marker.Ang.L.(char(names(i)))(1,:);
end
names = fieldnames(Marker.Ang.R);
for i = 1:length(names)
    Marker.Ang.R.(char(names(i))) = Marker.Ang.R.(char(names(i)))(1,:);
end
    
[Marker.JCs,Measures] = comp_JointCenters(Marker.Pos,Marker.Seg,Measures);

end

%% =================================================================================================
% Private function to determine type of data
function type = whichType(fileID)
switch fread(fileID,1,'int8')   % Length in bytes of each data element
    case -1, type='char'; % character data
    case  1, type='int8'; % 1 byte data
    case  2, type='int16';% integer data
    case  4, type='float';% real data (floating-point)
end
end

% Private function to construct a marker structure
% i.e. mkr.MARKERNAME.COORDINATE (e.g. mkr.RKNE.x)
function [mkr,name]=create_marker(labels,coord,jj,mkr)
factor=4;

if ~isempty(strfind(labels,':')), t1=strfind(':',labels)+1; else t1=1; end
if ~isempty(strfind(labels,'-')), labels(strfind('-',labels))='X'; end
if ~isempty(strfind(labels,' ')), name=(labels(t1:min(strfind(labels,' '))-1));
else name=labels(t1:end);
end

if strcmp(name(1),'L')
    if length(name) == 4
        if any([strcmp(name(end),'A'),strcmp(name(end),'L'),strcmp(name(end),'P'),strcmp(name(end),'O')]) ...
                & isempty(strfind(name,'SHO')) & isempty(strfind(name,'WRA'))
            mkr.Seg.L=setfield(mkr.Seg.L,name(2:end),coord((jj-1)*factor+1:(jj-1)*factor+4,:));
        else
            mkr.Pos.L=setfield(mkr.Pos.L,name(2:end),coord((jj-1)*factor+1:(jj-1)*factor+4,:));
        end
    elseif strfind(name,'Angle')
        mkr.Ang.L=setfield(mkr.Ang.L,name(2:end-6),coord((jj-1)*factor+1:(jj-1)*factor+4,:));
    end
elseif  strcmp(name(1),'R')
    if length(name) == 4
        if any([strcmp(name(end),'A'),strcmp(name(end),'L'),strcmp(name(end),'P'),strcmp(name(end),'O')]) ...
                & isempty(strfind(name,'SHO')) & isempty(strfind(name,'WRA'))
            mkr.Seg.R=setfield(mkr.Seg.R,name(2:end),coord((jj-1)*factor+1:(jj-1)*factor+4,:));
        else
            mkr.Pos.R=setfield(mkr.Pos.R,name(2:end),coord((jj-1)*factor+1:(jj-1)*factor+4,:));
        end
    elseif strfind(name,'Angle')
        mkr.Ang.R=setfield(mkr.Ang.R,name(2:end-6),coord((jj-1)*factor+1:(jj-1)*factor+4,:));
    end
end
end

% Private function to construct an analog structure
% i.e. analog.ANALOGNAME
function [analog,name]=create_analog(labels,analog_temp,jj,analog)
labels(strfind(labels,' '))='';
name=labels;

if strcmp(name(1),'R')
    analog.EMG.R.(name(2:end))=analog_temp(jj,:);
elseif strcmp(name(1),'L')
    analog.EMG.L.(name(2:end))=analog_temp(jj,:);
end

end


function [JCs,lg] = comp_JointCenters(Data,Seg,lg)
%% Estimate joint centers and calculate segment reference frames for lower body model
%
% Reference frames
% - Global  reference frame: GRF(Origin,X,Y,Z)
% - Body    reference frame: BRF(Origin,X,Y,Z)
% - Segment reference frame: SRF(JointCenter,postero-anterior,medio-lateral,distal-parietal)
%       Segments origins are given according to BRF
%
% Measurements
% - lg.RLeg, lg.LLeg -> leg length (cm)
% - lg.RKnee, lg.LKnee -> knee width (cm)
% - lg.RAnkle, lg.LAnkle -> ankle width (cm)
% - lg.MarkerRadius (cm)
%
% References 
% - Eric Desailly, Analyse biomecanique 3D de la marche de l'enfant deficient moteur -
%   Modelisation segmentaire et modelisation musculo-squelettique, 2008 Thesis, University of Poitier, p29-58
% - Gabriele Paolini, Interpreting PiG results: PiG biomechanical modelling, Presentation support material
% - Plug-In Gait manual
%
% Notes
%  - Ankle RF: rotate around dp axis to account for tibial torsion (subject measurement give the torsion angle)
%  - Foot RF: foot flat option where the dp axis is corrected to be horizontal when foot in on the
%          ground (the correction angle to apply is determined from recording during static calibration for example)
%
% Written by C.LeGoff, 03.March.2014
% modified by S.Anil, 25.Nov.2015
%

names = fieldnames(lg);
for i = 2:length(names)
    lg.(char(names(i))) =lg.(char(names(i)))./10; % convert mm to cm
end
if lg.LAsisTrocanter == 0, lg.LAsisTrocanter = 0.1288*lg.LLeg - 4.856; end
if lg.RAsisTrocanter == 0, lg.RAsisTrocanter = 0.1288*lg.RLeg - 4.856; end
lg.MarkerRadius = 0.8; % marker diameter = 16mm


%% Markers
PosL = Data.L; PosR = Data.R;
SegL = Seg.L; SegR = Seg.R;
LASI = PosL.ASI;
RASI = PosR.ASI;

%% Pelvis
% - origin: midpoint on ASIs line
Pelvis = (LASI+RASI)./2;

%% Hip joint centers
% - origin: based on Newington - Gage model
% create a direct orthonormal trihedron: origin at Pelvis
% computed by S.Anil

vGravity = repmat([0;0;-1],1,size(Pelvis,2));
vTan = [RASI(1,:)-LASI(1,:);RASI(2,:)-LASI(2,:);zeros(size(Pelvis(1,:)))];
vForward = cross(vTan,vGravity);

for i = 1:size(Pelvis,2)
    vTan(:,i) = vTan(:,i)/norm(vTan(:,i));
    vForward(:,i) = vForward(:,i)/norm(vForward(:,i));
end

for i = 1:size(Pelvis,2)
    BRF.(char(['frame' num2str(i)])) = [vTan(:,i),vForward(:,i),-vGravity(:,i)];
end %create the 3 orthonormal vectors attached to subject

lgAsis = mean(normvec(LASI-RASI));% inter ASIs distance
C = mean([lg.RLeg,lg.LLeg])*0.115 - 1.53;
theta = 0.5; beta = 0.314; % radians

OFFSET_Y = C*cos(theta)*sin(beta)-(lg.LAsisTrocanter+lg.MarkerRadius)*cos(beta);
OFFSET_Z = C*cos(theta)*cos(beta)-(lg.LAsisTrocanter+lg.MarkerRadius)*sin(beta);
LHJC_offset = [(C*sin(theta)-lgAsis/2); -OFFSET_Y; -OFFSET_Z];
RHJC_offset = [-(C*sin(theta)-lgAsis/2); -OFFSET_Y; -OFFSET_Z];

for i = 1:size(Pelvis,2)
    LHJC(:,i) = Pelvis(:,i) + BRF.(char(['frame' num2str(i)]))*LHJC_offset;
    RHJC(:,i) = Pelvis(:,i) + BRF.(char(['frame' num2str(i)]))*RHJC_offset;
end % in referential of the lab

%% Knee joint centers
% - origin: by CHORD function done by PlugIn Gait on VICON Nexus
if isfield(SegL,'FEO'), LKJC = SegL.FEO; 
else 
    LKJC = PosL.KNE;
    disp('[comp_JointCenters] - LKJC not computed --> assigned to LKNE marker')
end
if  isfield(SegR,'FEO'), RKJC = SegR.FEO; 
else 
    RKJC = PosR.KNE;
    disp('[comp_JointCenters] - RKJC not computed --> assigned to RKNE marker')
end

%% Ankle joint centers
% - origin: by CHORD function done by PlugIn Gait on VICON Nexus
if isfield(SegL,'TIO'), LAJC = SegL.TIO;
else
    LAJC = PosL.ANK;
    disp('[comp_JointCenters] - LAJC not computed --> assigned to LANK marker')
end
if isfield(SegR,'TIO'), RAJC = SegR.TIO;
else
    RAJC = PosR.ANK;
    disp('[comp_JointCenters] - RAJC not computed --> assigned to RANK marker')
end

%% Foot joint centers
% - origin: more or less on TOE marker, see PlugIn Gait on VICON Nexus
if isfield(SegL,'FOO'), LFJC = SegL.FOO;
else
    LFJC = PosL.TOE;
    disp('[comp_JointCenters] - LFJC not computed --> assigned to LTOE marker')
end
if isfield(SegR,'FOO'), RFJC = SegR.FOO;
else
    RFJC = PosR.TOE;
    disp('[comp_JointCenters] - RFJC not computed --> assigned to RTOE marker')
end

%% Create Output structure
%JCs.C.Pelvis = Pelvis;
JCs.L.HJC = LHJC; JCs.R.HJC = RHJC;
JCs.L.KJC = LKJC; JCs.R.KJC = RKJC;
JCs.L.AJC = LAJC; JCs.R.AJC = RAJC;
JCs.L.FJC = LFJC; JCs.R.FJC = RFJC;

end

function NORM = normvec(vector)
NORM = sqrt(sum(power(vector,2),1));
end


