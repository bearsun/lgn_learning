function lgn_search_behav_in_scanner
%LGN_SEARCH_BEHAV training code for lgn_learning project
%   4/14/16 by Liwei
% doing the behav session in scanner with button box

%% initialize everything
clc;
global ptb_RootPath %#ok<NUSED>

rng('shuffle');
sid = 0;
rect = [0 0 1024 768];

ntrial = 32;

monitorh=34.3; %30; %12;% in cm
distance=110.5; %55; %25;% in cm

subj=input('subject?','s');
% session = input('session?');

%% geometry
fixsi = 8;
righthalf = 180;
lefthalf = -180;
startAngle = 0;
fullAngle=360;

%% color
gray = [127 127 127];
white = [255 255 255];
black = [0 0 0];
% red = [175,87,87];
% green = [67,135,67];
truegreen = [0 255 0];
truered = [255 0 0];

% read iso colors
path_color = [pwd,'/subinfo/',subj,'_fmri.mat'];
rg = load(path_color);
red = rg.rg(1,:);
green = rg.rg(2,:);

fixcolor = white;
textcolor = black;
bgcolor = gray;

%% buttonbox setting
buttons = 49:52; %top, left, bottom, right
% trigger=53;
kn0 = buttons(2); %left
kn1 = buttons(1); %top
kn2 = buttons(4); %right
kn3 = buttons(3); %bottom

possiblekn = [kn0,kn1,kn2,kn3];

%% target & answer sequence
nrings=3;
stimPerRing=8;
nballs=nrings*stimPerRing;

targetoptions=[1:nballs,zeros(1,stimPerRing)]; %catch 25% 1,2,3 each 25%
targetindex = targetoptions(randperm(numel(targetoptions)));
% temp_ring = ceil(targetindex/8);
% corkeyarr=(temp_ring == 1) * kn1;%1st ring
% corkeyarr(temp_ring == 2) = kn2;
% corkeyarr(temp_ring == 3) = kn3;
% corkeyarr(temp_ring == 0) = kn0;

%% open scanner communication
IOPort('Closeall');
P4 = IOPort('OpenSerialPort', '/dev/ttyUSB0','BaudRate=115200'); %open port for receiving scanner pulse

%% open files
path_data = [pwd,'/data/pretrain-',subj,'-',num2str(session)];

outfile = fopen(path_data,'w');
fprintf(outfile,'%s %s %s %s %s %s %s\n','subject' ,'trial','targetindex' , 'keypressed' , 'cor','rt');

disp('data_file_opened');

%% initialize window
[mainwin,rect] = Screen('OpenWindow', sid, bgcolor,rect);

% open buffer
buffers = NaN(nframes,1);
for i = 1:nframes
    [buffers(i),~] = Screen('OpenOffscreenWindow', mainwin, bgcolor);
end


%% build up position arrays
% basic screen
center = [(rect(3)-rect(1))/2, (rect(4)-rect(2))/2];
fixRect = CenterRect([0 0 fixsi fixsi], rect);

% construct stimuli
corticalStimSize= 5.5;%in mm
proximalStimDist=.5;% in degrees!!
stimSeparation= .7; %in degrees, to be scaled
jitterDistance=.08;%in degrees; to be scaled
jitfreq = 5;

%empty variable
stimLocation = NaN(nrings,stimPerRing,4);
stimSize=NaN(nrings,4);
eccentricity=NaN(nrings,1);

separationAngle=360/stimPerRing;
compass = separationAngle:separationAngle:360; %rectangle
jpix=ang2pix(jitterDistance);

%enlarge for placeholder
jph=[-jpix,-jpix,jpix,jpix]';

%ecc in degree to size
for r=1:nrings
    ecc=proximalStimDist+stimSeparation*r^2;
    gratingSize=CorticalScaleFactor(corticalStimSize,ecc);
    stimSize(r,:)=[0 0 1 1] * ang2pix(gratingSize);
    eccentricity(r)=ang2pix(ecc);
end

%generate positions
for r=1:nrings
    for si=1:stimPerRing
        stimX=eccentricity(r)*cosd(compass(si));
        stimY=eccentricity(r)*sind(compass(si));
        stimLocation(r,si,:)=CenterRectOnPoint(stimSize(r,:),center(1)+stimX,center(2)+stimY);
    end
end

disp('pass_position_generation');


%% exp start
% accuracy
cor = NaN(ntrial,1);

DrawFormattedText(mainwin, ['Block No.', num2str(block)], 'center','center',textcolor);
Screen('Flip',mainwin);
IoWait;
WaitSecs(0.5);
IOPort('Read',P4);


for trial = 1:ntrial
    % prepare and wait to start
    t2 = GetSecs;
    ti=targetindex(block,trial);
    disp(ti);
    tring=ceil(ti/stimPerRing);
    disp(tring);
    tpos=ti-(tring-1)*stimPerRing;
    disp(tpos);
    jit=zeros(nrings,stimPerRing,4);
    for i = 1:nframes
        Screen('FillRect', buffers(i), bgcolor);
        Screen('FillRect', buffers(i), fixcolor, fixRect);
        for ring=1:nrings
            for stimIndex=1:stimPerRing
                if ~mod(i,floor(60/jitfreq))
                    xjit=jpix*(rand*2-1);
                    yjit=(-1)^randi(2)*sqrt(jpix^2-xjit^2);
                    jit(ring,stimIndex,:)=[xjit;yjit;xjit;yjit];
                end
                loc = squeeze(stimLocation(ring,stimIndex,:))+squeeze(jit(ring,stimIndex,:));
                if stimIndex==tpos&&ring==tring
                    Screen('FillArc',buffers(i),red,loc,startAngle,lefthalf);
                    Screen('FillArc',buffers(i),green,loc,startAngle,righthalf);
                else
                    Screen('FillArc',buffers(i),green,loc,startAngle,lefthalf);
                    Screen('FillArc',buffers(i),red,loc,startAngle,righthalf);
                end
                Screen('DrawArc',buffers(i),black,squeeze(stimLocation(ring,stimIndex,:))+jph,startAngle,fullAngle);
            end
        end
    end
    disp(t2-GetSecs);
    Screen('FillRect', mainwin, fixcolor, fixRect);
    Screen('Flip', mainwin);
    
    keypressed =NaN;
    rt = NaN;
    
    IoWait;
    IOPort('Read',P4);
    
    % show
    Screen('DrawTexture', mainwin, buffers(1));
    Screen('Flip', mainwin);
    t1 = GetSecs;
    for i = 2:nframes
        [data,when] = IOPort('Read',P4);
        if ~isempty(data)
            if numel(data) == 1
                if ismember(data,possiblekn)
                    keypressed=data;
                    rt = when - t1;
                    break
                end
            end
        end
        Screen('DrawTexture', mainwin, buffers(i));
        Screen('Flip', mainwin);
    end
    
    respring = find(ismember(possiblekn,keypressed)) - 1;
    if respring == tring
        cor(trial) = 1;
        Screen('FillRect', mainwin, truegreen, fixRect);
    else
        Screen('FillRect', mainwin, truered, fixRect);
        cor(trial) = 0;
    end
    fprintf(outfile,'%s %d %d %d %d %d\n',subj, trial, ti , keypressed , cor(trial), rt);
    Screen('Flip', mainwin);
    if trial == ntrialsperb
        WaitSecs(1);
    end
end
session_end;

    function session_end
        IOPort('Closeall');
        fclose(outfile);
        sca;
        return
    end

    function pixels=ang2pix(ang)
        pixpercm=rect(4)/monitorh;
        pixels=tand(ang/2)*distance*2*pixpercm;
    end

    function stimSize=CorticalScaleFactor(corticalSize,eccentricity)
        M=.065*eccentricity+.054;
        stimSize=M*corticalSize;
    end

    function IoWait
        while 1
            data = IOPort('Read',P4);
            if ~isempty(data)
                break
            end
            WaitSecs(0.001);
        end
    end

end