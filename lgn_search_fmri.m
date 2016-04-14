function lgn_search_fmri( debug )
%LGN_SEARCH_FMRI Code for pre and post scan for lgn_learning project
%4/13/16 by Liwei
% use 3 rings since only have 4 buttons
% 5 runs per session
% 32 trials per run
% TR = 2s, use 2mm resolution to cover LGN and V1
% ISI could be 2/3/4/5 TRs
% response period is included in ISI period (1TR to response)
% each sti presentation for 2 TRs

% 183 TRs, 6 min, actual images: 5 + 32 *(2 + 3.5) = 181

%% initialize everything
clc;
global ptb_RootPath %#ok<NUSED>
global monitorh
global distance
global rect

%AssertOpenGL;
%Priority(1);
rng('shuffle');

monitorh=30; %12;% in cm
distance=55; %25;% in cm
sid = 0;
srect = [0 0 1024 768];
rng('shuffle');
ntrial = 32;

%% geometry
fixsi = 8;
righthalf = 180;
lefthalf = -180;
startAngle = 0;
fullAngle=360;

%% mri param
tr = 0; %TR counter
pretr = 5; % so the 6th volume is actually the 1st volume of the 1st trial
%posttr = 2;
pool_isi = 2:5; % jitter 2 3 4 5 TRs
temp_isi = repmat(pool_isi, 1, ntrial/numel(pool_isi));
isi = temp_isi(randperm(ntrial));
searchtime = 2;% in trs
nframes = searchtime * 2 * 60;

tstart_tr = cumsum([1,isi(1:end-1)+searchtime])+pretr;

%% color
gray = [127 127 127];
white = [255 255 255];
black = [0 0 0];
red = [175,87,87];
green = [69,138,69];
truegreen = [0 255 0];
truered = [255 0 0];

fixcolor = white;
%textcolor = black;
bgcolor = gray;

% button
if debug
    kn0 = KbName('KP_Insert');
    kn1 = KbName('KP_End');
    kn2 = KbName('KP_Down');
    kn3 = KbName('KP_Next');
    BUFFER = [];
    BUTTONS = [kn1,kn0,kn2,kn3];
    buttons = BUTTONS;
    CODES = 49:52;
else
    %% buttonbox setting
    buttons = 49:52; %top, left, right, bottom
    trigger=53;
    kn0 = buttons(2); %left
    kn1 = buttons(1); %top
    kn2 = buttons(3); %right
    kn3 = buttons(4); %bottom
end

possiblekn = [kn0,kn1,kn2,kn3];

%% target & answer sequence
nrings=3;
stimPerRing=8;
nballs=nrings*stimPerRing;

targetoptions=[1:nballs,zeros(1,stimPerRing)]; %catch 25% 1,2,3 each 25%
targetindex = targetoptions(randperm(numel(targetoptions)));
temp_ring = ceil(targetindex/8);
corkeyarr=(temp_ring == 1) * kn1;%1st ring
corkeyarr(temp_ring == 2) = kn2;
corkeyarr(temp_ring == 3) = kn3;
corkeyarr(temp_ring == 0) = kn0;

%% open scanner communication
if ~debug
    IOPort('Closeall');
    P4 = IOPort('OpenSerialPort', '/dev/ttyUSB0','BaudRate=115200'); %open port for receiving scanner pulse
    fRead = @() ReadScanner(P4);
else
    fRead = @() ReadFakeTrigger;
    tr_tmr = timer('TimerFcn',@SetTrigger,'Period',2,'ExecutionMode','fixedDelay','Name','tr_timer');
end

%% open files
subj=input('subject?','s');
session = input('session? (pre/post)','s');
run=input('run? ');
path_data = [pwd,'/data/data-',subj,'-',session,'-run',num2str(run)];
log_data = [pwd,'/data/data-',subj,'-',session,'-run',num2str(run),'-log'];
design_data = [pwd,'/data/data-',subj,'-',session,'-run',num2str(run),'-design'];

outdesign = fopen(design_data,'w');
fprintf(outdesign,'%s\t %s\t %s\n','tstart','targetindex','corkey');
fprintf(outdesign,'%d %d %d\n', [tstart_tr;targetindex;corkeyarr]);
fclose(outdesign);

outlog = fopen(log_data,'w');
fprintf(outlog,'%s\t %s\n','when','what');

outfile = fopen(path_data,'w');
fprintf(outfile,'%s\t %s\t %s\t %s\t %s\t %s\n','subject' ,'session' ,'trial','targetindex' , 'keypressed' , 'cor');

%% initialize window
[mainwin,rect] = Screen('OpenWindow', sid, bgcolor, srect);

% open buffers
buffers = NaN(nframes,1);
for fr = 1:nframes
    [buffers(fr),~] = Screen('OpenOffscreenWindow', mainwin, bgcolor);
end

%% build up position arrays
% basic screen
center = [(rect(3)-rect(1))/2, (rect(4)-rect(2))/2];
fixRect = CenterRect([0 0 fixsi fixsi], rect);

% construct stimuli
corticalStimSize= 6.5;%in mm
proximalStimDist=.8;% in degrees!!
stimSeparation= 1.2; %in degrees, to be scaled
jitterDistance=.12;%in degrees; to be scaled
jitfreq = 5;

%empty variable
stimLocation = NaN(nrings,stimPerRing,4);
stimSize=NaN(nrings,4);
eccentricity=NaN(nrings);

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

% accuracy
cor = NaN(ntrial,1);

% for debug
if debug
    start(tr_tmr);
end

%% exp start
for trial = 1:ntrial
    % prep buffer
    [preptime,tring,ti] = PrepareTrial(trial);
    fprintf('Prep time: %d\n', preptime);
    
    bresponded = 0; %responded flag
    keypressed = NaN;
    
    % wait until certain TR to start
    TRWait(tstart_tr(trial));
    
    % present stimuli and wait for resp
    Screen('DrawTexture', mainwin, buffers(1));
    tstart = Screen('Flip', mainwin);
    for i = 2:nframes
        if ~bresponded
            [data, when] = fRead();
            if ~isempty(data)
                if ismember(data,buttons)
                    keypressed = buttons(ismember(buttons,data));
                    rt = when-tstart;
                    bresponded = 1;
                    fprintf('BUTTON RECIEVED: %d @ %d\n',keypressed,rt*1000);
                    fprintf(outlog,'BUTTON RECIEVED: %d @ %d\n',keypressed,rt);
                end
            end
        end
        Screen('DrawTexture', mainwin, buffers(i));
        Screen('Flip', mainwin);
        
        if tr > tstart_tr(trial) + 1
            break
        end
    end
    
    % if not responded yet, wait 2s with blank screen for response
    while ~bresponded && tr < tstart_tr(trial) + 3
        Screen('FillRect', mainwin, fixcolor, fixRect);
        Screen('Flip', mainwin);
        [data, when] = fRead();
        if ~isempty(data)
            if ismember(data,buttons)
                keypressed = buttons(ismember(buttons,data));
                rt = when-tstart;
                bresponded = 1;
                fprintf('BUTTON RECIEVED: %d @ %d\n',keypressed,rt*1000);
                fprintf(outlog,'BUTTON RECIEVED: %d @ %d\n',keypressed,rt);
            end
        end
    end
    
    % feedback and save
    respring = find(ismember(possiblekn,keypressed)) - 1;
    if respring == tring
        cor(trial) = 1;
        Screen('FillRect', mainwin, truegreen, fixRect);
    else
        Screen('FillRect', mainwin, truered, fixRect);
        cor(trial) = 0;
    end
    fprintf(outfile,'%s\t %s\t %d\t %d\t %d\t %d\n', ...
        subj, session, trial, ti, keypressed, cor(trial));
    Screen('Flip', mainwin);
    WaitSecs(.5);
    
    Screen('FillRect', mainwin, fixcolor, fixRect);
    Screen('Flip', mainwin);
end

%% save and close everything
fRead();
fclose(outfile);
fclose(outlog);
if debug
    StopTimer;
else
    IOPort('Closeall');
end

sca;

fprintf('Accuracy: %d\n',mean(cor));

%% function to prepare each stimuli presentation (4s, 2TRs)
    function [preptime,tring,ti] = PrepareTrial(trial)
        % prepare and wait to start
        t_prep_start = GetSecs;
        
        ti=targetindex(trial);
        fprintf('Target Index: %d\n', ti);
        
        tring=ceil(ti/stimPerRing);
        fprintf('Ring: %d\n', tring);
        
        tpos=ti-(tring-1)*stimPerRing;
        fprintf('Position: %d\n',tpos);
        
        jit=zeros(nrings,stimPerRing,4);
        
        for frame = 1:nframes
            Screen('FillRect', buffers(frame), bgcolor);
            Screen('FillRect', buffers(frame), white, fixRect);
            for ring=1:nrings
                for stimIndex=1:stimPerRing
                    if ~mod(frame,floor(60/jitfreq))
                        xjit=jpix*(rand*2-1);
                        yjit=(-1)^randi(2)*sqrt(jpix^2-xjit^2);
                        jit(ring,stimIndex,:)=[xjit;yjit;xjit;yjit];
                    end
                    loc = squeeze(stimLocation(ring,stimIndex,:))+squeeze(jit(ring,stimIndex,:));
                    %                 r = (loc(3)-loc(1))/2;
                    if stimIndex==tpos&&ring==tring
                        % Draw Target
                        Screen('FillArc',buffers(frame),red,loc,startAngle,lefthalf);
                        Screen('FillArc',buffers(frame),green,loc,startAngle,righthalf);
                    else
                        % Draw Distractors
                        Screen('FillArc',buffers(frame),green,loc,startAngle,lefthalf);
                        Screen('FillArc',buffers(frame),red,loc,startAngle,righthalf);
                    end
                    Screen('DrawArc',buffers(frame),black,squeeze(stimLocation(ring,stimIndex,:))+jph,startAngle,fullAngle);
                    %                 Screen('DrawLine', buffers(i), black, loc(1), loc(2) + r, loc(3), loc(2) + r);
                    %                 Screen('DrawLine', buffers(i), black, loc(1) + r, loc(2), loc(1) + r, loc(4));
                end
            end
        end
        Screen('FillRect', mainwin, fixcolor, fixRect);
        preptime = GetSecs - t_prep_start;
    end

%% function wrapper for IOPort('Read'),also counting the total TRs
    function [data, when] = ReadScanner(Port)
        [data, when] = IOPort('Read',Port);
        
        if ~empty(data)
            tr=tr+sum(data==trigger);
            fprintf('%d\t %d\n',when,tr);
            fprintf(outlog, '%d\t %d\n',when,tr);
        end
    end

%% function to wait until certain trs
    function TRWait(t)
        while t > tr
            fRead();
            WaitSecs(.001);
        end
    end

%% function to read the faked triggers from buffer
    function [data,when] = ReadFakeTrigger
        data = BUFFER;
        BUFFER = [];
        [~,~,kDown] = KbCheck;
        b = logical(kDown(BUTTONS));
        BUFFER = [BUFFER CODES(b)];
        when = GetSecs;
    end
%-----------------------------------------------------------------------------%
%% function to simulate triggers in the buffer, handle for timer
    function SetTrigger(varargin)
        tr=tr+1;
        fprintf('TR TRIGGER %d\n',tr);
        BUFFER = [BUFFER 53];
    end
%-----------------------------------------------------------------------------%
%% function to stop the simulation
    function StopTimer
        if isobject(tr_tmr) && isvalid(tr_tmr)
            if strcmpi(tr_tmr.Running,'on')
                stop(tr_tmr);
            end
            delete(tr_tmr);
        end
    end
%-----------------------------------------------------------------------------%
    function pixels=ang2pix(ang)
        pixpercm=rect(4)/monitorh;
        pixels=tand(ang/2)*distance*2*pixpercm;
    end
%-----------------------------------------------------------------------------%
    function stimSize=CorticalScaleFactor(corticalSize,eccentricity)
        M=.065*eccentricity+.054;
        stimSize=M*corticalSize;
    end
end

