function lgn_search_behav( debug )
%LGN_SEARCH_BEHAV training code for lgn_learning project
%   4/14/16 by Liwei
% in fmri, we have 5 runs of training with 5 X 32 = 160 trials in total
% here, just to keep them to be the same, we also have 5 blocks with 32
% trials per block, so 160 trials in total

%% initialize everything
clc;
global ptb_RootPath %#ok<NUSED>

rng('shuffle');
sid = 0;
rect = [0 0 1024 768];

monitorh=30; %12;% in cm
distance=57; %25;% in cm

subj=input('subject?','s');
session = input('session?');

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
path_color = [pwd,'/subinfo/',subj,'.mat'];
rg = load(path_color);
red = rg.rg(1,:);
green = rg.rg(2,:);

fixcolor = white;
textcolor = black;
bgcolor = gray;

%% keys
kesc = KbName('Escape');
kspace = KbName('space');
kreturn = KbName('Return');
kback = KbName('BackSpace');
kn0 = KbName('KP_Insert');
kn1 = KbName('KP_End');
kn2 = KbName('KP_Down');
kn3 = KbName('KP_Next');

possiblekn = [kn0,kn1,kn2,kn3];

% parameters
nblocks = 5;
ntrialsperb = 32;
ntrials = nblocks * ntrialsperb;
searchtime = 4; % 4s to search
nframes = searchtime * 60;

nrings=3;
stimPerRing=8;
nballs=nrings*stimPerRing;

%% random target position with no repeat
targetoptions=[1:nballs,zeros(1,stimPerRing)]; %catch trials 25%, each ring 25%
repeat=1;
while repeat
    targetindex=BalanceTrials(ntrials,1,targetoptions);
    targetindex=(reshape(targetindex,ntrialsperb,nblocks))' ;
    repeat=0;
    for i=2:ntrialsperb
        w=targetindex(:,i-1)==targetindex(:,i);
        if any(w)&&any(targetindex(w,i-1)~=0)
            repeat=1;
            break;
        end
    end
end


%% open files
path_data = [pwd,'/data/training-',subj,'-',num2str(session)];

outfile = fopen(path_data,'w');
fprintf(outfile,'%s\t %s\t %s\t %s\t %s\t %s\t %s\n','subject' ,'session' ,'trial','targetindex' , 'keypressed' , 'cor','rt');

abbreviatedFilename=[subj,'_',datestr(now,'mmdd')];
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Calibration
if ~debug
    Eyelink('Shutdown');
    Eyelink('Initialize');
    HideCursor;
    
    Eyelink('StartSetup')
    pause(2)
    
    
    whichKey=0;
    
    keysWanted=[kspace kreturn kback];
    FlushEvents('KeyDown');
    while 1
        pressed = 0;
        while pressed == 0
            [pressed, ~, kbData] = KbCheck;
        end;
        
        for keysToCheck = 1:length(keysWanted)
            if kbData(keysWanted(keysToCheck)) == 1
                
                keyPressed = keysWanted(keysToCheck);
                if keyPressed == kback
                    whichKey=9;
                    FlushEvents('KeyDown');
                    WaitSecs(.1)
                elseif keyPressed == kspace
                    whichKey=1;
                    FlushEvents('KeyDown');
                    WaitSecs(.1)
                elseif keyPressed == kreturn
                    whichKey=5;
                    FlushEvents('KeyDown');
                    WaitSecs(.1)
                else
                end
                FlushEvents('KeyDown');
                
            end;
        end;
        
        if whichKey == 1
            whichKey=0;
            [~, tx, ty] = Eyelink('TargetCheck');
            tx=tx*.64;
            ty=ty*.64;
            Screen('FillRect', mainwin ,black, [tx-20 ty-5 tx+20 ty+5]);
            Screen('FillRect', mainwin ,black, [tx-5 ty-20 tx+5 ty+20]);
            Screen('Flip', mainwin);
        elseif whichKey == 5
            whichKey=0;
            Eyelink('AcceptTrigger');
        elseif whichKey == 9
            break;
        end
    end;
    status = Eyelink('OpenFile',abbreviatedFilename);
    if status
        error(['openfile error, status: ', num2str(status)]);
    end
    Eyelink('StartRecording');
end

%% exp start
if ~debug
    Eyelink('Message','session_start');
end

for block = 1:nblocks
    DrawFormattedText(mainwin, ['Block No.', num2str(block)], 'center','center',textcolor);
    Screen('Flip',mainwin);
    KbStrokeWait;
    
    for trial = 1:ntrialsperb
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
        
        while 1 %wait to start
            [keyIsDown, ~, keyCode] = KbCheck;
            FlushEvents('keyDown');
            if keyIsDown
                nKeys = sum(keyCode);
                if nKeys == 1
                    if keyCode(kesc)
                        session_end;return
                    elseif keyCode(kspace)
                        break;
                    end
                end
            end
        end
        
        keypressed =NaN;
        rt = NaN;
        if ~debug
            Eyelink('Message','trial_start');
        end
        
        % show
        Screen('DrawTexture', mainwin, buffers(1));
        Screen('Flip', mainwin);
        t1 = GetSecs;
        for i = 2:nframes
            [keyIsDown, secs, keyCode] = KbCheck;
            FlushEvents('keyDown');
            if keyIsDown
                nKeys = sum(keyCode);
                if nKeys == 1
                    if keyCode(kesc)
                        session_end;return
                    elseif any(keyCode(possiblekn))
                        keypressed=find(keyCode);
                        rt = secs - t1;
                        break;
                    end
                end
            end
            Screen('DrawTexture', mainwin, buffers(i));
            Screen('Flip', mainwin);
        end
        
        if ~debug
            Eyelink('Message','trial_end');
        end
        
        respring = find(ismember(possiblekn,keypressed)) - 1;
        if respring == tring
            cor = 1;
            Screen('FillRect', mainwin, truegreen, fixRect);
        else
            Screen('FillRect', mainwin, truered, fixRect);
            cor = 0;
        end
        fprintf(outfile,'%s\t %d\t %d\t %d\t %d\t %d\t %d\n',subj, session ,trial, targetindex , keypressed , cor, rt);
        Screen('Flip', mainwin);
        if trial == ntrialsperb
            WaitSecs(1);
        end
    end
end
session_end;

    function session_end
        if ~debug
            Eyelink('Message','session_end');
            Eyelink('Stoprecording');
            Eyelink('CloseFile');
            Eyelink('ReceiveFile');
        end
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

end