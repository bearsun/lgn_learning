function meridian_loc( debug )
%MERIDIAN_LOC Meridian localizer for quick retinotopy
%   4/14/16 Liwei Sun

% 2 Meridians to localize, repeat 8 times
% each block: 12s-12s-12s-12s
% 2s TRs
% 200 TRs, 7.5min, actual image: 5+32*6 = 197 TRs

%% initialize everything
clc;
global ptb_RootPath %#ok<NUSED>
% global monitorh
% global distance
% global rect

%AssertOpenGL;
%Priority(1);
rng('shuffle');

% monitorh=30; %12;% in cm
% distance=55; %25;% in cm

postscreenwait = 6; % Screen wait after IOport closed
tr = 0;
tbeginning = NaN;
pretr = 5;
blocktr = 6;
nrepeat = 8;
blocktypes = [1 2]; %ring vertical/horizontal
block_pool = repmat(blocktypes,1,nrepeat);
block_seq = reshape([block_pool;zeros(size(block_pool))],1,[]);

block_start = cumsum([1,ones(1,numel(block_seq)-1)*blocktr])+pretr;

nblock = numel(block_seq);


% fixation task
fixcor = 0;
fixfreq = 4; %fixation task happens every 4 TR on average (randomized) >2
permtr = randperm(floor(blocktr*nblock/2))*2;
fixtr = sort(permtr(1:round(blocktr*nblock/fixfreq)))+pretr;
nfix = numel(fixtr);

% geometry
sid = 0;
srect = [0 0 1024 768];
fixsi = 8;
checkfreq = 8; %flickering at 8 Hz
nmeridians = 2;
radius_mask = srect(4); % mask to create meridians

coverangle = 22.5;
StartAngle = [coverangle/2 , coverangle/2-180;
              coverangle/2-90, coverangle/2+90];
arcAngle = 180-coverangle;

% colors
gray = [127 127 127];
white = [255 255 255 255];
% black = [0 0 0];
% green = [0 255 0];
red = [255 0 0 255];

bgcolor = gray;
fixcolor = white;
blinkcolor = red;

% button
if debug
    kn0 = KbName('KP_Insert');
    kn1 = KbName('KP_End');
    kn2 = KbName('KP_Down');
    kn3 = KbName('KP_Next');
    BUFFER = [];
    BUTTONS = [kn1,kn0,kn2,kn3];
    CODES = 49:52;
else
    %% buttonbox setting
    buttons = 49:52; %top, left, right, bottom
    trigger=53;
    kn0 = buttons(2); %left
    kn1 = buttons(1); %top
    kn2 = buttons(4); %right
    kn3 = buttons(3); %bottom
end

possiblekn = [kn0,kn1,kn2,kn3];

%% open scanner communication
if ~debug
    IOPort('Closeall');
    P4 = IOPort('OpenSerialPort', '/dev/ttyUSB0','BaudRate=115200'); %open port for receiving scanner pulse
    fRead = @() ReadScanner(P4);
else
    fRead = @() ReadFakeTrigger;
    tr_tmr = timer('TimerFcn',@SetTrigger,'Period',2,'ExecutionMode','fixedDelay','Name','tr_timer');
end

subj=input('subject?','s');

log_data = [pwd,'/data/meridian-loc-',subj,'-log'];
design_data = [pwd,'/data/meridian-loc-',subj,'-design'];

outdesign = fopen(design_data,'w');
fprintf(outdesign,'%s\t %s\n','block_start','block_seq');
fprintf(outdesign,'%d\t %d\n', [block_start;block_seq]);
fclose(outdesign);

outlog = fopen(log_data,'w');
fprintf(outlog,'%s\t %s\n','when','what');

%% initialize window
[mainwin,rect] = Screen('OpenWindow', sid, bgcolor, srect);

% basic screen
ifi = Screen('GetFlipInterval', mainwin);
%center = [(rect(3)-rect(1))/2, (rect(4)-rect(2))/2];
fixRect = CenterRect([0 0 fixsi fixsi], rect);
arcrect = CenterRect([-radius_mask,-radius_mask,radius_mask,radius_mask],rect);
%% frequency 8 Hz
checkt = 1/checkfreq;
waitframes = round(checkt/ifi);

%% build textures
npics = 8;
orig_pics = cell(8,1);
for p = 1:npics
    orig_pics{p} = imread(['./checkerboards/bildAll',num2str(p),'.png']);
end

% construct stimuli
merhandle = NaN(nmeridians,npics);
for m = 1:nmeridians
    merhandle(m,:) = domeridian(StartAngle(m,:),arcAngle);
end


tested = 1; %flag for response

% for debug
if debug
    start(tr_tmr);
end

TRWait(block_start(1));
vbl = Screen('Flip', mainwin);
for block = 1:nblock
    
    btype = block_seq(block);
    if btype
        dp = randi(npics);
    end
    
    while tr < block_start(block) + blocktr
        if btype
            Screen('DrawTexture', mainwin, merhandle(btype,dp));
            dp = mod(dp,8)+1;
        end
        Screen('FillRect', mainwin, fixcolor, fixRect);
        vbl = Screen('Flip', mainwin, vbl + (waitframes - 0.5) * ifi);
        data = fRead();
        
        % test fixation
        if ~tested
            if ismember(data,possiblekn)
                fixcor = fixcor + 1;
                tested = 1;
            end
            
            if tr > test_end
                tested = 1;
            end
        end
        
        % start test fixation
        if ismember(tr, fixtr)
            fixcolor = blinkcolor;
            fixtr(ismember(fixtr,tr))=[];
            
            tested = 0;
            test_end = tr+1;
        else
            fixcolor = white;
        end
    end
end

%% save and close everything
fRead();
fclose(outlog);
if debug
    StopTimer;
else
    IOPort('Closeall');
end

WaitSecs(postscreenwait);
sca;

fprintf('Accuracy: %d\n',fixcor/nfix);

%% function to create rings based on original textures
    function mhandle = domeridian(Start,Arc)
        mhandle = NaN(1,npics);
        for pic = 1:npics
            mhandle(pic) = Screen('MakeTexture',mainwin,orig_pics{pic});
            Screen('FillArc', mhandle(pic), gray,arcrect,Start(1),Arc);
            Screen('FillArc', mhandle(pic), gray,arcrect,Start(2),Arc);
        end
    end

%% function wrapper for IOPort('Read'),also counting the total TRs
    function [data, when] = ReadScanner(Port)
        [data, when] = IOPort('Read',Port);
        
        if ~isempty(data)
            tr=tr+sum(data==trigger);
            if tr == 1
                tbeginning = when;
            end
            fprintf('%d\t %d\n',when-tbeginning,tr);
            fprintf(outlog, '%d\t %d\n',when-tbeginning,tr);
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
end

