function ring_loc( debug )
%RING_LOC localizer for rings in lgn study
%   4/14/16 Liwei
% 3 rings to localize, repeat 8 times
% each block: 12s-12s-12s-12s
% 2s TRs
% 296TRs, 10.5 min, actual images: 5 + 48 * 6 = 293 TRs

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

tr = 0;
pretr = 5;
blocktr = 6;
nrepeat = 8;
blocktypes = [1 2 3]; %ring 1,2,3
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
checkfreq = 16; %flickering at 8 Hz
nrings = 3;


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
    kn2 = buttons(3); %right
    kn3 = buttons(4); %bottom
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

log_data = [pwd,'/data/ring-loc-',subj,'-log'];
design_data = [pwd,'/data/ring-loc-',subj,'-design'];

outdesign = fopen(design_data,'w');
fprintf(outdesign,'%s\t %s\n','block_start','block_seq');
fprintf(outdesign,'%d\t %d\n', [block_start;block_seq]);
fclose(outdesign);

outlog = fopen(log_data,'w');
fprintf(outlog,'%s\t %s\n','when','what');

%% initialize window
[mainwin,rect] = Screen('OpenWindow', sid, bgcolor, srect);
Screen('BlendFunction', mainwin, GL_SRC_ALPHA,GL_ONE_MINUS_SRC_ALPHA);

% basic screen
ifi = Screen('GetFlipInterval', mainwin);
%center = [(rect(3)-rect(1))/2, (rect(4)-rect(2))/2];
fixRect = CenterRect([0 0 fixsi fixsi], rect);

%% frequency 8 Hz
checkt = 1/checkfreq;
waitframes = round(checkt/ifi);

%% build textures
npics = 8;
orig_pics = cell(8,1);
for p = 1:npics
    orig_pics{p} = imread(['./checkerboards/bildAll',num2str(p),'.png']);
end

% construct stimuli %% same as in lgn_search!!!!!
corticalStimSize= 6.5;%in mm
proximalStimDist=.8;% in degrees!!
stimSeparation= 1.2; %in degrees, to be scaled
jitterDistance=.12;%in degrees; to be scaled

%generate radius
ringRadius = NaN(nrings,2); %need an inner radius and an outer radius
jpix=ang2pix(jitterDistance);

for r=1:nrings
    ecc=proximalStimDist+stimSeparation*r^2;
    gratingSize=CorticalScaleFactor(corticalStimSize,ecc);
    eccentricity=ang2pix(ecc);
    stimRadius = ang2pix(gratingSize)/2;
    ringRadius(r,1) = eccentricity - stimRadius - jpix;
    ringRadius(r,2) = eccentricity + stimRadius + jpix;
end

ringhandle = NaN(nrings,npics);
for ring = 1:nrings
    ringhandle(ring,:) = doring(ringRadius(ring,1),ringRadius(ring,2));
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
            Screen('DrawTexture', mainwin, ringhandle(btype,dp));
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

sca;

fprintf('Accuracy: %d\n',fixcor/nfix);

%% function to create rings based on original textures
    function ringhandle = doring(inner,outer)
        % Create circular aperture for the alpha-channel:
        [x,y]=meshgrid(1:rect(3), 1:rect(4));
        dist = (x-mean(mean(x))).^2 + (y-mean(mean(y))).^2;
        ringalpha = ((dist <= (outer)^2) & (dist >= (inner)^2)).*255;
        
        ringhandle = NaN(1,npics);
        for pic = 1:npics
            rgbapic = cat(3,orig_pics{pic},ringalpha);
            ringhandle(pic) = Screen('MakeTexture',mainwin,rgbapic);
        end
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

