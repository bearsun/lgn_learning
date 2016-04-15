function rg=IsoRG_fmri
% function for isolumiance between red/green
monitor=0; %primary monitor for display
rect = [0 0 1024 768];

subj = input('subject?: ','s');

IOPort('Closeall');
P4 = IOPort('OpenSerialPort', '/dev/ttyUSB0','BaudRate=115200');

rightKey = 52;
leftKey = 50;
upKey = 49;
keysWanted=[rightKey leftKey upKey]; %NB upkeyremoved for this block
nRepeats=10;
maxRGB=255;

incrementI=0.01;%increment for manipulating I
inc=[incrementI,-incrementI,0];
HSVred=0/360;
HSVgreen=120/360;
% HSVblue=240/360;
% HSVyellow=60/360;
HSVmiddle=180/360;
referenceColor=maxRGB*hsv2rgb([HSVmiddle .5 .5]); %blu
satu=.5;

reds=[ones(nRepeats,1)*[HSVred,satu],rand(nRepeats,1)];
greens=[ones(nRepeats,1)*[HSVgreen,satu],rand(nRepeats,1)];
% yellows=[ones(nRepeats,1)*[HSVyellow,satu],rand(nRepeats,1)];
% blues=[ones(nRepeats,1)*[HSVblue,satu],rand(nRepeats,1)];

newreds=NaN(nRepeats,3);
newgreens=NaN(nRepeats,3);
% newyellows=NaN(nRepeats,3);
% newblues=NaN(nRepeats,3);

fixColor=[255 255 255];
bgColor=[127 127 127];

fixRect=[0 0 4 4];
colorRect=[0 0 100 100];
frameRect=[0 0 110 110];

flickerRate = 30; %in Hz
[w,screenSize]=Screen('OpenWindow',monitor,bgColor,rect);
HideCursor;

refreshRate=Screen('FrameRate', w);
flipsToWait=refreshRate/flickerRate;

fixRect=CenterRect(fixRect, screenSize);
colorRect=CenterRect(colorRect, screenSize);
frameRect=CenterRect(frameRect, screenSize);

%draw starting screen
Screen('FrameRect',w,fixColor,frameRect);
Screen('FillRect',w,fixColor,fixRect);
Screen('Flip',w,[],0);
IoWait;

for repeat=1:nRepeats
    red=reds(repeat,:);
    green=greens(repeat,:);
%    yellow=yellows(repeat,:);
%    blue=blues(repeat,:);
    
    red=testone(red);
    green=testone(green);
%    yellow=testone(yellow);
%    blue=testone(blue);

    newreds(repeat,:)=red;
    newgreens(repeat,:)=green;
%    newyellows(repeat,:)=yellow;
%    newblues(repeat,:)=blue;
end
rg=round(hsv2ptb([mean(newreds(4:10,:),1);mean(newgreens(4:10,:),1)]));
% save;
disp('rg:');
disp(rg);
Screen('CloseAll');
IOPort('Closeall');
save([pwd,'/subinfo/',subj,'_fmri.mat'],'rg');

    function color=testone(color)
        while 1
            %draw color
            Screen('FrameRect',w,fixColor,frameRect);
            Screen('FillRect',w,hsv2ptb(color),colorRect);
            Screen('FillRect',w,fixColor,fixRect);
            for h=1:flipsToWait-1
                Screen('Flip',w,[],2); %wait for desired # refreshes-1
            end
            Screen('Flip',w,[],0);%flip one more time and clear the buffer
            %draw reference
            Screen('FrameRect',w,fixColor,frameRect);
            Screen('FillRect',w,referenceColor,colorRect);
            Screen('FillRect',w,fixColor,fixRect);
            for h=1:flipsToWait-1
                Screen('Flip',w,[],2); %wait for desired # refreshes-1
            end
            Screen('Flip',w,[],0);%flip one more time and clear the buffer
            
            data = IOPort('Read',P4);
            if ~isempty(data)
                if numel(data) == 1
                    if ismember(data,keysWanted)
                        color(3)=color(3)+inc*ismember(keysWanted,data)';
                        if color(3)>1
                            color(3)=1;
                        elseif color(3)<0
                            color(3)=0;
                        end
                        if data == keysWanted(3)
                            break;
                        end
                    end
                end
            end
            
        end
        Screen('FrameRect',w,fixColor,frameRect);
        Screen('FillRect',w,fixColor,fixRect);
        Screen('Flip',w,[],0);
        IoWait;
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
    function ptb=hsv2ptb(hsv)
        ptb=maxRGB*hsv2rgb(hsv);
    end
end
