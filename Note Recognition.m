clear;
clc;
close all;
% Note Recognition
%% Note Initialization
mainNames = char('C', 'C#', 'D', 'D#', 'E', 'F', 'F#', 'G', 'G#', 'A', 'A#', 'B');
names = char('A0', 'A#0' ,'B0');
A0 = 27.5;
ind = 4;
for i = 0:87
    data(i+1) = A0 * (2^(1/12))^i;
    if i>2
         a = [mainNames( rem((ind - 4), 12)+1,:)  num2str(fix((ind - 4)/12)+1)];
         names = char(names, a);
         ind = ind + 1;
    end
end

% Band Initialization
bands(1) = 20;
for i = 1:87
    bands(i+1) = +7+(data(i) + data(i+1))/2;
end
bands(89) = 4500;

%%
fileName = 'Records/Plug in baby.wav';  % File name
[y, Fs] = wavread(fileName);            % Read audio file
y = (y(:,1) + y(:,2))*4;                % Decrease 2 channels to 1 
%y = (y(:,1));
y(1:2:end) = 0;                         % Do decimation

frameLength = 4410*2; % 2 Güzel oldu

endPart = frameLength*ceil(length(y)/frameLength);  % complete the last frame
y(length(y)+1 : endPart) = 0;

f = linspace(1,Fs,frameLength);

%%
harmonics = 0;
for i = 1:round(length(y)/frameLength)      % For each frame
    % Divide audio into frames
    frames(i, 1:frameLength) = y( (frameLength*(i-1)+1):(frameLength*i) )';
    
    frame = y( (frameLength*(i-1)+1):(frameLength*i) )';
    frame = frame .* Hamming(length(frame))';   % Hamming Window
    fframe = abs(fft(frame));                   % FFT
    
    fp = sum(fframe);
    p = sum(abs(frame));
    b = true;
    if(p < 200 || fp < 1000) % Put a threshold for processing
       b = false; 
    end
    
    % Bands
    for i=1:88
        freqBand(i) = mean(fframe( round(bands(i)/(Fs/frameLength) ):round(bands(i+1)/(Fs/frameLength))))^2;
    end
    
    % Plotting
    subplot(3,1,1)
    stem(freqBand)
    subplot(3,1,2)
    plot(fframe)
    xlim([0,500])
    subplot(3,1,3)
    plot(frame)
    ylim([-1 1])
    
    hold off
    pause(0.1)
    wavplay(frame,Fs)
    
    
    % Desicion
    m = find(freqBand == max(freqBand(:)));
    if(b) disp(names(m,:));     % Print the result
    else disp('.'); end
    
    if(b)
        index = 1;
        for i = 1:88
            if(freqBand(i) > 2000)
                harmonics(index) = i;
                index = index+1;
            end 
        end
    
    end

end
