clear;
clc;
close all;
% Chord Recognition 
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
    bands(i+1) = +0+(data(i) + data(i+1))/2;
end
bands(89) = 4500;

%%
fileName = 'Records/Es.wav';    % File name 
[y, Fs] = wavread(fileName);    % Read audio file
wavplay(y,Fs)                   % Play audio file
y = (y(:,1) + y(:,2))*2;        % Decrease 2 channels to 1 

fy = abs(fft(y.*Hamming(length(y))));   % Hamming + FFT

frameLength = length(y);
for i=1:88                      % Create Bands
    freqBand(i) = mean(fy( round(bands(i)/(Fs/frameLength) ):round(bands(i+1)/(Fs/frameLength))))^2;
end


stem(freqBand)
xlabel('88-Keys')
ylabel('Amplitude')
title('Bands')

f = linspace(1,Fs,length(fy));

m = find(freqBand == max(freqBand(:)));
index = 1;
for i = 1:88                    % Find Harmonics
    if(freqBand(i) > 2000)
        harmonics(index) = i;
        index = index+1;
    end
end
    
rate = freqBand(m) ./ freqBand(harmonics);  % Put a rate
index = 1;                                  % According to rate
for i = 1:length(harmonics)                 % determine the harmonics 
    if(rate(i) < 5)
        newHarmonics(index) = harmonics(i);
        index = index +1;
    end
end
harmonics = newHarmonics;
fharmonics = freqBand(harmonics);
 
disp(names(harmonics,:));

% Check octaves
index = 1;
remove = 0;
for i = 1:length(harmonics)
    for j = i:length(harmonics)
        if(i ~= j)
            a1 = names(harmonics(i),1:2);
            a2 = names(harmonics(j),1:2);
            
            if(a1 == a2)
                if(fharmonics(j)> fharmonics(i))
                    % Do nothing
                else
                    remove(index) = j;
                    index = index +1;
                end 
                
            end
        end
    end
end

if(remove(1) ~= 0)
    harmonics(remove) = [];     % Remove if it is a harmonic of a note
end

disp('*****')               % Print
for i = 1:length(harmonics)
    if(harmonics(i) > 0)
        disp(names(harmonics(i),:));
    end
end


