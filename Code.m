% Clearing workspace and command window; and closing all opened figures
clc;clear all;close all;

% Initializing vectors
loge_test = [];
zcr_test  = [];
framedur  = 30;                       % Frame duration
shift     = 0.5;                      % Shift rate

%   Loading Training Data
load Dataset

%   Reading Testing file
[X,fs] = audioread('FSA1.wav');
X      = X - mean(X);                   % Remove the mean value to normalize
Ns     = length(X);
[m,n]  = size(X);

framelength = framedur * fs/1000;       % Frame length
shiftlength = shift*framedur * fs/1000; % Shift length
frameNb     = floor(Ns - framelength)/shiftlength; % Number of short-time frames

% Loop to get audio signal from one ear 
if (n==2)
    X1 = X(:,1);
    clear X
    X = X1;
end

%   Extracting Feature vectors
for i=1:16
    Feature(i,:) = [loge_Mean(i,:), zcr_Mean(i,:), loge_Std(i,:), zcr_Std(i,:)];
    Label (i,:) = Threshold(i,:);
end

% Taking Transpose of Feature and Label
Feature = Feature.';
Label = Label.';

% Defining Neural Network
net = feedforwardnet(1);

% Training Neural Network
[nett,tr] = train(net,Feature,Label);

% Processing Testing signal
for k=1:frameNb
    frame       = X((k-1)*shiftlength+1:(k-1)*shiftlength+framelength); % frame
    frametime(k)=((k-1)*shiftlength+1)*1000/fs; % Time index for each frame
    loge_test   = [loge_test 10*log10(sum(frame.^2))/framelength];
    zcr_test    = [zcr_test sum(abs(diff(sign(frame))))/(2*framelength)];
end
     
% Extracting Feature vectors from testing signal 
% Taking Mean and standard deviation  of log of energy
loge_test_Mean = sum(loge_test)/length(loge_test);
loge_test_Std = sqrt(sum(loge_test.^2)/length(loge_test)-loge_test_Mean^2);

% Taking Mean and standard deviation  of zero crossings rate
zcr_test_Mean = sum(zcr_test)/length(zcr_test);
zcr_test_Std = sqrt(sum(zcr_test.^2)/length(zcr_test)-zcr_test_Mean^2);    

% Combining all features
test = [loge_test_Mean, zcr_test_Mean, loge_test_Std, zcr_test_Std];
test = test.';

% predicting output
threshold = floor(nett(test));

% Getting VUV using threshold
for k=1:frameNb
    VUV(k) = zcr_test(k)<(threshold/framelength) & loge_test(k)> -0.02; % if short-time ZCR is lower
end

% Plotting Logarithm of Energy
subplot(311),plot(frametime/1000,loge_test,'r','LineWidth',1.5);
title('Logarithm of Energy')
hold on
         
% Plotting Zero crossing rate
subplot(312),plot(frametime/1000,zcr_test,'g','LineWidth',1.5);
title('Zero Crossing Rate')
hold on

% Plotting speech signal and VUV
t = (0 : Ns - 1)/fs;
subplot(313),plot(t,X,'r','LineWidth',1.5);
hold on 
plot(frametime/1000,VUV*max(X))
xlabel('Time(sec)');
legend('speech signal','VUV')
title('Voicing activity in speech signal')