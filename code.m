%EEA Hackathon, IIT Madras.
%Team of Joseph Samuel, Green Rosh KS, Govind K and Vipin S.
%Problem Statement 1
%March 2, 2014

%Based on the algorithm by Rabiner and Sambur(1975). 
%More references in the README file.

Fs = 44100; %Sampling Frequency which will be used throughout.

disp('Recording starts(for 5s)')
%Making the recorder object at 8 bits per sample and mono
recc = audiorecorder(Fs,8,1);
%Recording from mic
recordblocking(recc,5);
disp('End of Recording.')

%Getting samples from the recorded audio
x=getaudiodata(recc);

%Processing starts now.
Ini = 0.1; %Initial silence duration in seconds
Ts = 0.01; %Frame width in seconds
Tsh = 0.005; %Frame shift in seconds


%Initialising counters to be used later
counter1 = 0; 
counter2 = 0;
counter3 = 0;
counter4 = 0;
ZCRCountf = 0; %Stores forward count of crossing rate > IZCT
ZCRCountb = 0; %As above, for backward count
ZTh = 40; %Zero crossing comparison rate for threshold

w_sam = fix(Ts*Fs); %No of Samples/window 
o_sam = fix(Tsh*Fs); %No of samples/overlap
lengthX = length(x);
segs = fix((lengthX-w_sam)/o_sam)+1; %Number of segments in speech signal
sil = fix((Ini-Ts)/Tsh)+1; %Number of segments in silent period
win = hamming(w_sam);

Limit = o_sam*(segs-1)+1; %Start index of last segment
FrmIndex = 1:o_sam:Limit; %Vector containing starting index 
                                      %for each segment
ZCR_Vector = zeros(1,segs); %Vector to hold zero crossing rate 
                                      %for all segments
                                     
%Below code computes and returns zero crossing rates for all segments in
%speech sample 
for t = 1:segs
    ZCRCounter = 0; 
    nextIndex = (t-1)*o_sam+1;
    for r = nextIndex+1:(nextIndex+w_sam-1)
        if (x(r) >= 0) && (x(r-1) >= 0)
         
        elseif (x(r) >= 0) && (x(r-1) < 0)
         ZCRCounter = ZCRCounter + 1;
        elseif (x(r) < 0) && (x(r-1) < 0)
         
        elseif (x(r) < 0) && (x(r-1) >= 0)
         ZCRCounter = ZCRCounter + 1;
        end
    end
    ZCR_Vector(t) = ZCRCounter;
end


%Below code computes and returns frame energy for all segments in speech
%sample
Erg_Vector = zeros(1,segs);
for u = 1:segs
    nextIndex = (u-1)*o_sam+1;
    Energy = x(nextIndex:nextIndex+w_sam-1).*win;
    Erg_Vector(u) = sum(abs(Energy));
end

%Energy calculations
IMN = mean(Erg_Vector(1:sil)); %Mean silence energy (noise energy)
IMX = max(Erg_Vector); %Maximum energy for entire utterance
I1 = 0.03 * (IMX-IMN) + IMN; %I1 & I2 are Initial thresholds
I2 = 4 * IMN;
ITL = min(I1,I2); %Lower energy threshold
ITU = 5 * ITL; %Upper energy threshold

IZC = mean(ZCR_Vector(1:sil)); %mean zero crossing rate for silence region
stdev = std(ZCR_Vector(1:sil)); %standard deviation of crossing rate for 
                                %silence region
IZCT = min(ZTh,IZC+2*stdev); %Zero crossing rate threshold
indexi = zeros(1,lengthX); %Four single-row vectors are created
indexj = indexi; %in these lines to facilitate computation below
indexk = indexi;
indexl = indexi;

%Search forward for frame with energy greater than ITU
for i = 1:length(Erg_Vector)
    if (Erg_Vector(i) > ITU)
        counter1 = counter1 + 1;
        indexi(counter1) = i;
    end
end
ITUs = indexi(1);

if (ITUs==0)
    subplot(2,1,1), plot(x);
    subplot(2,1,2);
    disp('There was no voice');
    return;
end

%Search further forward for frame with energy greater than ITL
for j = ITUs:-1:1
    if (Erg_Vector(j) < ITL)
        counter2 = counter2 + 1;
        indexj(counter2) = j;
    end
end
start = indexj(1)+1;

Erg_Vectorf = fliplr(Erg_Vector);%Flips round the energy vector 
%Search forward for frame with energy greater than ITU
%This is equivalent to searching backward from last sample for energy > ITU
for k = 1:length(Erg_Vectorf)
    if (Erg_Vectorf(k) > ITU)
        counter3 = counter3 + 1;
        indexk(counter3) = k;
    end
end
ITUf = indexk(1);

%Search further forward for frame with energy greater than ITL
for l = ITUf:-1:1
    if (Erg_Vectorf(l) < ITL)
        counter4 = counter4 + 1;
        indexl(counter4) = l;
    end
end

finish = length(Erg_Vector)-indexl(1)+1;%Tentative finish index
    
%Search back from start index for crossing rates higher than IZCT
   
BackSearch = min(start,25);
for m = start:-1:start-BackSearch+1
    rate = ZCR_Vector(m);
    if rate > IZCT
        ZCRCountb = ZCRCountb + 1;
        realstart = m;
    end
end
if ZCRCountb > 3
    start = realstart; %If IZCT is exceeded in more than 3 frames
                                %set start to last index where IZCT is exceeded
end

%Search forward from finish index for crossing rates higher than IZCT
FwdSearch = min(length(Erg_Vector)-finish,25);
for n = finish+1:finish+FwdSearch
    rate = ZCR_Vector(n);
    if rate > IZCT
        ZCRCountf = ZCRCountf + 1;
        realfinish = n;
    end
end
if ZCRCountf > 3
    finish = realfinish; %If IZCT is exceeded in more than 3 frames
                                %set finish to last index where IZCT is 
                                %exceeded
end

x_start = FrmIndex(start); %actual sample index for frame 'start'
x_finish = FrmIndex(finish-1); %actual sample index for frame 'finish'
trimmedX = x(x_start:x_finish); %Trim speech sample by start and finish 
                                %indices
subplot(2,1,1), plot(x)
subplot(2,1,2), plot(trimmedX)
a=length(x)-length(trimmedX);
soundsc(trimmedX,Fs) %Playing speech alone
voicelen=length(trimmedX)/Fs; %Length of voice in seconds
%Displaying the final output
disp(sprintf('Speech is there for %f s',voicelen));
