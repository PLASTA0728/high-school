%Subsample audio signal in the time domain, try to recover it by exploiting sparsity in the cosine frequency domain.

clear;
close all;

%read audio
[y,Fs] = audioread('1.m4a');

%Clip off dead air
%x=y(21625:end);
x=y(30000:end);

%Select a 1 second segment
x = x(1:16000);

%What does it sound like?
sound(x,Fs)

%Length, time T, array of sample times
n = length(x);
T = n/Fs; %Sample period
ts = ([0:n-1]+0.5)/Fs; %Sample times

%Display the selected signal
figure(5);
plot(ts,x);
title('One Second of Audio');

%Generate about m random integers in the range 1 to n
m = 2000; %Number of time domain samples we'll take, roughly.
thing = randi(n,1,m); %Choose m random times to sample
thing = unique(thing);
m = length(thing);

%Select m subsamples from the data from x (not truly random, 
%constrained by original sampling interval)
d = x(thing);
ts = ((thing-0.5)/16000); %These are the times we sampled x.

%The sensing matrix is m X N, where N is the number of ck we will
%try to estimate in ck*cos(k*pi*t). So highest frequency is (N-1)/2.
N = 8000; %Up to 4000 Hz, should be high enough.
kf = [0:N-1]; %k values in cos(k*pi*t)
A = zeros(m,N);
for j=1:m
    A(j,:) = cos(ts(j)*pi*kf);
end

%Now try to recover c (DCT coefficients) from A*c = d, take c to be
%p-sparse where p = m/4.
p = floor(m/4);

%Let's call OMP
opts.omp = 1;
opts.maxits = p;
[alpha, rhist] = omp(A, d, opts);

%How does the estimate of sparse c look, compare to the DCT
%of x on the same frequency range?
figure(1);
plot(kf/2,alpha,'r');
title('OMP Recovered Spectrum');
X = kdct(x);

figure(2)
plot(kf/2,X(1:N),'b');
title('DCT Computed Spectrum');

%Now form the time domain reconstruction, by using the ck in
%sum(ck*cos((k-1)*pi*t),k=0..N) evaluated at 0<=t<=1.
ts2 = ((1:16000)-0.5)/16000;
xrecon = zeros(1,16000);
for k=1:N
    xrecon = xrecon + alpha(k)*cos(ts2*pi*kf(k));
end;
xrecon = xrecon';

%Simple error estimate.
norm(x-xrecon,2)/norm(x,2)

figure(3)
plot(1:n,x,'r',1:n,xrecon,'b')
title('OMP Recovered and Actual Gong Signal');

%Now play the signal as reconstructed from a p-sparse estimate.
sound(xrecon,Fs)