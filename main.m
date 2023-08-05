%% Adaptive Filtering - 02-03-2023
% Robin Abrahamse - robinab@kth.se
% Shashang Murali - shashang@kth.se

%% Setup
clear;close all;clc;
[mic1,fs] = audioread('EQ2401project2data2023_mic1.wav');
[mic2,fs] = audioread('EQ2401project2data2023_mic2.wav');
n = size(mic1,1);

%% Parameters
order = 12;
mu_lms = 0.08;
mu_nlms = 0.05;
l = 0.999;

%% LMS
lms = dsp.LMSFilter(order,'Method','LMS','StepSize',mu_lms);
[~,out_lms] = lms(mic2,mic1);

% Own implementation
% x_lms = zeros(order,1);
% theta_lms = zeros(order,1);
% out_lms2 = zeros(n,1);
% 
% for i = 1:n
%     x_lms = [mic2(i); x_lms(1:end-1)];
%     out_lms2(i) = mic1(i) - theta_lms'*x_lms;
%     theta_lms = theta_lms + mu_lms*out_lms2(i)*x_lms;
% end

%% LMS parameter estimation
% k = 10;
% X = zeros(k,k);
% Y = zeros(k,k);
% for i = 1:k
%     X(i,:) = i*0.1/k;
%     Y(:,i) = 6+i*30/k;
% end
% X = X(:);
% Y = Y(:);
% Z = zeros(k*k,1);
% for i = 1:k*k
%     lms = dsp.LMSFilter(Y(i),'Method','LMS','StepSize',X(i));
%     [~,out_lms] = lms(mic2,mic1);
%     Z(i) = sum(msesim(lms,mic2,mic1));
% end
% X = X(1:k);
% Y = Y(1:k:k*k);
% Z = reshape(Z,k,k)';
% surf(X,Y,Z)
% title('Estimated LMS MSE residue')
% xlabel('\mu')
% ylabel('N')

%% NLMS
nlms = dsp.LMSFilter(order,'Method','Normalized LMS','StepSize',mu_nlms);
[~,out_nlms] = nlms(mic2,mic1);

% Own implementation
% x_nlms = zeros(order,1);
% theta_nlms = zeros(order,1);
% out_nlms2 = zeros(n,1);
% 
% for i = 1:n
%     x_nlms = [mic2(i); x_nlms(1:end-1)];
%     out_nlms2(i) = mic1(i) - theta_nlms'*x_nlms;
%     power = x_nlms'*x_nlms;
%     theta_nlms = theta_nlms + mu_nlms/power*out_nlms2(i)*x_nlms;
% end

%% NLMS parameter estimation
% k = 10;
% X = zeros(k,k);
% Y = zeros(k,k);
% for i = 1:k
%     X(i,:) = i*0.1/k;
%     Y(:,i) = 6+i*30/k;
% end
% X = X(:);
% Y = Y(:);
% Z = zeros(k*k,1);
% for i = 1:k*k
%     nlms = dsp.LMSFilter(Y(i),'Method','Normalized LMS','StepSize',X(i));
%     [~,out_nlms] = nlms(mic2,mic1);
%     Z(i) = sum(msesim(nlms,mic2,mic1));
% end
% X = X(1:k);
% Y = Y(1:k:k*k);
% Z = reshape(Z,k,k)';
% surf(X,Y,Z)
% title('Estimated NLMS MSE residue')
% xlabel('\mu')
% ylabel('N')

%% RLS
rls = dsp.RLSFilter(order,'Method','Conventional RLS','ForgettingFactor',l);
[~,out_rls] = rls(mic2,mic1);
rlssw = dsp.RLSFilter(order,'Method','Sliding-window RLS','ForgettingFactor',l,'SlidingWindowBlockLength',1400);
[~,out_rlssw] = rlssw(mic2,mic1);

% Own implementation
% P0 = 1000;
% theta_rls = zeros(order,1);
% P = P0*eye(order);
% x_rls = zeros(order,1);  
% out_rls2 = zeros(n,1);
% for i = 1:n
%     x_rls = [mic2(i); x_rls(1:end-1)];
%     out_rls2(i) = mic1(i) - theta_rls'*x_rls;
%     g = P*x_rls/(l + x_rls'*P*x_rls);
%     P = (P - g*x_rls'*P)/l;
%     theta_rls = theta_rls + P*x_rls*out_rls2(i);
% end

%% RLS parameter estimation
% order & lambda
% k = 10;
% X = zeros(k,k);
% Y = zeros(k,k);
% for i = 1:k
%     X(i,:) = 0.99+i*0.01/k;
%     Y(:,i) = 6+i*30/k;
% end
% X = X(:);
% Y = Y(:);
% Z = zeros(k*k,1);
% for i = 1:k*k
%     rls = dsp.RLSFilter(Y(i),'Method','Conventional RLS','ForgettingFactor',X(i));
%     [~,out_rls] = rls(mic2,mic1);
%     Z(i) = sum(msesim(rls,mic2,mic1));
% end
% X = X(1:k);
% Y = Y(1:k:k*k);
% Z = reshape(Z,k,k)';
% surf(X,Y,Z)
% title('Estimated RLS MSE residue')
% xlabel('\lambda')
% ylabel('N')

% Sliding window
% k = 21;
% X = (800:1000/(k-1):1800)';
% Z = zeros(k,1);
% for i = 1:k
%     rls = dsp.RLSFilter(order,'Method','Sliding-window RLS','ForgettingFactor',l,'SlidingWindowBlockLength',X(i));
%     [~,out_rls] = rls(mic2,mic1);
%     Z(i) = sum(msesim(rls,mic2,mic1));
% end
% plot(X,Y)


%% Frequency responses & Spectra
% close all
[r_mic1,c] = xcov(mic1);
r_mic2 = xcov(mic2);
r_lms = xcov(out_lms);
r_nlms = xcov(out_nlms);
r_rls = xcov(out_rls);
r_rlssw = xcov(out_rlssw);
fft_mic1 = fft(r_mic1(n:2*n-1));
fft_mic2 = fft(r_mic2(n:2*n-1));
fft_lms = fft(r_lms(n:2*n-1));
fft_nlms = fft(r_nlms(n:2*n-1));
fft_rls = fft(r_rls(n:2*n-1));
fft_rlssw = fft(r_rlssw(n:2*n-1));

figure('Name','Spectra')
hold on
plot(0:fs/n:fs/2-fs/n, 20*log10(abs(fft_mic2(1:n/2))), 'DisplayName', 'noise sequence');
plot(0:fs/n:fs/2-fs/n, 20*log10(abs(fft_mic1(1:n/2))), 'DisplayName', 'input sequence');
plot(0:fs/n:fs/2-fs/n, 20*log10(abs(fft_lms(1:n/2))), 'DisplayName', 'LMS filter');
plot(0:fs/n:fs/2-fs/n, 20*log10(abs(fft_nlms(1:n/2))), 'DisplayName', 'NLMS filter');
plot(0:fs/n:fs/2-fs/n, 20*log10(abs(fft_rls(1:n/2))), 'DisplayName', 'RLS filter');
legend('Location','northwest')
xlabel('frequency (Hz)')
ylabel('signal magnitude (dB)')
hold off

figure('Name','Filter responses')
hold on
plot(0:fs/n:fs/2-fs/n, 20*log10(abs(fft_mic2(1:n/2))), 'DisplayName', 'noise sequence');
plot(0:fs/n:fs/2-fs/n, 20*log10(abs(fft_lms(1:n/2))./abs(fft_mic1(1:n/2))), 'DisplayName', 'LMS filter');
plot(0:fs/n:fs/2-fs/n, 20*log10(abs(fft_nlms(1:n/2))./abs(fft_mic1(1:n/2))), 'DisplayName', 'NLMS filter');
plot(0:fs/n:fs/2-fs/n, 20*log10(abs(fft_rls(1:n/2))./abs(fft_mic1(1:n/2))), 'DisplayName', 'RLS filter');
legend('Location','northwest')
xlabel('frequency (Hz)')
ylabel('signal magnitude (dB)')
hold off

figure('Name','Sliding window vs conventional RLS')
hold on
plot(0:fs/n:fs/2-fs/n, 20*log10(abs(fft_rls(1:n/2))./abs(fft_mic1(1:n/2))), 'DisplayName', 'RLS filter');
plot(0:fs/n:fs/2-fs/n, 20*log10(abs(fft_rlssw(1:n/2))./abs(fft_mic1(1:n/2))), 'DisplayName', 'Sliding window RLS filter');
legend('Location','northwest')
xlabel('frequency (Hz)')
ylabel('signal magnitude (dB)')
hold off

%% Outputs
% soundsc(out_lms,fs)
% soundsc(out_nlms,fs)
% soundsc(out_rls,fs)
% soundsc(out_rlssw,fs)

% audiowrite('out_lms.wav',out_lms,fs)
% audiowrite('out_nlms.wav',out_nlms,fs)
% audiowrite('out_rls.wav',out_rls,fs)
% audiowrite('out_rlssw.wav',out_rlssw,fs)