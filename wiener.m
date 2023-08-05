%% Wiener Filtering - 05-02-2023
% Robin Abrahamse - robinab@kth.se
% Shashang Murali - shashang@kth.se

%% Setup
clear;close;clc;
% import audio
[mic1,fs] = audioread('EQ2401project2data2023_mic1.wav');
[mic2,fs] = audioread('EQ2401project2data2023_mic2.wav');
n = size(mic1,1);
% select sections with noise only
noise = mic2;
n_noise = size(noise,1);

%% FIR Wiener filter
n_fir = 30;
Ree = xcorr(noise,n_fir,'unbiased');
Ree = Ree(n_fir+1:2*n_fir+1);
Ryy = xcorr(mic1,n_fir,'unbiased');
Ryy = Ryy(n_fir+1:2*n_fir+1);
Syy = toeplitz(Ryy);
Rxx = Ryy - Ree;
Syx = Rxx;
theta = Syy\Syx;
out_fir = filter(theta,1,mic1);

%% NC Wiener filter
order = 30;
[A_noise,v_sigma] = arcov(noise,order);
[A,w_sigma] = arcov(mic1,order);

% H = phi_xx/phi_yy = (phi_yy-phi_ee)/phi_yy
[phi_yy_num,phi_yy_den] = filtspec(1,A,w_sigma);                        % phi_yy = B/A
[phi_ee_num,phi_ee_den] = filtspec(1,A_noise,v_sigma);                  % phi_ee = D/C
H_num = conv(phi_yy_num,phi_ee_den) - conv(phi_ee_num,phi_yy_den);      % H_num = BC-DA
H_den = conv(phi_yy_num,phi_ee_den);                                    % H_den = BC
out_nc = ncfilt(H_num, H_den, mic1);

%% C Wiener filter
% order = 30;
% [A_noise,v_sigma] = arcov(noise,order);
% [A,w_sigma] = arcov(src,order);

% note: phi_xx = phi_yy-phi_ee
[phi_yy_num,phi_yy_den] = filtspec(1,A,w_sigma);                        % phi_yy = B/A
[phi_ee_num,phi_ee_den] = filtspec(1,A_noise,v_sigma);                  % phi_ee = D/C
phi_xx_num = conv(phi_yy_num,phi_ee_den) - conv(phi_ee_num,phi_yy_den); % H_num = BC-DA
phi_xx_den = conv(phi_yy_den,phi_ee_den);                               % H_den = AC
[out_c,num,den]=cw(mic1,phi_xx_num,phi_xx_den,phi_yy_num,phi_yy_den,0);


%% Plot frequency responses
n_pts = 1000;
[h,w] = freqz(theta,1,n_pts,'whole',fs);
% figure
plot(w(1:n_pts/2+1),20*log10(abs(h(1:n_pts/2+1))), DisplayName='FIR')

hold on
[h,w] = freqz(H_num,H_den,n_pts,'whole',fs);
plot(w(1:n_pts/2+1),20*log10(abs(h(1:n_pts/2+1))), DisplayName='NC')

[h,w] = freqz(num,den,n_pts,'whole',fs);
plot(w(1:n_pts/2+1),20*log10(abs(h(1:n_pts/2+1))), DisplayName='C')
xlabel('Frequency (Hz)')
ylabel('Magnitude (dB)')
title('Frequency responses of filters')
legend()
hold off

%% Outputs
% soundsc(src(1:n/4),fs)
% soundsc(src,fs)

% soundsc(out_fir(1:n/4),fs)
% soundsc(out_fir,fs)

% soundsc(out_nc(1:n/4),fs)
% soundsc(out_nc,fs)

% soundsc(out_c(1:n/4),fs)
% soundsc(out_c,fs)
% audiowrite('out_fir.wav',out_fir,fs)
% audiowrite('out_nc.wav',out_nc,fs)
% audiowrite('out_c.wav',out_c,fs)