%% Characteristics & Transmitter
c =[1 1];
N = length(c);
L=512; % Subcarriers
P0 = 1000; % Transmit power
Lx = 100; % simulate Lx blocks
G = abs(fft(c,L)); % L point FFT of channel
SNR = (G.^2); % SNR for a channel, with noise variance 1.
D = waterfilling(SNR,P0);
Diag = diag(D);
T = dftmtx(L)/sqrt(L); % L-point DFT matrix, scaled to be unitary
x = randn(1,L*Lx); XX = zeros(L,Lx); % signal
XX(:) = x; % demultiplex
XX2 = Diag*XX;
X = T'*XX2; % IDFT
S = [X(end-N+2:end,:); X]; % cyclic prefix
s = S(:); % multiplex
%% Noisy and Dispersive Channel
v = 0.0 * (randn(Lx*(L+N-1),1) + 1i*randn(Lx*(L+N-1),1))/sqrt(2);
r = filter(c,1,s) + v; % transmission over channel
%% OFDM Receiver and Equaliser
R= zeros(L+N-1,Lx);
R(:) = r; % demultiplex
R2 = T*R(N:end,:); % prefix removal and DFT
S_hat = pinv(Diag)*pinv(diag(fft(c,L)))*R2; % ZF equalisation
%% MSE
for i= 1:Lx,
MSE(:,i) = sum((S_hat(:,i) - XX(:,i)).^2)/L;
end
