function D = waterfilling(SNR,P0);
% D = waterfilling(SNR,P0);
%
% Perform waterfilling to maximise the capacity of a system with
% subchannels with given SNRs under a constraint P0 for the transmit
% power.
%
% Input parameters:
%   SNR     column vector of subchannel SNRs  
%   P0      total transmit power budget
%
% Output parameter:
%   D       column vector containing optimised subchannel gains
%           applicable at the transmitter's precoder.
%  
% S. Weiss, UoS, 2/5/2007
  
L = length(SNR);      % number of subchannels  
epsilon = 10^(-14);   % avoid division by zero
  
% order the inverse subchannel SNRs from smallest to largest 
[SNRinv,Index] = sort(1./(SNR+epsilon)); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% part 1:
% determine how many subchannels K<=L can be supported 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
K=L;   % first consider taking all (K=L) subchannels
mu = SNRinv(K);
P = sum(mu-SNRinv);
while (P>P0),
  K=K-1;
  mu = SNRinv(K);
  P = sum((mu-SNRinv).*((mu-SNRinv)>0));
end;    
% K now contains the number of subchannels that can be "covered by
% water" (i.e. utilised for transmission)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% part 2:
% P is the power required to service the K strongest subchannels,
% i.e. water just reaches up to the Kth subchannel.
% The remaining (P0-P) is uniformly poored across these K subchannels.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mu = mu + (P0-P)/K;
SNRinv = SNRinv(:);   % enforce column vector 
D2 = max(zeros(L,1), mu-SNRinv);

% turn power into gain and re-order subchannels
D(Index) = sqrt(D2);
