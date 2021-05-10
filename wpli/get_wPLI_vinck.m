function [wpli,wpli_biased] = get_wPLI_vinck(X, Y, f, Fs)
% Rewrite of get_wPLI.mat (feels like there's something wrong?)
% Based on Vinck et al. 2011
% Computes the debiased WPLI-square estimator (Eq 31)
% 
% Here X and Y are N_trials x N_samples matrices from 2 electrodes/sources


% init
if size(X) ~= size(Y)
    error('Inputs are not the same size, halting.')
end    
    
Ns = width(X);
Ntrials = height(X);
hz = linspace(0, Fs, Ns);

FX = fft(X, [], 2)*(2/Ns);  % This scaling is incorrect for DC and Nyquist but we don't care
FY = fft(Y, [], 2)*(2/Ns);

xspectrum = FX .* conj(FY);
xspectrum = reshape(xspectrum, Ntrials, 1, Ns);

f_ix = dsearchn(hz', f);

wpli = ft_connectivity_wpli(xspectrum, 'dojack', 0, 'debias', 1);
wpli_biased = ft_connectivity_wpli(xspectrum, 'dojack', 0, 'debias', 0);

wpli = abs(wpli(f_ix));
wpli_biased = abs(wpli_biased(f_ix));
