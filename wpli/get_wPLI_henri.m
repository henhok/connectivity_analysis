function wpli = get_wPLI_henri(X, Y, f, Fs)
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

f_ix = dsearchn(hz', f);
FX = squeeze(FX(:,f_ix));
FY = squeeze(FY(:,f_ix));

xspectrum = FX .* conj(FY);

% Compute the pairwise product of all imaginary components
xmatrix = imag(xspectrum) * imag(xspectrum)';
xmatrix = xmatrix - diag(diag(xmatrix));

% Sum them up and divide by the sum of the magnitudes of the same imag
% components (Equation 31, Vinck et al. 2011)
wpli = abs(sum(xmatrix, 'all'))/sum(abs(xmatrix), 'all');

% UNCOMMENT TO GET PLOTS
% figure;
% subplot(121)
% polarplot(xspectrum,'o');
% sgtitle(['dWPLI^2 = ' num2str(wpli)]);
%  
% subplot(122)
% timevec = linspace(0, Ntrials*Ns/Fs, Ns*Ntrials);
% plot(timevec, squeeze(reshape(X', 1, Ntrials*Ns))); hold on
% plot(timevec, squeeze(reshape(Y', 1, Ntrials*Ns))); hold off