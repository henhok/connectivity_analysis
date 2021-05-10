%% Some parameters
len_s = 1;
srate_hz = 500;
N = len_s * srate_hz;
xvec = linspace(0,len_s,N);

f = 4;
Ntrials = 60; 
spreadvec = [0.02 0.1 0.5]*pi;
phasediffvec = linspace(0,pi/2,10);

%% Demo plot
figure;
tiledlayout(length(spreadvec),length(phasediffvec), 'TileSpacing', 'compact', 'Padding', 'none');
for j = 1:length(spreadvec)
    for i = 1:length(phasediffvec)
        ph1 = 0 + spreadvec(j)*randn(1,Ntrials);
        ph2 = phasediffvec(i) + spreadvec(j)*randn(1,Ntrials);
        nexttile
        polarscatter(ph1-ph2, ones(Ntrials,1), '.');
        title({
            ['\theta = ' num2str(phasediffvec(i)/pi,2) ' x \pi']
            ['sd = ' num2str(spreadvec(j)/pi,2) 'x \pi']
            });
        ax = gca;
%         ax.ThetaGrid = 'off';
        ax.RGrid = 'off';
        ax.RTickLabel = []; 
        ax.ThetaTickLabel = [];
    end
end

%% Main computation 
M = zeros(length(phasediffvec), length(spreadvec), 4);
Nexp = 1;  % Number of 60 trial experiments

for i = 1:length(phasediffvec)
    for j = 1:length(spreadvec)
        
        phasediff = phasediffvec(i);
        spread = spreadvec(j);
        
        for n=1:Nexp
            phasevec1 = 0 + spread*randn(1,Ntrials);
            phasevec2 = phasediff + spread*randn(1,Ntrials);

            y1 = zeros(Ntrials, N);
            for k=1:Ntrials
                y1(k,:) = sin(2*pi*f*xvec + phasevec1(k));
            end

            y2 = zeros(Ntrials, N);
            for k=1:Ntrials
                y2(k,:) = sin(2*pi*f*xvec + phasevec2(k));
            end
            
            % HH method
            hh_dwpli(n) = get_wPLI_henri(y1,y2,f,srate_hz);
            
            % Vinck method
            [wpli,wpli_biased] = get_wPLI_vinck(y1,y2,f,srate_hz);
            vinck_dwpli(n) = wpli;
            vinck_wpli_biased(n) = wpli_biased;
            
            % AT method
            y1_supertrial = reshape(y1', 1, Ntrials*N);
            y2_supertrial = reshape(y2', 1, Ntrials*N);
            baba_dwpli(n) = abs(get_wPLI_baba(y1_supertrial, y2_supertrial, 1));
        end

        M(i,j,1) = mean(hh_dwpli);
        M(i,j,2) = mean(baba_dwpli);
        M(i,j,3) = mean(vinck_dwpli);
        M(i,j,4) = mean(vinck_wpli_biased);
    end
end

figure;
ax1 = plot(phasediffvec/pi, M(:,:,1), 'k'); hold on
ax2 = plot(phasediffvec/pi, M(:,:,2), 'b');
ax3 = plot(phasediffvec/pi, M(:,:,3), 'ro'); hold off  % Change to M(:,:,4) for biased wPLI
set(ax1, {'LineWidth'}, {3;2;1});
set(ax2, {'LineWidth'}, {3;2;1});
set(ax3, {'LineWidth'}, {3;2;1});
xlabel('Phase difference (x \pi)');
ylabel('dwPLI-squared');
