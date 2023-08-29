%*************************************************************************
%   FUNCTION:      get_eeg_features.m
%   =========      ==================
%
%   DESCRIPTION:   Get features using linear fits to eeg bands 
%   ============   and ratio of psd in bands.
%                  Add to above cross spectral density for 
%                  opposite sides of the brain
%                  
%                  Entry 05
%                  
%
%   BY:            Allan Moser and Jackie Le
%   ===
%
%   DATE CREATED:  02-03-2023
%   =============
%
%   LAST CHANGED:  07-27-2023
%   =============
%
%**************************************************************************

function eeg_features=get_eeg_features(dat,sampling_frequency,patient_id,hour_first_last)
    fs = sampling_frequency; % Sampling frequency

%     % Wavelet denoise
%     for i = 1:18
%         temp = wdenoise(dat(:,i));
%         dat(:,i) = temp;
%     end

    % Calculate mean, std
    mn_all_sig  = mean(dat);
    std_all_sig = std(dat);

    % Label individual signals to help interpretation 
    Fp1F7 = dat(:,1);  F7T3 = dat(:,2);  T3T5 = dat(:,3);  T5O1 = dat(:,4);
    Fp2F8  = dat(:,5);  F8T4 = dat(:,6);  T4T6 = dat(:,7);  T6O2 = dat(:,8);
    Fp1F3 = dat(:,9);  F3C3 = dat(:,10); C3P3 = dat(:,11); P3O1 = dat(:,12);
    Fp2F4 = dat(:,13); F4C4 = dat(:,14); C4P4 = dat(:,15); P4O2 = dat(:,16);
    FzCz  = dat(:,17); CzPz = dat(:,18);
   
    % plot signal data
    % plt_eeg_chnnls(hour_first_last, ...
                   % Fp1F7,F7T3,T3T5,F7T3, ...
                   % F8T4,T4T6,T6O2,Fp1F3,F3C3,C3P3,P3O1, ...
                   % Fp2F4,F4C4,C4P4,P4O2, ...
                   % FzCz,CzPz);
    % plt_wve(hour_first_last,Fp1F7,Fp1F3,F7T3,F3C3);
    % plt_wve(hour_first_last,T3T5,C3P3,T5O1,P3O1);
    % plt_wve(hour_first_last,FzCz,CzPz);
    % plt_wve(hour_first_last,Fp2F4,Fp2F8,F4C4,F8T4);
    % plt_wve(hour_first_last,C4P4,T4T6,P4O2,T6O2);

    % geographic breakdown
    % Fp1F7,Fp1F3        Fp2F4,Fp2F8
    %   F7T3,F3C3  FzCz  F4C4,F8T4
    %   T3T5,C3P3  CzPz  C4P4,T4T6
    %   T5O1,P3O1        P4O2,T6O2


    % Get power for signal, delta, theta, alpha, beta, theta-alpha-beta for all EEG lead
    bp_total = bandpower(dat,fs,[0,26]);
    bp_delta = bandpower(dat,fs,[1,3]);
    bp_theta = bandpower(dat,fs,[3,6]);
    bp_alpha = bandpower(dat,fs,[6,10]);
    bp_beta  = bandpower(dat,fs,[10,26]);
    bp_tab   = bandpower(dat,fs,[3,26]);

    % For opposite side of brain calculations the order is:
    % 1: Fp1 - F7 / Fp2 - F8
    % 2: F7  - T3 / F8  - T4
    % 3: T3  - T5 / T4  - T6
    % 4: T5  - O1 / T6  - O2
    % 5: Fp1 - F3 / Fp2 - F4
    % 6: F3  - C3 / F4  - C4
    % 7: C3  - P3 / C4  - P4
    % 8: P3  - O1 / P4  - O2

    % Cross ratio of bandpower for opposite sides of brain
%     for i = 1:4
%         cross_rat_tot(i)   = 2*(bp_total(i) - bp_total(i+4))/(bp_total(i) + bp_total(i+4));
%         cross_rat_delta(i) = 2*(bp_delta(i) - bp_delta(i+4))/(bp_delta(i) + bp_delta(i+4));
%         cross_rat_theta(i) = 2*(bp_theta(i) - bp_theta(i+4))/(bp_theta(i) + bp_theta(i+4));
%         cross_rat_alpha(i) = 2*(bp_alpha(i) - bp_alpha(i+4))/(bp_alpha(i) + bp_alpha(i+4));
%         cross_rat_beta(i) =  2*(bp_beta(i)  - bp_beta(i+4)) /(bp_beta(i)  + bp_beta(i+4));
%         cross_rat_tab(i)  =  2*(bp_tab(i)   - bp_tab(i+4))  /(bp_tab(i)   + bp_tab(i+4));
%         k = i+4; m = i+8;
%         cross_rat_tot(k)   = 2*(bp_total(m) - bp_total(m+4))/(bp_total(m) + bp_total(m+4));
%         cross_rat_delta(k) = 2*(bp_delta(m) - bp_delta(m+4))/(bp_delta(m) + bp_delta(m+4));
%         cross_rat_theta(k) = 2*(bp_theta(m) - bp_theta(m+4))/(bp_theta(m) + bp_theta(m+4));
%         cross_rat_alpha(k) = 2*(bp_alpha(m) - bp_alpha(m+4))/(bp_alpha(m) + bp_alpha(m+4));
%         cross_rat_beta(k) =  2*(bp_beta(m)  - bp_beta(m+4)) /(bp_beta(m)  + bp_beta(m+4));
%         cross_rat_tab(k)  =  2*(bp_tab(m)   - bp_tab(m+4))  /(bp_tab(m)   + bp_tab(m+4));
%     end
   
    % Find the power spectrum for all EEG leads
    [psd,pf] = pspectrum(dat,fs);
    % Find index into psd for individual bands
    pindx_delta = find(pf >  1 & pf <  3); % delta band
    pindx_theta = find(pf >  3 & pf <  6); % theta band
    pindx_alpha = find(pf >  6 & pf < 10); % alpha band
    pindx_beta  = find(pf > 10 & pf < 26); % beta band
    pindx_tab   = find(pf >  3 & pf < 26); % theta through beta bands

    % Cross-spectral coherence for opposite sides of brain
    if length(dat) > 1024
        for i = 1:4 
            [cpsd(:,i),cf] = mscohere(dat(:,i),dat(:,i+4),1024,[],[],fs);
            k=i+4; m = i+8;
            [cpsd(:,k),cf] = mscohere(dat(:,m),dat(:,m+4),1024,[],[],fs);
        end
    else
        for i = 1:8
            cpsd(:,i) = 0;
        end
    end
    cindx_delta = find(cf >  1 & cf <  3); % delta band
    cindx_theta = find(cf >  3 & cf <  6); % theta band
    cindx_alpha = find(cf >  6 & cf < 10); % alpha band
    cindx_beta  = find(cf > 10 & cf < 26); % beta band
    cindx_tab   = find(cf >  3 & cf < 26); % theta through beta bands  

    % Process delta band
    freq = pf(pindx_delta);
    dbps  = 10*log10(psd(pindx_delta,:));
    std_delta = std(dbps);
    for i = 1:18
        [coefs,parm] = polyfit(freq,dbps(:,i),1);
        rsq_delta(i) = 1 - parm.normr^2 / norm(dbps(:,i) - mean(dbps(:,i)))^2;
        slp_delta(i) = coefs(1);
        rat_delta_total(i)  = bp_delta(i)/bp_total(i);
        rat_delta_theta(i)  = bp_delta(i)/bp_theta(i);
        rat_delta_alpha(i)  = bp_delta(i)/bp_alpha(i);
        rat_delta_beta(i)   = bp_delta(i)/bp_beta(i);
    end
    cohere_delta = mean(cpsd(cindx_delta,:));
    
    % Process theta band
    freq = pf(pindx_theta);
    dbps  = 10*log10(psd(pindx_theta,:));
    std_theta = std(dbps);
    for i = 1:18
        [coefs,parm] = polyfit(freq,dbps(:,i),1);
        rsq_theta(i) = 1 - parm.normr^2 / norm(dbps(:,i) - mean(dbps(:,i)))^2;
        slp_theta(i) = coefs(1);
        rat_theta_total(i)  = bp_theta(i)/bp_total(i);
        rat_theta_alpha(i)  = bp_theta(i)/bp_alpha(i);
        rat_theta_beta(i)   = bp_theta(i)/bp_beta(i);
    end
    cohere_theta = mean(cpsd(cindx_theta,:));
    
    % Process alpha band
    freq = pf(pindx_alpha);
    dbps  = 10*log10(psd(pindx_alpha,:));
    std_alpha = std(dbps);
    for i = 1:18
        [coefs,parm] = polyfit(freq,dbps(:,i),1);
        rsq_alpha(i) = 1 - parm.normr^2 / norm(dbps(:,i) - mean(dbps(:,i)))^2;
        slp_alpha(i) = coefs(1);
        rat_alpha_total(i)  = bp_alpha(i)/bp_total(i);
        rat_alpha_beta(i)   = bp_alpha(i)/bp_beta(i);
    end
    cohere_alpha = mean(cpsd(cindx_alpha,:));

    % Process beta band
    freq = pf(pindx_beta);
    dbps  = 10*log10(psd(pindx_beta,:));
    std_beta = std(dbps);
    for i = 1:18
        [coefs,parm] = polyfit(freq,dbps(:,i),1);
        rsq_beta(i) = 1 - parm.normr^2 / norm(dbps(:,i) - mean(dbps(:,i)))^2;
        slp_beta(i) = coefs(1);
        rat_beta_total(i)  = bp_beta(i)/bp_total(i);
    end
    cohere_beta = mean(cpsd(cindx_beta,:));

    % Process theta_alpha_beta bands
    freq = pf(pindx_tab);
    dbps  = 10*log10(psd(pindx_tab,:));
    std_tab = std(dbps);
    for i = 1:18
        [coefs,parm] = polyfit(freq,dbps(:,i),1);
        rsq_tab(i) = 1 - parm.normr^2 / norm(dbps(:,i) - mean(dbps(:,i)))^2;
        slp_tab(i) = coefs(1);
        rat_tab_total(i)  = bp_tab(i)/bp_total(i);
    end
    cohere_tab = mean(cpsd(cindx_tab,:));

%     % Get average fractal dimension
%     kmax = 100;
%     for i = 1:18
%         hgfrac(i) = higuchi_frac_dim(dat(:,i),kmax);
%     end
%     mn_hgfrac = mean(hgfrac);
%     sd_hgfrac = std(hgfrac);
% 
%      % Wavelet decomposition
%     for i = 1:18
%         [coef,lev] = wavedec(dat(:,i),12,'sym8');
%         [cd1,cd2,cd3,cd4,cd5,cd6,cd7,cd8,cd9,cd10,cd11,cd12] = detcoef(coef,lev,[1 2 3 4 5 6 7 8 9 10 11 12]);
%         wstd1(i) = std(cd1); wstd2(i) = std(cd2); wstd3(i) = std(cd3); wstd4(i) = std(cd4); 
%         wstd5(i) = std(cd5); wstd6(i) = std(cd6); wstd7(i) = std(cd7); wstd8(i) = std(cd8);
%         wstd9(i) = std(cd9); wstd10(i) = std(cd10); wstd11(i) = std(cd11); wstd12(i) = std(cd12); 
%     end

    % channel indices
    vchan = [3 7 2 6 17 18]; % use all the channels
    vchan = [1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18];
    xchan = [3 2];
    xchan = [1 2 3 4 5 6 7 8];
    eeg_features = [std_all_sig, ...
                    bp_total(vchan),bp_delta(vchan),bp_theta(vchan),bp_alpha(vchan),bp_beta(vchan), ... 
                    std_delta(vchan),slp_delta(vchan),rsq_delta(vchan),rat_delta_total(vchan),rat_delta_theta(vchan),rat_delta_alpha(vchan),rat_delta_beta(vchan), ...
                    std_theta(vchan),slp_theta(vchan),rsq_theta(vchan),rat_theta_total(vchan),rat_theta_alpha(vchan),rat_theta_beta(vchan), ...
                    std_alpha(vchan),slp_alpha(vchan),rsq_alpha(vchan),rat_alpha_total(vchan),rat_alpha_beta(vchan), ...
                    std_beta(vchan),slp_beta(vchan),rsq_beta(vchan),rat_beta_total(vchan), ...
                    cohere_delta(xchan),cohere_theta(xchan),cohere_alpha(xchan),cohere_beta(xchan)];

% additional features

% correlation coefficients (-1<cc<1)
% organize by frequency range & brain region (latter in the works)
bp_cc = [bp_total(vchan);bp_delta(vchan);bp_theta(vchan);bp_alpha(vchan);bp_beta(vchan)].'; % length 5
std_cc = [std_delta(vchan);std_theta(vchan);std_alpha(vchan);std_beta(vchan)].'; % 4
slp_cc = [slp_delta(vchan);slp_theta(vchan);slp_alpha(vchan);slp_beta(vchan)].'; % 4
rsq_cc = [rsq_delta(vchan);rsq_theta(vchan);rsq_alpha(vchan);rsq_beta(vchan)].'; % 4
rat_cc = [rat_delta_total(vchan);rat_theta_total(vchan);rat_alpha_total(vchan);rat_beta_total(vchan); ...
          rat_delta_theta(vchan);rat_delta_alpha(vchan);rat_delta_beta(vchan); ...
          rat_theta_alpha(vchan);rat_theta_beta(vchan);rat_alpha_beta(vchan)].'; % 10
co_cc = [cohere_delta(xchan);cohere_theta(xchan);cohere_alpha(xchan);cohere_beta(xchan)].'; % 4

pca_threshold = 0.3;

[bp_cc_list,bp_cc_vals,bp_cc_percent] = cc_check(bp_cc,pca_threshold);
[std_cc_list,std_cc_vals,std_cc_percent] = cc_check(std_cc,pca_threshold);
[slp_cc_list,slp_cc_vals,slp_cc_percent] = cc_check(slp_cc,pca_threshold);
[rsq_cc_list,rsq_cc_vals,rsq_cc_percent] = cc_check(rsq_cc,pca_threshold);
[rat_cc_list,rat_cc_vals,rat_cc_percent] = cc_check(rat_cc,pca_threshold);
[co_cc_list,co_cc_vals,co_cc_percent] = cc_check(co_cc,pca_threshold);

% subplot correlated nodes
% dir_grph(patient_id,bp_cc_list,std_cc_list,slp_cc_list,rsq_cc_list,rat_cc_list,co_cc_list);

% scree plot

% inspect graph formed by nodes of correlated features
% figure(10);
bp_pcmp = pca_cllpse(bp_cc,[1,2,3,4],5);
std_pcmp = pca_cllpse(std_cc,[1,2,4],3);
slp_pcmp = pca_cllpse(slp_cc,[1,2],[3,4]);
rsq_pcmp = pca_cllpse(rsq_cc,[1,2],3,4);
rat_pcmp = pca_cllpse(rat_cc,[1,2,3,4],[5,6,7,8,9,10]);
co_pcmp = pca_cllpse(co_cc,[1,2,3],4);
% hold off

% reshape to singular row of features
bp_pcmp = reshape(bp_pcmp,1,[]);
std_pcmp = reshape(std_pcmp,1,[]);
slp_pcmp = reshape(slp_pcmp,1,[]);
rsq_pcmp = reshape(rsq_pcmp,1,[]);
rat_pcmp = reshape(rat_pcmp,1,[]);
co_pcmp = reshape(co_pcmp,1,[]);
eeg_features = [std_all_sig,bp_pcmp,std_pcmp,slp_pcmp,rsq_pcmp,rat_pcmp,co_pcmp];

end
%% add'l fcns
function hfd = higuchi_frac_dim(sig,kmax)
    % Higuchi Fractal Dimension
    % Inputs:  sig - input signal
    %          kmax - maximum value of k (typically 8)
    n = length(sig);
    % Preallocate arrrays for efficiency
    Ldif = zeros(1,kmax);
    x    = zeros(1,kmax);
    y    = zeros(1,kmax);
    for k = 1:kmax
        for m = 1:k
            nfac    = (n-1)/(round((n-m)/k)*k); 
            da      = sum(abs(diff(sig(m:k:n)))); 
            Ldif(m) = da*nfac/k; 
        end
        y(k) = log(sum(Ldif)/k); 
        x(k) = log(1/k); 
    end
    parms = polyfit(x,y,1); % Fit line
    hfd = parms(1);         % slope
end
%
function [cc_list,cc_vals,cc_percent] = cc_check(features,pca_min)
    cc_matrix = corrcoef(features); % assumes tranpose to apply corrcoef()
    [nrows,mcols] = size(cc_matrix); % pre-allocate for speed
    cc_list = zeros(nrows*mcols,2); % 2 for correlated features
    % scan matrix
    for i = 1:nrows
        for j = 1:mcols
            if cc_matrix(i,j) > pca_min
                cc_list((i-1)*mcols+j,:) = [i,j]; % sub in features
            end
        end
    end
    cc_list = cc_list(~all(cc_list==0,2),:); % remove zero rows
    % sort from most to least correlated
    [nrows,~] = size(cc_list);
    cc_vals = zeros(nrows,1);
    for k = 1:nrows
            cc_vals(k) = cc_matrix(cc_list(k,1),cc_list(k,2)); % index cc
    end
    [cc_vals,idx] = sort(cc_vals,'descend');
    cc_list = cc_list(idx,:);
    % evaluate remaining cc percent including duplicates
    cc_percent = numel(cc_vals)/numel(cc_matrix);
    % remove duplicate feature combinations, from symmetric matrix
    % may exist non-unique cc
    [cc_list,uniq_idx,~] = unique(sort(cc_list,2),'rows','stable');
    cc_vals = cc_vals(uniq_idx);
end
%
% inputs: 1st is array of features, 1st+x is grouped features
function pcmp_list = pca_cllpse(varargin)
    [pca_mrow,num_feat] = size(varargin{1}); % assumes varargin{1}(2,2+x) mxn
    feat_arg = varargin{1};
    uniq_tot_chck = 0;
    for q = 2:nargin
        uniq_elts = sum(numel(unique(varargin{q})));
        uniq_tot_chck = uniq_tot_chck + uniq_elts;
    end
    assert(uniq_tot_chck == num_feat,'Not all nodes defined.');
    pcmp_list = zeros(pca_mrow,num_feat); % pre-allocate for speed
    sing_grp = zeros(1,num_feat);
    for i = 2:nargin % perform pca per group of correlated features
        arg = varargin{i};
        if numel(arg) == 1
            frst_zro = find(sing_grp == 0,1,'first');
            sing_grp(frst_zro) = arg;
            continue
        end
        feat_arg_grp = feat_arg(:,arg);
        [pcmp,~,lam,~,~] = pca_arm(feat_arg_grp); % input = mxn = output
        % tot_lam = sum(lam); % calculate variance accounted for per lambda
        % num_lam = numel(lam);
        % lam_per = zeros(num_lam,1);
        pcmp_num = 2; % non-zero for constant thresholding instead

        % for j = 1:num_lam % adaptive thresholding
        %     lam_per(j) = lam(j)/tot_lam;
        %     % select principal components accounting for minimum threshold
        %     if pcmp_num ~= 0
        %         continue
        %     end
        %     lam_var_tot = sum(lam_per(1:j));
        %     if lam_var_tot > 0.9
        %         pcmp_num = j;
        %     end
        % end

        %lam_var = (cumsum(lam)/tot_lam)*100; % scree plot
        % plot(1:numel(lam_var),lam_var,'-o','LineWidth',1);
        % legend('show');
        % title('Scree plot')
        % xlabel('Principal Component');
        % ylabel('Variance accounted for by principal component');
        % grid on
        % hold on

        pcmp = pcmp(:,1:pcmp_num); % truncate list of principal components
        [~,pcmp_col] = size(pcmp); % combine principal component features
        zro_col = all(pcmp_list == 0,1);
        frst_zro_col = find(zro_col,1,'first');
        pcmp_list(:,frst_zro_col:frst_zro_col+pcmp_col-1) = pcmp;
    end
    sing_grp = sing_grp(~all(sing_grp==0,1));
    [~,sing_grp_col] = size(sing_grp);
    zro_col = all(pcmp_list == 0,1);
    frst_zro_col = find(zro_col,1,'first'); % redefined
    pcmp_list(:,frst_zro_col:frst_zro_col+sing_grp_col-1) = varargin{1}(:,sing_grp);
    pcmp_list = pcmp_list(:,~all(pcmp_list == 0,1)); % remove zero cols
    % have yet to perform pca second time on the reduced features
end
%
function dir_grph(varargin)
    figure(8);
    for i = 2:nargin
        subplot(ceil((nargin-1)/3),3,i-1);
        node1 = varargin{i}(:,1);
        node2 = varargin{i}(:,2);
        plot(graph(node1,node2,'omitselfloops'));
        title(sprintf('Directed Graph of Patient %s:\n%s', ...
              varargin{1},inputname(i)));
    end
end
% plot eeg signals of all channels
function plt_eeg_chnnls(varargin) % first input is hour_first_last
    figure(1);
    title('EEG signal data of each bipolar channel');
    xlabel('Hour stamp');
    ylabel('Amplitude (mV)');
    hr_frst = varargin{1}(1);
    hr_lst = varargin{1}(2);
    for i = 2:nargin
        [mrows,~] = size(varargin{i});
        hr_mrk = linspace(hr_frst,hr_lst,mrows);
        plot(hr_mrk,varargin{i}.');
        hold on
    end
end
% plot wavelet scalograms
function plt_wve(varargin)
    for i = 2:nargin
        figure(i+100);
        scales = 1:128; % experiment with
        Fs = 100;
        cwt(varargin{i},Fs);
        title(sprintf('Magnitude Scalogrom for %s',inputname(i)));
        xlabel('Time (hr)');
    end
end
