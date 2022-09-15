% Initial processing of NKI-RS physio data. Align with fMRI scan, extract
% respiratory and cardiac measures. All variables in OUT_p are aligned to fMRI.


function physio_preproc_nki(sub_id, my_path)
% function physio_preproc_nki(sub_id, my_path)
% sub_id: NKI-RS subject ID (sub-A*)
% my_path: path to save results

    fs_phys = 62.5; % nki sampling rate physio
    TR = 1.4; % fMRI interval between volumes
    nframes = 404; % number of fMRI time frames
        
    % directory for saving
    dir = append(my_path, id, "/");
    if ~isempty(dir)
        if ~exist(dir,'dir')
            mkdir(dir);
        end
    end
    dir = append(dir, '/');
    
    % read physio file (nki)
    tsv_path = append('/home/songrw/neurdy/tsvfiles/', id, '_ses-BAS1_task-rest_acq-1400_physio.tsv');
    disp(tsv_path)
    phys = load(tsv_path);
    [ro, co] = size(phys); 
    
    % ------------------------------------------ %
    if co == 2
        TABLE(:,1) = phys(:, 1);%mri trigs;% change to read col of tsv file
        TABLE(:,2) = phys(:, 1);%cardiac;
        TABLE(:,3) = phys(:, 2);%respiratory;
    elseif co == 3
        TABLE(:,1) = phys(:, 1);%mri trigs;% change to read col of tsv file
        TABLE(:,2) = phys(:, 2);%cardiac;
        TABLE(:,3) = phys(:, 3);%respiratory;
    else 
        TABLE(:,1) = phys(:, 4);%mri trigs;% change to read col of tsv file
        TABLE(:,2) = phys(:, 2);%cardiac;
        TABLE(:,3) = phys(:, 3);%respiratory;
    end 
    dt_phys = 1/fs_phys; 
    tax = [0:dt_phys:size(TABLE,1)*dt_phys-dt_phys];
    
    % check raw data
    % ------------------------------------------ %
    figure; clf; set(gcf,'color','w');
    g1 = subplot(311);
    if co == 4
        plot(TABLE(:,1)); title('MR triggers');
    end 
    g2 = subplot(312);
    plot(TABLE(:,2)); title('cardiac (raw)');
    g3 = subplot(313);
    plot(TABLE(:,3)); title('respiration (raw)');
    xlabel('time (samples)');
    linkaxes([g1 g2 g3],'x'); %% what does this do
    % saving first figure
    save_path_1 = append(dir, id, '_ses-BAS1_task-rest_acq-1400_physio_fig1.png');
    saveas(gcf, save_path_1); 
    
    % extract respiration data 
    % ------------------------------------------ %
    resp_dat = TABLE(:,3); % aligned with fMRI
    resp.dt = dt_phys;
    resp.wave = double(resp_dat); %resp is 1x1 struct with dt & wave
        
    % extract cardiac data & run peak detection
    % after some initial filtering
    % ------------------------------------------ %
    card_dat = TABLE(:,2);  % raw cardiac trace, *aligned with fMRI scan*
    % band-pass filter to help peak detection
    fcut_BPF = [0.5,2];
    Fn = fs_phys/2; % sampling rate/2
    Wn = fcut_BPF/Fn; Nb = 2; % ??
    [B, A] = butter(Nb,Wn);
    card_bpf =  filtfilt(B,A,double(card_dat));
    % peak detection
    card_rng = iqr(card_bpf); 
    %[maxtab, mintab] = peakdet(card_bpf, 0.05*card_rng);
    % card_trig_samples = maxtab(:,1);
    minHeight = 0.05*card_rng;
    %minDist = 100+(fs_phys/2); % 04.01.17
    [pks,locs] = findpeaks(card_bpf,'minpeakheight',minHeight); %,'minpeakdistance',minDist);
    maxtab_c(:,1) = locs; maxtab_c(:,2) = pks;
    % extract cardiac trigger times and approximate heart rate
    card_trig_samples = locs;
    card_trig_times = card_trig_samples*dt_phys;
    % ibi & "instantaneous" hr
    IBI = (diff(card_trig_times));
    HR = (1./diff(card_trig_times))*60; %bpm
    
    
    % check peak detection 
    % ------------------------------------------ %
    figure; 
    g1 = subplot(3,1,1); hold on;
    plot(card_bpf); 
    plot(maxtab_c(:,1),maxtab_c(:,2),'r.','markersize',20); % does well!
    legend('cardiac data (band-pass filt)','detected beats');
    xlabel('physio sample #');
    % check heart rate as a function of time
    g2 = subplot(3,1,2); hold on;
    plot(card_trig_samples(1:end-1),HR,'b'); ylabel('beats per min'); title('cardiac rate');
    hold on;
    plot(card_trig_samples(1:end-1),HR,'c.'); % +dots
    xlabel('physio sample #');
    linkaxes([g1,g2],'x');
    % IBI time series (by sample number - flag outliers)
    g3 = subplot(3,1,3); hold on;
    plot(IBI); 
    xlabel('index'); ylabel('IBI');
    drawnow;
    % save second figure
    save_path_2 = append(dir, id, '_ses-BAS1_task-rest_acq-1400_physio_fig2.png');
    saveas(gcf, save_path_2); 
    
    % make output struct 
    % ------------------------------------------ %
    % note, all are aligned with fMRI
    OUT_p.IBI_raw = IBI;
    OUT_p.HR_raw = HR;
    OUT_p.card_trig_times_s = card_trig_times; %sec
    OUT_p.card_trig_samples = card_trig_samples;
    OUT_p.card_dat = card_dat;
    OUT_p.card_bpf = card_bpf;
    OUT_p.card_pks = maxtab_c;
    OUT_p.dt_phys = dt_phys;
    OUT_p.resp = resp;
    
    % remove outlier beats?
    % ------------------------------------------ %
    [OUT_p.IBI_clean,OUT_p.HR_clean,outliers] = despike_hr(OUT_p);
    OUT_p.outlier_beats = outliers;
    if ~isempty(outliers)
        save_path_4 = append(dir, id, '_ses-BAS1_task-rest_acq-1400_physio_fig4.png');
        saveas(gcf, save_path_4);       
    end
    
    % generate fMRI regressors
    % ------------------------------------------ %
    REGS = build_fmri_regs(OUT_p,TR,nframes);
    % save third figure
    save_path_3 = append(dir, id, '_ses-BAS1_task-rest_acq-1400_physio_fig3.png');
    saveas(gcf, save_path_3); 
    
    % save - MODIFY PATH AND FILENAME TO SAVE AS
    % ------------------------------------------ %
    savePath = append(dir, id, '_ses-BAS1_task-rest_acq-1400_physio','_physOUT.mat');
    save(savePath,'OUT_p','TABLE','REGS');
    
end 

function [IBI_clean,HR_clean,outliers] = despike_hr(IN_p)
% for manually despiking 
    outliers = [];
    figure;
    subplot(211);
    IBI_clean = interp_ts(IN_p.IBI_raw, outliers, 1);
    subplot(212);
    HR_clean = interp_ts(IN_p.HR_raw, outliers, 1);
    
end
    
    
function REGS = build_fmri_regs(IN_p,TR,nframes)
% make some regressors sampled in TR bins
% hr: need smaller bins to see hrv?? check this.
    
    dt_phys = IN_p.dt_phys;
    resp = IN_p.resp;
    card_dat = IN_p.card_dat;
    card_bpf = IN_p.card_bpf;
    IBI_clean = IN_p.IBI_clean;
    
    % lookup table of IBI (denoised) v. cardiac trig
    % time (assigned to halfway between the respective beats)
    t_ibi = 0.5*(IN_p.card_trig_times_s(2:end) + IN_p.card_trig_times_s(1:end-1));
    assert(length(t_ibi)==length(IBI_clean))
    
    % sampling to match fMRI tr (center of each tr)
    Twin = 6; % sec windows (3s on either side)
    TR_s = TR;
    t_fmri = (TR_s/2)+[0:TR_s:TR_s*nframes-TR_s];
    
    % make RV & HR regressors, as well as pulseox amplitude (stdev)
    rv = [];
    pa = [];
    hr = [];
    pa_bpf = [];
    for kk=1:nframes
        t = t_fmri(kk);
        
        % heart rate
        % ---------------------- %
        % get time bin centered at this TR
        t1 = max(0,t-Twin*0.5);
        t2 = min(TR_s*nframes,t+Twin*0.5);
        % find IBI's falling within this interval
        inds = intersect(find(t_ibi<=t2),find(t_ibi>=t1));
        hr(kk) = (60./median(IBI_clean(inds)));
        hrv(kk) = sqrt(mean(diff(IBI_clean(inds)).^2)); %rmssd
        
        % pulse amplitude
        % ---------------------- %
        if length(resp.wave)~=length(card_dat)
            error('resp & card sampled at different rates');
        else
            np = length(resp.wave);
        end
        % window (in samples)
        i1 = max(1,floor((t - Twin*0.5)/dt_phys)); 
        i2 = min(np, floor((t + Twin*0.5)/dt_phys));
        pa(kk) = std(card_dat(i1:i2));
 
        % respiration variation
        % ---------------------- %
        rv(kk) = std(resp.wave(i1:i2));
    end
    
    % regressors for fmri
    % ---------------------- %
    REGS.rv = rv(:);
    REGS.pa = pa(:);
    REGS.hr = hr(:);
    REGS.hrv = hrv(:);
    
    figure;
    subplot(411); 
    plot(hr); title('heart rate');
    subplot(412);
    plot(hrv); title('heart rate variability (rmssd)');
    subplot(413);
    plot(pa); title('pulse wave amplitude');
    subplot(414);
    plot(rv); title('respiratory volume');
    xlabel('fmri TR');
    
end