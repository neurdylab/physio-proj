function [OUT_p,REGS] = make_hcp_physio(po_wave,resp_wave,fs_phys,TR_s,nframes,DISP)
%
% function REGS = make_hcp_physio(po_wave,resp_wave,fs_phys,TR_sec,nframes);  
% make regressors for hcp physio data
%
% Inputs:
% ========
% po_wave: PPG signal
% resp_wave: respiration signal
% fs_phys: physio sampling rate
% TR_s: fMRI TR (in sec)
% nframes: length of fMRI scan (frames)
% DISP: 1 to display figures, 0 otherwise
%
% Outputs:
% ========
% REGS: struct containing respiratory volume (rv), heart rate (hr),
% and pulse wave amplitude signals, sampled at the fMRI TR
% OUT_p: struct containing intermediate, derived outputs from raw
% physiological signals

  
dt_phys = 1/fs_phys;

%% resp:
resp.wave = resp_wave - mean(resp_wave);

%% detect cardiac peaks
card_dat = po_wave - mean(po_wave);
fcut_BPF = [0.5,4];
Fn = fs_phys/2;
Wn = fcut_BPF/Fn; Nb = 2;
[B, A] = butter(Nb,Wn);
card_bpf =  filtfilt(B,A,double(card_dat));
card_rng = iqr(card_bpf);


minHeight = 0.05*card_rng;
minDist = (fs_phys/2); 

[pks,locs] = findpeaks(card_bpf,'minpeakheight',minHeight,'minpeakdistance',minDist);
clear maxtab_c; maxtab_c(:,1) = locs; maxtab_c(:,2) = pks;
% extract cardiac trigger times and approximate heart rate
card_trig_samples = locs;
card_trig_times = card_trig_samples*dt_phys;
% ibi & "instantaneous" hr
IBI = (diff(card_trig_times));
HR = (1./diff(card_trig_times))*60; %bpm


% check peak detection 
% ------------------------------------------ %
if (DISP)
  figure(101); clf; set(gcf,'color','w'); 
  g1 = subplot(3,1,1); hold on;
  plot(card_bpf); 
  plot(card_dat-mean(card_dat),'color',0.8*[1 1 1]);
  plot(maxtab_c(:,1),maxtab_c(:,2),'r.','markersize',20); 
  xlabel('physio sample #');
  % check heart rate as a function of time
  figure(101);
  g2 = subplot(3,1,2); hold on;
  plot(card_trig_samples(1:end-1),HR,'b'); ylabel('beats per min'); title('cardiac rate');
  hold on;
  plot(card_trig_samples(1:end-1),HR,'c.'); % +dots
  xlabel('physio sample #');
  linkaxes([g1,g2],'x');
  %ylim([40 95]);
  % IBI time series (by sample number - flag outliers)
  figure(101);
  g3 = subplot(3,1,3); hold on;
  plot(IBI); 
  xlabel('index'); ylabel('IBI');
  drawnow;
end

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

% generate fMRI regressors
% ------------------------------------------ %
REGS = build_fmri_regs(OUT_p,TR_s,nframes,DISP); 



function REGS = build_fmri_regs(IN_p,TR_s,nframes,DISP)
% make some regressors sampled in TR bins
  
  dt_phys = IN_p.dt_phys;
  resp = IN_p.resp;
  card_dat = IN_p.card_dat;
  card_bpf = IN_p.card_bpf;
  IBI_clean = IN_p.IBI_raw;
    
  % lookup table of IBI (denoised) v. cardiac trig
  % time (assigned to halfway between the respective beats)
  t_ibi = 0.5*(IN_p.card_trig_times_s(2:end) + IN_p.card_trig_times_s(1:end-1));
  assert(length(t_ibi)==length(IBI_clean))
  
  % sampling to match fMRI tr (center of each tr)
  Twin = 6; % sec windows (3s on either side)
  t_fmri = (TR_s/2)+[0:TR_s:TR_s*nframes-TR_s];
  
  % make RV, HR, and pulse wave amplitude regressors
  rv = [];
  pa = [];
  hr = [];
  
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
    pa_bpf(kk) = std(card_bpf(i1:i2));
    % respiration variation
    % ---------------------- %
    rv(kk) = std(resp.wave(i1:i2));
  end
    

  % regressors for fmri
  % ---------------------- %
  REGS.rv = rv(:);
  REGS.hr = hr(:);
  REGS.hrv = hrv(:);
  REGS.pa = pa(:);
  
  