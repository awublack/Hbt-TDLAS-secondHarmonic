% Due to spectral leakage, a higher sampling rate is required to
% numerically simulate the harmonic signals of TDLAS at low concentrations.
% To verify this, you can change the sample rate, run the code, and look at the spectrum of the signal.
% Therefore, the settings of sampling rate, modulation frequency, 
% and scanning frequency in the simulation are different from those of the experiments.
clear
%% Initialization
fs=168e6;           % sampling frequency(Hz)
L=fs*0.125+1;       % sampling number
t=(0:L-1)'/fs;      % time series
f=(0:L-1)'*fs/L;    % frequency series
fre=fs/100;         % modulation frequency(Hz)
fsaw=4;             % scanning frequency(Hz)
phi=pi/6;           % phase delay between wavelength modulation and intensity modulation
dc=70;              % dc bias of injection current(mA)
a_saw=20;           % amplitude of sawtooth(mA)
a_sine= 2.3491;     % amplitude of sinewave(mA), corresponding to the modulation index of 2.2
%% Laser parameters
laser=@(a) 0.1*(a-5);                       % intensity of laser(mW),a(mA)
wn=@(a) 1e7./(0.026*a+2001.178);            % wavenumber of laser(cm-1),a(mA)
current_lag=@(t) dc+a_saw*sawtooth(2*pi*fsaw*t,0.5)+...
               a_sine*sin(2*pi*fre*t-phi);  % injection current(mA)
%% Emitted light of laser
saw=a_saw*sawtooth(2*pi*fsaw*t,0.5);    % sawtooth
sine=a_sine*sin(2*pi*fre*t);            % sinewave
current=saw+dc+sine;                    % injection current(mA)
Io=laser(current_lag(t));               % laser emitting intensity(mW)
wavenumber=wn(current);                 % wavenumber(cm-1)
twavenumber=wn(saw+dc);                 % wavenumber of scanning signal
%% Lorentian profile
S=0.0306;                   % the spectral line intensity(cm-2/atm) at normal atmosphere,296K(Tref=296K)
v0=4992.516;                % the wavenumber of the spectral line(cm-1),wavelength 2002.998nm
L_path=1.1;                 % path length(cm)            
gamma=0.0692;               % HMHW(cm-1)

N=15;                   % number of data points, depending on the performance of the computer. 
                        % if your computer is good, you can you can set a higher value

C=zeros(N,1);           % concentration of CO2
alpha_v0=zeros(N,1);    % absorbance in the center of spectral line
error_rel=zeros(N,1);   % relative error in the center of spectral line

for i=1:N
    alpha_v0(i)=0.05/N*i;
    C(i)=alpha_v0(i)*pi*gamma/S/L_path;
    if C(i)>=1
        C(i)=[];
        alpha_v0(i)=[];
        break
    end
    alpha=@(v) S*C(i)*L_path*gamma./(pi*(gamma^2+(v-v0).^2));% absorbance
    absorb=exp(-alpha(wavenumber)); % absorption function
    It=Io.*absorb;                  % transmitted light intensity(mW)
    ave=round(fs/fre);
    I=It-movmean(movmean(It,4*ave),4*ave);
    fft_buffer=fft(I);
    fft_buffer((0.8*fre>f | f>1.2*fre)&(1.9997*fre>f | f>2.0003*fre)& ...
        ((fs-1.2*fre)>f | f>(fs-0.8*fre))&((fs-2.0003*fre)>f | f>(fs-1.9997*fre)))=0;
    I=ifft(fft_buffer);             % band-pass filtered signal of It, Only the 1f,2f component is retained
    z1=abs(hilbert(I));             % envelope of I
    fft_buffer=fft(z1-movmean(z1,ave));
    fft_buffer((0.9997*fre>f | f>1.0003*fre)&((fs-1.0003*fre)>f | f>(fs-0.9997*fre)))=0;
    z1f=ifft(fft_buffer);           % 1f component of z1
    Ht2=abs(hilbert(z1f));          % Second harmonic demodulation based on Hilbert transform
    
    H2=LIA(It,fs,2*fre,fre);        % Second harmonic demodulated by lock-in amplification

    error_rel(i)=(Ht2(round(L/2))-H2(round(L/2)))./H2(round(L/2));
    fprintf('%d ',round(i/N*100));
end

% figure;plot(alpha_v0,error_rel);

xx=linspace(0,100*alpha_v0(end),200);
yy=spline(100*alpha_v0,100*error_rel,xx);   % data interpolation to make the curve look smoother

%% Plot
figure;
    plot(xx,yy,'k',LineWidth=1);ylim([-4,0]);
    ax1=gca;ax1.Box='off';
    set(ax1,'FontSize',11,'FontName','Times New Roman','FontWeight','bold')
    xlabel('Absorbance (%)');ylabel('Relative error (%)')
    ax1.XTick=linspace(ax1.XTick(1),ax1.XTick(end),6);
    ax2=axes('Position',ax1.Position,'Color','none', ...
        'XAxisLocation','top','YAxisLocation','right',...
        'Box','off');
    ax2.XTick=round(100*ax1.XTick*pi*gamma/S/L_path)/100;
    ax2.XLim=ax1.XLim*pi*gamma/S/L_path;
    ax2.YTick=ax1.YTick;ax2.YTickLabel=[];
    ax2.YLim=ax1.YLim;
    xlabel('Concentration of CO_2 (%)')
    set(ax2,'FontSize',11,'FontName','Times New Roman','FontWeight','bold') 
    ax1.Position=[ax1.Position*0.95];
    ax2.Position=ax1.Position;
    axes(ax1);
    grid on;
    axes(ax2);

function [out]=LIA(fcn,fs,fre,filter) %lock-in amplification
    % [output]=LIA(input,sampling_frequency,reference_frequency,cutoff_frequency)
    t=(1:size(fcn))'/fs;
    ave=round(fs/filter);
    sinw=sin(2*pi*fre*t);
    cosw=cos(2*pi*fre*t);
    mixs=sinw.*fcn;
    mixc=cosw.*fcn;
    outdcs=movmean(movmean(mixs,ave),ave);
    outdcc=movmean(movmean(mixc,ave),ave);
    out=2*sqrt(outdcc.^2+outdcs.^2);
    out(1:ave)=out(ave+1);
    out(end-ave+1:end)=out(end-ave);
end