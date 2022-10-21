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
C=400/1e6;                  % concentration of CO2,400 ppm
v0=4992.516;                % the wavenumber of the spectral line(cm-1),wavelength 2002.998nm
L_path=1.1;                 % path length(cm)
gamma=0.0692;               % HMHW(cm-1)
alpha=@(v) S*C*L_path*gamma./(pi*(gamma^2+(v-v0).^2));% absorbance
%% Transmitted light
absorb=exp(-alpha(wavenumber)); % absorption function
It=Io.*absorb;                  % transmitted light intensity(mW)
%    plot(t,It);
%    figure; hold on;
%    fplot(alpha,[v0-1,v0+1]);
%    plot(wn(dc-a_saw:0.01:dc+a_saw),alpha(wn(dc-a_saw:0.01:dc+a_saw)),'Color','r','LineWidth',3);% Range of scan
%    hold off
%% Second harmonic demodulation based on Hilbert transform
ave=round(fs/fre);
I=It-movmean(movmean(It,4*ave),4*ave);
fft_buffer=fft(I);
fft_buffer((0.8*fre>f | f>1.2*fre)&(1.9997*fre>f | f>2.0003*fre)& ...
    ((fs-1.2*fre)>f | f>(fs-0.8*fre))&((fs-2.0003*fre)>f | f>(fs-1.9997*fre)))=0;
I=ifft(fft_buffer);     % band-pass filtered signal of It, Only the 1f,2f component is retained
z1=abs(hilbert(I));     % envelope of I
fft_buffer=fft(z1-movmean(z1,ave));
fft_buffer((0.8*fre>f | f>1.2*fre)&((fs-1.2*fre)>f | f>(fs-0.8*fre)))=0;
z1f=ifft(fft_buffer);   % 1f component of z1

Ht2=abs(hilbert(z1f));  % Second harmonic demodulation based on Hilbert transform
H2=LIA(It,fs,2*fre,fre);% Second harmonic demodulated by lock-in amplification

error=Ht2-H2;
error_rel=error./H2;    % relative error
%% plot

Iwant=find(abs(twavenumber-v0)<7.492*gamma);    % Intercept the signal near the center of the spectral
twant=linspace(0,t(Iwant(end))-t(Iwant(1)),length(Iwant));% time series

N1=7001;    %There are too many data points, 7001 of them are plotted
t1=linspace(0,t(Iwant(end))-t(Iwant(1)),N1);% time series
L1=round(linspace(Iwant(1),Iwant(end),N1)');% The data points that are plotted


%% Signal processing process
figure('Name',"Signal processing process");
set(gcf,"Position",[100 100 1120-28 840]);
% transmitted light
    ax1=axes("Position",[0.143/0.975,0.16,0.76/0.975,0.765]/2+[0,.5,0,0]);
    plot(twant,It(Iwant),'Color','#2486b9');
    xlim([0,0.05]);ylim([4.5,7.6]);
    ylabel('Intensity (a.u.)');xlabel("Time (s)");
    set(gca,'FontSize',18,'FontName','Times New Roman','FontWeight','bold')
    xLim=ax1.XLim;yLim=ax1.YLim;
    title('(a)','Position',[1.1*xLim(1)-0.1*xLim(2),yLim(2)]);

    ax1_1=axes("Position",[0.37/0.975,0.280,0.49/0.975,0.14]/2+[0,.5,0,0]);
    plot(twant,It(Iwant),'Color','#2486b9');
    ax1_1.YTick=[5.9,6.4];axis([0.01446,0.014475,5.7,6.6]);
    annotation('arrow',[0.3643/0.975,0.3821/0.975]/2,[0.5048,0.4286]/2+[0.5,.5]);
% band-pass filtered
    axes('Position',ax1.Position+[0.487,0,0,0]);
    plot(twant,I(Iwant),'Color','#2486b9');
    xlim([0,0.05]);ylim([-.6,.6]);
    ylabel('Intensity (a.u.)');xlabel("Time (s)");
    set(gca,'FontSize',18,'FontName','Times New Roman','FontWeight','bold');
    xLim=get(gca,'XLim');yLim=get(gca,'YLim');
    title('(b)','Position',[1.1*xLim(1)-0.1*xLim(2),yLim(2)]);
 % envelope
    ax3=axes('Position',[0.143/0.975,0.16,0.76/0.975,0.765]/2);
    plot(twant,z1(Iwant),'Color','#2486b9');
    xlim([0,0.05]); 
    ylim([0,0.299]);
    xlabel("Time (s)");ylabel('Intensity (a.u.)');
    set(gca,'FontSize',18,'FontName','Times New Roman','FontWeight','bold');
    xLim=get(gca,'XLim');yLim=get(gca,'YLim');
    title('(c)','Position',[1.1*xLim(1)-0.1*xLim(2),yLim(2)]);
    
    ax3_1=axes("Position",[0.42/0.975,0.236,0.42/0.975,0.42]/2);
    plot(twant,z1(Iwant),'Color','#2486b9')
    set(gca,'FontSize',11,'FontName','Times New Roman');
    xlim([0,0.05]);
    annotation('arrow',[.545,.546]/2/0.975,[.745,.664]/2,'Color',[0.30,0.75,0.93]);
    annotation("rectangle",[.1446/0.975,.7167,.7571/0.975,.0857]/2, ...
        'Color',[0.30,0.75,0.93],'LineStyle','--');
% 1f component
    axes('Position',ax3.Position+[0.487,0,0,0]);
    plot(twant,z1f(Iwant),'Color','#2486b9');
    ylabel('Intensity (a.u.)');xlabel("Time (s)");xlim([0,0.05]);ylim([-1.6e-4,1.6e-4]);
    set(gca,'FontSize',18,'FontName','Times New Roman','FontWeight','bold');
    xLim=get(gca,'XLim');yLim=get(gca,'YLim');
    title('(d)','Position',[1.1*xLim(1)-0.1*xLim(2),yLim(2)]);
%% Result
figure("Name",'The demodulated signal')
    plot(twant,Ht2(Iwant),'Color','#2486b9',LineWidth=2);xlim([0,0.05]);
    ylabel('Intensity (a.u.)');xlabel("Time (s)");
    set(gca,'FontSize',16,'FontName','Times New Roman','FontWeight','bold');
%% Comparison of the two methods
figure('Name',"对比图");
set(gcf,"Position",[100 100 1120-60 420]);
% plot together
   axes("Position",[0.143/0.946/2,0.16,0.76/0.946/2,0.765]);
    plot(t1,H2(L1),'b',t1,Ht2(L1),'r',LineWidth=1);
    xlim([0,0.05]);ylim([0,1.08*max(Ht2(L1))]);
    legend({'lock-in amplifier','HT based method'}, ...
        Position=[0.573/2,0.79,0.4/2,0.117]);legend('boxoff');
    xLim=get(gca,'XLim');yLim=get(gca,'YLim');
    title('(a)','Position',[1.1*xLim(1)-0.1*xLim(2),yLim(2)]);
    xlabel('Time (s)');ylabel('Intensity (a.u.)')
    set(gca,'FontSize',16,'FontName','Times New Roman','FontWeight','bold')

    axes('position',[0.176/2/0.946,0.664,0.22/2/0.946,0.20]);
    plot(t1,H2(L1),'b',t1,Ht2(L1),'r',LineWidth=1);
    axis([0.02493,0.02517,1.3807e-4,1.3817e-4]);
    set(gca,'xticklabel',[],'yticklabel',[]);
    annotation('rectangle',[0.49/2/0.946,0.84,0.064/2/0.946,0.057]);
    annotation('arrow',[.498/2/0.946,.384/2/0.946],[.86,.78]);
% relative error
    axes("Position",[0.143/0.946/2,0.16,0.76/0.946/2,0.765]+[.48,0,0,0]);
    plot(t1,100*error_rel(L1),'k',LineWidth=1);
    xlim([0,0.05]);
    xLim=get(gca,'XLim');yLim=get(gca,'YLim');
    title('(b)','Position',[1.1*xLim(1)-0.1*xLim(2),yLim(2)]);
    xlabel('Time (s)');ylabel('Relative error (%)')
    set(gca,'FontSize',16,'FontName','Times New Roman','FontWeight','bold')
    hold on;scatter(0.025,100*error_rel(round(L/2)),'*',MarkerEdgeColor=([0 0.4470 0.7410]));hold off
    text(0.023,100*error_rel(round(L/2)+1)-0.05*yLim(2)+0.05*yLim(1),[num2str(round(1000000*error_rel(round(L/2)+1))/10000),'%'], ...
        'FontSize',12,'FontName','Times New Roman','FontWeight','bold');
    xlim(xLim);ylim(yLim);

figure();
plot(t1,error(L1),'LineWidth',1,'Color','#2486b9');
xlim([0,0.05]);ylim([-6e-8,7e-8]);
xlabel("Time (s)");ylabel('Error (a.u.)');
% title('Error between HT based method and LIA');
set(gca,'FontSize',16,'FontName','Times New Roman','FontWeight','bold');
