clear

It=readmatrix('experiment_data.csv','Range','C3:C6000003');% transmitted light

fs=25e6;                    % sampling frequency
fre=21e3;                   % modulation frequency
L=length(It);               % data number 4e6+1
t=(0:L-1)'/fs;              % time
f=(0:L-1)'*fs/L;            % frequency



H2=LIAlp(It,fs,2*fre);  % Second harmonic demodulated by lock-in amplification

ave=round(fs/fre);

I=It-movmean(movmean(movmean(It,1.1*ave),1.1*ave),1.1*ave);
fft_tem=fft(I);
fft_tem((0.5*fre>f | f>2.5*fre)& ...
    ((f>(fs-0.5*fre))|((fs-2.5*fre)>f)))=0;
I=ifft(fft_tem);     % band-pass filtered signal of It, Only the 1f,2f component is retained
z1=abs(hilbert(I));  % the envelope of I
fft_tem=fft(z1);
fft_tem((0.8*fre>f | f>1.2*fre)&((fs-1.2*fre)>f | f>(fs-0.8*fre)))=0;
z1f=ifft(fft_tem);   % the 1f component of z1
Ht2=abs(hilbert(z1f));  
fft_tem=fft(Ht2);
fft_tem((f>1000)&((fs-1000)>f))=0;
Ht2=ifft(fft_tem);   % second harmonic demodulation based on Hilbert transform

%% Plot
t=t-0.0518554;

It=It(t>=0 & t<=0.15);
Ht2=Ht2(t>=0 & t<=0.15);
H2=H2(t>=0 & t<=0.15);
t=t(t>=0 & t<=0.15);
tcenter=((abs(t-0.02865)==min(abs(t-0.02865)))| ...
    (abs(t-0.0746)==min(abs(t-0.0746)))| ...
    (abs(t-0.129005)==min(abs(t-0.129005)))...
    );

error=Ht2-H2;
error_rel=error./H2; %relative error

N1=7001;             %There are too many data points, 7001 of them are plotted
L1=round(linspace(1,size(t,1),N1)');


%% Transmitted light
figure("Name","the It");
    plot(t(L1),It(L1),'Color','#2486b9');
    ylabel('Intensity (a.u.)');xlabel("Time (s)");
    set(gca,'FontSize',16,'FontName','Times New Roman','FontWeight','bold')


figure('Name','plot together');
subplot(2,1,1);
  plot(t(L1),H2(L1),'b',t(L1),Ht2(L1),'r');
    legend('lock-in amplification','HT based method');legend('boxoff');
    xlabel('Time (s)');ylabel('Intensity (a.u.)');ylim([0,0.035]);
    yLim=get(gca,'YLim');
    ax=gca;
    set(gca,'FontSize',11,'FontName','Times New Roman','FontWeight','bold')
    subplot(2,1,1);title('(a)','Position',[ax.YLabel.Position(1),yLim(2)]);
subplot(2,1,2);
  plot(t(L1),error_rel(L1)*100,'k');
  ylim([-10,10])
    xlabel('Time (s)');ylabel('Relative error (%)');
    ax=gca;
    set(gca,'FontSize',11,'FontName','Times New Roman','FontWeight','bold');
    xLim=get(gca,'XLim');yLim=get(gca,'YLim');
    subplot(2,1,2);title('(b)','Position',[ax.YLabel.Position(1),yLim(2)]);
    hold on;scatter(t(tcenter),100*error_rel(tcenter),'*',MarkerEdgeColor=([0 0.4470 0.7410]));hold off
    error_rel0=error_rel(tcenter);
    text(0.0239,-0.1377*yLim(2)+0.1377*yLim(1),[num2str(round(10000*error_rel0(1))/100),'%'], ...
        'FontSize',8,'FontName','Times New Roman','FontWeight','bold');
    text(0.0675,-0.1377*yLim(2)+0.1377*yLim(1),[num2str(round(10000*error_rel0(2))/100),'%'], ...
        'FontSize',8,'FontName','Times New Roman','FontWeight','bold');
    text(0.1245,-0.1377*yLim(2)+0.1377*yLim(1),[num2str(round(10000*error_rel0(3))/100),'%'], ...
        'FontSize',8,'FontName','Times New Roman','FontWeight','bold');
    grid on;set(gca,'GridLineStyle','--');

function [out]=LIAlp(fcn,fs,fre) %lock-in amplification
% [output]=LIA(input,sampling_frequency,reference_frequency)
    L=size(fcn,1);
    t=(0:L-1)'/fs;
    f=(0:L-1)'*fs/L;
    sinw=sin(2*pi*fre*t);
    cosw=cos(2*pi*fre*t);
    mixs=sinw.*fcn;
    mixc=cosw.*fcn;
    
    fft_tap=fft(mixs);
    fft_tap(1000<f & f<fs-1000)=0;
    outdcs=ifft(fft_tap);
    
    fft_tap=fft(mixc);
    fft_tap(1000<f & f<fs-1000)=0;
    outdcc=ifft(fft_tap);
    
    out=2*sqrt(outdcc.^2+outdcs.^2);
    
    fft_tap=fft(out);
    fft_tap(1000<f & f<fs-1000)=0;
    out=ifft(fft_tap);
end
