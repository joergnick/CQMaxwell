clear all
load('data/cond.mat')
n=length(freqCond(1,2:end));
figure(1)
 subplot(1,2,1)
 title('Wirklich gelöste Frequenzen')
 subplot(1,2,2)
 title('Alle Frequenzen nach oben gespiegelt')
 
writerObj = VideoWriter('Frequencies.avi');
writerObj.FrameRate=40;
open(writerObj);
for j=2:n
 subplot(1,2,1)
  title('Wirklich gelöste s')
s=freqCond(1,j);
scatter(real(s),imag(s),'b')
ylim([-400,400])
xlim([0,600])
hold on
subplot(1,2,2)
title('Alle s nach oben gespiegelt')
if imag(s)<0
    s=conj(s);
end

scatter(real(s),imag(s),'b')
ylim([-400,400])
xlim([0,600])


frame=getframe(gcf);
drawnow
hold on
writeVideo(writerObj,frame);
end
for j=1:100
   writeVideo(writerObj,frame);
end
close(writerObj);
