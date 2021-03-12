clear all
load('data/cond.mat')
set(groot, 'defaultAxesTickLabelInterpreter','latex');
 set(groot, 'defaultLegendInterpreter','latex');
freqs=freqCond(1,2:end);
conds=freqCond(2,2:end);
norms=freqCond(3,2:end);
norminvs=freqCond(4,2:end);

for j=1:length(freqs)
   if imag(freqs(j))<0
       freqs(j)=conj(freqs(j));
   end
end


[~, indices]= sort(real(freqs));
freqscomplete=[freqs(indices),flip(conj(freqs(indices)))];
condscomplete=[conds(indices),flip(conds(indices))];
normscomplete=[norms(indices),flip(norms(indices))];
norminvscomplete=[norminvs(indices),flip(norminvs(indices))];


h=figure(1)
marksize=100;
subplot(1,2,2)
semilogy(condscomplete,'k')

xlim([0, 606])
hold on

semilogy(normscomplete, '-.k')
semilogy(norminvscomplete, '--k')
scatter([3,604],condscomplete([3,604]),marksize,'d','k')
scatter([60,607-60],condscomplete([60,607-60]),marksize,'s','k')
scatter([141,607-141],condscomplete([141,607-141]),marksize,'o','k')

legend('cond$_2(A_h(s))$', '$||A_h(s)||_2$', ' $||A_h(s)^{-1}||_2$')
% freqs=[freqCond(1,2:end) ,conj(freqCond(1,2:end))];
% conds=[freqCond(1,2:end) ,freqCond(1,2:end)];
% norms=[freqCond(1,2:end) ,freqCond(1,2:end)];
% norminvs=[freqCond(1,2:end) ,freqCond(1,2:end)];

subplot(1,2,1)
smmarksize=8;
scatter(real(freqscomplete),imag(freqscomplete),smmarksize,'k.')
hold on
scatter(real(freqscomplete([3,604])),imag(freqscomplete([3,604])),marksize,'d','k')
scatter(real(freqscomplete([60,607-60])),imag(freqscomplete([60,607-60])),marksize,'s','k')
scatter(real(freqscomplete([141,607-141])),imag(freqscomplete([141,607-141])),marksize,'o','k')
xlabel('Re $s$','Interpreter','Latex')
ylabel('Im $s$','Interpreter','Latex')

set(h,'Position',[10 10 900 420])
%saveas(gcf,'Contour_Conditions','epsc')  
