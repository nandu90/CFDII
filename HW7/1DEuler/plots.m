clear all
B = importdata('plot.dat',' ')

C = importdata('out2.txt',' ')
A = importdata('out1.txt',' ') 

figure
set(gcf,'color','w');
set(gcf,'Position', [50, 50, 720, 500]);
plot(B(:,1),B(:,3),'-b','LineWidth',2)
hold on
plot(A(:,1),A(:,2),'*-r')

plot(C(:,1),C(:,2),'*-k')
legend('Exact','Minmod Limiter','MC Limiter')
title('Comparison of Limiter Schemes (Conservative Variables), 1-Stage R-K','FontSize',12,'FontWeight','bold')
ylabel('rho','FontSize',12,'FontWeight','bold')
xlabel('x','FontSize',12,'FontWeight','bold')
ylim([0 1.1])
hold off

figure
set(gcf,'color','w');
set(gcf,'Position', [50, 50, 720, 500]);
plot(B(:,1),B(:,2),'-b','LineWidth',2)
hold on
plot(A(:,1),A(:,4),'*-r')

plot(C(:,1),C(:,4),'*-k')
legend('Exact','Minmod Limiter','MC Limiter')
title('Comparison of Limiter Schemes (Conservative Variables), 1-Stage R-K','FontSize',12,'FontWeight','bold')
ylabel('u','FontSize',12,'FontWeight','bold')
xlabel('x','FontSize',12,'FontWeight','bold')


hold off

figure
set(gcf,'color','w');
set(gcf,'Position', [50, 50, 720, 500]);
plot(B(:,1),B(:,4),'-b','LineWidth',2)
hold on
plot(A(:,1),A(:,6),'*-r')

plot(C(:,1),C(:,6),'*-k')
legend('Exact','Minmod Limiter','MC Limiter')
title('Comparison of Limiter Schemes (Conservative Variables), 1-Stage R-K','FontSize',12,'FontWeight','bold')
ylabel('p','FontSize',12,'FontWeight','bold')
xlabel('x','FontSize',12,'FontWeight','bold')
ylim([0 1.1])

hold off

figure
set(gcf,'color','w');
set(gcf,'Position', [50, 50, 720, 500]);
for i=1:length(B(:,2))
    M(i)=B(i,2)/B(i,5);
end
plot(B(:,1),M,'-b','LineWidth',2);
hold on
plot(A(:,1),A(:,8),'*-r');
plot(C(:,1),C(:,8),'*-k');
legend('Exact','Minmod Limiter','MC Limiter')
title('Comparison of Limiter Schemes (Conservative Variables), 1-Stage R-K','FontSize',12,'FontWeight','bold')
ylabel('Mach','FontSize',12,'FontWeight','bold')
xlabel('x','FontSize',12,'FontWeight','bold')
ylim([0 1.1])
hold off

%Steger-Warming
%Van Leer