clear all
A = importdata('out3.txt',' ') 
B = importdata('out4.txt',' ') 
set(gcf,'color','w');
set(gcf,'Position', [50, 50, 720, 500]);
plot(A(:,1),A(:,3),'-b','LineWidth',2)
hold on
plot(A(:,1),A(:,2),'-*r')
title('Comparison of Minmod & MC limiter, t=5','FontSize',12,'FontWeight','bold')
ylabel('u(x,t)','FontSize',12,'FontWeight','bold')
xlabel('x','FontSize',12,'FontWeight','bold')
plot(B(:,1),B(:,2),'-*k')
%plot(A(:,1),A(:,4),'ok')

legend('Exact','Minmod Limiter','MC Limiter')