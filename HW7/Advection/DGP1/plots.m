clear all
A = importdata('out.txt',',') 
plot(A(:,1),A(:,2),'-*r')
title('DG(P_1) Method, 3-Stage R-K, (100 cells, CFL=0.1)','FontSize',12,'FontWeight','bold')
ylabel('u(x,t)','FontSize',12,'FontWeight','bold')
xlabel('x','FontSize',12,'FontWeight','bold')
hold on
grid on
plot(A(:,1),A(:,3),'-b')
plot(A(:,1),A(:,4),'-ok')

legend('t=1','t=0','t=5')