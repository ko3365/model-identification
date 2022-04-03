clc; clear; close all;

load data/u1_impulse.mat
y11 = u1_impulse.Y(3).Data;
y21 = u1_impulse.Y(4).Data;
u1 = u1_impulse.Y(1).Data; %%% note that the pulse magnitude is 5
[m,mi] = max(u1>0); %%% find index where pulse occurs
load data/u2_impulse.mat
y12 = u2_impulse.Y(3).Data;
y22 = u2_impulse.Y(4).Data;
u2 = u2_impulse.Y(2).Data;
%%% remove any offsets in output data using data prior to pulse application
y11 = y11 - mean(y11([1:mi-1]));
y12 = y12 - mean(y12([1:mi-1]));
y21 = y21 - mean(y21([1:mi-1]));
y22 = y22 - mean(y22([1:mi-1]));
%%% rescale IO data so that impulse input has magnitude 1
y11 = y11/max(u1);
y12 = y12/max(u2);
y21 = y21/max(u1);
y22 = y22/max(u2);
u1 = u1/max(u1);
u2 = u2/max(u2);
ts = 1/40; %%%% sample period
N = length(y11); %%%% length of data sets
t = [0:N-1]*ts - 1;

% =========== H100 ============
number = 100;
H = zeros(2*number, 2*number);
for i = 1:number
    for j = 1:number
        k = i+j-1;
        H(2*i-1,2*j-1) = y11(k+mi);
        H(2*i,2*j-1) = y21(k+mi);
    
        H(2*i-1,2*j) = y12(k+mi);
        H(2*i,2*j) = y22(k+mi);
    end
end

[U, S, V] = svd(H);

singular_value = diag(S);
x = 1:number*2;

figure(3);
semilogy(x,singular_value(1:number*2),'ro','MarkerFaceColor','r','Markersize',4)
grid on
ylabel('Hankel singular value of $H_{100}$','Fontsize',14,'Interpreter','Latex');
xlabel('singular value index','Fontsize',14,'Interpreter','Latex')
% legend('H_2_0', 'H_4_0', 'H_8_0', 'H_1_0_0')

H_tilt = zeros(2*number, 2*number);
for i = 1:number
    for j = 1:number
        k = i+j;
        H_tilt(2*i-1,2*j-1) = y11(k+mi);
        H_tilt(2*i,2*j-1) = y21(k+mi);
    
        H_tilt(2*i-1,2*j) = y12(k+mi);
        H_tilt(2*i,2*j) = y22(k+mi);
    end
end

%% H_100, ns = [6,7,10,20]

% build ns = [6,7,10,20] model
ns = 6;
H100_6 = U(:,1:ns)*S(1:ns,1:ns)*V(:,1:ns)';
ns = 7;
H100_7 = U(:,1:ns)*S(1:ns,1:ns)*V(:,1:ns)';
ns = 10;
H100_10 = U(:,1:ns)*S(1:ns,1:ns)*V(:,1:ns)';
ns = 20;
H100_20 = U(:,1:ns)*S(1:ns,1:ns)*V(:,1:ns)';

% ========== H100_6 impulse response ===========
y11_6 = zeros(1, 201);
y21_6 = zeros(1, 201);
y12_6 = zeros(1, 201);
y22_6 = zeros(1, 201);

for i = 1:number
    for j = 1:number
        k = i+j-1;
        y11_6(k) = H100_6(2*i-1,2*j-1);
        y21_6(k) = H100_6(2*i,2*j-1);
    
        y12_6(k) = H100_6(2*i-1,2*j);
        y22_6(k) = H100_6(2*i,2*j);
    end
end

% ========== H100_7 impulse response ===========
y11_7 = zeros(1, 201);
y21_7 = zeros(1, 201);
y12_7 = zeros(1, 201);
y22_7 = zeros(1, 201);

for i = 1:number
    for j = 1:number
        k = i+j-1;
        y11_7(k) = H100_7(2*i-1,2*j-1);
        y21_7(k) = H100_7(2*i,2*j-1);
    
        y12_7(k) = H100_7(2*i-1,2*j);
        y22_7(k) = H100_7(2*i,2*j);
    end
end

% ========== H100_10 impulse response ===========
y11_10 = zeros(1, 201);
y21_10 = zeros(1, 201);
y12_10 = zeros(1, 201);
y22_10 = zeros(1, 201);

for i = 1:number
    for j = 1:number
        k = i+j-1;
        y11_10(k) = H100_10(2*i-1,2*j-1);
        y21_10(k) = H100_10(2*i,2*j-1);
    
        y12_10(k) = H100_10(2*i-1,2*j);
        y22_10(k) = H100_10(2*i,2*j);
    end
end

% ========== H100_20 impulse response ===========
y11_20 = zeros(1, 201);
y21_20 = zeros(1, 201);
y12_20 = zeros(1, 201);
y22_20 = zeros(1, 201);

for i = 1:number
    for j = 1:number
        k = i+j-1;
        y11_20(k) = H100_20(2*i-1,2*j-1);
        y21_20(k) = H100_20(2*i,2*j-1);
    
        y12_20(k) = H100_20(2*i-1,2*j);
        y22_20(k) = H100_20(2*i,2*j);
    end
end

t = [0:ts:5];

figure(4); 
hold on
title('Original vs Simulated Impulse Response for $y_{11}$','FontSize',14,'Interpreter','Latex')
plot(t,y11_6,'go','MarkerFaceColor','g','Markersize', 6)
plot(t,y11_7,'bo','MarkerFaceColor','b','Markersize', 6)
plot(t,y11_10,'co','MarkerFaceColor','c','Markersize', 6)
plot(t,y11_20,'mo','MarkerFaceColor','m','Markersize', 6)
plot(t,y11(mi+1:mi+1+200),'ko','LineWidth',0.8)
grid on
axis([0 2 -0.1 0.1])
ylabel('$y_{11}$ (volts)','FontSize',14,'Interpreter','Latex');
xlabel('second','FontSize',14,'Interpreter','Latex')
legend('ns=6', 'ns=7', 'ns=10', 'ns=20', 'original')
hold off

figure(5); 
hold on
title('Original vs Simulated Impulse Response for $y_{21}$','FontSize',14,'Interpreter','Latex')
plot(t,y21_6,'go','MarkerFaceColor','g','Markersize', 6)
plot(t,y21_7,'bo','MarkerFaceColor','b','Markersize', 6)
plot(t,y21_10,'co','MarkerFaceColor','c','Markersize', 6)
plot(t,y21_20,'mo','MarkerFaceColor','m','Markersize', 6)
plot(t,y21(mi+1:mi+1+200),'ko','LineWidth',0.8)
grid on
axis([0 2 -0.1 0.1])
ylabel('$y_{21}$ (volts)','FontSize',14,'Interpreter','Latex');
xlabel('second','FontSize',14,'Interpreter','Latex')
legend('ns=6', 'ns=7', 'ns=10', 'ns=20', 'original')
hold off

figure(6); 
hold on
title('Original vs Simulated Impulse Response for $y_{12}$','FontSize',14,'Interpreter','Latex')
plot(t,y12_6,'go','MarkerFaceColor','g','Markersize', 6)
plot(t,y12_7,'bo','MarkerFaceColor','b','Markersize', 6)
plot(t,y12_10,'co','MarkerFaceColor','c','Markersize', 6)
plot(t,y12_20,'mo','MarkerFaceColor','m','Markersize', 6)
plot(t,y12(mi+1:mi+1+200),'ko','LineWidth',0.8)
grid on
axis([0 2 -0.1 0.1])
ylabel('$y_{12}$ (volts)','FontSize',14,'Interpreter','Latex');
xlabel('second','FontSize',14,'Interpreter','Latex')
legend('ns=6', 'ns=7', 'ns=10', 'ns=20', 'original')
hold off

figure(7); 
hold on
title('Original vs Simulated Impulse Response for $y_{22}$','FontSize',14,'Interpreter','Latex')
plot(t,y22_6,'go','MarkerFaceColor','g','Markersize', 6)
plot(t,y22_7,'bo','MarkerFaceColor','b','Markersize', 6)
plot(t,y22_10,'co','MarkerFaceColor','c','Markersize', 6)
plot(t,y22_20,'mo','MarkerFaceColor','m','Markersize', 6)
plot(t,y22(mi+1:mi+1+200),'ko','LineWidth',0.8)
grid on
axis([0 2 -0.1 0.1])
ylabel('$y_{22}$ (volts)','FontSize',14,'Interpreter','Latex');
xlabel('second','FontSize',14,'Interpreter','Latex')
legend('ns=6', 'ns=7', 'ns=10', 'ns=20', 'original')
hold off

%% (A,B,C) model of H_100, ns = [6,7,10,20] 

ns = 6;
O_6 = U(:,1:ns)*sqrt(S(1:ns,1:ns));
C_6 = sqrt(S(1:ns,1:ns))*V(:,1:ns)';
A6 = pinv(O_6) * H_tilt * pinv(C_6);
C6 = O_6(1:2, 1:ns);
B6 = C_6(1:ns, 1:2);

ns = 7;
O_7 = U(:,1:ns)*sqrt(S(1:ns,1:ns));
C_7 = sqrt(S(1:ns,1:ns))*V(:,1:ns)';
A7 = pinv(O_7) * H_tilt * pinv(C_7);
C7 = O_7(1:2, 1:ns);
B7 = C_7(1:ns, 1:2);

ns = 10;
O_10 = U(:,1:ns)*sqrt(S(1:ns,1:ns));
C_10 = sqrt(S(1:ns,1:ns))*V(:,1:ns)';
A10 = pinv(O_10) * H_tilt * pinv(C_10);
C10 = O_10(1:2, 1:ns);
B10 = C_10(1:ns, 1:2);

ns = 20;
O_20 = U(:,1:ns)*sqrt(S(1:ns,1:ns));
C_20 = sqrt(S(1:ns,1:ns))*V(:,1:ns)';
A20 = pinv(O_20) * H_tilt * pinv(C_20);
C20 = O_20(1:2, 1:ns);
B20 = C_20(1:ns, 1:2);
