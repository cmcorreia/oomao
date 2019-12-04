%% Impact of dL0 on r0 estimation.
L0     = 30;
r0     = 0.15;
k2     = 1.005;
nL0    = 101;
dL0    = linspace(-L0/2,L0,nL0);
r0     = 0.2;
Darray = [4 10 40];
nD     = length(Darray);
wr0p   = zeros(nD,nL0);
wr0p2  = zeros(nD,nL0);

for iD=1:nD
    D = Darray(iD);
    for i=1:nL0
        L0p = L0+dL0(i);
        r0p = @(x,y) r0.*(L0/L0p).* ...
            ((k2 - (2*pi.*hypot(x,y)./L0).*besselk(5./6,2*pi.*hypot(x,y)./L0) ) ...
            ./(k2 - (2*pi.*hypot(x,y)./L0p).*besselk(5./6,2*pi.*hypot(x,y)./L0p) ) ).^(3./5);
        
        wr0p(iD,i) = integral2(r0p,-D/2,D/2,-D/2,D/2)./D^2;
        
        f1 = @(x,y) (k2 - (2*pi.*hypot(x,y)./L0).*besselk(5./6,2*pi.*hypot(x,y)./L0) ).^(3./5);
        f2 = @(x,y) (k2 - (2*pi.*hypot(x,y)./L0p).*besselk(5./6,2*pi.*hypot(x,y)./L0p) ).^(3./5);
        wr0p2(iD,i) = r0.*(L0/L0p).*integral2(f1,-D/2,D/2,-D/2,D/2)./integral2(f2,-D/2,D/2,-D/2,D/2);
    end
end

figure;
xlabel('Signed relative error on L_0 [%]');
ylabel('Signed relative error on r_0 [%]');
hold on;
plot( 100.*dL0./L0,100.*(r0-wr0p(1,:))./r0);
plot( 100.*dL0./L0,100.*(r0-wr0p(2,:))./r0);
plot( 100.*dL0./L0,100.*(r0-wr0p(3,:))./r0);
legend('4 m telescope','10 m telescope','40 m telescope');
grid on;

figure;
xlabel('Signed relative error on L_0 [%]');
ylabel('Signed relative error on r_0 [%]');
hold on;
plot( 100.*dL0./L0,100.*(r0-wr0p2(1,:))./r0);
plot( 100.*dL0./L0,100.*(r0-wr0p2(2,:))./r0);
plot( 100.*dL0./L0,100.*(r0-wr0p2(3,:))./r0);
legend('4 m telescope','10 m telescope','40 m telescope');
grid on;

% FUN = @(x,xdata) x(1).*log((xdata + x(2))) + x(3);
% X1 = lsqcurvefit(FUN,[1,1,1],100.*dL0./L0,100.*(r0-wr0p(1,:))./r0);
% M1 = FUN(X1,100.*dL0./L0);
% X2 = lsqcurvefit(FUN,[1,1,1],100.*dL0./L0,100.*(r0-wr0p(2,:))./r0);
% M2 = FUN(X2,100.*dL0./L0);
% X3 = lsqcurvefit(FUN,[1,1,1],100.*dL0./L0,100.*(r0-wr0p(3,:))./r0);
% M3 = FUN(X3,100.*dL0./L0);