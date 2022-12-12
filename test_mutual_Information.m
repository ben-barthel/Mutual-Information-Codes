clear all
close all
clc
%% Inputs

%rho = 0.9*j/100;
a = 4; b = 25;
n = 5*10^4
% NBV = [0,3,4,5,6,7,8,9,10];
NBV = [4,5,6,7]
for jb = 1:length(NBV)
    if NBV == 0
        nbins = 0;
    else
        nbins = 2^NBV(jb);
    end
    tic
    for j = 1:20

        %Testing Mutual Information Function
        rho = 0.9*(j/20);
        Sigma = [a^2,rho*a*b;rho*a*b,b^2];mu =[0,0];
        R = mvnrnd(mu,Sigma,n);
        q1 = R(:,1)';
        q2 = R(:,2)';
        Rt = cov(q1,q2);
        rho2 = Rt(1,2)/sqrt(Rt(1,1)*Rt(2,2));

        test1(j) = -.5*log(1-rho^2);
        test2(j) = -.5*log(1-rho2^2);
        testH12(j) = (2/2)*(1+log(2*pi)) +(1/2)*log(det(Sigma));
        testH1(j) = 1/2 + (1/2)*log(2*pi*a^2);
        testH2(j) = 1/2 + (1/2)*log(2*pi*b^2);

        H1T = testH12(j);

        I(j) = mutualInfo_MC(q1,q2,nbins);
        Ikde(j) = mutualInfo_KDE(q1,q2,nbins);

    end
    toc



    disp('Mutual Information Computation Complete')






    %% Plot


    jv = linspace(1,100,100);
    figure(111);
    subplot(2,2,jb)
    hold on;
    test3 = testH1 + testH2 -testH12;
    plot(jv(1:length(I)),I,'k','linewidth',3)
    plot(jv(1:length(I)),Ikde,'ok','linewidth',3)
    plot(jv(1:length(I)),test1,'<r','linewidth',6)
    plot(jv(1:length(I)),test2,'og','linewidth',3)
    plot(jv(1:length(I)),test3,'xb','linewidth',4)
    legend('MC','KDE','Truth 1.0','Truth 2.0','Truth 3.0')
    set(gca,'fontsize',20,'fontname','times')
    xlabel('$s$','interpreter','latex')
    ylabel('$I$','interpreter','latex')
    if nbins == 0
        title('$Optimal N_{bins}$','Interpreter','latex')
    else
        title(['$N_{bins} = $',num2str(nbins) ],'Interpreter','latex')
    end
    figure(1111);
    subplot(2,2,jb)
    hold on;
    test3 = testH1 + testH2 -testH12;
        plot(jv(1:length(I)),test1,'or','linewidth',6)
   
    plot(jv(1:length(I)),I,'k','linewidth',3)
    plot(jv(1:length(I)),Ikde,'<k','linewidth',3)

    legend('Truth','MC','KDE')
    set(gca,'fontsize',20,'fontname','times')
    xlabel('$s$','interpreter','latex')
    ylabel('$I$','interpreter','latex')
    if nbins == 0
        title('$Optimal N_{bins}$','Interpreter','latex')
    else
        title(['$N_{bins} = $',num2str(nbins) ],'Interpreter','latex')
    end
end