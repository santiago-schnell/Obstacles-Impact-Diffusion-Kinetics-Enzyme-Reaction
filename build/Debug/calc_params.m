function [avgs,stds]=calc_params(folder,nsimuls,niterations,obst,tam,conc_S,plts)

if isempty(plts)
    params=zeros(nsimuls,5);
    for i=1:nsimuls
        tt_matrix1=load(sprintf('%s/Results1-%d-%d-%d-%i.txt',folder,obst,tam,conc_S,i));
        tt_matrix2=load(sprintf('%s/Results2-%d-%d-%d-%i.txt',folder,obst,tam,conc_S,i));
        sz=size(tt_matrix1);
        [A,DS,C,H,DO,~,~,~,~,~,~]=calc_vars(tt_matrix1,tt_matrix2,sz(1));
        params(i,:)=[A DS C H DO/DS];
        clear tt_matrix1;
        clear tt_matrix2;
    end
    avgs=mean(params,1);
    stds=std(params,0,1);
else
    % Get minimum number of iterations from the simulations
    min_size=niterations;
    for i=1:nsimuls
        matrix1=load(sprintf('%s/Results1-%d-%d-%d-%i.txt',folder,obst,tam,conc_S,i));
        [x,~]=size(matrix1);
        if x<min_size
            min_size=x;
        end
    end
    
    %min_size
    
    % Get averages
    tt_matrix1=zeros(min_size,10);
    %tt_matrix2=zeros(min_size,19);
    slopes=zeros(1,nsimuls);
    for i=1:nsimuls
        matrix1=load(sprintf('%s/Results1-%d-%d-%d-%i.txt',folder,obst,tam,conc_S,i));
        %matrix2=load(sprintf('%s/Results2-%d-%d-%d-%i.txt',folder,obst,tam,conc_S,i));
        %size(matrix1)
   
        p=polyfit((1:30)',matrix1(1:30,5),1); % Obtain slope of initial product formation
        slopes(i)=p(1);
        
        tt_matrix1=tt_matrix1+matrix1(1:min_size,:);
        %tt_matrix2=tt_matrix2+matrix2(1:min_size,:);
        clear matrix1;
        %clear matrix2;
    end
    tt_matrix1=tt_matrix1./nsimuls;
    %tt_matrix2=tt_matrix2./nsimuls;
    tt_matrix2=[];
    %[A,DS,C,H,DO,Mk1,Sk1,Mk2,Sk2,Mk3,Sk3]=calc_vars(tt_matrix1,tt_matrix2,min_size);
    %Mk3,Sk3
    %avgs=[A DS C H DO/DS Mk1 Mk2 Mk3];
    %stds=[0 0 0 0 0 Sk1 Sk2 Sk3];
    
    %max_time=floor(min_size*0.1); % Maximum number of time steps to obtain slope
    %figure(12);
    %hold on;
    %plot(tt_matrix1(1:max_time,5));
    %pause(2);
    
    %p=polyfit((15:30)',tt_matrix1(15:30,5),1); % Obtain slope of initial product formation
    %IPS=p(1); % Sets IPS
    avgs=mean(slopes);
    stds=std(slopes);
    
    figure(10);
    subplot(1,2,1);
    hold on;
    plot(1:30,tt_matrix1(1:30,5),'k');
    xlabel('t','FontWeight','bold','FontSize',16);
    ylabel('[P]','FontWeight','bold','FontSize',16);
    set(gca,'FontWeight','bold','FontSize',16);
    grid on;
    
    figure(10);
    subplot(1,2,2);
    hold on;
    max=1000;
    if (min_size<max)
        max=min_size;
    end
    plot(1:max,tt_matrix1(1:max,5),'k');
    xlabel('t','FontWeight','bold','FontSize',16);
    ylabel('[P]','FontWeight','bold','FontSize',16);
    set(gca,'FontWeight','bold','FontSize',16);
    grid on;
    
    % Plot concentrations
    if (sum(plts==1)>0)
        color='m';
        figure(1);
        subplot(2,2,1);
        hold on;
        plot(1:min_size,tt_matrix1(1:min_size,2),color);
        title('Enzyme','FontWeight','bold','FontSize',16);
        xlabel('Time','FontWeight','bold','FontSize',16);
        ylabel('Conc','FontWeight','bold','FontSize',16);
        grid on;
        subplot(2,2,2);
        hold on;
        plot(1:min_size,tt_matrix1(1:min_size,3),color);
        title('Substrate','FontWeight','bold','FontSize',16);
        xlabel('Time','FontWeight','bold','FontSize',16);
        ylabel('Conc','FontWeight','bold','FontSize',16);
        grid on;
        subplot(2,2,3);
        hold on;
        plot(1:min_size,tt_matrix1(1:min_size,4),color);
        title('Complex','FontWeight','bold','FontSize',16);
        xlabel('Time','FontWeight','bold','FontSize',16);
        ylabel('Conc','FontWeight','bold','FontSize',16);
        grid on;
        subplot(2,2,4);
        hold on;
        plot(1:min_size,tt_matrix1(1:min_size,5),color);
        title('Product','FontWeight','bold','FontSize',16);
        xlabel('Time','FontWeight','bold','FontSize',16);
        ylabel('Conc','FontWeight','bold','FontSize',16);
        grid on;
        xlim([0 25]);
    end
    
    % Plot diffusions
    if (sum(plts==2)>0)
        figure(2);
        subplot(2,2,1);
        hold on;
        plot(log(1:min_size),log(diffusions(:,1)),'b');
        title('Enzyme','FontWeight','bold','FontSize',16);
        xlabel('Log(t)','FontWeight','bold','FontSize',16);
        ylabel('Log(<r^2> / t)','FontWeight','bold','FontSize',16);
        axis([0 8 -1 0.2]);
        grid on;
        subplot(2,2,2);
        hold on;
        plot(log(1:min_size),log(diffusions(:,2)),'b');
        yfit=polyval(p,log(1:20)');
        hold on;
        plot(log(1:20),yfit,'r--');
        hold on;
        plot(log(1:min_size-5+1),ones(1,min_size-5+1)*log(DS),'r--');
        title('Substrate','FontWeight','bold','FontSize',16);
        xlabel('Log(t)','FontWeight','bold','FontSize',16);
        ylabel('Log(<r^2> / t)','FontWeight','bold','FontSize',16);
        axis([0 8 -1 0.2]);
        grid on;
        subplot(2,2,3);
        hold on;
        plot(log(1:min_size),log(diffusions(:,3)),'b');
        title('Complex','FontWeight','bold','FontSize',16);
        xlabel('Log(t)','FontWeight','bold','FontSize',16);
        ylabel('Log(<r^2> / t)','FontWeight','bold','FontSize',16);
        axis([0 8 -1 0.2]);
        grid on;
        subplot(2,2,4);
        hold on;
        plot(log(1:min_size),log(diffusions(:,4)),'b');
        title('Product','FontWeight','bold','FontSize',16);
        xlabel('Log(t)','FontWeight','bold','FontSize',16);
        ylabel('Log(<r^2> / t)','FontWeight','bold','FontSize',16);
        axis([0 8 -1 0.2]);
        grid on;
    end
    
    % Plot mean square displacements over time
    if (sum(plts==3)>0)
        figure(3);
        hold on;
        plot(1:min_size,tt_matrix2(1:min_size,8:11));
        xlabel('t','FontWeight','bold','FontSize',16);
        ylabel('<r^2>','FontWeight','bold','FontSize',16);
        legend('Enzyme','Substrate','Complex','Product');
        grid on;
    end
    
    % Plot k1
    if (sum(plts==4)>0)
        figure(7);
%         subplot(1,2,1);
%         hold on;
%         plot(1:min_size-1,res(1:min_size-1,1));
%         xlabel('t','FontWeight','bold','FontSize',16);
%         ylabel('k1','FontWeight','bold','FontSize',16);
%         set(gca,'FontWeight','bold','FontSize',16);
%         grid on;
        subplot(2,1,1);
        hold on;
        plot(log(1:min_size-1),log(res(1:min_size-1,1)),'color','k');
        %x=log(min_size-floor(min_size/2):min_size-1)';
        %y=log(res(min_size-floor(min_size/2):min_size-1,1));
%         p=polyfit(x,y,1)
%         hold on;
%         yfit=polyval(p,x);
%         plot(x,yfit,'r--');
        xlabel('Log(t)','FontWeight','bold','FontSize',16);
        ylabel('Log(k1)','FontWeight','bold','FontSize',16);
        set(gca,'FontWeight','bold','FontSize',16);
        grid on;
        %yresid=y-yfit;
        %SSresid=sum(yresid.^2);
        %SStotal=(length(y)-1)*var(y);
        %rsq=1-SSresid/SStotal
    end
    
    % Plot k2
    if (sum(plts==5)>0)
        figure(5);
        hold on;
        plot(2:100,res(:,2));
        xlabel('t','FontWeight','bold','FontSize',16);
        ylabel('k_{-1} (reversability rate)','FontWeight','bold','FontSize',16);
        set(gca,'XTickLabel',cell_obst,'FontWeight','bold','FontSize',16);
        axis([2 100 0 0.03]);
        grid on;
    end
    
    % Plot k3
    if (sum(plts==6)>0)
        figure(6);
        hold on;
        plot(2:100,res(:,3));
        xlabel('t','FontWeight','bold','FontSize',16);
        ylabel('k_2 (catalytic rate)','FontWeight','bold','FontSize',16);
        set(gca,'XTickLabel',cell_obst,'FontWeight','bold','FontSize',16);
        axis([2 100 0 0.03]);
        grid on;
    end
   
end

    % This function obtains rates and diffusion parameters
    function [A,DS,C,H,DO,Mk1,Sk1,Mk2,Sk2,Mk3,Sk3]=calc_vars(tt_matrix1,tt_matrix2,min_size)
        % Obtain kinetic rates
        tt_matrix_deriv=diff(tt_matrix1)./repmat(diff(tt_matrix1(:,1)),1,10);
        k1=1./(tt_matrix1(2:min_size,2).*tt_matrix1(2:min_size,3)).*tt_matrix_deriv(:,10);
        k2=1./tt_matrix1(2:min_size,4).*(tt_matrix_deriv(:,10)+tt_matrix_deriv(:,3));
        k3=1./tt_matrix1(2:min_size,4).*(tt_matrix_deriv(:,5));
        res=[k1 k2 k3];
        
%         % Obtains Kopelman parameter
%         x=log(min_size-floor(min_size/2):min_size-1)';
%         y=log(res(min_size-floor(min_size/2):min_size-1,1));
%         p=polyfit(x,y,1);
%         H=-p(1);
%         
%         % Find alpha and diffusion value using the first max_alpha_time steps
%         diffusions=tt_matrix2(1:min_size,8:12)./tt_matrix2(1:min_size,2:6);
%         max_alpha_time=5;
%         p=polyfit(log(1:max_alpha_time)',log(diffusions(1:max_alpha_time,2)),1);
%         A=p(1)+1; % Sets alpha
%         max_diff_time=250;
%         DS=mean(diffusions(min_size-max_diff_time:min_size,2)); % Substrate diffusion
%         DO=mean(diffusions(min_size-max_diff_time:min_size,5)); % Obstacle diffusion
%         
%         % Find crossover time
%         C=exp((log(DS)-p(2))/p(1));
        
        A=0;
        DS=0;
        C=0;
        H=0;
        DO=0;

        % Find mean and std's of reaction rates
        max_time=floor(min_size/10);
        %Mk1=mean(res(min_size-100:min_size-1,1));
        %Sk1=std(res(min_size-100:min_size-1,1));
        Mk1=mean(res(1:max_time,1));
        Sk1=std(res(1:max_time,1));
        Mk2=mean(res(:,2));
        Sk2=std(res(:,2));
        Mk3=mean(res(:,3));
        Sk3=std(res(:,3));
    end

end
