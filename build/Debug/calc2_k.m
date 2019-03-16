function calc2_k(nsimuls,niterations,wfi,obst,tam,spin,folder,plts,color)

%clc;
%close all;

part1();
%part2();
%part3();

    function part1()
        
        % Get minimum number of iterations from the simulations
        min_size=niterations;
        for i=1:nsimuls
            matrix1=load(sprintf('%s/Results1-%d-%d-%d-%d.txt',folder,obst,tam,spin,i));
            [x,~]=size(matrix1);
            if x<min_size
                min_size=x;
            end
        end
        
        % Get averages
        tt_matrix1=zeros(min_size,10);
        tt_matrix2=zeros(min_size,19);
        for i=1:nsimuls
            matrix1=load(sprintf('%s/Results1-%d-%d-%d-%d.txt',folder,obst,tam,spin,i));
            matrix2=load(sprintf('%s/Results2-%d-%d-%d-%d.txt',folder,obst,tam,spin,i));
            tt_matrix1=tt_matrix1+matrix1(1:min_size,:);
            tt_matrix2=tt_matrix2+matrix2(1:min_size,:);
        end
        tt_matrix1=tt_matrix1./nsimuls;
        tt_matrix2=tt_matrix2./nsimuls;
        
        % Obtain kinetic rates
        tt_matrix_deriv=diff(tt_matrix1)./repmat(diff(tt_matrix1(:,1)),1,10);
        k1=1./(tt_matrix1(2:min_size,2).*tt_matrix1(2:min_size,3)).*tt_matrix_deriv(:,10);
        k2=1./tt_matrix1(2:min_size,4).*(tt_matrix_deriv(:,10)+tt_matrix_deriv(:,3));
        k3=1./tt_matrix1(2:min_size,4).*(tt_matrix_deriv(:,5));
        %beta=nlinfit(2:min_size,k1',@kopelman_model,[1 0]) % Obtains parameters from Kopelman fit
        %newp=mean(tt_matrix2(ceil((min_size)/1.2):min_size,16)); % Obtains newp from diffusion
        res=[k1 k2 k3];
      
%         if (wfi)
%             TR1=zeros(min_size-1,1);
%             TR2=zeros(min_size-1,1);
%             TR3=zeros(min_size-1,1);
%             YR1=zeros(min_size-1,4);
%             YR2=zeros(min_size-1,4);
%             YR3=zeros(min_size-1,4);
%         else
%             [TR1,YR1]=solveDiffEquations(1,0,1,min_size-1); % Law of mass action
%             [TR2,YR2]=solveDiffEquations(beta(1),beta(2),1,min_size-1); % Fractal Kinetics
%             [TR3,YR3]=solveDiffEquations(beta(1),beta(2),newp,min_size-1); % Modified Fractal Kinetics
%         end
        
        % Plot concentrations
        if (sum(plts==1)>0)
            %figure(10);
            subplot(2,2,1);
            max=100;
            if (min_size<max)
                max=min_size;
            end
            hold on;
            plot(1:max,tt_matrix1(1:max,2),color);
%             hold on;
%             plot(TR1,YR1(:,1),'k--o');
%             hold on;
%             plot(TR2,YR2(:,1),'k--s','MarkerFaceColor','k');
%             hold on;
%             plot(TR3,YR3(:,1),'k--');
%             legend('MC','CK','FK','Modified FK');
            %title('Enzyme (E)','FontWeight','bold','FontSize',16);
            xlabel('t','FontWeight','bold','FontSize',16);
            ylabel('[E]','FontWeight','bold','FontSize',16);
            set(gca,'FontWeight','bold','FontSize',16);
            %axis([1 60 0 8e-3]);
            %grid on;
            figure(1);
            subplot(2,2,2);
            hold on;
            plot(1:max,tt_matrix1(1:max,3),color);
%             hold on;
%             plot(TR1,YR1(:,2),'k--o');
%             hold on;
%             plot(TR2,YR2(:,2),'k--s','MarkerFaceColor','k');
%             hold on;
%             plot(TR3,YR3(:,2),'k--');
%             legend('MC','CK','FK','Modified FK');
            %title('Substrate (S)','FontWeight','bold','FontSize',16);
            xlabel('t','FontWeight','bold','FontSize',16);
            ylabel('[S]','FontWeight','bold','FontSize',16);
            set(gca,'FontWeight','bold','FontSize',16);
            %axis([1 60 0.18 0.2]);
            %grid on;
            subplot(2,2,3);
            %figure(10);
            %subplot(1,2,1);
            hold on;
            plot(1:max,tt_matrix1(1:max,4),color);
%             hold on;
%             plot(TR1,YR1(:,3),'k--o');
%             hold on;
%             plot(TR2,YR2(:,3),'k--s','MarkerFaceColor','k');
%             hold on;
%             plot(TR3,YR3(:,3),'k--');
%             legend('MC','CK','FK','Modified FK');
            %title('Complex (C)','FontWeight','bold','FontSize',16);
            xlabel('t','FontWeight','bold','FontSize',16);
            ylabel('[C]','FontWeight','bold','FontSize',16);
            set(gca,'FontWeight','bold','FontSize',16);
            %axis([1 100 0 10e-3]);
            %grid on;
            figure(1);
            subplot(2,2,4);
            hold on;
            plot(1:max,tt_matrix1(1:max,5),color);
%             hold on;
%             plot(TR1,YR1(:,4),'k--o');
%             hold on;
%             plot(TR2,YR2(:,4),'k--s','MarkerFaceColor','k');
%             hold on;
%             plot(TR3,YR3(:,4),'k--');
%             legend('MC','CK','FK','Modified FK');
            %title('Product (P)','FontWeight','bold','FontSize',16);
            xlabel('t','FontWeight','bold','FontSize',16);
            ylabel('[P]','FontWeight','bold','FontSize',16);
            set(gca,'FontWeight','bold','FontSize',16);
            %axis([1 100 0 0.04]);
            %grid on;
            
            %figure(10);
            %hold on;
            %plot(1:min_size-1,diff(tt_matrix1(1:min_size,4)),color);
        end
        
        % Plot diffusions
        diffusions=tt_matrix2(1:min_size,8:11)./tt_matrix2(1:min_size,2:5);
        
        if (sum(plts==2)>0)
            figure(10);
            %subplot(1,2,1);
            hold on;
            max=100;
            if (min_size<100)
                max=min_size;
            end    
            plot(log(1:max),log(diffusions(1:max,2)),color);
            %title('Enzyme','FontWeight','bold','FontSize',16);
            xlabel('Log(t)','FontWeight','bold','FontSize',16);
            ylabel('Log(<r^2> / \tau)','FontWeight','bold','FontSize',16);
            axis([-1 log(max) -1 0.5]);
            set(gca,'FontWeight','bold','FontSize',16);
            %legend('\theta : 0','\theta : 0.4; \delta : 1','\theta : 0.4; \delta : 100');
            grid on;
%             subplot(2,2,2);
%             figure(10);
%             subplot(1,2,1);
%             hold on;
%             plot(log(1:min_size),log(diffusions(:,2)),color);
%             %title('Substrate','FontWeight','bold','FontSize',16);
%             xlabel('Log( t )','FontWeight','bold','FontSize',16);
%             ylabel('Log( <r^2> / \tau )','FontWeight','bold','FontSize',16);
%             axis([-1 log(min_size) -1 0.5]);
%             set(gca,'FontWeight','bold','FontSize',16);
%             grid on;
%             legend('\theta - 0%; \delta - 1','\theta - 40%; \delta - 1','\theta - 40%; \delta - 100');
%             subplot(2,2,3);
%             hold on;
%             plot(log(1:min_size),log(diffusions(:,3)),color);
%             title('Complex','FontWeight','bold','FontSize',16);
%             xlabel('Log(t)','FontWeight','bold','FontSize',16);
%             ylabel('Log(<r^2> / t)','FontWeight','bold','FontSize',16);
%             axis([-1 log(min_size) -1 0.5]);
%             set(gca,'FontWeight','bold','FontSize',16);
%             grid on;
%             subplot(2,2,4);
%             hold on;
%             plot(log(1:min_size),log(diffusions(:,4)),color);
%             title('Product','FontWeight','bold','FontSize',16);
%             xlabel('Log(t)','FontWeight','bold','FontSize',16);
%             ylabel('Log(<r^2> / t)','FontWeight','bold','FontSize',16);
%             axis([-1 log(min_size) -1 0.5]);
%             set(gca,'FontWeight','bold','FontSize',16);
%             grid on;
        end
        
        % Plot mean square displacements over time
        if (sum(plts==3)>0)
            figure(3);
            hold on;
            plot(1:min_size,(tt_matrix2(1:min_size,8:11)));
            xlabel('t','FontWeight','bold','FontSize',16);
            ylabel('<r^2>','FontWeight','bold','FontSize',16);
            legend('Enzyme','Substrate','Complex','Product');
            set(gca,'FontWeight','bold','FontSize',16);
            grid on;
        end
        
        % Plot k1
        if (sum(plts==4)>0)
            figure(4);
            %subplot(1,2,1);
            %hold on;
            %plot(1:min_size-1,res(1:min_size-1,1),color);
            %yfit=beta(1).*(2:min_size).^(-beta(2));
            %hold on;
            %plot(2:min_size,yfit,'r--');
            %xlabel('t','FontWeight','bold','FontSize',16);
            %ylabel('k_1','FontWeight','bold','FontSize',16);
            %set(gca,'FontWeight','bold','FontSize',16);
            %grid on;
            max=1000;
            if (min_size<1000)
                max=min_size;
            end
            %figure(7);
            %subplot(2,1,1);
            hold on;
            plot(log(1:max-1),log(res(1:max-1,1)),color);
            %plot(1:max-1,res(1:max-1,1),color);
            %mean(res(max-200:max-1,1))
            %x=log(min_size-floor(min_size/2):min_size-1)';
            %y=log(res(min_size-floor(min_size/2):min_size-1,1));
            %p=polyfit(x,y,1)
            %hold on;
            %yfit=polyval(p,x);
            %plot(x,yfit,'r--');
            xlabel('log(t)','FontWeight','bold','FontSize',16);
            ylabel('log(f)','FontWeight','bold','FontSize',16);
            set(gca,'FontWeight','bold','FontSize',16);
            %legend('\theta : 0','\theta : 0.4; \delta : 1','\theta : 0.4; \delta : 100');
            grid on;
            %yresid=y-yfit;
            %SSresid=sum(yresid.^2);
            %SStotal=(length(y)-1)*var(y);
            %rsq=1-SSresid/SStotal
        end
        
        % Plot k-1
        if (sum(plts==5)>0)
            figure(5);
            subplot(1,2,1);
            hold on;
            plot(2:100,res(2:100,2),color);
            xlabel('t','FontWeight','bold','FontSize',16);
            ylabel('k_{-1} (reversability rate)','FontWeight','bold','FontSize',16);
            set(gca,'FontWeight','bold','FontSize',16);
            axis([2 100 0 0.03]);
            grid on;
        end
        
        % Plot k2
        if (sum(plts==6)>0)
            figure(5);
            subplot(1,2,2);
            hold on;
            plot(2:100,res(2:100,3),color);
            xlabel('t','FontWeight','bold','FontSize',16);
            ylabel('k_2 (catalytic rate)','FontWeight','bold','FontSize',16);
            set(gca,'FontWeight','bold','FontSize',16);
            axis([2 100 0 0.05]);
            grid on;
        end
        
    end

% This function solves differential equations
    function [TR,YR] = solveDiffEquations(k1,h,newp,min_size)
        function dx_dt=solODE(t,x)
            k2=0.02;
            k3=0.04;
            
            dx_dt(1)=-k1*t^(-h)*x(1)*x(2)+(k2+k3)*newp*x(3);
            dx_dt(2)=-k1*t^(-h)*x(1)*x(2)+k2*newp*x(3);
            dx_dt(3)=+k1*t^(-h)*x(1)*x(2)-(k2+k3)*newp*x(3);
            dx_dt(4)=+k3*newp*x(3);
            dx_dt=dx_dt';
        end
        
        initConds=[0.01 0.1 0 0];
        [TR,YR]=ode15s(@solODE,[1 min_size],initConds);
    end

% This function extracts parameters from the Kopelman model
    function yhat=kopelman_model(beta,t)
          yhat=beta(1).*t.^(-(beta(2)));
    end

% This function builds an animation of the simulation
    function part2()
        matrix3=load(sprintf('%s/Results3-%d-%d.txt',folder,obst,1));
        vidObj=VideoWriter('~/Desktop/movie.avi');
        vidObj.FrameRate=8;
        open(vidObj);
        for i=1:max(matrix3(:,1))
            aux=find(matrix3(:,1)==i);
            scatter(matrix3(aux,3),matrix3(aux,4),7,matrix3(aux,2));
            axis([1 100 1 100]);
            currFrame=getframe;
            writeVideo(vidObj,currFrame);
            pause(0.1);
        end
        close(vidObj);
    end

% This function plots distance until reaction occurs
    function part3()
        maximum=[0 0];
        parts=[378 632];
        titles={'Enzyme','Substrate'};
        XX1=zeros(2,parts(1));
        XX2=zeros(2,parts(2));
        for i=1:nsimuls
            matrix4=load(sprintf('%s/Results4-%d-%d.txt','MM_obst0_mob0_size0',0,i));
            matrix5=load(sprintf('%s/Results4-%d-%d.txt','MM_obst40_mob0_size1',40,i));
            for j=0:1
                aux1=matrix4(matrix4(:,1)==j,2);
                aux2=matrix5(matrix5(:,1)==j,2);
                aux3=max([aux1;aux2]); % Obtains maximum distance
                if (aux3>maximum(j+1))
                    maximum(j+1)=aux3;
                end
                [X1,~]=hist(aux1,parts(j+1));
                [X2,~]=hist(aux2,parts(j+1));
                if (j==0)
                    XX1(1,:)=XX1(1,:)+X1;
                    XX1(2,:)=XX1(2,:)+X2;
                else
                    XX2(1,:)=XX2(1,:)+X1;
                    XX2(2,:)=XX2(2,:)+X2;
                end
            end
        end
        %maximum
        % Obtain averages
        XX1=XX1./nsimuls;
        XX2=XX2./nsimuls;
        aux={'XX1','XX2'};
        for j=0:1
            res=eval(aux{j+1});
            figure(7);
            subplot(1,2,j+1);
            bar(res(1,:)./sum(res(1,:)),'b');
            hold on;
            bar(res(2,:)./sum(res(2,:)),'r');
            axis([1 100 0 0.15]);
            title(titles{j+1},'FontWeight','bold','FontSize',16);
            legend('0% Obstacles','40% Obstacles');
            xlabel('Distance','FontWeight','bold','FontSize',16);
            ylabel('Fraction of reactions','FontWeight','bold','FontSize',16);
            grid on;
        end
    end

end


