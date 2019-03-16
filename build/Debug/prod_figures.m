function prod_figures(folder,nsimuls,niterations,figs)

for tt_fig=1:8
    if (sum(figs==tt_fig)>0)
        fn=strcat('fig',num2str(tt_fig));
        eval(fn);
        break;
    end
end

    % Produces figure 1
    function fig1 %#ok<DEFNU>
        cell_obst={'0','0.10','0.20','0.30','0.40'};
        cell_tam={'1','16','49','100'};
        obst=str2double(cell_obst)*100;
        tam=str2double(cell_tam);
        f_rate=100;
        spin=0;
        wfi=0;
        
        avg_tt1=zeros(length(obst),length(tam),4);
        std_tt1=zeros(length(obst),length(tam),4);
        avg_tt2=zeros(length(obst),length(tam),4);
        std_tt2=zeros(length(obst),length(tam),4);
        for i=1:length(obst)
            for j=1:length(tam)
                [avgs1,stds1]=calc_params(folder,nsimuls,niterations,obst(i),tam(j),spin,f_rate,wfi,[]);
                [avgs2,stds2]=calc_params(folder,nsimuls,niterations,obst(i),tam(j),spin,f_rate,wfi,10);
                avg_tt1(i,j,:)=avgs1(1:4);
                std_tt1(i,j,:)=stds1(1:4);
                avg_tt2(i,j,:)=avgs2(1:4);
                std_tt2(i,j,:)=stds2(1:4);
            end
        end
        
        figure(1);
        cell_type={'\alpha','D*','t_{CR}','h'};
        line_type={'ks-','ko-','ko-','ks-'};
        line_color={'w','k','w','k'};
        for i=1:4
            subplot(2,2,i);
            if (i==4)
                for j=1:4
                    hold on;
                    plot(avg_tt2(:,j,i),line_type{j},'MarkerFaceColor',line_color{j});
                end
            else
                for j=1:4
                    hold on;
                    errorbar(avg_tt1(:,j,i),std_tt1(:,j,i),line_type{j},'MarkerFaceColor',line_color{j});
                end
            end
            xlabel('\theta','FontWeight','bold','FontSize',16);
            ylabel(cell_type{i},'FontWeight','bold','FontSize',16);
            legend(cell_tam{1},cell_tam{2},cell_tam{3},cell_tam{4});
            xlim([0.5 5.5]);
            set(gca,'XTickLabel',cell_obst,'FontWeight','bold','FontSize',16);
            grid on;
        end
    end

    % Produces figure 2
    function fig2 %#ok<DEFNU>
        cell_obst={'0.10','0.20','0.30','0.40'};
        cell_tam={'1','16','49','100'};
        obst=str2double(cell_obst)*100;
        tam=str2double(cell_tam);
        f_rate=100;
        spin=0;
        wfi=0;
        
        avg_tt=zeros(length(obst),length(tam));
        std_tt=zeros(length(obst),length(tam));
        for i=1:length(obst)
            for j=1:length(tam)
                [avgs1,stds1]=calc_params(folder,nsimuls,niterations,obst(i),tam(j),spin,f_rate,wfi,[]);
                avg_tt(i,j)=avgs1(5);
                std_tt(i,j)=stds1(5);
            end
        end
        
        figure(2);
        line_type={'ks-','ko-','ko-','ks-'};
        line_color={'w','k','w','k'};
        for j=1:4
            hold on;
            errorbar(avg_tt(:,j),std_tt(:,j),line_type{j},'MarkerFaceColor',line_color{j});
        end
        xlabel('\theta','FontWeight','bold','FontSize',16);
        ylabel('D*(\theta) / D*','FontWeight','bold','FontSize',16);
        legend(cell_tam{1},cell_tam{2},cell_tam{3},cell_tam{4});
        xlim([0.5 4.5]);
        ylim([0 1.1]);
        set(gca,'XTickLabel',cell_obst,'FontWeight','bold','FontSize',16);
        grid on;
    end

    % Produces figure 3
    function fig3 %#ok<DEFNU>
        cell_obst={'0'};
        cell_tam={'1'};
        cell_spin={'0','0.20','0.40','0.60','0.80','1'};
        obst=str2double(cell_obst)*100;
        tam=str2double(cell_tam);
        spin=str2double(cell_spin)*100;
        f_rate=100;
        wfi=0;
        
        avg_tt=zeros(1,length(spin));
        for i=1:length(spin)
            [avgs,~]=calc_params(folder,nsimuls,niterations,obst(1),tam(1),spin(i),f_rate,wfi,10);
            avg_tt(i)=avgs(4);
        end
        
        figure(3);
        plot(avg_tt,'color','k');
        xlabel('Spin','FontWeight','bold','FontSize',16);
        ylabel('h','FontWeight','bold','FontSize',16);
        xlim([1 6]);
        set(gca,'XTickLabel',cell_spin(1:6),'FontWeight','bold','FontSize',16);
        grid on;
    end

    % Produces figure 4
    function fig4 %#ok<DEFNU>
        cell_obst={'0.10','0.20','0.30','0.40'};
        cell_tam={'1','16','49','100'};
        cell_spin={'0','0.20','0.40','0.60','0.80','1'};
        obst=str2double(cell_obst)*100;
        tam=str2double(cell_tam);
        spin=str2double(cell_spin)*100;
        f_rate=100;
        wfi=0;
        
        avg_tt=zeros(length(spin),length(tam));
        count=length(spin)*length(tam)*length(obst);
        for i=1:length(spin)
            for j=1:length(tam)
                for k=1:length(obst)
                    [avgs,~]=calc_params(folder,nsimuls,niterations,obst(k),tam(j),spin(i),f_rate,wfi,10);                    
                    avg_tt(i,j,k)=avgs(4);
                    count=count-1
                end
            end
        end
        
        figure(4);
        line_type={'ks-','ko-','ko-','ks-'};
        line_color={'w','k','w','k'};
        for i=1:4
            subplot(2,2,i);
            for j=1:4
                hold on;
                plot(avg_tt(:,j,i),line_type{j},'MarkerFaceColor',line_color{j});
            end
            xlabel('Spin','FontWeight','bold','FontSize',16);
            ylabel('h','FontWeight','bold','FontSize',16);
            str=sprintf('\\theta : %s',cell_obst{i});
            title(str,'FontWeight','bold','FontSize',16);
            xlim([1 6]);
            set(gca,'XTickLabel',cell_spin(1:6),'FontWeight','bold','FontSize',16);
            legend(cell_tam{1},cell_tam{2},cell_tam{3},cell_tam{4});
            grid on;
        end
    end

    % Produces figure 5
    function fig5 %#ok<DEFNU>
        cell_obst={'0','0.20','0.40'};
        cell_tam={'1','16','49','100'};
        cell_f_rate={'0.20','0.40','0.60','0.80','1'};
        obst=str2double(cell_obst)*100;
        tam=str2double(cell_tam);
        f_rate=str2double(cell_f_rate)*100;
        spin=0;
        wfi=0;
        
        avg_tt=zeros(length(f_rate),length(tam));
        count=length(f_rate)*length(tam)*length(obst);
        for i=1:length(f_rate)
            for j=1:length(tam)
                for k=1:length(obst)
                    [avgs,~]=calc_params(folder,nsimuls,niterations,obst(k),tam(j),spin,f_rate(i),wfi,10);
                    avg_tt(i,j,k)=avgs(4);
                    count=count-1
                end
            end
        end
        
        figure(5);
        line_type={'ks-','ko-','ko-','ks-'};
        line_color={'w','k','w','k'};
        for i=1:3
            subplot(3,1,i);
            for j=1:4
                hold on;
                plot(avg_tt(:,j,i),line_type{j},'MarkerFaceColor',line_color{j});
            end
            xlabel('f','FontWeight','bold','FontSize',16);
            ylabel('h','FontWeight','bold','FontSize',16);
            str=sprintf('\\theta : %s',cell_obst{i});
            title(str,'FontWeight','bold','FontSize',16);
            xlim([1 5]);
            set(gca,'XTickLabel',cell_f_rate(1:5),'FontWeight','bold','FontSize',16);
            legend(cell_tam{1},cell_tam{2},cell_tam{3},cell_tam{4});
            grid on;
        end
    end

    % Produces figure 6
    function fig6 %#ok<DEFNU>
        cell_obst={'0','0.10','0.20','0.30','0.40'};
        cell_tam={'1'};%,'16','49','100'};
        obst=str2double(cell_obst)*100;
        tam=str2double(cell_tam);
        
        avg_tt=zeros(length(obst),length(tam),3);
        std_tt=zeros(length(obst),length(tam),3);
        count=length(obst)*length(tam);
        for i=1:length(obst)
            for j=1:length(tam)
                [avgs1,stds1]=calc_params(folder,nsimuls,niterations,obst(i),tam(j),20,10);
                avg_tt(i,j,:)=avgs1(6:8);
                std_tt(i,j,:)=stds1(6:8);
                count=count-1
            end
        end
        
        figure(6);
        line_type={'ks-','ko-','ko-','ks-'};
        line_color={'w','k','w','k'};
        var={'k1','k2','k3'};
        for i=1:3
            subplot(3,1,i);
            for j=1:1
                hold on;
                errorbar(avg_tt(:,j,i),std_tt(:,j,i),line_type{j},'MarkerFaceColor',line_color{j});
            end
            xlabel('\theta','FontWeight','bold','FontSize',16);
            ylabel(var(i),'FontWeight','bold','FontSize',16);
            legend(cell_tam{1});%,cell_tam{2},cell_tam{3},cell_tam{4});
            xlim([0.5 5.5]);
            set(gca,'XTickLabel',cell_obst,'FontWeight','bold','FontSize',16);
            grid on;
        end
        
        figure(10);
        subplot(1,2,1);
        plot(1:length(obst),avg_tt(:,1,3).*0.01,'kx');
        xlim([0.5 5.5]);
        set(gca,'XTickLabel',cell_obst,'FontWeight','bold','FontSize',16);
        grid on;
        subplot(1,2,2);
        plot(1:length(obst),(avg_tt(:,1,2)+avg_tt(:,1,3))./avg_tt(:,1,1),'kx');
        xlim([0.5 5.5]);
        set(gca,'XTickLabel',cell_obst,'FontWeight','bold','FontSize',16);
        grid on;
    end

    % Produces figure 7
    function fig7 %#ok<DEFNU>
        cell_wfi={'0','0.2','0.4','0.6','0.8','1'};
        cell_wf_rate={'0.01','0.01','0.1'};
        cell_wr_rate={'0.1','0.01','0.01'};
        wfi=str2double(cell_wfi)*100;
        wf_rate=str2double(cell_wf_rate)*100;
        wr_rate=str2double(cell_wr_rate)*100;
        
        avg_tt=zeros(length(wfi),length(wf_rate));
        for i=1:length(wfi)
            for j=1:length(wf_rate)
                [avgs,~]=calc_params(folder,nsimuls,niterations,40,1,0,100,wfi(i),wf_rate(j),wr_rate(j),10);
                avg_tt(i,j)=avgs(4);
            end
        end
        
        figure(7);
        subplot(2,1,2);
        line_type={'ks-','ko-','kd-'};
        line_color={'w','k','w'};
        for j=1:3
            hold on;
            plot(avg_tt(:,j),line_type{j},'MarkerFaceColor',line_color{j});
        end
        xlabel('wfi','FontWeight','bold','FontSize',16);
        ylabel('h','FontWeight','bold','FontSize',16);
        xlim([1 6]);
        set(gca,'XTickLabel',cell_wfi(1:6),'FontWeight','bold','FontSize',16);
        legend('wf/wr : 0.01 / 0.1','wf/wr : 0.01 / 0.01','wf/wr : 0.1 / 0.01');
        grid on;
    end

    % Produces figure 8
    function fig8 %#ok<DEFNU>
        cell_obst={'0.10','0.20','0.30','0.40'};
        cell_tam={'1'};
        cell_conc_S={'0.01','0.02','0.03','0.04','0.05','0.06','0.07','0.08','0.09','0.10','0.11,','0.12','0.13','0.14','0.15','0.16','0.17','0.18','0.19','0.20'};
        %,'0.21','0.22','0.23','0.24','0.25','0.26','0.27','0.28','0.29','0.30'};
        obst=str2double(cell_obst)*100;
        tam=str2double(cell_tam);
        conc_S=str2double(cell_conc_S)*100;
        
        avg_tt=zeros(length(obst),length(tam),length(conc_S));
        std_tt=zeros(length(obst),length(tam),length(conc_S));
        count=length(obst)*length(tam)*length(conc_S);
        for i=1:length(obst)
            for j=1:length(tam)
                for k=1:length(conc_S)
                    %obst(i)
                    %tam(j)
                    %fix(conc_S(k))
                    [avgs1,stds1]=calc_params(folder,nsimuls,niterations,obst(i),tam(j),fix(conc_S(k)),10);
                    avg_tt(i,j,k)=avgs1;
                    std_tt(i,j,k)=stds1;
                    count=count-1
                end
            end
        end
        
        % Calculate Vmax and Km
        count=1;
        CI_L=zeros(length(obst),length(tam),2);
        CI_U=zeros(length(obst),length(tam),2);
        final_res=zeros(length(obst),length(tam),2);
        line_type={'ks','ko','ko','ks'};
        line_color={'w','k','w','k'};
        %color=['k','g','b','m','y'];
        figure(15);
        %subplot(3,1,3);
        for i=1:length(obst)
            for j=1:length(tam)
                x=str2double(cell_conc_S);
                y(:,1)=avg_tt(i,j,:);
                z(:,1)=std_tt(i,j,:)./sqrt(nsimuls); % Standard Error
                hold on;
                %plot(x,y','o','color',color(count));
                errorbar(x,y',z',line_type{i},'MarkerFaceColor',line_color{i});
                [MM,R,~,S]=nlinfit(x',y,@MM_model,[0.00030 0.04]); % Obtains parameters from MM fit
                CI=nlparci(MM,R,'covar',S);
                CI_L(i,j,1)=abs(MM(1)-CI(1,1));
                CI_U(i,j,1)=abs(MM(1)-CI(1,2));
                CI_L(i,j,2)=abs(MM(2)-CI(2,1));
                CI_U(i,j,2)=abs(MM(2)-CI(2,2));
                yfit=(MM(1)*x)./(MM(2)+x);
                C=corrcoef(y,yfit);
                rsq1=C(1,2)^2
                rsq2=1-sum(R.^2)/sum((y-mean(y)).^2)
                hold on;
                plot(x,yfit,'color','k');
                grid on;
                xlabel('[S]','FontWeight','bold','FontSize',16);
                ylabel('v','FontWeight','bold','FontSize',16);
                %title('k_3 : 0.50; k_{-3} : 0.01','FontWeight','bold','FontSize',16);
                set(gca,'FontWeight','bold','FontSize',16);
                final_res(i,j,1)=MM(1);
                final_res(i,j,2)=MM(2);
                count=count+1;
                xlim([0.01 0.2]);
            end
        end
        
        line_type={'ks-','ko-','ko-','ks-'};
        var={'Vmax','Km'};
        for i=1:2
            for j=1:1
                figure(9);
                subplot(1,2,i);
                hold on;
                final_res(:,j,i)
                plot(1:length(obst),final_res(:,j,i),line_type{2},'MarkerFaceColor',line_color{2});
                %figure(9+i);
                %subplot(2,2,j);
                %hold on;
                %errorbar(1:length(obst),final_res(:,j,i),CI_L(:,j,i),CI_U(:,j,i),line_type{j},'MarkerFaceColor',line_color{j});
                %xlabel('\theta','FontWeight','bold','FontSize',16);
                %ylabel(var(i),'FontWeight','bold','FontSize',16);
                %legend(cell_tam{1},cell_tam{2},cell_tam{3},cell_tam{4});
                %legend(cell_tam{j});
                %xlim([0.5 4.5]);
                %set(gca,'XTickLabel',cell_obst,'FontWeight','bold','FontSize',16);
                %grid on;
            end
            figure(9);
            xlabel('\theta','FontWeight','bold','FontSize',16);
            ylabel(var(i),'FontWeight','bold','FontSize',16);
            %legend(cell_tam{1},cell_tam{2},cell_tam{3},cell_tam{4});
            xlim([0.5 4.5]);
            set(gca,'XTickLabel',cell_obst,'FontWeight','bold','FontSize',16);
            grid on;
        end
    end

    % This function extracts parameters from the Michaelis Menten model
    function yhat=MM_model(beta,S)
          yhat=(beta(1).*S)./(beta(2)+S);
    end

end