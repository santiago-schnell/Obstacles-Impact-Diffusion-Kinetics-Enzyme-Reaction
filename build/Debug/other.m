% Evaluate k1
%x=1:min_size-2;
%y=res(x,1);
%p=polyfit(x',y,1);
%yfit=polyval(p,x)';
%display(p);
%yresid=y-yfit;
%SSresid=sum(yresid.^2);
%SStotal=(length(y)-1)*var(y);
%rsq_adj=1-SSresid/SStotal*(length(y)-1)/(length(y)-length(p)-1);
%display(rsq_adj);

% Produces figure 4
    function fig4
        cell_obst={'0.10','0.20','0.30','0.40'};
        cell_tam={'1','16','49','100'};
        cell_spin={'0','0.20','0.40','0.60','0.80','1'};
        obst=str2double(cell_obst)*100;
        tam=str2double(cell_tam);
        spin=str2double(cell_spin)*100;
        f_rate=100;
        wfi=0;
        
        avg_tt=zeros(length(spin)-1,length(tam),4);
        count=length(spin)*length(tam)*length(obst);
        H0=[];
        for i=1:length(spin)
            for j=1:length(tam)
                for k=1:length(obst)
                    [avgs,~]=calc_params(folder,nsimuls,niterations,obst(k),tam(j),spin(i),f_rate,wfi,10);
                    if (i==1)
                        H0=[H0 avgs(4)];
                    else
                        avg_tt(i-1,j,k)=avgs(4);
                    end
                    count=count-1
                end
            end
        end
        
        figure(4);
        line_type={'ks-','ko-','ko-','ks-'};
        line_color={'w','k','w','k'};
        for i=1:4
            subplot(2,2,i);
            vec=i:length(obst):length(tam)*length(obst);
            norm_h=repmat(H0(vec),length(spin)-1,1);
            for j=1:4
                hold on;
                plot(avg_tt(:,j,i)./norm_h(:,j),line_type{j},'MarkerFaceColor',line_color{j});
            end
            xlabel('Spin','FontWeight','bold','FontSize',16);
            ylabel('h/h0','FontWeight','bold','FontSize',16);
            str=sprintf('\\theta : %s',cell_obst{i});
            title(str,'FontWeight','bold','FontSize',16);
            xlim([1 5]);
            set(gca,'XTickLabel',cell_spin(2:6),'FontWeight','bold','FontSize',16);
            legend(cell_tam{1},cell_tam{2},cell_tam{3},cell_tam{4});
            ylim([0 1]);
            grid on;
        end
    end