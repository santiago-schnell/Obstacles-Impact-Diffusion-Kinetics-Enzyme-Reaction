function test_files(folder,nsimuls)

cell_obst={'0.1','0.2','0.3','0.4'};
%cell_obst={'0'};
cell_tam={'100'};
%cell_tam={'1'};
cell_conc_S={'0.02','0.03','0.04','0.05','0.06','0.07','0.08','0.09','0.10','0.11','0.12','0.13','0.14','0.15','0.16','0.17','0.18','0.19','0.20'};
%cell_conc_S={'0.55'};

obst=str2double(cell_obst)*100;
tam=str2double(cell_tam);
conc_S=str2double(cell_conc_S)*100;

count=length(obst)*length(tam)*length(conc_S);
for i=1:length(obst)
    for j=1:length(tam)
        for k=1:length(conc_S)
            for l=1:nsimuls
                try
                    tt_matrix1=load(sprintf('%s/Results1-%d-%d-%d-%i.txt',folder,obst(i),tam(j),fix(conc_S(k)),l));
                    %tt_matrix2=load(sprintf('%s/Results2-%d-%d-%d-%i.txt',folder,obst(i),tam(j),conc_S(k),l));
                catch
                    msg=sprintf('ERROR1: %d,%d,%d,%d\n',obst(i),tam(j),fix(conc_S(k)),l);
                    display(msg);
                    continue;
                end
                sz1=size(tt_matrix1);
                %sz2=size(tt_matrix2);
                if ((sz1(1)==0) || (sz1(2)==0))% || (sz2(1)==0) || (sz2(2)==0))
                    msg=sprintf('ERROR2: %d,%d,%d,%d\n',obst(i),tam(j),fix(conc_S(k)),l);
                    display(msg);
                end
                clear tt_matrix1;
                clear tt_matrix2;
            end
            count=count-1
        end
    end
end

end

