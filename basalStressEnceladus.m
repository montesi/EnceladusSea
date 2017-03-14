% BasalStress.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Similar to stressProfile but for the base of the shell (Figure 9)
% Doesn't loop: use currently loaded model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Uses LiveLink to matlab and Comsol Multiphysics 5.2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Used in Johnston and Montesi, Journal of Geophysical Research, 2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for ThicknessShell=20:20:80;
    ModelAll=[1:2];%[1:4]; %note: boundary numbering is off for model 4
    Lb=2e3; % distance for the sea edge for stress determination
    ThickAll=ThicknessShell/2; %sea thickness to consider
    nT=numel(ThickAll);
    for ModelType=ModelAll
        %Produce label
        switch ModelType
            case 1; %No slip boundary between core and ice shell
                LabelType='Fixed';
            case 2; %Free slip boundary between core and ice shell
                LabelType='Roller';
            case 3; %Constant pressure ocean
                LabelType='Ocean';
            case 4; %Constant pressure ocean with indentations at both poles
                LabelType='North';
        end
        Label=sprintf('EnceladusT%g%s',ThicknessShell,LabelType);
        
        figure(1); clf; hold on; set(gca,'fontsize',12); box on;
        
        xlabel('Sea angle (\circ)','fontsize',18);
        ylabel('Vertical stress (kPa)','fontsize',18);
        title(sprintf('T%g%s',ThicknessShell,LabelType));
        
        for iT=1:nT;
            D=ThickAll(iT);
            Fname=sprintf('EnceladusT%gD%g%s',ThicknessShell,D,LabelType)
            
            model=mphload(sprintf('T%g/%s',ThicknessShell,Fname));
            
            % Extract Model parameters
            Rsurf=model.param().evaluate('Rsurface');
            Rcore=model.param().evaluate('Rcore');
            dA=Lb/Rcore; %angular distance for sampling point
            
            % what models do we have?
            dataset='dset2';
            INFO2=mphsolinfo(model,'dataset',dataset);
            nmod=numel(INFO2.batch.sol);
            
            figure(2); clf;
            hold on; set(gca,'fontsize',12,'linewidth',1); box on;
            Slift=NaN(nmod,1);
            Aall=NaN(nmod,1);
            Akeep=[];ik=0;
            for is=2:1:nmod
                Inow=mphsolinfo(model,'soltag',INFO2.batch.sol(is));
                Runall(is).info=Inow;
                Runall(is).OceanAngle=Inow.paramsweepvals;
                As=Runall(is).OceanAngle;
                Aall(is)=As*180/pi;
                
                %% extract data at all boundaries
                data=mpheval(model,{'solid.sr','solid.sz','solid.srz','TH'},...
                    'dataset','dset2','outersolnum',is,...
                    'edim','boundary');
                %%
                xp=data.p(1,:);yp=data.p(2,:);
                Rp=(xp.^2+yp.^2).^(1/2);
                Ip=find(abs(Rp-Rcore)<10); 
                Tp=acos(-yp(Ip)./Rp(Ip));
                [Ts,Ic]=sort(Tp);Is=Ip(Ic);
                % rotate stress to vertical
                Tnn=(data.d1(Is)+data.d2(Is))/2-(data.d1(Is)-data.d2(Is))/2.*cos(2*Ts)...
                    -data.d3(Is).*sin(2*Ts);
                
                iu=min(find(diff(Ts(2:end-1))==0)); 
                if isempty(iu); iu=numel(Ts)-2; end                    
                if and(Ts(iu+1)>=As+dA,iu>=3);
                Slift(is)=interp1(Ts(2:iu+1),Tnn(2:iu+1),As+dA);
                else
                    Slift(is)=NaN;
                end
                if mod(is-1,5)==0
                    plot((Ts(2:end-1)-As)*180/pi,Tnn(2:end-1)/1e3,'linewidth',1);
                    ik=ik+1;
                    Akeep(ik)=As*180/pi;
                end
                %%
                
            end
            % plot profile of stress (not used in paper)
            xlabel('Angular distance from the sea edge (\circ)','fontsize',18);
            ylabel('Vertical stress (kPa)','fontsize',18);
            title(sprintf('T%gD%g%s',ThicknessShell,D,LabelType));
            axis([0,45,-20,30])
            legend(num2str(Akeep'),'location','Northeast');
            plot([0,45],[0,0],'k');
            set(gca,'DataAspectRatio',[2,5,1])
            print(2,'-dpdf',sprintf('BaseProfileT%gD%g%s',ThicknessShell,D,LabelType));
            % plot stress close to sea edge as a function of sea angle
            % (Figure 9)
            figure(1);
            plot(Aall,Slift/1000);           
            
            
        end
        plot([0,90],[0,0],'k');
        if ModelType==2;
        axis([0,90,-25,125]);
        set(gca,'dataaspectratio',[1,3,1])
        elseif ModelType==3
        axis([0,90,-100/4,50/4]);
        set(gca,'dataaspectratio',[4,3,1])
        else
        axis([0,90,-25/4,125/4]);
        set(gca,'dataaspectratio',[4,3,1])
        end
        print(1,'-dpdf',sprintf('BaseT%g%s',ThicknessShell,LabelType));
    end
end
% return