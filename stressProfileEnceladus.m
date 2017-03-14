function [Tectonics, MidFail]=stressProfileEnceladus(model,Yshear,Ycrack,Label);
% [Tectonics, MidFail]=stressProfileEnceladus(model,Yshear,Ycrack,Label)%%%
% Master file used to automatically generate many models of the stress
% field generated in the ice shell of Enceladus upon pressurization of an
% internal ocean or regional sea.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Uses LiveLink to matlab and Comsol Multiphysics 5.2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Used in Johnston and Montesi, Journal of Geophysical Research, 2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Extract Model parameters
Rsurf=model.param().evaluate('Rsurface');
Rcore=model.param().evaluate('Rcore');
Tshell=model.param().evaluate('ThicknessShell');
Tocean=model.param().evaluate('ThicknessOcean');
RoceanTop=model.param().evaluate('RoceanTop');
OceanPressure=model.param().evaluate('OceanPressure');

Ycenter=@(OceanAngle)(RoceanTop^2-Rcore^2)./(RoceanTop-Rcore*cos(OceanAngle))/2;
OceanCurvature=@(OceanAngle)RoceanTop-Ycenter(OceanAngle);

% Define coordinates of points at the surface
nA=181;
Alin=linspace(0,pi,nA);
Xsurf=Rsurf*sin(Alin);
Ysurf=Rsurf*cos(Alin); %0 at the south

% what models do we have?
dataset='dset2';
INFO2=mphsolinfo(model,'dataset',dataset);
nmod=numel(INFO2.batch.sol);

% initialize stress arrays
Srr=NaN(nmod,nA); Szz=Srr; Srz=Srr; Spp=Srr; Sp1=Srr;
Trn=Srr; Trt=Trn; Ttn=Trn; Ttt=Trn; Tzn=Trn; Tzt=Trn; Tnn=Trn; Tnt=Trn; Tll=Trn;
Norms=NaN(3,nA);
Norms(1,:)=sin(Alin); Norms(2,:)=cos(Alin); Norms(3,:)=0;
Alind=180-Alin*180/pi;
%%
for is=1:nmod
    % Extract stress
    [Smises(is,:),Sp1(is,:),Srr(is,:),Szz(is,:),Spp(is,:),Srz(is,:)]=mphinterp(model,...
        {'solid.mises','solid.sp1','solid.sr','solid.sz','solid.sphi','solid.srz'},...
        'coord',[Xsurf(:)'; Ysurf(:)'],'dataset',dataset,'outersolnum',is);
    % tractions on normal and tangential surface
    Trn(is,:)=[Srr(is,:).*Norms(1,:)+Srz(is,:).*Norms(2,:)];
    Tzn(is,:)=[Srz(is,:).*Norms(1,:)+Szz(is,:).*Norms(2,:)];
    Trt(is,:)=[Srr(is,:).*Norms(2,:)-Srz(is,:).*Norms(1,:)];
    Tzt(is,:)=[Srz(is,:).*Norms(2,:)-Szz(is,:).*Norms(1,:)];
    % rotate to normal/tangential frame => gives stress in
    % normal/tangential frame
    Tnn(is,:)=[Trn(is,:).*Norms(1,:)+Tzn(is,:).*Norms(2,:)];
    Ttn(is,:)=[Trn(is,:).*Norms(2,:)-Tzn(is,:).*Norms(1,:)];
    Ttt(is,:)=[Trt(is,:).*Norms(2,:)-Tzt(is,:).*Norms(1,:)];
    Tnt(is,:)=[Trt(is,:).*Norms(1,:)+Tzt(is,:).*Norms(2,:)];
    Tll(is,:)=Spp(is,:);
    % store ocean angle
    Inow=mphsolinfo(model,'soltag',INFO2.batch.sol(is));
    Runall(is).info=Inow;
    Runall(is).OceanAngle=Inow.paramsweepvals;
    
end
%% Detect failure
Crack=double(Sp1>=Ycrack); %cracking criterion
Yield=double(Smises>=Yshear); %faulting criterion
Fmode=Crack-Yield.*(1-Crack); %combined failure mode index
Regime=(Tnn>Ttt)+2*(Ttt>Tll)+2*(Tnn>Tll)+1; %stress regime
Regime=Regime.*Fmode; %stress regim at failure
OceanAll=[Runall.OceanAngle]*180/pi; %needed for plots
%% Classify Failure
[Tectonics,Asouth,Anorth,MidFail]=classifyTectonics(Crack, Fmode, Regime, Alind);

%% figures reporting stress regime and failure
%First plot: kind of failure, just grey for cracks, red for faults
figure(1); clf; 
hold on; box on; set(gca','fontsize',12);
image(Alind,OceanAll,Fmode,'CDataMapping','scaled');
colormap([1,0,0;1,1,1;0.5,0.5,0.5]); set(gca,'Clim',[-1,1]);
axis([0,180,0,90])
box on;
xlabel('Colatitude (\circ)','fontsize',18);
ylabel('Ocean angle (\circ)','fontsize',18);
print(1,'-dpdf',sprintf('%s_Failure.pdf',Label));

%Second plot: Stress regime at failure; makes panel for ocean angle
colpositive=[1,0,0;0,1,0;1,0,1;0,0,1;1,1,0;0,1,1];
colall=[flipud(colpositive)*0.5;[1,1,1];colpositive];

figure(2); clf; 
hold on; box on; set(gca','fontsize',12)
image(Alind,OceanAll,Regime,'CDataMapping','scaled');% shading flat;
colormap(colall);
set(gca,'Clim',[-6.5,6.5]);
axis([0,180,0,90])
box on;
xlabel('Colatitude (\circ)','fontsize',18);
ylabel('Ocean angle (\circ)','fontsize',18);
colorbar
print(2,'-dpdf',sprintf('%s_Regime.pdf',Label));

%% Stress profiles (Figure 5, 7)
YL=[-25,50];
figure(5); clf; hold on; set(gca','fontsize',12)
    axis([0,180,YL(1),YL(2)]);    box on;
    xlabel('Colatitude (\circ)','fontsize',18);
    ylabel('First Principal Stress (kPa)','fontsize',18);
    title(Label)
    plot([0,180],Ycrack*[1,1]/1e3,'k','linewidth',0.5)
figure(6); clf; hold on; set(gca','fontsize',12)
    axis([0,180,YL(1),YL(2)]);    box on;
    xlabel('Colatitude (\circ)','fontsize',18);
    ylabel('von Mises Stress (kPa)','fontsize',18);
    title(Label)
    plot([0,180],Yshear*[1,1]/1e3,'k','linewidth',0.5)

for is=1:nmod
    figure(5);     plot(Alind,Sp1(is,:)/1e3);
    figure(6);     plot(Alind,Smises(is,:)/1e3);
    
    figure(6+is); clf; hold on; set(gca','fontsize',12)
    fill([0,Asouth(is),Asouth(is),0],YL([1,1,2,2]),[1,1,1]*.75,'edgecolor','none','FaceAlpha',0.5);
    fill([180,Anorth(is),Anorth(is),180],YL([1,1,2,2]),[1,1,1]*.75,'edgecolor','none','FaceAlpha',0.5);
    faulting=find(Fmode(is,:)==-1);
    if ~isempty(faulting);
        fmin=min(faulting);fmax=max(faulting);
        fill(Alind([fmin,fmax,fmax,fmin]),YL([1,1,2,2]),[1,0.75,0.75],'edgecolor','none','FaceAlpha',0.5);
    end
    H(1)=plot(Alind,Tll(is,:)/1e3,'r','linewidth',2);
    H(2)=plot(Alind,Tnn(is,:)/1e3,'g','linewidth',2);
    H(3)=plot(Alind,Ttt(is,:)/1e3,'b','linewidth',2);
    H(4)=plot(Alind,Smises(is,:)/1e3,'k','linewidth',2);
    legend(H,'Longitudinal','Vertical','Latitudinal','von Mises');
    plot([0,180],Yshear*[1,1]/1e3,'k','linewidth',1);
    plot([0,180],Ycrack*[1,1]/1e3,'k','linewidth',1);
    plot([1,1]*Runall(is).OceanAngle*180/pi,YL,'k--')
    axis([0,180,YL(1),YL(2)]);
    box on;
    % legend('nn','tt','ll','mises','first prp','nt','tn')
    xlabel('Colatitude (\circ)','fontsize',18);ylabel('Stress (kPa)','fontsize',18);
    OA=OceanAll(is); C=OceanCurvature(OA*pi/180);
    title(sprintf('%sA%gC%g  Regime %d',Label,OA,C/1000,Tectonics(is)));
%     print(7,'-dpdf',sprintf('%s%gA_profile.pdf',Label,OA))
end
% set(5,'ylim',[0,YL(2)]);
print(5,'-dpdf',sprintf('%s_Sp1.pdf',Label))

% set(5,'ylim',[0,YL(2)]);
print(6,'-dpdf',sprintf('%s_Smises.pdf',Label))

return