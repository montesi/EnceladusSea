% EnceladusBatch.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Master file used to automatically generate many models of the stress
% field generated in the ice shell of Enceladus upon pressurization of an
% internal ocean or regional sea.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Uses LiveLink to matlab and Comsol Multiphysics 5.2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Used in Johnston and Montesi, Journal of Geophysical Research, 2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Range of parameters to be considered.
ThicknessShell=40; %km
ModelAll=[1:4]; %codes for model configurations
nA=46; % number of sea angles to consider
dT=2; % increment of sea thickness

% generate model configurations
OceanAngle=linspace(0,90,nA);
ThickAll=[dT:dT:ThicknessShell-1];nT=numel(ThickAll);
% initialize matrices to store global tectonic regimes
MidFail=NaN(nA,nT);
Tectonics=NaN(nA,nT);

% Loop thought the model configurations 
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

    for iT=1:nT;
        % solve model
        [ModelOut,LabelMod]=EnceladusBuild(ThicknessShell,ThickAll(iT),ModelType,nA);        
        % makefigureEnceladus(ModelOut,Label) %uncomment to make map of stress field
        % [Tn,Mn]=analyseStressEnceladus(ModelOut,20e3,2e3,Label);
        % [Tectonics(:,iT), MidFail(:,iT)]=analyseStressEnceladus(ModelOut,15000,11000,LabelMod);
        % Determine tectonic regime
        [Tectonics(:,iT), MidFail(:,iT)]=stressProfileEnceladus(ModelOut,15000,11000,LabelMod);
    end
    %% Plot global regime (for Figure 9 or S1 in paper)
    Rcore=252-ThicknessShell;
    Rtop=Rcore+repmat(ThickAll,[nA,1]);
    OC=repmat(OceanAngle',[1,nT]);
    Yc=(Rtop.^2-Rcore.^2)./(Rtop-Rcore.*cosd(OC))/2;
    C=Rtop-Yc; %indentation curvature
    % Each color corresponds to a code and global tectonic regime
    ColorSequence=[0,0,1;0,1,0;1,0,1;1,0,0;1,1,0];
    ColorMapCustom=[flipud(ColorSequence)*0.5;[1,1,1];(ColorSequence)];
    
    % Plot as a function of ocean angle (Figure 9)
    figure(3); clf;
    hold on; box on; set(gca','fontsize',12)
    imagesc(OceanAngle,ThickAll,Tectonics');% shading flat;
    colormap(ColorMapCustom);
    set(gca,'Clim',[-5.5,5.5]);
    axis([0,90,0,ThicknessShell])
    hold on;
    contour(OceanAngle,ThickAll,MidFail',0.5*[1,1],'k')
    box on;
    xlabel('Ocean angle (\circ)','fontsize',18);
    ylabel('Sea Thickness (km)','fontsize',18);
    colorbar
    print(3,'-dpdf',sprintf('%s_RegimeAngle.pdf',Label));
    
    % Plot as a function of curvature (not used in the paper)
    figure(4); clf;
    hold on; box on; set(gca','fontsize',12)
    pcolor(C,Rtop-Rcore,Tectonics); shading flat;
    colormap(ColorMapCustom);
    set(gca,'Clim',[-5.5,5.5]);
    axis([0,Rcore,0,ThicknessShell]);
    contour(C,Rtop-Rcore,MidFail,0.5*[1,1],'k')
    box on;
    xlabel('Sea Curvature (km)','fontsize',18);
    ylabel('Sea Thickness (km)','fontsize',18);
    colorbar
    print(4,'-dpdf',sprintf('%s_RegimeCurvature.pdf',Label));
end
