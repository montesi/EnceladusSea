function makefigureEnceladus(model,Label)
% makefigureEnceladus(model,Label) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Visualizes the stress field in model "model"; Uses label "Label" for
% printout name a
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Uses LiveLink to matlab and Comsol Multiphysics 5.2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Used in Johnston and Montesi, Journal of Geophysical Research, 2017;
% Especially, makes Figure 4
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('Making figures: model section\n')
% close all
FigList={'pg1','pg2'};nFig=numel(FigList);
T=model.result.tags;

for iFig=1:nFig; % may need to reset figure structure
    if ~isempty(strmatch(FigList(iFig),T))
        model.result.remove(FigList(iFig))
    end
end
% Setup figure with von Mises stress
model.result.create('pg1', 2);
model.result('pg1').set('data', 'dset2');

model.result('pg1').create('surf1', 'Surface');
model.result('pg1').feature('surf1').set('expr', 'solid.mises');
model.result('pg1').feature('surf1').set('descr', 'von Mises stress');
model.result('pg1').feature('surf1').set('unit', 'kPa');
model.result('pg1').feature('surf1').set('rangecoloractive', 'on');
model.result('pg1').feature('surf1').set('rangecolormin', '0');
model.result('pg1').feature('surf1').set('rangecolormax', '40');
model.result('pg1').feature('surf1').set('colortable', 'Thermal');

model.result('pg1').create('con1', 'Contour');
model.result('pg1').feature('con1').set('expr', 'solid.sp1');
model.result('pg1').feature('con1').set('descr', 'First principal stress');
model.result('pg1').feature('con1').set('unit', 'kPa');
model.result('pg1').feature('con1').set('levelmethod', 'levels');
model.result('pg1').feature('con1').set('levels', 'range(0,20,100)');
model.result('pg1').feature('con1').set('coloring', 'uniform');
model.result('pg1').feature('con1').set('color', 'black');
model.result('pg1').feature('con1').set('colorlegend', 'off');


% Setup figure with first principal stress
model.result.create('pg2', 2);
model.result('pg2').set('data', 'dset2');

model.result('pg2').create('surf2', 'Surface');
model.result('pg2').feature('surf2').set('expr', 'solid.sp1');
model.result('pg2').feature('surf2').set('descr', 'First principal stress');
model.result('pg2').feature('surf2').set('unit', 'kPa');
model.result('pg2').feature('surf2').set('rangecoloractive', 'on');
model.result('pg2').feature('surf2').set('rangecolormin', '0');
model.result('pg2').feature('surf2').set('rangecolormax', '40');
model.result('pg2').feature('surf2').setIndex('const', '0', 0, 1);
model.result('pg2').feature('surf2').set('colortablesym', 'off');
model.result('pg2').feature('surf2').set('colortable', 'AuroraBorealis');
model.result('pg2').feature('surf2').set('colortablerev', 'on');

model.result('pg2').create('con2', 'Contour');
model.result('pg2').feature('con2').set('levelmethod', 'levels');
model.result('pg2').feature('con2').set('expr', 'solid.sp1');
model.result('pg2').feature('con2').set('descr', 'First principal stress');
model.result('pg2').feature('con2').set('unit', 'kPa');
model.result('pg2').feature('con2').set('levels', 'range(0,20,100)');
model.result('pg2').feature('con2').set('coloring', 'uniform');
model.result('pg2').feature('con2').set('colorlegend', 'off');
model.result('pg2').feature('con2').set('color', 'black');

%% access solution set
INFO=mphsolutioninfo(model);
MAP=INFO.sol2.map;
% how many models are there (i.e. values of sea angle in parameter sweep
nmod=size(MAP,1); 

for imod=1:nmod; %Repeat for every sea angle value
    A=rad2deg(MAP(imod,1));
    LabelA=sprintf('%sA%g',Label,A);
    
    model.result('pg1').setIndex('looplevel', imod, 0);
    model.result('pg1').run;
    model.result('pg2').setIndex('looplevel', imod, 0);
    model.result('pg2').run;
    
    % von Mises figure
    figure(10); mphplot(model,'pg1');
    Rs=model.param.evaluate('Rsurface');
    axis equal;axis tight;
    set(gca,'Visible','off');
    colormap(colortable('Thermal')); set(gca,'clim',[0,40])
    colorbar('southoutside')
    print(10,'-dpng',sprintf('%sVonMises.png',LabelA));
    
    % first principal stress figure
    figure(20); mphplot(model,'pg2');
    Rs=model.param.evaluate('Rsurface');
    axis equal;axis tight;
    set(gca,'Visible','off');
    colormap(flipud(colortable('AuroraBorealis'))); set(gca,'clim',[0,40])
    colorbar('southoutside');
    print(20,'-dpng',sprintf('%sTension.png',LabelA));
    
end
