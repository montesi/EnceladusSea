function [out, Label] = EnceladusBuild(ThicknessShell,SeaThick,ModelType,nA);
% [out, Label] = EnceladusBuild(ThicknessShell,SeaThick,ModelType,nA); %%%%
% Setup and solve COMSOL multiphysics model to calculate the stress field
% in the ice shell of Enceladus upon pressurization of an internal global
% ocean or regional sea
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Uses LiveLink to matlab and Comsol Multiphysics 5.2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Used in Johnston and Montesi, Journal of Geophysical Research, 2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Model exported on Jan 17 2017, 10:25 by COMSOL 5.2.1.152.

import com.comsol.model.*
import com.comsol.model.util.*

model = ModelUtil.create('Model');

model.modelPath('/Users/laurentmontesi/Documents/COMSOL/Enceladus');

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
fprintf('Model Type: %s, with T=%gkm and D=%gkm\n',LabelType,ThicknessShell, SeaThick);
Label=sprintf('EnceladusT%gD%g%s.mph',ThicknessShell,SeaThick,LabelType);

% Prepare for parametric study varying systematically sea angle
Amin=0;Amax=90; dA=(Amax-Amin)/(nA-1);
Range=sprintf('range(%g,%g,%g)',Amin,dA,Amax);
model.label(Label);
model.comments(Label);

%% setup parameters and variables
% Specific to this run
model.param.set('ThicknessShell', sprintf('%g [km]',ThicknessShell));%'20 [km]');
model.param.set('ThicknessOcean', sprintf('%g [km]',SeaThick));%'10 [km]');
% Rest is default
model.param.set('Rsurface', '252 [km]');
model.param.set('Rcore', 'Rsurface-ThicknessShell');
model.param.set('RoceanTop', 'Rcore+ThicknessOcean');
model.param.set('OceanAngle', '40 [deg]'); %overwriten during parameter sweep
model.param.set('Ycenter', '(RoceanTop^2-Rcore^2)/(RoceanTop-Rcore*cos(OceanAngle))/2');
model.param.set('OceanCurvature', 'RoceanTop-Ycenter');
if le(ModelType,2);
    model.param.set('OceanPressure', '-10 [kPa]');
else
    model.param.set('OceanPressure','-10[kPa]/(3/2/((Rsurface/Rcore)^3-1))');
end
model.param.set('Gravity', '0.113 [m/s^2]');
model.param.set('RhoIce', '930[kg/m^3]');
model.param.set('GIce', '9 [GPa]');
if ModelType==4; %Need a cavity at North Pole
    model.param.set('ThicknessNorth', 'ThicknessOcean/2');%'10 [km]');
    model.param.set('RNorthTop', 'Rcore+ThicknessNorth');
    model.param.set('NorthCurvature', 'OceanCurvature');
    model.param.set('YNorth', 'RNorthTop-NorthCurvature');
%     model.param.set('NorthAngle', '(RNorthTop/ThicknessNorth)-ThicknessNorth*(Rcore+RNorthTop)/(2*Rcore*YNorth)');
end

%% define geometrical elements
model.modelNode.create('comp1');

model.geom.create('geom1', 2);

model.modelNode('comp1').defineLocalCoord(false);

model.geom('geom1').axisymmetric(true);

model.mesh.create('mesh1', 'geom1');

% Base of the ice shell
model.geom('geom1').repairTolType('relative');
model.geom('geom1').create('c1', 'Circle');
model.geom('geom1').feature('c1').set('r', 'Rcore');
model.geom('geom1').feature('c1').set('rot', '-90');
model.geom('geom1').feature('c1').set('angle', '180');
% Surface of the satellite
model.geom('geom1').create('c2', 'Circle');
model.geom('geom1').feature('c2').set('r', 'Rsurface');
model.geom('geom1').feature('c2').set('rot', '-90');
model.geom('geom1').feature('c2').set('angle', '180');
% Ice shell
model.geom('geom1').create('dif1', 'Difference');
model.geom('geom1').feature('dif1').selection('input2').set({'c1'});
model.geom('geom1').feature('dif1').selection('input').set({'c2'});
% south polar indentation
model.geom('geom1').create('c3', 'Circle');
model.geom('geom1').feature('c3').set('r', 'OceanCurvature');
model.geom('geom1').feature('c3').set('rot', '-90');
model.geom('geom1').feature('c3').set('angle', '180');
model.geom('geom1').feature('c3').set('pos', {'0' '-Ycenter'});
model.geom('geom1').create('dif2', 'Difference');
model.geom('geom1').feature('dif2').selection('input2').set({'c3'});
model.geom('geom1').feature('dif2').selection('input').set({'dif1'});
if ModelType==4; %Need north polar indentation
    model.geom('geom1').create('c4', 'Circle');
    model.geom('geom1').feature('c4').set('r', 'NorthCurvature');
    model.geom('geom1').feature('c4').set('rot', '-90');
    model.geom('geom1').feature('c4').set('angle', '180');
    model.geom('geom1').feature('c4').set('pos', {'0' 'YNorth'});
    model.geom('geom1').create('dif3', 'Difference');
    model.geom('geom1').feature('dif3').selection('input2').set({'c4'});
    model.geom('geom1').feature('dif3').selection('input').set({'dif2'});
end
model.geom('geom1').run;
% model.geom('geom1').run('fin');

%% Define variables
model.variable.create('var1');
model.variable('var1').model('comp1');
model.variable('var1').set('Dc', 'sqrt(r^2+z^2)', 'Distance to core');
model.variable('var1').set('TH', 'acos(z/Dc)', 'colatitude');
model.variable('var1').set('un', 'u*sin(TH)+w*cos(TH)', 'radial displacement');

%% Define materials
model.material.create('mat1', 'Common', 'comp1');

%% Define boundary conditions
model.physics.create('solid', 'SolidMechanics', 'geom1');
% Ocean pressure
model.physics('solid').create('bndl1', 'BoundaryLoad', 1);
model.physics('solid').feature('bndl1').set('FperArea', {'0'; '0'; 'OceanPressure'});
model.physics('solid').feature('bndl1').set('coordinateSystem', 'sys1');
% Surface
model.physics('solid').create('bndl2', 'BoundaryLoad', 1);
model.physics('solid').feature('bndl2').set('FperArea', {'-un*RhoIce*Gravity*sin(TH)'; '0'; '-un*RhoIce*Gravity*cos(TH)'});
% Base 
model.physics('solid').create('fix1', 'Fixed', 1);
model.physics('solid').create('roll1', 'Roller', 1);
model.physics('solid').create('disp1', 'Displacement0', 0);
model.physics('solid').feature('disp1').set('Direction', {'0'; '0'; '1'});
% Assign boundary conditions
switch ModelType
    case 1; % Fixed boundary
        model.physics('solid').feature('bndl1').selection.set([4]);
        model.physics('solid').feature('bndl2').selection.set([3 6]);
        model.physics('solid').feature('fix1').selection.set([5 7]);
        model.physics('solid').feature('roll1').active(false);
        model.physics('solid').feature('disp1').active(false);
    case 2; % roller boundary
        model.physics('solid').feature('bndl1').selection.set([4]);
        model.physics('solid').feature('bndl2').selection.set([3 6]);
        model.physics('solid').feature('roll1').selection.set([5 7]);
        model.physics('solid').feature('fix1').active(false);
        model.physics('solid').feature('disp1').active(false);
    case 3; % Ocean
        model.physics('solid').feature('bndl1').selection.set([4 5 7]);
        model.physics('solid').feature('bndl2').selection.set([3 6]);
        model.physics('solid').feature('disp1').selection.set([6]);
        model.physics('solid').feature('fix1').active(false);
        model.physics('solid').feature('roll1').active(false);
    case 4; % North
        model.physics('solid').feature('bndl1').selection.set([4 5 7 8]);
        model.physics('solid').feature('bndl2').selection.set([3 6]);
        model.physics('solid').feature('disp1').selection.set([5]);
        model.physics('solid').feature('fix1').active(false);
        model.physics('solid').feature('roll1').active(false);
    otherwise
        error(sprintf('Not ready for option %g',ModelType));
end
%% Mesh
model.mesh('mesh1').autoMeshSize(1);

% model.view('view1').axis.set('abstractviewrratio', '0.049999941140413284');
% model.view('view1').axis.set('abstractviewlratio', '-0.05000000447034836');
% model.view('view1').axis.set('abstractviewxscale', '1046.0311279296875');
% model.view('view1').axis.set('abstractviewbratio', '-0.10990232229232788');
% model.view('view1').axis.set('xmax', '265122.125');
% model.view('view1').axis.set('xmin', '-13122.1416015625');
% model.view('view1').axis.set('abstractviewyscale', '1046.0311279296875');
% model.view('view1').axis.set('ymax', '325838.6875');
% model.view('view1').axis.set('ymin', '-325838.6875');
% model.view('view1').axis.set('abstractviewtratio', '0.10990235209465027');

model.material('mat1').propertyGroup('def').set('density', 'RhoIce');
model.material('mat1').propertyGroup('def').set('youngsmodulus', 'GIce');
model.material('mat1').propertyGroup('def').set('poissonsratio', '0.3');

model.physics('solid').prop('ShapeProperty').set('order_displacement', '2');
%% Solve model!
model.mesh('mesh1').run;

model.study.create('std1');
model.study('std1').create('param', 'Parametric');
model.study('std1').create('stat', 'Stationary');
model.study('std1').feature('param').set('pname', {'OceanAngle'});
model.study('std1').feature('param').set('punit', {'deg'});
model.study('std1').feature('param').set('sweeptype', 'filled');
model.study('std1').feature('param').set('plistarr', {Range});

model.sol.create('sol1');
model.sol('sol1').study('std1');

model.study('std1').feature('stat').set('notlistsolnum', 1);
model.study('std1').feature('stat').set('notsolnum', '1');
model.study('std1').feature('stat').set('listsolnum', 1);
model.study('std1').feature('stat').set('solnum', '1');

model.sol('sol1').create('st1', 'StudyStep');
model.sol('sol1').feature('st1').set('study', 'std1');
model.sol('sol1').feature('st1').set('studystep', 'stat');
model.sol('sol1').create('v1', 'Variables');
model.sol('sol1').feature('v1').set('control', 'stat');
model.sol('sol1').create('s1', 'Stationary');
model.sol('sol1').feature('s1').create('fc1', 'FullyCoupled');
model.sol('sol1').feature('s1').feature('fc1').set('termonres', 'auto');
model.sol('sol1').feature('s1').feature('fc1').set('reserrfact', 1000);
model.sol('sol1').feature('s1').feature('fc1').set('linsolver', 'dDef');
model.sol('sol1').feature('s1').feature('fc1').set('termonres', 'auto');
model.sol('sol1').feature('s1').feature('fc1').set('reserrfact', 1000);
model.sol('sol1').feature('s1').feature.remove('fcDef');
model.sol('sol1').attach('std1');

model.batch.create('p1', 'Parametric');
model.batch('p1').study('std1');
model.batch('p1').create('so1', 'Solutionseq');
model.batch('p1').feature('so1').set('seq', 'sol1');
model.batch('p1').feature('so1').set('store', 'on');
model.batch('p1').feature('so1').set('clear', 'on');
model.batch('p1').feature('so1').set('psol', 'none');
model.batch('p1').set('pname', {'OceanAngle'});
model.batch('p1').set('plistarr', {Range});
model.batch('p1').set('sweeptype', 'filled');
model.batch('p1').set('probesel', 'all');
model.batch('p1').set('probes', {});
model.batch('p1').set('plot', 'off');
model.batch('p1').set('err', 'on');
model.batch('p1').attach('std1');
model.batch('p1').set('control', 'param');

model.sol.create('sol2');
model.sol('sol2').study('std1');
model.sol('sol2').label('Parametric Solutions 1');

model.batch('p1').feature('so1').set('psol', 'sol2');

model.batch('p1').run;

%% Save! 

mphsave(model,Label);
% Uncomment if you want to see the stress field (Figure 4 in paper)
% makefigureEnceladus(model,Label) 

out = model;
