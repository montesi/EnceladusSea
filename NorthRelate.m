%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot sea angle for north identation
% Used in Johnston and Montesi, Journal of Geophysical Research, 2017,
% Figure 6
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Rs=252; % surface radius
for T=[10,20,40,60,80];
    for Ds=[0.1:0.1:0.9]*T; % which sea thickness
    
    %% find curvature in the South
    Rc=Rs-T; %core/base of shell radius 
    Rt=Rc+Ds; % radius at the top of the South Polar indentation
    As=linspace(0,90,100); % set of south polar indentation angles
    Ys=(Rt.^2-Rc.^2)./(Rt-Rc.*cosd(As))/2; % center of the south polar indentation
    Cs=Rt-Ys; %curvature of the south polar indentation
    
    %%
    Cn=Cs; Dn=Ds/2; %curvature and thickness of the north polar indentation 
    Yn=Rc+Dn-Cn; % Center of the north polar indentation
    An=acosd((Rc.^2-Cn.^2+Yn.^2)./(2.*Yn.*Rc)); %north indentation angle
    
    %%
    figure(1); clf;
    hold on; box on; set(gca','fontsize',12)
    plot(180-An,As,'k');
    plot([0,90],[0,90],'k');
    xlabel('Colatitude (\circ)','fontsize',18);
    ylabel('Ocean angle (\circ)','fontsize',18);
    
    print(1,'-dpdf',sprintf('NorthRelateT%gD%g.pdf',T,Ds));
    end
end
%% Example geometry (last set of parameters)
i=100;
Th=linspace(0,180,100);xc=sind(Th);yc=cosd(Th);
figure (2); clf; hold on;
plot(Rs*xc,Rs*yc,'k'); %surface
plot(Rc*xc,Rc*yc,'k'); %ice shell
plot(Cs(i)*xc,Cs(i)*yc-Ys(i),'r'); %south polar indentation
plot([0,Rc*sind(As(i))],[0,-Rc*cosd(As(i))],'r');
plot(Cn(i)*xc,Cn(i)*yc+Yn(i),'b'); %north polar indentation
plot([0,Rc*sind(An(i))],[0,Rc*cosd(An(i))],'b');
axis equal;
axis tight;

