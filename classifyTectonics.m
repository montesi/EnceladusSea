function [Tectonics,Asouth,Anorth,MidFail]=classifyTectonics(Crack, Fmode, Regime, Alind)
% [Tectonics,Asouth,Anorth,MidFail]=classifyTectonics(Crack, Fmode, Regime, Alind) %
% Analyses global stress profile and assign regime code, according to
% Figure 7 in the paper.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Uses LiveLink to matlab and Comsol Multiphysics 5.2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Used in Johnston and Montesi, Journal of Geophysical Research, 2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nmod=size(Crack,1);nA=numel(Alind);
for is=1:nmod;
    % Bounds to the region of tensile failure.
    NoCrack=find(Crack(is,:)==0);
    if ~isempty(NoCrack);
        NoCrackMin=min(NoCrack);
        if NoCrackMin~=1;
            An=(Alind(NoCrackMin)+Alind(NoCrackMin-1))/2;
        else
            An=Alind(NoCrackMin);
        end
        NoCrackMax=max(NoCrack);
        if NoCrackMax~=nA;
            As=(Alind(NoCrackMax)+Alind(NoCrackMax-1))/2;
        else
            As=Alind(NoCrackMax);
        end
    else; %cracks everywhere
        [A,NoCrackMin]=min(abs(Alind-90));
        As=90+A;
        An=90+A;NoCrackMax=NoCrackMin;
    end
    MidFail(is)=max(abs(Fmode(is,[NoCrackMin:NoCrackMax])));
    RegCg=diff(Regime(is,[nA-1:-1:NoCrackMax]));
    % classify tectonic regime;
    iA=nA-1; %avoid pole
    Regnow=Regime(is,iA);
    switch Regnow; %Failure mode at the South Pole
        case 1; % NS cracks: A, E, or C
            Tn=1; %assume no change but we'll look into it
            iCg=min(find(RegCg~=0)); % First change
            if RegCg(iCg)==2; %switch to regime 3 EW cracks: E or C
                RegCg(iCg)=0;
                iCg=min(find(RegCg~=0)); % Next change
                if RegCg(iCg)==-2; %switch back to regime 1 NS cracks
                    Tn=3; %C
                else
                    Tn=5; %E
                end
            end
        case 3; %EW Cracks: B or D
            iCg=min(find(RegCg~=0)); % change
            if RegCg(iCg)==-2; %switch to regime 1 NS cracks
                Tn=2; %B
            else
                Tn=4; %D
            end
        otherwise
            Tn=0;
    end
    Tectonics(is)=Tn;
    Asouth(is)=As;
    Anorth(is)=An;
end
Tectonics=Tectonics'.*(-1).^Fmode(:,2);