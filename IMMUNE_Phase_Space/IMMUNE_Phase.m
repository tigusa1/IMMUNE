function IMMUNE_Phase(handles)
%----------------------------------------------------------------------------------------------
% 2D and 3D phase space analysis
%----------------------------------------------------------------------------------------------
% Initialize
%----------------------------------------------------------------------------------------------
% global Slim0 Sslope Ilim0 Islope Clim0 Cslope SIClast SScCoeff
% global Imult ICmult Smult
% global S0 I0 C0 
global St Ct It t Ssaved Csaved Isaved tsaved Slast Clast Ilast

set(0,'DefaultTextFontSize',12), set(0,'DefaultAxesFontSize',12)

if nargin
    flag_3D = false;
else
    flag_3D = false; % option for 3D
    handles = [];
end

[ Slim0,Clim0,Ilim0,Sslope,Cslope,Islope,S0,C0,I0,doSave ] = setParam(nargin,handles);

if nargin
    set(handles.C0_txt,     'String',sprintf('C0: %.2f',      C0    ))
    set(handles.Islope_txt, 'String',sprintf('Islope: %.2f',  Islope))
    set(handles.Cslope_txt, 'String',sprintf('Cslope: %.2f',  Cslope))
    set(handles.Sslope_txt, 'String',sprintf('Sslope: %.2f',  Sslope))
end

if doSave, return, end

[ tnew,X ] = ode45( @(t,X) fd3(t,X, Slim0,Clim0,Ilim0,Sslope,Cslope,Islope), ...
    [0 100], [S0;C0;I0] ); % solve ODE

Jt = find( max(abs(X')) < 1.2 );                 % only keep values less than 1.2

if isempty(tsaved)
    tlast = 0;
else
    tlast = tsaved(end);
end
t  = [tsaved;tnew(Jt)+tlast];                    % truncate
St = [Ssaved;X(Jt,1)]; Ct = [Csaved;X(Jt,2)]; It = [Isaved;X(Jt,3)];
Slast = X(1,1);        Clast = X(1,2);        Ilast = X(1,3);

plotResults( Slim0,Clim0,Ilim0,Sslope,Cslope,Islope,nargin,handles )

if flag_3D
    plotResults3D(nargin)
end


function plotResults3D(Nargin)
%----------------------------------------------------------------------------------------------
% Plot animated series of points in 3D
%----------------------------------------------------------------------------------------------
global St Ct It t

os = [1 1];

drawnow % finish drawing old figures before animation
if Nargin, fprintf('\n***** Paused for 3D plot *****\n'), pause, end

fig = figure(250); fig.Name = '3D fields';
% quiver3(SS,CS,IS,SIx,CIy,CSIz), axis tight
% plot3(St,Ct,It,'ro')

for j=1:2:round(length(St)*0.7)
%   plot3(St,Ct,It,'ro')
    plot3(St(1:j),Ct(1:j),It(1:j),'ko') % redraw plot
    hold on
    plot3(St(j),Ct(j),It(j),'ko')
    axis([0 1 0 1 0 0.5])
    xlabel('SSc'), ylabel('\color{red}C'), zlabel('I')
    plot3([0 St(j)],Ct(j)*os,It(j)*os,'k-')
    plot3(St(j)*os,[0 Ct(j)],It(j)*os,'r-')
    plot3(St(j)*os,Ct(j)*os,[0 It(j)],'g-'), grid on
    ax = gca; ax.YColor = 'red';
    pause(0.025)
    hold off
end


function plotResults( Slim0,Clim0,Ilim0,Sslope,Cslope,Islope,Nargin,handles )
%----------------------------------------------------------------------------------------------
% 2D plotting: I, C, SSc
%----------------------------------------------------------------------------------------------
global St Ct It t

n   = 26;                                   % number of points in each coordinate

Imax = 0.8;                                 % limit for I-axis
Is  = linspace(0,Imax,n);                   % 1D coordinates
Ss  = linspace(0,1,n);
Cs  = linspace(0,1,n);
[ SS    ] = meshgrid( Ss,Is );              % 2D coordinates
[ CS,IS ] = meshgrid( Cs,Is );
SIx = fSIx( SS, IS, Slim0, Sslope );                       % 2D vector fields
CIx = fCIx( CS, IS, Clim0, Cslope );
SICend= fSIy( SS, Ct(end), IS, Ilim0, Islope );
SISend= fSIy( St(end), CS, IS, Ilim0, Islope );
[ ~,Sx ] = fSIx( 0,  Is, Slim0, Sslope );                  % boundaries
[ ~,Cx ] = fCIx( 0,  Is, Slim0, Sslope );
[ ~,IyCend] = fSIy( Ss, Ct(end), 0, Ilim0, Islope  );
[ ~,IySend] = fSIy( St(end), Cs, 0, Ilim0, Islope  );

if ~Nargin, figure(110), clf, set(gcf,'Name','basic 2D fields','color','white'), end

if Nargin, axes(handles.axes1), else, subplot(2,2,1), end
hold off, p1a = plot(Ss,IyCend,'k-',Sx,Is,'b-'); hold on
axis tight, xlabel('SSc'), ylabel('Immune response'), axis([0 1 0 Imax])

if Nargin, axes(handles.axes2), else, subplot(2,2,2), end
hold off, p2a = plot(Cs,IySend,'k-',Cx,Is,'r-'); hold on
ax = gca; ax.XDir = 'reverse';
axis tight, xlabel('Cancer'), ylabel('Immune response'), axis([0 1 0 Imax])

if Nargin, axes(handles.axes3), else, subplot(2,2,[3 4]), end
cla

js = 1:3:round(length(St)*0.99);

if Nargin, axes(handles.axes1), else, subplot(2,2,1), end
p1b = plot(St(js),It(js),'g.','MarkerSize',9); plot(St(1),It(1),'g^','MarkerSize',9)
p1c = quiver(SS,IS,SIx,SICend);
legend([p1a;p1b;p1c],'I equilibrium','SSc equilibrium','disease trajectory','flow field',...
    'Location','NW')

if Nargin, axes(handles.axes2), else, subplot(2,2,2), end
p2b = plot(Ct(js),It(js),'g.','MarkerSize',9); plot(Ct(1),It(1),'g^','MarkerSize',9)
p2c = quiver(CS,IS,CIx*8,SISend);
legend([p2a;p2b;p2c],'I equilibrium','C equilibrium','disease trajectory','flow field')

if Nargin, axes(handles.axes3), else, subplot(2,2,[3 4]), end
plot(t(js),St(js),'b-',t(js),Ct(js),'r-',t(js),It(js),'k-','LineWidth',1)
ylabel('SSc, Cancer, Immune response'), xlabel('time'), hold on
axis([0 max(t) 0 0.8]), legend('Scleroderma','Cancer','Immune response')

%----------------------------------------------------------------------------------------------
% 3D plotting
%----------------------------------------------------------------------------------------------
[ SS,CS,IS ] = meshgrid( Ss,Cs,Is );        % coordinate change: SS wrt 2, CS wrt 1, IS wrt 3
SIx  = fSIx( SS, IS, Slim0, Sslope );                      % 3D vector fields
CIy  = fCIx( CS, IS, Clim0, Cslope );
CSIz = fSIy( SS, CS, IS, Ilim0, Islope );                  % simple sum of S and C used to determine change in I
[ ~,SIs ] = fSIx( 0, IS, Slim0, Sslope );
[ ~,CIs ] = fCIx( 0, IS, Clim0, Cslope );
[ ~,cSIs] = fSIy( SIs, CS,  0, Ilim0, Islope );
[ ~,CsIs] = fSIy( SS,  CIs, 0, Ilim0, Islope );
[ ~,CSIs] = fSIy( SS,  CS,  0, Ilim0, Islope );
Js = IS>cSIs;
Jc = IS<CsIs;

ss   = [SS(  :,[1 end],1)   SS([1 end],:,1)'];
cs   = [CS(  :,[1 end],1)   CS([1 end],:,1)'];
csis = [CSIs(:,[1 end],1) CSIs([1 end],:,1)'];
os   = ones(size(t));

axs = [0 1 0 1 0 0.8]; axs2 = [0 1 0 1];

if ~Nargin, fig = figure(220); fig.Name = 'overhead'; end

if Nargin, axes(handles.axes5), else, subplot(3,2,1), end
surf(squeeze(SIs(:,1,:)), squeeze(CS(:,1,:)), squeeze(IS(:,1,:)), 'FaceColor', [0.5 1 0.5], ...
    'FaceAlpha', 0.5, 'EdgeAlpha', 0.5 )
% hold on
surf(SS(:,:,1), CS(:,:,1), CSIs(:,:,1), 'FaceColor', [0.5 1 0.5], ...
    'FaceAlpha', 0.5, 'EdgeAlpha', 0.0 ), hold on
plot3(ss, cs, csis, 'k-', 'LineWidth', 0.3)
hold off
axis(axs), xlabel('SSc'), ylabel('C'), zlabel('I')
hold on, plot3(St,Ct,It,'r.',St(1),Ct(1),It(1),'r^'), hold off

if Nargin, axes(handles.axes6), else, subplot(3,2,3), end
surf(squeeze(SS(1,:,:)), squeeze(CIs(1,:,:)), squeeze(IS(1,:,:)), 'FaceColor', [0.5 0.5 1], ...
    'FaceAlpha', 0.5, 'EdgeAlpha', 0.5 )
hold on
surf(SS(:,:,1), CS(:,:,1), CSIs(:,:,1), 'FaceColor', [1 1 1]*0.8, ...
    'FaceAlpha', 0.5, 'EdgeAlpha', 0.0 )
plot3(ss, cs, csis, 'k-', 'LineWidth', 0.3)
hold off
axis(axs), xlabel('SSc'), ylabel('C'), zlabel('I')
hold on, plot3(St,Ct,It,'g.',St(1),Ct(1),It(1),'g^'), hold off

if Nargin, axes(handles.axes4), else, subplot(3,2,6), end
surf(squeeze(SS(1,:,:)), squeeze(CIs(1,:,:)), squeeze(IS(1,:,:)), 'FaceColor', [0.5 0.5 1], ...
    'FaceAlpha', 0.5, 'EdgeAlpha', 0.0 )
hold on
surf(squeeze(SIs(:,1,:)), squeeze(CS(:,1,:)), squeeze(IS(:,1,:)), 'FaceColor', [0.5 1 0.5], ...
    'FaceAlpha', 0.5, 'EdgeAlpha', 0.0 )
surf(SS(:,:,1), CS(:,:,1), CSIs(:,:,1), 'FaceColor', [1 1 1]*1.0, ...
    'FaceAlpha', 1.0, 'EdgeAlpha', 0.0 )
hold off
axis(axs2), xlabel('SSc'), ylabel('C')
hold on, plot3(St,Ct,os,'g.',St(1),Ct(1),os,'g^','MarkerSize',9), hold off, view(2)

if ~Nargin
    subplot(3,2,2)
    surf(squeeze(SIs(:,1,:)), squeeze(CS(:,1,:)), squeeze(IS(:,1,:)), 'FaceColor', [0.5 1 0.5], ...
        'FaceAlpha', 0.5, 'EdgeAlpha', 0.0 )
    hold on
    surf(SS(:,:,1), CS(:,:,1), CSIs(:,:,1), 'FaceColor', [1 1 1]*1.0, ...
        'FaceAlpha', 1.0, 'EdgeAlpha', 0.0 )
    hold off
    axis(axs2), xlabel('SSc'), ylabel('C')
    hold on, plot3(St,Ct,os,'g.',St(1),Ct(1),os,'g^','MarkerSize',9), hold off, view(2)
    
    subplot(3,2,4)
    surf(squeeze(SS(1,:,:)), squeeze(CIs(1,:,:)), squeeze(IS(1,:,:)), 'FaceColor', [0.5 0.5 1], ...
        'FaceAlpha', 0.5, 'EdgeAlpha', 0.0 )
    hold on
    surf(SS(:,:,1), CS(:,:,1), CSIs(:,:,1), 'FaceColor', [1 1 1]*1.0, ...
        'FaceAlpha', 1.0, 'EdgeAlpha', 0.0 )
    hold off
    axis(axs2), xlabel('SSc'), ylabel('C')
    hold on, plot3(St,Ct,os,'g.',St(1),Ct(1),os,'g^','MarkerSize',9), hold off, view(2)
    
    subplot(3,2,5)
    surf(squeeze(SS(1,:,:)), squeeze(CIs(1,:,:)), squeeze(IS(1,:,:)), 'FaceColor', [0.5 0.5 1], ...
        'FaceAlpha', 0.5, 'EdgeAlpha', 0.5 )
    hold on
    surf(squeeze(SIs(:,1,:)), squeeze(CS(:,1,:)), squeeze(IS(:,1,:)), 'FaceColor', [0.5 1 0.5], ...
        'FaceAlpha', 0.5, 'EdgeAlpha', 0.5 )
    surf(SS(:,:,1), CS(:,:,1), CSIs(:,:,1), 'FaceColor', [1 1 1]*0.8, ...
        'FaceAlpha', 0.5, 'EdgeAlpha', 0.0 )
    plot3(ss, cs, csis, 'k-', 'LineWidth', 0.3)
    hold off
    axis(axs), xlabel('SSc'), ylabel('C'), zlabel('I')
    hold on, plot3(St,Ct,It,'g.',St(1),Ct(1),It(1),'g^'), hold off
end


function dX = fd3( t, X, Slim0, Clim0, Ilim0, Sslope, Cslope, Islope )
%----------------------------------------------------------------------------------------------
% 3D vector field
%----------------------------------------------------------------------------------------------
% Basic shape of equilibrium curve for x-direction vectors (boundary between -> and <- vectors)
%
%    ^  |       Slim0 is the upper bound in S
%    |  |       
%  I | /        Sslope is the initial slope
%    |/
%    -------->
%           S
%----------------------------------------------------------------------------------------------
% Basic shape of equilibrium curve for y-direction vectors (boundary between ^ and v vectors)
%
%    ^  _____   Ilim0 is the upper bound in I
%  I | /        Islope is the initial slope
%    |/
%    -------->
%           S
%----------------------------------------------------------------------------------------------
Sx = X(1); Cy = X(2); Iz = X(3);              % vector components

dX    = 0*X;                                  % initialize
dX(1) = fSIx( Sx,     Iz, Slim0, Sslope );                   % S-component
dX(2) = fCIx( Cy,     Iz, Clim0, Cslope );                   % C-component
dX(3) = fSIy( Sx, Cy, Iz, Ilim0, Islope );                   % I-component

return


function [ SIx,Sx ] = fSIx( S, I, Slim0, Sslope )
%----------------------------------------------------------------------------------------------
% Vector field in the S (horizontal) direction
%----------------------------------------------------------------------------------------------
% Slims = a*(1 - exp(-b*I)), initial slope = a*b = Sslope*SScCoeff*Smult
% dS/dx = Slims - S
%
%       Slim0
%    ^-> | <-
%    |  /       
%  I | /        
%    |/
%    -------->
%           S
%----------------------------------------------------------------------------------------------
% Slims  = @(y) Slim0*(1-exp(-1/Sslope/Slim0*SScCoeff*y))*Smult;
Slims  = @(y) Slim0*(1-exp(-1/Sslope/Slim0*y));
Sx     = Slims(I);                           % equilibrum curve in S,I space
SIx    = Sx - S;                             % difference from Sequil and S: Sx - S
return


function [ CIx,Cx ] = fCIx( C, I, Clim0, Cslope )
%----------------------------------------------------------------------------------------------
% Vector field in the C (horizontal) direction
%
%   ^            ^  <- / -> |
%   |  ____    I |    /     |
%   | /          | __/y0    |
%   |/           |/         |
%   ------->     ------------->
%                0      C   1
%----------------------------------------------------------------------------------------------
% Clims = y0 + (I/0.8 - y0^3)^(1/3)
% 
%----------------------------------------------------------------------------------------------
% i = (c-1/2)^3 + b,  b = (1/2)^3
% (i-b)^1/3 = c - 1/2
% c = (i-b)^1/3 + 1/2
%   = (i-(1/2)^3)^1/3 + 1/2
%----------------------------------------------------------------------------------------------
% I3     = (10/ICmult)*I.^4 + I*0.0;           % change from 10 to 100
% y0     = 0.6;
% Clims  = @(y) (abs(y-y0^3).^(1/3).*sign(y-y0^3) + y0);
% Cx     = Clims(I3);                          % equilibrium curve in C,I space
% Cx     = Clims(I/0.8)*ICmult;
Clims  = @(y) Clim0*(1-exp(-1/Cslope/Clim0*y));
Cx     = Clims(I);
CIx    = C .* (C - Cx);                        % make it zero on vertical axis
% CIx  = C .* (C - Cx).*(1 - C);               % make it zero on vertical axes

return


function [ SIy,Iy ] = fSIy( S, C, I, Ilim0, Islope )
%----------------------------------------------------------------------------------------------
% Vector field in the I (vertical) direction
%----------------------------------------------------------------------------------------------
% Ilims = a*(1-exp(-b*(S+C)))
%----------------------------------------------------------------------------------------------
% Ilims  = @(x) Ilim0*(1-exp(-Islope/Ilim0*x))*Imult;
Ilims  = @(x) Ilim0*(1-exp(-Islope/Ilim0*x));
SC     = (S+C)/1;
Iy     = Ilims(SC);                          % equilibrium curve in S+C,I space
SIy    = Iy - I;                             % here both S and C have additive effect on I

return


function [ Slim0,Clim0,Ilim0,Sslope,Cslope,Islope,S0,C0,I0,doSave ] = setParam( Nargin,handles )
%----------------------------------------------------------------------------------------------
% Set parameters
%----------------------------------------------------------------------------------------------
% doReset -> set Ssaved = [] and return
% doSave  -> save current S into Ssaved, save last point and return
% else    -> run and append to Ssaved
%----------------------------------------------------------------------------------------------
global St Ct It t Ssaved Csaved Isaved tsaved Slast Clast Ilast
global iRuns

if Nargin
    h = handles.initialValues;
    doSave  = get(handles.save_button, 'Value');
    doReset = get(handles.reset_button,'Value');
    doClear = get(handles.clear_button,'Value');      % Clear also sets Reset
    doInitialValues = h.UserData;                     % user values in initialValues textbox
else
    doSave  = false; doReset = false; doClear = false; doInitialValues = false;
end

if isempty(iRuns),  doInitial = true; else, doInitial = false; end

initial0 = 2;
initialS = 2;                                         % 1-initialParam, 2-sliders, 3-textbox, 4-saved
if ~isempty(tsaved), initial0 = 4; end

if (doInitial || doClear || ~Nargin), initial0 = 1; initialS = 1; end
if (doInitialValues),                 initial0 = 1; initialS = 3; end
if (doReset),                         initial0 = 1;               end

if (doInitial || doClear), iRuns = 0; end

[ Slim0,Sslope,Ilim0,Islope,Clim0,Cslope,S0,I0,C0 ] = initialParam;

if doInitialValues                                   % user values in initialValues textbox
    hValues = str2num(h.String(h.Value,:));
    cValues = num2cell(hValues);
    [ S0,I0,C0 ] = cValues{2:4};        % get saved slopes
    [ Slim0,Ilim0,Clim0 ] = cValues{5:7};            % get saved limits
end
if initial0==2
    C0 = get(handles.C0, 'Value');
elseif initial0==4
    C0 = Csaved(end); I0 = Isaved(end); S0 = Ssaved(end);
end
if initialS==2
    Cslope      = get(handles.Cslope, 'Value');
    Sslope      = get(handles.Sslope, 'Value');
    Islope      = get(handles.Islope, 'Value');
elseif Nargin
    set(handles.Cslope,'Value',Cslope);
    set(handles.Sslope,'Value',Sslope);
    set(handles.Islope,'Value',Islope);
end
%----------------------------------------------------------------------------------------------
% Handle values
%----------------------------------------------------------------------------------------------
% Initial run, Clear, ~Nargin: initialParam          for 0,             '
% Initial parameters:          textbox               for 0,             '
% Other:                       initialParam or saved for 0, sliders for '
% Reset:                       initialParam          for 0, sliders for ', clear saved arrays [return]
% Save:                        save into saved arrays [return]
%----------------------------------------------------------------------------------------------
% Delete (local, deletes one line)
% Read   (local, populate textbox)
%----------------------------------------------------------------------------------------------
if doSave
    iRuns  = iRuns + 1;
    params = [Slast Ilast Clast Slim0 Ilim0 Clim0];
    paramsTxt = sprintf('%2i %.2f %.2f %.2f %.2f %.2f %.2f',iRuns,params);
    if isempty(h.String)
        h.String = paramsTxt;
    else
        set(h,'String',char(h.String,paramsTxt))            % append new row in the character vector
    end
    for i=1:size(h.String,1)
        fprintf('\n''%s'',...',h.String(i,:))
    end

    tsaved = t; Ssaved = St; Csaved = Ct; Isaved = It;
elseif doReset || doInitialValues
    tsaved = []; Ssaved = []; Csaved = []; Isaved = [];
    handles.C0.Value = C0;
end


function [ Slim00,Sslope,Ilim00,Islope,Clim00,Cslope,S0,I0,C0 ] = initialParam
%----------------------------------------------------------------------------------------------
% Initial values for the parameters
%----------------------------------------------------------------------------------------------
Krun     =  1;                      % 1(01), 2(00), 3(11), 4(11), 5(10)

Irun     =  1;                      % index for vector of parameters describing 2D equilibria
I0run    =  1;                      % index for initial value

if Krun==1
    [Irun,I0run] = deal(1,1);      % basic run
elseif Krun==2
    [Irun,I0run] = deal(2,1);
end

Slim0s = [ 0.8 0.8 ];  % set of values for Slim0 (see above) indexed by Irun
Sslopes= [ 0.5 0.5 ];  %                   Sslope

Ilim0s = [ 0.5 0.5 ];   % set of values for Ilim0 (see above) indexed by Irun
Islopes= [ 0.5 0.5 ];   %                   Islope

Clim0s = [ 0.6 0.6 ];   % set of values for Clim0 (see above) indexed by Irun
Cslopes= [ 0.5 0.5 ];   %                   Cslope

S0s    = [ 0.0 0.0 0.0 ];      % initial value for S indexed by I0run
C0s    = [ 0.2 0.3 0.4 ];      %                   C
I0s    = [ 0.0 0.0 0.0 ];      %                   I

SICs   = [ Slim0s;Sslopes;Ilim0s;Islopes;Clim0s;Cslopes ];     % matrix of parameters
SICc   = num2cell(SICs (:,Irun )');                            % make it into a cell

[ Slim00,Sslope,Ilim00,Islope,Clim00,Cslope ] = deal(SICc{:}); % deal the parameters

S0 = S0s(I0run); I0 = I0s(I0run); C0 = C0s(I0run);
