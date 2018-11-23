function IMMUNE_Phase(handles)
%----------------------------------------------------------------------------------------------
% 2D and 3D phase space analysis
%----------------------------------------------------------------------------------------------
% Initialize
%----------------------------------------------------------------------------------------------
global Slim0 Sslope Ilim0 Islope Clim0 Cslope SIClast SScCoeff
global Imult ICmult Smult
global S0 I0 C0 doContinue doSave doReset St Ct It t Ssaved Csaved Isaved tsaved
global iRuns tlast Slast Clast Ilast

set(0,'DefaultTextFontSize',12), set(0,'DefaultAxesFontSize',12)
lbl = {'SSc','C','I'};                        % plot labels

n   = 26;                                     % number of points in each coordinate

if ~nargin
    handles = [];
end

[ Imult,ICmult,Smult,SScCoeff ] = setParam(nargin,handles);

if nargin
    set(handles.C0_txt,    'String',sprintf('C0: %.2f',     C0    ))
    set(handles.Imult_txt, 'String',sprintf('Imult: %.2f',  Imult ))
    set(handles.ICmult_txt,'String',sprintf('ICmult: %.2f', ICmult))
    set(handles.Smult_txt, 'String',sprintf('Smult: %.2f',  Smult ))
end

if doSave, return, end

[ tnew,X ] = ode45( @fd3, [0 100], [S0;C0;I0] ); % solve ODE

Jt = find( max(abs(X')) < 1.2 );                 % only keep values less than 1.2

if isempty(tsaved)
    tlast = 0;
else
    tlast = tsaved(end);
end
t  = [tsaved;tnew(Jt)+tlast];              % truncate
St = [Ssaved;X(Jt,1)]; Ct = [Csaved;X(Jt,2)]; It = [Isaved;X(Jt,3)];
Slast = X(1,1);        Clast = X(1,2);        Ilast = X(1,3);

%----------------------------------------------------------------------------------------------
% 2D plotting: I, C, SSc
%----------------------------------------------------------------------------------------------
Imax = 0.8;
Is  = linspace(0,Imax,n);                   % 1D coordinates
Ss  = linspace(0,1,n);
Cs  = linspace(0,1,n);
[ SS    ] = meshgrid( Ss,Is );              % 2D coordinates
[ CS,IS ] = meshgrid( Cs,Is );
SIx = fSIx( SS, IS );                       % 2D vector fields
CIx = fCIx( CS, IS );
SICend= fSIy( SS, Ct(end), IS );
SISend= fSIy( St(end), CS, IS );
[ ~,Sx ] = fSIx( 0,  Is );                  % boundaries
[ ~,Cx ] = fCIx( 0,  Is );
[ ~,IyCend] = fSIy( Ss, Ct(end), 0  );
[ ~,IySend] = fSIy( St(end), Cs, 0  );

if ~nargin, figure(110), clf, set(gcf,'Name','basic 2D fields','color','white'), end

if nargin, axes(handles.axes1), else, subplot(2,2,1), end
hold off, p1a = plot(Ss,IyCend,'k-',Sx,Is,'b-'); hold on
axis tight, xlabel('SSc'), ylabel('Immune response'), axis([0 1 0 Imax])

if nargin, axes(handles.axes2), else, subplot(2,2,2), end
hold off, p2a = plot(Cs,IySend,'k-',Cx,Is,'r-'); hold on
ax = gca; ax.XDir = 'reverse';
axis tight, xlabel('Cancer'), ylabel('Immune response'), axis([0 1 0 Imax])

if nargin, axes(handles.axes3), else, subplot(2,2,[3 4]), end
cla

js = 1:3:round(length(St)*0.99);

if nargin, axes(handles.axes1), else, subplot(2,2,1), end
p1b = plot(St(js),It(js),'g.','MarkerSize',9); plot(St(1),It(1),'g^','MarkerSize',9)
p1c = quiver(SS,IS,SIx,SICend);
legend([p1a;p1b;p1c],'I equilibrium','SSc equilibrium','disease trajectory','flow field',...
    'Location','NW')

if nargin, axes(handles.axes2), else, subplot(2,2,2), end
p2b = plot(Ct(js),It(js),'g.','MarkerSize',9); plot(Ct(1),It(1),'g^','MarkerSize',9)
p2c = quiver(CS,IS,CIx*8,SISend);
legend([p2a;p2b;p2c],'I equilibrium','C equilibrium','disease trajectory','flow field')

if nargin, axes(handles.axes3), else, subplot(2,2,[3 4]), end
plot(t(js),St(js),'b-',t(js),Ct(js),'r-',t(js),It(js),'k-','LineWidth',1)
ylabel('SSc, Cancer, Immune response'), xlabel('time'), hold on
axis([0 max(t) 0 0.8]), legend('Scleroderma','Cancer','Immune response')

% for j=1:9:round(length(St)*0.99)                                  % point-by-point plotting
%     if nargin, axes(handles.axes1), else subplot(2,2,1), end
%     plot(St(j),It(j),'ko',St(1),It(1),'ro'), hold on
%     if nargin, axes(handles.axes2), else subplot(2,2,2), end
%     plot(Ct(j),It(j),'ro',Ct(1),It(1),'ro'), hold on
%     if nargin, axes(handles.axes3), else subplot(2,2,[3 4]), end
%     plot(t(j),X(Jt(j),1),'ko',t(j),X(Jt(j),2),'ro',t(j),X(Jt(j),3),'go','LineWidth',1)
%     ylabel('\color{red}Cancer,\color{black}SSc,\color{green}Immune response'), xlabel('time'), hold on
%     axis([0 max(t)/2 0 0.8])
%     pause(0.00)
% end

%----------------------------------------------------------------------------------------------
% 3D plotting
%----------------------------------------------------------------------------------------------
[ SS,CS,IS ] = meshgrid( Ss,Cs,Is );        % coordinate change: SS wrt 2, CS wrt 1, IS wrt 3
SIx  = fSIx( SS, IS );                      % 3D vector fields
CIy  = fCIx( CS, IS );
CSIz = fSIy( SS, CS, IS );                  % simple sum of S and C used to determine change in I
[ ~,SIs ] = fSIx( 0, IS );
[ ~,CIs ] = fCIx( 0, IS );
[ ~,cSIs] = fSIy( SIs, CS,  0 );
[ ~,CsIs] = fSIy( SS,  CIs, 0 );
[ ~,CSIs] = fSIy( SS,  CS,  0 );
Js = IS>cSIs;
Jc = IS<CsIs;

ss   = [SS(  :,[1 end],1)   SS([1 end],:,1)'];
cs   = [CS(  :,[1 end],1)   CS([1 end],:,1)'];
csis = [CSIs(:,[1 end],1) CSIs([1 end],:,1)'];
os   = ones(size(t));

axs = [0 1 0 1 0 0.8]; axs2 = [0 1 0 1];

if ~nargin, fig = figure(220); fig.Name = 'overhead'; end

if nargin, axes(handles.axes5), else, subplot(3,2,1), end
surf(squeeze(SIs(:,1,:)), squeeze(CS(:,1,:)), squeeze(IS(:,1,:)), 'FaceColor', [0.5 1 0.5], ...
    'FaceAlpha', 0.5, 'EdgeAlpha', 0.5 )
% hold on
surf(SS(:,:,1), CS(:,:,1), CSIs(:,:,1), 'FaceColor', [0.5 1 0.5], ...
    'FaceAlpha', 0.5, 'EdgeAlpha', 0.0 ), hold on
plot3(ss, cs, csis, 'k-', 'LineWidth', 0.3)
hold off
axis(axs), xlabel('SSc'), ylabel('C'), zlabel('I')
hold on, plot3(St,Ct,It,'r.',St(1),Ct(1),It(1),'r^'), hold off

if nargin, axes(handles.axes6), else, subplot(3,2,3), end
surf(squeeze(SS(1,:,:)), squeeze(CIs(1,:,:)), squeeze(IS(1,:,:)), 'FaceColor', [0.5 0.5 1], ...
    'FaceAlpha', 0.5, 'EdgeAlpha', 0.5 )
hold on
surf(SS(:,:,1), CS(:,:,1), CSIs(:,:,1), 'FaceColor', [1 1 1]*0.8, ...
    'FaceAlpha', 0.5, 'EdgeAlpha', 0.0 )
plot3(ss, cs, csis, 'k-', 'LineWidth', 0.3)
hold off
axis(axs), xlabel('SSc'), ylabel('C'), zlabel('I')
hold on, plot3(St,Ct,It,'g.',St(1),Ct(1),It(1),'g^'), hold off

if nargin, axes(handles.axes4), else, subplot(3,2,6), end
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

if ~nargin
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

return

%----------------------------------------------------------------------------------------------
% Plot 3D
%----------------------------------------------------------------------------------------------
os = [1 1];
figure(210), clf, set(gcf,'Name','3D fields')    % plot 3D fields and trajectories
% quiver3(SS,CS,IS,SIx,CIy,CSIz), axis tight
% plot3(St,Ct,It,'ro')
axis([0 1 0 1 0 0.5]), grid on
xlabel('SSc'), ylabel('\color{red}C'), zlabel('I')
set(gca,'YColor','red')

for j=1:1:round(length(St)*0.7)
%   plot3(St,Ct,It,'ro')
    plot3(St(1:j),Ct(1:j),It(1:j),'ko')
    hold on
    plot3(St(j),Ct(j),It(j),'ko')
    axis([0 1 0 1 0 0.5])
    xlabel('SSc'), ylabel('\color{red}C'), zlabel('I')
    plot3([0 St(j)],Ct(j)*os,It(j)*os,'k-')
    plot3(St(j)*os,[0 Ct(j)],It(j)*os,'r-')
    plot3(St(j)*os,Ct(j)*os,[0 It(j)],'g-'), grid on
    set(gca,'YColor','red')
    pause(0.05)
    hold off
end

return


function dX = fd3( t, X )
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
dX(1) = fSIx( Sx,     Iz );                   % S-component
dX(2) = fCIx( Cy,     Iz );                   % C-component
dX(3) = fSIy( Sx, Cy, Iz );                   % I-component

return


function [ SIx,Sx ] = fSIx( S, I )
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
global Slim0 Sslope SScCoeff Smult

Slims  = @(y) Slim0*(1-exp(-1/Sslope/Slim0*SScCoeff*y))*Smult;
Sx     = Slims(I);                           % equilibrum curve in S,I space
SIx    = Sx - S;                             % difference from Sequil and S: Sx - S
return


function [ CIx,Cx ] = fCIx( C, I )
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
global Clim0 Cslope ICmult

% I3     = (10/ICmult)*I.^4 + I*0.0;           % change from 10 to 100
I3     = I;           % change from 10 to 100
Clims  = @(y) Clim0*(1-exp(-1/Cslope/Clim0*y));
Cx     = Clims(I3);                          % equilibrium curve in C,I space
% y0     = 0.6;
% Clims  = @(y) (abs(y-y0^3).^(1/3).*sign(y-y0^3) + y0);
% Cx     = Clims(I/0.8)*ICmult;
CIx    = C .* (C - Cx).*(1 - C);             % make it zero on vertical axes

return


function [ SIy,Iy ] = fSIy( S, C, I )
%----------------------------------------------------------------------------------------------
% Vector field in the I (vertical) direction
%----------------------------------------------------------------------------------------------
% Ilims = a*(1-exp(-b*(S+C)))
%----------------------------------------------------------------------------------------------
global Ilim0 Islope SScCoeff Imult

Ilims  = @(x) Ilim0*(1-exp(-Islope/Ilim0*x))*Imult;
SC     = (S+C)/1;
Iy     = Ilims(SC);                          % equilibrium curve in S+C,I space
SIy    = Iy - I;                             % here both S and C have additive effect on I

return


function [ Imult,ICmult,Smult,SScCoeff ] = setParam( Nargin,handles )
%----------------------------------------------------------------------------------------------
% Set parameters
%----------------------------------------------------------------------------------------------
% doReset -> set Ssaved = [] and return
% doSave  -> save current S into Ssaved, save last point and return
% else    -> run and append to Ssaved
%----------------------------------------------------------------------------------------------
global S0 I0 C0 SIClast doContinue doSave doReset St Ct It t Ssaved Csaved Isaved tsaved iRuns
global Slim0 Sslope Ilim0 Islope Clim0 Cslope Slast Clast Ilast

Krun     =  1;                      % 1(01), 2(00), 3(11), 4(11), 5(10)

Irun     =  2;                      % index for vector of parameters describing 2D equilibria
I0run    =  3;                      % index for initial value

if Krun==1
    [Irun,I0run] = deal(2,3);      % basic run
elseif Krun==2
    [Irun,I0run] = deal(4,3);
elseif Krun==3
    [Irun,I0run] = deal(4,4);
elseif Krun==4
    [Irun,I0run] = deal(3,3);
elseif Krun==5
    [Irun,I0run] = deal(5,3);
end

Slim0s = [ 0.8 0.8 0.8 0.8 0.01 ];  % set of values for Slim0 (see above) indexed by Irun
Sslopes= [ 0.9 0.3 0.3 0.3 0.9  ];  %                   Sslope

Ilim0s = [ 0.5 0.8 0.8 0.8 0.4 ];   % set of values for Ilim0 (see above) indexed by Irun
Islopes= [ 0.6 0.7 0.3 0.2 0.3 ];   %                   Islope

Clim0s = [ 0.6 0.8 0.4 0.8 0.4 ];   % set of values for Clim0 (see above) indexed by Irun
Cslopes= [ 0.2 0.4 0.4 0.3 0.8 ];   %                   Cslope

S0s    = [ 0.2 0.7 0.0 0.0  ];      % initial value for S indexed by I0run
C0s    = [ 0.5 0.6 0.4 0.20 ];      %                   C
I0s    = [ 0.1 0.3 0.0 0.0  ];      %                   I

SICs   = [ Slim0s;Sslopes;Ilim0s;Islopes;Clim0s;Cslopes ];  % matrix of parameters
SIC0s  = [ S0s;I0s;C0s ];                                   % matrix of initial conditions

SICc   = num2cell(SICs (:,Irun )');                         % make it into a cell
[ Slim0,Sslope,Ilim0,Islope,Clim0,Cslope ] = deal(SICc{:}); % deal the parameters

Slim0  = 0.8;
Sslope = 0.3;
Ilim0  = 0.8;
Islope = 0.7;
Clim0  = 0.8;
Cslope = 0.4;

ICmult0  = 0.60; Smult0   = 1.20; Imult0   = 0.70; 
%----------------------------------------------------------------------------------------------
% Handle values
%----------------------------------------------------------------------------------------------
if Nargin
    h = handles.initialValues;
    doSave     = get(handles.save_button,'Value');
    doReset    = get(handles.reset_button,'Value');
    doClear    = get(handles.clear_button,'Value');
    if isempty(iRuns) || doClear, iRuns = 0; end
    if h.UserData
        hValues = str2num(h.String(h.Value,:));
        SIClast = num2cell(hValues(2:4));
        cValues = num2cell(hValues);
        [ Smult,Imult,ICmult,SScCoeff ] = cValues{5:8};
        set(handles.ICmult,'Value',ICmult)
        set(handles.Smult, 'Value',Smult)
        set(handles.Imult, 'Value',Imult)
    else
        C0         = get(handles.C0,'Value');
        ICmult     = get(handles.ICmult,'Value');
        Smult      = get(handles.Smult, 'Value');
        Imult      = get(handles.Imult,'Value');
    end
else
    C0         = SIC0s(3,I0run);
    doSave     = false;
    doReset    = true;
    ICmult     = ICmult0; Smult = Smult0; Imult = Imult0;
end

SScCoeff =  1.00;                   % not used

if doSave
    iRuns  = iRuns + 1;
    params = [Slast Ilast Clast Smult Imult ICmult SScCoeff];
    paramsTxt = sprintf('%2i %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f',iRuns,params);
    set(handles.initialValues,'String',char(handles.initialValues.String,paramsTxt))
    for i=1:size(handles.initialValues.String,1)
        fprintf('\n''%s'',...',handles.initialValues.String(i,:))
    end

    tsaved = t; Ssaved = St; Csaved = Ct; Isaved = It; doContinue = true;
    SIClast = {St(end),It(end),Ct(end)};                    % save the end point
    return
elseif doReset
    tsaved = []; Ssaved = []; Csaved = []; Isaved = []; doContinue = false;
    SIClast = num2cell(SIC0s(:,I0run)');                     % initial values for new runs
    ICmult     = ICmult0; Smult = Smult0; Imult = Imult0;
    if Nargin
        handles.ICmult.Value = ICmult; handles.Smult.Value = Smult; handles.Imult.Value = Imult; 
    end
end

[ S0,I0 ] = deal(SIClast{1:2});                           % deal the initial conditions

if ~isempty(tsaved)
    C0 = deal(SIClast{3});
end

%----------------------------------------------------------------------------------------------
% .33 .71 .77  S .
%         .59  . .
%     .67      S C
%----------------------------------------------------------------------------------------------
