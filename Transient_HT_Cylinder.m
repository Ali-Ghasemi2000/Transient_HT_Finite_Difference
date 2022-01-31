clear, clc, close all

%   ----- BY : Mir Ali Ghasemi  ------

% This MATLAB script models the heat transfer from a cylinder exposed to a fluid.
% I used Finite Difference (Explicit) for cylindrical coordinates in order to derive formulas.
% Temperature matrix of the cylinder is plotted for all time steps.
% Three points are of interest: T(0,0,t), T(r0,0,t), T(0,L,t)
% Finally, a video of changing temp is generated.

%% Known

k = 17.4;               % conductivity
rho = 7900;             % density
c = 526;                % specific heat 
Ti = 600;               % Initial Temp
Tinf = 300;             % Infinite Temp
tend = 5*60;            % time of end of process (s)
r = 0.04;               % Radius of cylinder (m)
L = 0.3;                % Length of cylinder (m)
p = 7;                  % number of partitions
dt = 1;                 % time step      

%% Temp Matrix for h=500

% [TT, T00, T0L, Tr0]               = Temps(h, k, c, rho, Ti, Tinf, p, tend, r, L, dt)
[TT_500, T00_500, T0L_500, Tr0_500] = Temps(500, k, c, rho, Ti, Tinf, p, tend, r, L, dt);

figure
plot(T00_500, 'LineWidth',1)
hold on
plot(Tr0_500, 'LineWidth',1)
plot(T0L_500, 'LineWidth',1)
legend('T(0,0,t)', 'T(r0,0,t)','T(0,L,t)')
title("h = 500")
ylabel("Temperature (C)")
xlabel("Time (sec)")

%% T(0,0,t) plot for h=500 & h=1000
figure
plot(T00_500, 'LineWidth',1)
hold on

[TT_1000, T00_1000, T0L_1000, Tr0_1000] = Temps(1000, k, c, rho, Ti, Tinf, p, tend, r, L, dt);
plot(T00_1000, 'LineWidth',1)
title("T(0,0,t)")
ylabel("Temperature (C)")
xlabel("Time (sec)")
legend({'h=500', 'h=1000'})
hold off

%% Cylinder Temp Plotting
NR = 2*size(TT_1000{1},1);  % number of rows in the final temp matrix
NC = 2*size(TT_1000{1},2);  % number of columns in the final temp matrix
Tcyl = Ti*ones(NR, NC);

% Initialize video
myVideo = VideoWriter('Cylinder_Temp', 'MPEG-4'); % open video file
myVideo.FrameRate = 10;  % Video Frame Rate
open(myVideo)

% Plot in a loop and grab frames

hm = heatmap(Tcyl);
Ax = gca;
Ax.XDisplayLabels = nan(size(Ax.XDisplayData));
Ax.YDisplayLabels = nan(size(Ax.YDisplayData));
set(gcf,'Position',[600 50 500 1000])
hm.ColorLimits = [Tinf Ti];
colormap('jet')
ni = ceil(tend/dt);     % number of iterations                                                      

for i = 1:ni
    T = TT_500{i};
    Tcyl = [flip(T,2) T; flip([flip(T,2) T])];
%     pause(0.01) % Pause and grab frame
    frame = getframe(gcf); %get frame
    writeVideo(myVideo, frame);    
    hm.ColorData = Tcyl;
end
close(myVideo)

%% Temp Matrix Generator
function [TT, T00, T0L, Tr0] = Temps(h, k, c, rho, Ti, Tinf, p, tend, r, L, dt)

dx = min([L, r])/p;
nr = round(r/dx  + 1);  % number of radial elements
nl = round(L/dx  + 1); % number of vertical elements
Bi = h*dx/k;            % Biot number
Fo = k*dt/(rho*c*dx^2); % Fourier number
FB = Fo * Bi;         % Fourier x Biot
ni = ceil(tend/dt);     % number of iterations

T0 = Ti + zeros(nl, nr);% Initial Temp matrix    
Tnew = T0;              % new Temp matrix
TT = cell(1, ni);       % array of Temp matrices
TT{1} = T0;               
T_0_0_to_10min = zeros(1,ni);   % Temp at (0, 0) from time 0 to 10 min
T_0_0_to_10min(1) = Ti;         % Initialize Temp at (0, 0)
T_r0_0_to_10min = zeros(1,ni);  % Temp at (r0, 0) from time 0 to 10 min
T_r0_0_to_10min(1) = Ti;        % Initialize Temp at (r0, 0)
T_0_L_to_10min = zeros(1,ni);   % Temp at (0, L) from time 0 to 10 min
T_0_L_to_10min(1) = Ti;         % Initialize Temp at (0, L)

%             convection
%             -----------
%            | B   D   G |
% insulated  | A   E   H |  convection
%            | c   F   I |
%             ------------
%              insulated

% Temp Matrix Computation
for k = 2:ni
    for i=1:nl
        for j=1:nr
            n = j-1;
            if (i==1 && j==1) % B
                Fm = 2*Bi*Tinf + 4*T0(i,j+1) + 2*T0(i+1,j); % Fo Multiplier
                Tm = 1-2*FB-6*Fo;                         % T0 Multiplier
                
            elseif (i==1 && j==nr) % G
                Fm = n*Bi*Tinf/(n/2-1/8) + 2*Bi*Tinf + 2*T0(i+1,j) + ((n-0.5)/(n/2-1/8))*T0(i,j-1);
                Tm = 1-((n-0.5)/(n/2-1/8))*Fo-2*FB-(n/(n/2-1/8))*FB-2*Fo;
                
            elseif (i==nl && j==nr) % I
                Fm = 2*T0(i-1,j) + ((2*n-1)/(n-1/4))*T0(i,j-1) + 2*n*Bi*Tinf/(n-1/4);
                Tm = 1-((2*n-1)/(n-1/4))*Fo-2*Fo-FB*((2*n)/(n-1/4));
                
            elseif (i==nl && j==1) % C
                Fm = 2*T0(i-1,j) + 4*T0(i,j+1);
                Tm = 1-6*Fo;
                
            elseif i==1 % D
                Fm = 2*Bi*Tinf + 2*T0(i+1,j) + ((n+0.5)/n)*T0(i,j+1) + ((n-0.5)/n)*T0(i,j-1);
                Tm = 1-4*Fo-2*FB;
                
            elseif i==nl % F
                Fm = 2*T0(i-1,j) + ((n+0.5)/n)*T0(i,j+1) + ((n-0.5)/n)*T0(i,j-1);
                Tm = 1-4*Fo;
                
            elseif j==1 % A
                Fm = T0(i-1,j) + T0(i+1,j)+ 4*T0(i,j+1);
                Tm = 1-6*Fo;
                
            elseif j==nr % H
                Fm = Bi*n*Tinf + T0(i+1,j) + T0(i-1,j)+ ((n-0.5)/(n/2-1/8))*T0(i,j-1);
                Tm = 1-n*FB-2*Fo-Fo*((n-0.5)/(n/2-1/8));
                
            else % E (midlle)
                Fm = T0(i+1,j) + T0(i-1,j) + ((n-0.5)/n)*T0(i,j-1) + ((n+0.5)/n)*T0(i,j+1);
                Tm = 1-4*Fo;
                
            end
            if Tm <= 0
                Fo
                error('NOT Stable')
            end
            Tnew(i,j) = Fm*Fo + Tm*T0(i,j);
        end
    end
    T_0_0_to_10min(k) = Tnew(end,1);
    T_r0_0_to_10min(k) = Tnew(end,end);
    T_0_L_to_10min(k) = Tnew(1, 1);
    TT{k} = Tnew;
    T0 = Tnew;
end

T00 = T_0_0_to_10min;
T0L = T_0_L_to_10min;
Tr0 = T_r0_0_to_10min;

end