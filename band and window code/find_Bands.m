function find_Bands(a, b, arc)
%takes x-axis (a) and y-axis (b) half-lengths, and Let-99 band arclength 
%arc and finds the corresponding Cartesian angles, T1 and T2, of the Let-99 
%band end points

%Erin Angelini, 5.25.18

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%parameter marking 60:40 egg-length point
tL1 = acos(0.2);

syms t ts te;
%arc length integrand
gammaprime = sqrt((a*sin(t))^2+ (b*cos(t))^2);

%upper band starting point
t1 = double(vpasolve(0.5*arc-int(gammaprime,t,ts,tL1)==0,ts));
%upper band ending point
t2 = double(vpasolve(0.5*arc-int(gammaprime,t,tL1,te)==0,te));

%find Cartesian coordinates of start and end points
x1 = a*cos(t1); y1 = b*sin(t1);
x2 = a*cos(t2); y2 = b*sin(t2);

%use atan2 to get the Cartesian angles
T1 = atan2(y1,x1);
T2 = atan2(y2,x2);

%adjust them accordingly
if T1 < 0
    T1 = T1 + 2*pi;
elseif T1 >= 2*pi
    T1 = T1 - 2*pi;
end

if T2 < 0
    T2 = T2 + 2*pi;
elseif T2 >= 2*pi
    T2 = T2 - 2*pi;
end

fprintf(strcat('const float_T start=', num2str(T1), ';', '\n', ...
'const float_T end=', num2str(T2), ';', '\n'))
end
