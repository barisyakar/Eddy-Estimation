function vortex=getvortex(i,s,t,firstx,firsty,aa,bb,k,eventTime)
% The input columns contain arrival time, lifetime in minutes, initial center
% coordinates, amplitude, and radius. Lifetimes are multiples of 15 minutes.
% k indexes occupied timestamps; eventTime is elapsed time in 15-minute blocks.

if nargin < 9
    eventTime = k - 1;
end

vortex.s=s(i);
vortex.t=t(i);
vortex.firstx=firstx(i);
vortex.firsty=firsty(i);
vortex.aa=aa(i);
vortex.bb=bb(i);
vortex.kk=k;
vortex.time=eventTime;
