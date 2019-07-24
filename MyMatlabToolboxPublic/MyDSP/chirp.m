function y = chirp(t,f0,t1,f1,method,phi,quadtype)
%CHIRP  Swept-frequency cosine generator.
%   Y = CHIRP(T,F0,T1,F1) generates samples of a linear swept-frequency
%   signal at the time instances defined in array T.  The instantaneous
%   frequency at time 0 is F0 Hertz.  The instantaneous frequency F1
%   is achieved at time T1.  By default, F0=0, T1=1, and F1=100.
%
%   Y = CHIRP(T,F0,T1,F1,method) specifies alternate sweep methods.
%   Available methods are 'linear','quadratic', and 'logarithmic'; the
%   default is 'linear'.  Note that for a logarithmic-sweep, F0>=1e-6 is 
%   required and by default F0=1e-6.
%
%   Y = CHIRP(T,F0,T1,F1,method, PHI) allows an initial phase PHI to
%   be specified in degrees.  By default, PHI=0.
%
%   Y = CHIRP(T,FO,T1,F1,'quadratic',PHI,'concave') generates samples of
%   a quadratic swept-frequency signal whose spectrogram is a parabola with
%   its concavity in the positive frequency axis.
%
%   Y = CHIRP(T,FO,T1,F1,'quadratic',PHI,'convex') generates samples of
%   a quadratic swept-frequency signal whose spectrogram is a parabola with
%   its convexity in the positive frequency axis.
%
%   Default values are substituted for empty or omitted trailing input
%   arguments.
%
%   EXAMPLE 1: Compute the spectrogram of a linear chirp.
%     t=0:0.001:2;                    % 2 secs @ 1kHz sample rate
%     y=chirp(t,0,1,150);             % Start @ DC, cross 150Hz at t=1sec 
%     spectrogram(y,256,250,256,1E3); % Display the spectrogram
%
%   EXAMPLE 2: Compute the spectrogram of a quadratic chirp.
%     t=-2:0.001:2;                   % +/-2 secs @ 1kHz sample rate
%     y=chirp(t,100,1,200,'q');       % Start @ 100Hz, cross 200Hz at t=1sec 
%     spectrogram(y,128,120,128,1E3); % Display the spectrogram
%
%   EXAMPLE 3: Compute the spectrogram of a "convex" quadratic chirp
%     t= 0:0.001:1;                     % 1 second @ 1kHz sample rate
%     fo=25;f1=100;                     % Start at 25Hz, go up to 100Hz
%     y=chirp(t,fo,1,f1,'q',[],'convex');
%     spectrogram(y,256,200,256,1000);  % Display the spectrogram.
%
%   EXAMPLE 4: Compute the spectrogram of a "concave" quadratic chirp
%     t= 0:0.001:1;                      % 1 second @ 1kHz sample rate
%     fo=100;f1=25;                      % Start at 100Hz, go down to 25Hz
%     y=chirp(t,fo,1,f1,'q',[],'concave');
%     spectrogram(y,256,200,256,1000);   % Display the spectrogram.
%
%   EXAMPLE 5: Compute the spectrogram of a logarithmic chirp
%     t= 0:0.001:10;                   % 10 seconds @ 1kHz sample rate
%     fo=10;f1=400;                    % Start at 10Hz, go up to 400Hz
%     y=chirp(t,fo,10,f1,'logarithmic');
%     spectrogram(y,256,200,256,1000); % Display the spectrogram.
%
%   See also GAUSPULS, SAWTOOTH, SINC, SQUARE.

%   Author(s): D. Orofino, T. Krauss, 3/96
%   Copyright 1988-2009 The MathWorks, Inc.
%   $Revision: 1.1.6.1 $  $Date: 2009/04/21 04:28:19 $

%   References: 
%   [1] Agilent 33220A 20 MHz Function/Arbitrary Waveform Generator
%       Users Guide, Agilent Technologies, March 2002, pg 298

% Parse inputs, and substitute for defaults:
error(nargchk(1,7,nargin,'struct'));

if nargin<7, quadtype=[]; end
if nargin<6, phi=[]; end
if nargin<5, method=[]; end
if nargin<4, f1=[]; end
if nargin<3, t1=[]; end
if nargin<2, f0=[]; end
if isempty(phi), phi=0; end
if isempty(method), method='linear'; end
if isempty(f1), f1=100; end
if isempty(t1), t1=1; end
if isempty(f0), 
    if strncmpi('log',method,2)
        f0 = 1e-6; % Minimum output frequency in Agilent Function generator is 1 microHz
    else
        f0 = 0; 
    end
end

if (t1==0)
    error(generatemsgid('InvalidRange'),'End time T1 must be greater than 0');
end

% Parse the method string:
if length(f0)>1,
   %   Y = CHIRP(T,P) specifies a polynomial vector P for the
   %   instantaneous frequency trajectory of the chirp. 
   warnStr = 'Specifying a polynomial sweep vector as the second input argument will not be supported in future releases.';
   warning(generatemsgid('syntaxObselete'),warnStr)

   method='polynomial';
else
   % Set p=1 for linear, 2 for quadratic, 3 for logarithmic
   strs = {'linear','quadratic','logarithmic'};
   p = find(strncmpi(method,strs,length(method))); 
   if isempty(p),
      error(generatemsgid('InvalidParam'),'Unknown method selected.');
   elseif length(p)>1,
      error(generatemsgid('SignalErr'),'Ambiguous method selected.');
   end
   method = strs{p};
end

% Parse the  quadtype string and display an error message
% if quadtype is used with a sweep mode besides quadratic.
if ~isempty(quadtype) && ~strcmpi(method,'quadratic')
     error(generatemsgid('SigErr'),...
         'The "%s" mode is not defined for the %s sweep method.',...
        quadtype, method);
end


% Compute beta, phase and the output.

switch method
case 'polynomial'
    % Polynomial chirp
    y = cos( 2*pi * polyval(polyint(f0),t) );
    
case {'linear'},
    % Polynomial chirp: p is the polynomial order
    y = calculateChirp(f0,f1,t1,p,t,phi);
    
case {'quadratic'},
    
    % Determine the shape of the quadratic sweep - concave or convex
    if isempty(quadtype) && f1>f0
        quadtype = 'concave'; % Default for upsweep
    elseif isempty(quadtype) && f1<f0
        quadtype = 'convex'; % Default for downsweep.
    end
    
    % Polynomial chirp: p is the polynomial order
    % Compute the quadratic chirp output based on quadtype
    y = computequadchirp(f0,f1,t1,p,t,phi,quadtype);
    
    
case 'logarithmic',
    % Logarithmic chirp:
    if (f0 < 1e-6), error(generatemsgid('InvalidRange'),'F0 > 1e-6 is required for a log-sweep.'); end
    instPhi = t1/log(f1/f0)*(f0*(f1/f0).^(t/t1)-f0);
    y = cos(2*pi * (instPhi + phi/360));
    
end

%---------------------------------------------------------------------------
function yvalue = calculateChirp(f0,f1,t1,p,t,phi)
% General function to compute beta and y for both 
% linear and quadratic modes.
  
beta   = (f1-f0).*(t1.^(-p));
yvalue = cos(2*pi * ( beta./(1+p).*(t.^(1+p)) + f0.*t + phi/360));

%---------------------------------------------------------------------------
function y=computequadchirp(f0,f1,t1,p,t,phi,quadtype)
% Compute the quadratic chirp (upsweep or downsweep) for
% complex or concave modes.

% For the default 'concave-upsweep' and 'convex=downsweep' modes
% call calculateChirp without any changes to the input parameters.
% For the forced 'convex-upsweep' and 'concave-downsweep' call
% calculateChirp with f0 and f1 swapped and t= fliplr(-t)

% For 'convex-upsweep' and 'concave-downsweep' modes
if ((f0<f1) && strcmpi(quadtype,'convex')) || ((f0>f1) &&...
    strcmpi(quadtype,'concave'))
    t = fliplr(-t);
    ftemp=f0; f0 = f1; f1 = ftemp;
end

y = calculateChirp(f0,f1,t1,p,t,phi);

%----------------------------------------------------------------------------

% [EOF] chirp.m

