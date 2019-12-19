function [freq, centers, edges] = Saltykov(r, varargin)
%SALTYKOV Perform the Scheil-Schwartz-Saltykov method (called Saltykov 
% method for short) to unfold the distribution of apparent grain sizes in 
% 2D sections, giving the distribution of actual (3D) grain size, 
% considering grains as spheres.
%
%   [FREQ,CENTERS] = SALTYKOV(R) performs the Saltykov method on the
%   population of radii given in the array R and return the frequencies and
%   the centers of each bin of the unfolded distribution. By default, the
%   Saltykov method is applied using 15 bins and the upper limit of each
%   bin is used as the radius reference when computing the Wicksell's 
%   equation.
%
%   [...] = SALTYKOV(R, N) uses N bins.
%
%   [FREQ,CENTERS,EDGES] = SALTYKOV(...) also returns the bin edges.
%
%   [...] = SALTYKOV(...,'method',METHOD) sets the method to be used for 
%   choosing the reference radius. limit can be:
%       - 'lower'
%       - 'center'
%       - 'upper' (default)
%   
%   SALTYKOV(...) without output arguments produces a histogram bar plot of
%   the results.
%
% References:
%   Higgins (2000)                          doi:10.2138/am-2000-8-901
%   Lopez-Sanchez and Llana-Funez (2016)	doi:10.21105/joss.00863
%
% See also autoSaltykov

    if nargin>1 && isnumeric(varargin{1})
        nbins=varargin{1};
        varargin(1)=[];
    else
        nbins=15;
    end
    p = inputParser;
    addParameter(p,'method','upper')
    parse(p,varargin{:});

    [counts, centers]=hist(r,nbins);
    freq=counts/length(r);
    dr=centers(2)-centers(1);
    edges=[centers(1)-dr/2 centers+dr/2];
    edges(edges<0)=0;
    
    if strcmpi(p.Results.method,'lower')
        limits=edges(1:nbins);
    elseif strcmpi(p.Results.method,'center')
        limits=centers;
    elseif strcmpi(p.Results.method,'upper')
        limits=edges(2:nbins+1);
    else
        error('The ''method'' option can be either ''upper'', ''center'' or lower''.')
    end
    
    for i=1:nbins
        I=nbins+1-i;	% Start from upper classes
        R_I=edges(I+1);	% Use the maximum value in that class
        Pi=Wicksell(R_I, edges(I), edges(I+1));
        for j=1:I-1
            Pj=Wicksell(limits(I), edges(j), edges(j+1));
            Pnorm=Pj*freq(I)/Pi;
            freq(j)=max(freq(j)-Pnorm,0);   % Avoid negative values
            freq = freq/sum(freq);          % Keep the distribution normalized
        end
    end
    
    if nargout==0
        bar(centers,freq);
    end
end

function P=Wicksell(R,r1,r2)
    P=(sqrt(R^2-r1^2)-sqrt(R^2-r2^2))/R;
end

