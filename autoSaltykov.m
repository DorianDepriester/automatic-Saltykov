function [freq, centers, edges, param] = autoSaltykov(r, varargin)
%AUTOSALTYKOV Find the best parameters for performing the Saltykov method,
%with respect to the Cramer-von Mises (CvM) criterion.
%
%	[FREQ,CENTERS] = AUTOSALTYKOV(R) performs the Saltykov method on the
%	sample R, using different numbers of bins (ranging in 10 to 25) and 
%	different methods (namely 'lower', 'center' and 'upper') and returns 
%   the frequencies and bin centers of the best candidate, with respect to 
%   the CvM criterion.
%
%	[FREQ,CENTERS,EDGES] = AUTOSALTYKOV(...) also returns the bin edges.
%
%	[FREQ,CENTERS,EDGES,PARAM] = AUTOSALTYKOV(...) returns the parameters
%	used for minimization of the CvM criterion. PARAM is a structure array
%	comprising:
%       - the method (namely 'lower', 'center' or 'upper');
%       - the number of bins;
%       - the resulting value for the CvM criterion.
%
%	[...] = AUTOSALTYKOV(...,'range',RANGE) sets the range so that the
%	investigated values for the number of bins are comprised between
%	min(RANGE) and max(RANGE). Default value: [10 25].
%
%	[...] = AUTOSALTYKOV(...,'method',METHOD) constrains the minization to
%	the specified method (can be 'lower', 'center' or 'upper'). 
%	If METHOD=='all', all the aforementioned methods are investigated 
%	(default behaviour).
%
%	[...] = AUTOSALTYKOV(...,'integration_points',N) uses N points for
%	numerical integration when performing the CvM test (default: N=1000).
%
%	AUTOSALTYKOV(...) without output arguments produces a histogram bar 
%   plot of the results.
%
%	Reference:
%	  Depriester and Kubler (2019)    doi:10.5566/ias.2133
%
%	  See also Saltykov, WicksellHistogram

    p = inputParser;
    addParameter(p,'method','all')
    addParameter(p,'range',[10 25])
    addParameter(p,'integration_points',1000)
    parse(p,varargin{:});
    
    if ~strcmpi(p.Results.method,'all')
        methods={p.Results.method};
    else
        methods={'lower', 'center', 'upper'};
    end
    range=p.Results.range;
    nbins=min(range):max(range);
    npt_int=p.Results.integration_points;

    n=length(r);
    r_int=linspace(min(r),max(r),npt_int);  %	Interpolation points for numerical integration
    Fn=zeros(1,npt_int);                    %	Empirical CDF
    for j=1:npt_int
        Fn(j)=nnz(r<=r_int(j))/n;
    end
    
    w2_min=inf;
    for i=1:length(methods)
        for j=1:length(nbins)
            [f, c, e]=Saltykov(r, nbins(j), 'method', methods{i}); %	Perform the Saltykov method
            F=WicksellHistogram(r_int, f, e);    %	refold the distribution
            w2=n*trapz(F,(Fn-F).^2);                        %	Cramer-von Mises criterion
            if w2<w2_min
                w2_min=w2;
                method_min=methods{i};
                N_min=nbins(j);
                freq=f;
                centers=c;
                edges=e;
            end
        end
    end
    
    %%	Format output arguments
    if nargout==0
        bar(centers,freq);
    else
        param=struct('method',method_min,'N_bins',N_min,'Cramer_vonMises', w2_min);
    end
end

