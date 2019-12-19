function F = WicksellHistogram(r,freq,edges)
%WICKSELLHISTOGRAM Compute the Wicksell's transform of a finite histogram
%and return the folded Cumulative Density Function (CDF).
% 
%   F = WICKSELLHISTOGRAM(R,FREQ,EDGES) uses R as interpolation points to
%   compute the CDF of the folded distribution, considering that the
%   underlying distribution is given by an histogram with frequencies FREQ
%   and bins defined by EDGES.
%
% Reference:
%   Depriester and Kubler (2019)    doi:10.5566/ias.2133
%
% See also WicksellUniform, autoSaltykov
	freq=freq(:)';
	if length(freq)~=length(edges)-1
		error('The size of the freq must be N-1 with length(edges)=N')
	end
	N=length(freq);
	mid_points=(edges(2:N+1)+edges(1:N))/2;
	n_int=length(r);
	Fk=zeros(n_int,N);
	for k=1:N
		Fk(:,k)=WicksellUniform(r,edges(k),edges(k+1))*mid_points(k)*freq(k);
	end
	E=sum(mid_points.*freq);
	F=sum(Fk,2)'/E;
end