function F = WicksellUniform(r,Rmin,Rmax)
%WICKSELLUNIFORM Compute the Wicksell's transform of a uniform distribution
%and return the folded Cumulative Density Function (CDF).
%
%  F = WICKSELLUNIFORM(R,MIN,MAX) uses R as interpolation points to
%   compute the CDF of the folded distribution, considering that the
%   underlying distribution is homogeneous between bounds MIN and MAX.
%
% Reference:
%   Depriester and Kubler (2019)    doi:10.5566/ias.2133
%
% See also WicksellHistogram
	F=ones(size(r));
	if Rmax<=Rmin
		error('Rmax must be greater than Rmin')
	end
	gamma=Rmax*sqrt(Rmax^2-r.^2)-r.^2.*log(Rmax+sqrt(Rmax^2-r.^2));
	F(r<Rmin)=1-(gamma(r<Rmin)+r(r<Rmin).^2.*log(Rmin+sqrt(Rmin^2-r(r<Rmin).^2))-Rmin*sqrt(Rmin^2-r(r<Rmin).^2))/(Rmax^2-Rmin^2);
	F(Rmin<=r & r<Rmax)=1-(gamma(Rmin<=r & r<Rmax)+r(Rmin<=r & r<Rmax).^2.*log(r(Rmin<=r & r<Rmax)))/(Rmax^2-Rmin^2);
end

