%% anisotropic_ndgrid
%
% Function to build an anisotropic grid in dimension Ndim.
% Depending on the mode choice -'corners' / 'centre'- the size of cells
% increases -respectively decreases- as a power of Ndim.
%
% Author & support : nicolas.douillet (at) free.fr, 2017-2020.
%
%% Syntax
%
% [X, Y] = anisotropic_ndgrid(xmin, xmax);
%
% [X, Y] = anisotropic_ndgrid(xmin, xmax, nbsamples);
%
% [X, Y, Z,...] = anisotropic_ndgrid(xmin, xmax, nbsamples, Ndim);
%
% [X, Y, Z,...] = anisotropic_ndgrid(xmin, xmax, nbsamples, Ndim, mode);
%
% [X, Y, Z,...] = anisotropic_ndgrid(xmin, xmax, nbsamples, Ndim, mode, pow);
%
%% Description
% [X, Y] = anisotropic_ndgrid(xmin, xmax) returns 2 matrices of
% samples X, Y (Ndim = 2 by default), of size size(X) = [49 49]. Default
% mode is 'corners' and default power is pow = 2.
%
% [X, Y] = anisotropic_ndgrid(xmin, xmax, nbsamples) uses nbsamples
% samples. For grid construction matters, nbsamples must be an odd number.
%
% In case input nbsamples is even, it will be replaced int he code by
% nbsamples+1 (first greter odd number). Default nbsamples value is 49. In
% cases nbsamples<=3, anisotropic_ndgrid behaviour is strictly equivalent
% to ndgrid's (no sampling step difference between corners and centre samples)
%
% [X, Y, Z,...] = anisotropic_ndgrid(xmin, xmax, nbsamples, Ndim) performs
%  in dimensions Ndim. 
%
% [X, Y, Z,...] = anisotropic_ndgrid(xmin, xmax, nbsamples, Ndim, mode)
% uses mode ('corners' or 'centre') choice for the grid layout.
% Default mode is 'corners'.
%
% [X, Y, Z,...] = anisotropic_ndgrid(xmin, xmax, nbsamples, Ndim, mode, pow)
% uses pow choice for the power function. Default power is pow = 2. In case
% pow = 1, anisotropic_ndgrid behaviour is strictly equivalent to ndgrid's
% (no sampling step difference between corners and centre samples)
%
%% See also
%
% <https://fr.mathworks.com/help/matlab/ref/ndgrid.html?s_tid=srchtitle ndgrid> |
% <https://fr.mathworks.com/matlabcentral/fileexchange/73041-n-dimension-regular-triangular-grid ndtrigrid> |
% <https://fr.mathworks.com/help/matlab/ref/meshgrid.html?s_tid=srchtitle meshgrid> |
% <https://fr.mathworks.com/help/matlab/ref/linspace.html?s_tid=srchtitle linspace> |
% <https://fr.mathworks.com/help/matlab/ref/logspace.html logspace> |
%
%% Input Arguments
%
% - xmin : real scalar double, the minimum value of the grid in each dimension.
%
% - xmax : real scalar double, xmax > xmin, the maximum value of the grid in each dimension.
%
% - nbsamples : positive integer scalar, the number of samples. nbsamples >=3.
%
% - Ndim : positive integer scalar, the number of dimensions. dimensions. Ndim >= 1.
%
% - mode : character string in the set{'corners','centre','center'}. Case insensitive.
%
% - pow : positive real scalar double, the exponent of the power function.
%
%% Output Argument
%
% - varargout : [X, Y, Z, ...] the resulting sample matrices where
%               size(X) = size(Y) = size(Z) = nbsamples*ones(1,Ndim).
%
%% Example #1
% 2D example

[X, Y] = anisotropic_ndgrid(1,8);
figure;
for i = 1:size(X,1)
    line(X(i,:), Y(i,:)), hold on;
end

for j = 1:size(Y,2)
    line(X(:,j), Y(:,j)), hold on;
end

%% Example #2
% 3D example + 3D rescale

[X, Y, Z] = anisotropic_ndgrid(-pi,pi,32,3,'centre',2);
figure;
plot3(X(:),Y(:),Z(:),'b.');
axis equal, axis tight;