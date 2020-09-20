function varargout = anisotropic_ndgrid(xmin, xmax, nbsamples, Ndim, mode, pow)
%% anisotropic_ndgrid : function to build an anisotropic grid in dimension Ndim
% Depending on the mode choice -'corners' / 'centre'- the size of cells
% increases -respectively decreases- as a power of Ndim.
%
% Author & support : nicolas.douillet (at) free.fr, 2017-2020.
%
%
% Syntax
%
% [X, Y] = anisotropic_ndgrid(xmin, xmax);
% [X, Y] = anisotropic_ndgrid(xmin, xmax, nbsamples);
% [X, Y, Z,...] = anisotropic_ndgrid(xmin, xmax, nbsamples, Ndim);
% [X, Y, Z,...] = anisotropic_ndgrid(xmin, xmax, nbsamples, Ndim, mode);
% [X, Y, Z,...] = anisotropic_ndgrid(xmin, xmax, nbsamples, Ndim, mode, pow);
%
%
% Description
% [X, Y] = anisotropic_ndgrid(xmin, xmax) returns 2 matrices of
% samples X, Y (Ndim = 2 by default), of size size(X) = [49 49]. Default
% mode is 'corners' and default power is pow = 2.
%
% [X, Y] = anisotropic_ndgrid(xmin, xmax, nbsamples) uses nbsamples
% samples. For grid construction matters, nbsamples must be an odd number.
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
%
% See also : NDGRID, NDTRIGRID, MESHGRID, LINSPACE, LOGSPACE
%
%
% Input Arguments
%
% - xmin : real scalar double, the minimum value of the grid in each dimension.
% - xmax : real scalar double, xmax > xmin, the maximum value of the grid in each dimension.
% - nbsamples : positive integer scalar, the number of samples. nbsamples >=3.
% - Ndim : positive integer scalar, the number of dimensions. dimensions. Ndim >= 1.
% - mode : character string in the set{'corners','centre','center'}. Case insensitive.
% - pow : positive real scalar double, the exponent of the power function.
%
%
% Output Argument
%
% - varargout : [X, Y, Z, ...] the resulting sample matrices where
%               size(X) = size(Y) = size(Z) = nbsamples*ones(1,Ndim).
%
%
% Example #1 : 2D example
%
% [X, Y] = anisotropic_ndgrid(1, 8);
% figure;
% for i=1:size(X,1)
%   line(X(i,:), Y(i,:)), hold on;
% end
%
% for j=1:size(Y,2)
%   line(X(:,j), Y(:,j)), hold on;
% end
%
%
% Example #2 : 3D example + 3D rescale
% using different values for begin and end
% [X, Y, Z] = anisotropic_ndgrid(-pi, pi, 16, 3, 'corners', 2);
% zmin = min(min(min(Z)));
% Z = Z - sign(zmin)*abs(zmin);
% zmax = max(max(max(Z)));
% Z = Z/zmax;
% Z = 2*Z-1;
% figure;
% plot3(X(:),Y(:),Z(:), 'b.');
% axis equal, axis tight;


%% Input parsing
assert(nargin >= 2, 'Error : not enough input arguments.');

if nargin > 1
    
    assert(xmax > xmin, 'Error : maximum sample value must be greater than mimimum sample value. ');
    
    if nargin == 2
        
        nbsamples = 49;
        Ndim = 2;
        mode = 'corners';
        pow = 2;
        
    elseif nargin > 2
        
        assert(nbsamples > 0, 'Error : nbsamples must be a positive integer.');
        
        if nbsamples < 5
            
            warning('Using anisotropic_ndgrid with nbsamples < 5 is not relevant. Prefer the use of ndgrid in this case.');
            
        end
        
        if nargin == 3
            Ndim = 2;
            mode = 'corners';
            pow = 2;
            
        elseif nargin > 3
            
            assert(Ndim > 0, 'Error : Ndim must be a positive integer.');
            
            if nargin == 4
                
                mode = 'corners';
                pow = 2;
                
            elseif nargin > 4
                
                assert(strcmpi(mode, 'corners') == 1 || ...
                       strcmpi(mode, 'centre') == 1 || ...
                       strcmpi(mode, 'center') == 1, 'Error : mode must be either ''corners'', ''centre'', or ''center''.');
                   
                if nargin == 5
                    
                    pow = 2;
                    
                elseif nargin > 5
                    
                    assert(pow > 0, 'Error : pow must be a positive real number.');
                    
                end
            end
        end
    end
end


%% Body
% m-power spaced sampling vector
negativalues = 0;

if xmin < 0
    negativalues = 1;
    umin = xmin;
    xmax = xmax + abs(xmin);
    xmin = 0;    
end

if strcmpi(mode,'corners') == 1
    
    u = linspace(nthroot(xmin, pow),nthroot(xmax, pow),1+floor(nbsamples/2));
    sample_u = 0.5*(abs(xmin) + sign(u).*abs(u).^pow);    
    
elseif strcmpi(mode,'centre') == 1 || strcmpi(mode,'center') == 1
    
    u = linspace(xmin^pow,xmax^pow,1+floor(nbsamples/2));
    sample_u = 0.5*(abs(xmin)+sign(u).*abs(u).^(1/pow));
    
end

rsample_u = fliplr(2*sample_u(1,end)-sample_u);
sample_u  = [sample_u rsample_u(1,2:end)];

if negativalues
    sample_u = sample_u + umin;
end

% Grid creation
ndgrid_varargin = cell(1,Ndim);
varargout =       cell(1,Ndim);

for i = 1:Ndim
    ndgrid_varargin{i} = sample_u;
end

[varargout{:}] = ndgrid(ndgrid_varargin{:});


end % anisotropic_ndgrid