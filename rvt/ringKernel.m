%% This function makes a ring shaped mask for the radial symmetry kernels
% one way is to go via the circle equation, the other is via the angle theta.
% this one takes the circle equation approach, its slower but more accurate
%
% written by Andre Gemeinhardt
% -- June 2020
%
% inputs: - r:     radius of the ring
%         - rmax:  size of the image
%                  if not given, its equal to r
% outputs: - img:  image with the circle

function [img] = ringKernel(r,rmax)

if nargin < 1
    r = 5;
end

if nargin < 2
    rmax = r;
end

imgsize = 2*rmax+1;
img = zeros(imgsize);

for i = 1:imgsize
    for j = 1:imgsize
        rij=sqrt((i-rmax-1)^2+(j-rmax-1)^2);
        if floor(rij)==r
            img(i,j) = 1;
        end
    end
end

img = img./sum(img(:)); % do not forget to normalize :-)