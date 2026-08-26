function V = bleb_volume(d_core, d_bleb, neck_overlap)
% BLEB_VOLUME  Union volume of a blebbed LNP's two intersecting spheres,
% NOT the naive sum (which would double-count the overlapping
% lens-shaped region).
%
% Standard sphere-sphere intersection volume formula (Wolfram MathWorld):
%   V_overlap = (pi/(12*d)) * (r1+r2-d)^2 * (d^2 + 2*d*(r1+r2) - 3*(r1-r2)^2)
% where d = center-to-center distance. Since center_dist was built as
% r1+r2-neck_overlap when the mesh was generated, (r1+r2-d) here is
% exactly neck_overlap.

r1 = d_core / 2;
r2 = d_bleb / 2;
d  = r1 + r2 - neck_overlap;   % center-to-center distance

V1 = (4/3) * pi * r1^3;
V2 = (4/3) * pi * r2^3;

if d >= r1 + r2
    V_overlap = 0;   % spheres don't actually touch
elseif d <= abs(r1 - r2)
    V_overlap = (4/3) * pi * min(r1, r2)^3;   % one fully contains the other
else
    V_overlap = (pi / (12*d)) * (r1 + r2 - d)^2 * ...
        (d^2 + 2*d*(r1 + r2) - 3*(r1 - r2)^2);
end

V = V1 + V2 - V_overlap;

end