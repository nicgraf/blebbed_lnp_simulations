function set_cbrt_colorbar(ax, cmax_percent)
% SET_CBRT_COLORBAR  Configure a colorbar for cube-root-transformed
% contrast data (i.e. data plotted as nthroot(contrast, 3)), labeling
% the ticks with the ORIGINAL percentage values rather than the
% transformed ones -- so the display gets the dynamic-range compression
% (an extreme outlier pixel won't wash out the rest of the image) while
% staying physically interpretable at a glance.
%
%   ax           : axes handle whose image was plotted as nthroot(contrast,3)
%   cmax_percent : the max |contrast| percentage actually present in the
%                   data being shown -- used to pick a sensible tick range
%
% Usage:
%   imagesc(ax, x, x, nthroot(contrast, 3), nthroot([-cmax cmax], 3));
%   set_cbrt_colorbar(ax, cmax);

decade_ticks = [1 2 5];
all_ticks = [];
for e = 0:4
    all_ticks = [all_ticks, decade_ticks * 10^e]; %#ok<AGROW>
end
all_ticks = unique(all_ticks);
all_ticks = all_ticks(all_ticks <= cmax_percent);

if isempty(all_ticks)
    all_ticks = cmax_percent;   % fallback for very small ranges
end

tick_vals_percent = [-fliplr(all_ticks), 0, all_ticks];

cb = colorbar(ax);
cb.Ticks = nthroot(tick_vals_percent, 3);
cb.TickLabels = arrayfun(@(v) sprintf('%g%%', v), tick_vals_percent, 'UniformOutput', false);

end