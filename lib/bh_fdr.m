function q = bh_fdr(p)
%BH_FDR Benjamini-Hochberg adjusted P values (NaN ignored).
q = nan(size(p));
ok = ~isnan(p);
[ps, order] = sort(p(ok));
m = numel(ps);
qs = min(ps(:) * m ./ (1:m)', 1);
qs = flipud(cummin(flipud(qs)));
tmp = nan(m, 1);
tmp(order) = qs;
q(ok) = tmp;
end
