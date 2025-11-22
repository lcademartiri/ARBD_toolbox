function isFlat = tost_flat_slope(y, W, delta, alpha)
  % Two‐One‐Sided‐Test that slope ∈ [–delta, +delta] at level alpha
  y_win = y(end-W+1:end);
  t_win = (1:W)';
  X = [ones(W,1), t_win];
  [b, ~, ~, ~, stats] = regress(y_win, X, alpha);
  slope = b(2);
  s2    = stats(4);
  XtX_inv = inv(X' * X);
  se_b    = sqrt(s2 * XtX_inv(2,2));
  df    = W - 2;

  t1 = (slope + delta) / se_b;       % slope > -delta
  p1 = 1 - tcdf(t1, df);
  t2 = (delta - slope) / se_b;       % slope < +delta
  p2 = 1 - tcdf(t2, df);

%   fprintf('Slope = %.3g; testing in [%.3g,%.3g]\n', slope, -delta, delta);
%   fprintf('p1 (>-δ): %.4f, p2 (<+δ): %.4f\n', p1, p2);

  isFlat = (p1 < alpha) && (p2 < alpha);
%   if isFlat
%     fprintf('✅ Slope is practically zero (within ±%.3g) at %.1f%% confidence.\n', delta, 100*(1-alpha));
%   else
%     fprintf('❌ Cannot declare flat: slope outside ±%.3g at %.1f%% level.\n', delta, 100*(1-alpha));
%   end
end
