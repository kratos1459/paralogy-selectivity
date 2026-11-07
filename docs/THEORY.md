# Theory: noise ceiling and empirical-Bayes anchor calibration

Notation. For a compound i in publication d of a pair, observed selectivity y_id = s_i + o_d + e_id, with compound-specific selectivity s_i, publication offset o_d ~ (0, tau^2) and residual error e_id ~ (0, sigma^2), all independent.

## 1. Noise ceiling
For a compound measured in two publications, y_id - y_id' = (o_d - o_d') + (e_id - e_id'), so s_i cancels and Var(y_id - y_id')/2 = tau^2 + sigma^2 =: V_x (per measurement). A model that predicts s_i perfectly still leaves o_d + e_id unexplained on a new publication, so
R2_max = 1 - V_x / V_tot, with V_tot = Var(s) + V_x.
Scaffold-split scores above this bound require information about o_d, i.e. leakage of publication offsets.

## 2. Plain m-anchor offset
Let the model's residual on compound i of publication d be r_id = o_d + e_id (assuming the model captures s_i). With m anchors, the plain offset estimate is o_hat = mean of r over anchors = o_d + e_bar, Var(e_bar) = sigma^2/m. Error on a non-anchor compound: o_d - o_hat + e = -e_bar + e. MSE_plain = sigma^2 (1 + 1/m).

## 3. Empirical-Bayes offset
Shrink toward the mean offset mu: o_EB = mu + w (r_bar - mu), w = tau^2 / (tau^2 + sigma^2/m) (the posterior-mean weight). Then Var(o_d - o_EB) = (1-w)^2 tau^2 + w^2 sigma^2/m = w sigma^2/m, so
MSE_EB = sigma^2 (1 + w/m), and the MSE reduction over the plain offset is (1 - w) sigma^2 / m.
Because w < 1, EB always improves on the plain offset in expectation; the gain is largest when m is small and tau^2 is small relative to sigma^2/m, and vanishes as m grows (w -> 1). With no anchors the best estimate is mu (w = 0).

## 4. Empirical check (code/12)
tau^2 and sigma^2 are estimated from publication-grouped out-of-fold residuals on training data; sigma^2 for the test set from within-publication residual variance. Predicted vs observed MSE across 47 pairs: Pearson 0.95 to 0.99 (m = 1..5); predicted gains modestly exceed observed (e.g. 0.137 vs 0.117 at m = 1), consistent with imperfect model residual structure and estimation error in tau^2.

## 5. Why out-of-bag residuals fail
Random-forest out-of-bag residuals still use same-publication compounds in other trees, so they understate o_d and hence tau^2; w is underestimated and the offset is over-shrunk, which explains the negative result in code/07.
