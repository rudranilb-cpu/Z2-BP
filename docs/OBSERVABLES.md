# Gauge/dual observable dictionary and QPU priority

All principal observables require only global Z-basis shots, global X-basis shots, or
post-processing of those same bitstrings. No ancilla or Hadamard test is required.

## Exact operator mappings

| Gauge-theory quantity | Dual-spin operator | QPU basis/cost |
|---|---|---|
| Electric field on an internal direct link | \(\sigma^x_\ell=Z_rZ_{r'}\) | Z basis; two-bit parity |
| Electric-flux occupation | \(n_\ell=(1-Z_rZ_{r'})/2\) | same Z shots |
| Boundary flux in `fixed_exterior` | \((1-Z_r)/2\) | same Z shots; boundary-dependent |
| Total flux | \(\Phi=\sum_\ell n_\ell\) | sum of link estimates |
| Magnetic plaquette | \(B_p=X_r\) | X basis; one-bit mean |
| Wilson contour around plaquette set \(A\) | \(W(\partial A)=\prod_{p\in A}B_p=\prod_{r\in A}X_r\) | X basis; bitstring parity |
| Electric string along a dual path from \(r\) to \(s\) | \(\prod_{\ell\perp\gamma}\sigma^x_\ell=Z_rZ_s\) | Z basis; endpoint parity |
| Magnetic two-plaquette correlator | \(\langle X_rX_s\rangle\) | X basis |
| Connected flux correlation | \(\langle n_bn_{b'}\rangle-\langle n_b\rangle\langle n_{b'}\rangle\) | Z basis; bitstring post-processing |

In the charge-free bulk, different paths for the electric string differ by Gauss
operators and are equivalent. A single \(Z_r\) instead requires a path to an exterior
spin and is boundary-convention dependent. `sites.csv` retains one-site Z as a useful
dual-model diagnostic, but it is **not** a principal gauge observable for `open_dual`.

## Derived spatial observables

`links.csv` preserves horizontal and vertical bonds separately. A site-centered map is
also provided by assigning half of each internal link occupation to either endpoint
(and all of a fixed-exterior boundary link to its simulated endpoint), so that its sum
is exactly \(\Phi\).

For bond midpoint \(\mathbf r_b\) relative to the initial defect, the analysis uses

\[
\Delta n_b(t)=\langle n_b\rangle_{\rm loop}-\langle n_b\rangle_{\rm SCV}.
\]

Radial/angular profiles are binned means. Three front coordinates are stored separately:

\[
R_{\rm mean}=\frac{\sum_b|\Delta n_b|r_b}{\sum_b|\Delta n_b|},\qquad
R_{\rm rms}=\sqrt{\frac{\sum_b|\Delta n_b|r_b^2}{\sum_b|\Delta n_b|}},
\]

and the configured weighted radial quantile (default 90%). Effective bond
participation, \((\sum_b w_b)^2/\sum_b w_b^2\), and the eigenvalue anisotropy of the
weighted spatial second-moment tensor quantify spreading and directional distortion.
Absolute weights make the
front well-defined even when background subtraction generates a signed wake. Velocity
is a fit parameter over an explicitly reported pre-boundary time window. The plotting
driver fits both RMS and weighted-quantile fronts and records them as separate rows;
changing the front definition/window is a systematic check, not a statistical error bar.

## Loop persistence and spreading

The central elementary-loop survival operator is

\[
S_\square=\left\langle\prod_{b\ni r_c}n_b\right\rangle.
\]

For a bulk defect it is a diagonal five-dual-spin projector, equals one initially in
the loop state, and is computable from the same Z-basis shots. This is preferred over
the state-return probability, which would require an inverse circuit or overlap
protocol. Spreading is quantified by the radial profiles, front radii, total absolute
differential weight, and the full link map.

## Prioritized QPU set

1. **Every internal \(Z_rZ_{r'}\), hence every \(n_b\)**, from one global Z-basis
   setting. This produces local maps, \(\Phi_{\rm SCV}\), \(\Phi_{\rm loop}\),
   \(\Delta\Phi\), radial/angular profiles, fronts, and velocities.
2. **Central-loop survival and connected link-flux correlations**, extracted from the
   same Z bitstrings. No additional circuits are needed.
3. **Short electric strings \(Z_rZ_s\)** along axial/radial cuts, again from the same
   Z shots.
4. **Every local magnetic plaquette \(X_r\)** from one global X-basis setting; retain
   local maps and the SCV-subtracted signal.
5. **Small Wilson contours and connected magnetic correlations**, obtained as X-basis
   parities. Start with areas 1, 2, and 4; large-area Wilson products will become
   shot/noise limited quickly.

One-site Z, state fidelity, large arbitrary Pauli strings, and observables requiring
controlled/ancilla circuits are diagnostic or secondary and should not drive the main
Nighthawk allocation.
