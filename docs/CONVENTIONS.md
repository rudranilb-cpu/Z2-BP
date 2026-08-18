# Frozen physics and circuit conventions

## Gauge and dual Hamiltonians

On direct-lattice links,

\[
H_{\rm gauge}=-\Gamma\sum_\ell\sigma^x_\ell-K\sum_p B_p,
\qquad B_p=\prod_{\ell\in\partial p}\sigma^z_\ell .
\]

In the charge-free bulk the duality dictionary is

\[
X_r=B_p,\qquad \sigma^x_\ell=Z_rZ_{r'}
\]

for the direct link crossed by the internal dual bond \(\langle r,r'\rangle\).
Consequently,

\[
H_{\rm dual}=-K\sum_rX_r-\Gamma\sum_{\langle r,r'\rangle}Z_rZ_{r'},
\qquad n_\ell=\frac{1-Z_rZ_{r'}}{2}.
\]

The code uses `Up` as \(Z=+1\) and `Dn` as \(Z=-1\). SCV is all `Up`; the loop
state flips the configured dual site. A bulk flip frustrates exactly four bonds, so

\[
\Phi_{\rm SCV}(0)=0,\quad \Phi_{\rm loop}(0)=4,\quad
\Delta\Phi(0)\equiv\Phi_{\rm loop}(0)-\Phi_{\rm SCV}(0)=4.
\]

For an edge/corner defect in `open_dual`, the expected value is its internal degree;
the metadata records this instead of falsely demanding four.

## Open-boundary alternatives

`open_dual` is the model in the public scripts. An \(N_x\times N_y\) grid has

\[
N_b=(N_x-1)N_y+N_x(N_y-1)
\]

measured links and no boundary fields. It is an exact bulk duality with a rough/open
boundary interpretation, but it is not the full open direct lattice with all physical
boundary links.

`fixed_exterior` fixes exterior dual spins to \(+1\). Each missing neighbor contributes
\(-\Gamma Z_r\) to the Hamiltonian and \((1-Z_r)/2\) to the flux. Corners have two such
terms. The link count becomes \(N_b+2N_x+2N_y\), exactly the number of links in the
corresponding \((N_x+1)\times(N_y+1)\) direct vertex lattice.

## Gate signs and order

Write \(H=H_X+H_D\), where \(H_D\) contains internal ZZ and, when requested,
fixed-exterior Z terms. Since

\[
e^{-i(-KX)\,dt}=e^{+iKdtX}=R_x(-2Kdt),
\]

and similarly \(e^{+i\Gamma dtZZ}=R_{zz}(-2\Gamma dt)\), the negative Qiskit
rotation angles in the code are correct.

The frozen first-order reference is

\[
U_1(dt)=U_D(dt)U_X(dt)+O(dt^2).
\]

Because a gate list acts sequentially on the state, it must list X first and diagonal
gates second. The available orderings are:

| Configuration | Sequential list | Mathematical product | Role |
|---|---|---|---|
| `lie_x_zz` | X, diagonal | \(U_DU_X\) | frozen QPU/reference circuit |
| `lie_zz_x` | diagonal, X | \(U_XU_D\) | legacy mismatch reproduction only |
| `strang` | X/2, diagonal, X/2 | \(U_X(dt/2)U_D(dt)U_X(dt/2)\) | second-order reference |

All X terms commute with one another. All diagonal terms commute with one another.
Edge coloring changes simple-update scheduling and may change finite-D truncation
details, but not the ideal unitary.

## Phase-regime labels

For the square-lattice TFIM written as above, \(K\) is the transverse field and
\(\Gamma\) the Ising coupling. The thermodynamic critical ratio is approximately
\((K/\Gamma)_c=3.044330(6)\) ([Huang et al.](https://arxiv.org/abs/2005.10066)).
Thus `K=3, Gamma=1` is a near-critical point on the ordered/confined side, not a
generic deep strong-coupling Hamiltonian. “Strong-coupling vacuum” names the prepared
initial state; it does not describe every post-quench Hamiltonian.
