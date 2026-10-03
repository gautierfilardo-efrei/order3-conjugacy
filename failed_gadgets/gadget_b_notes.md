# Gadget b track — hardness of P3 (x a x^2 = c, a an involution) through the conjugate b = x a x^{-1}

Conventions: permutations as image tuples, (pq)(j) = p(q(j)); t := x^3 = b c, b := x a x^{-1}.
Files: `gadget_b_track.py` (generator, intended-family constructor, complete enumeration, shape classifier),
`gadget_b_results.json` (all instances, solution counts, unintended shapes).

## 1. Structural lemmas (proved)

**L1 (bookkeeping).** For any solution x: the t-cycles are permuted by x (x commutes with t = x^3); for every
t-cycle B, #(a-endpoints in B) = #(b-endpoints in x(B)), because b = x a x^{-1} maps supp(a) onto supp(b)
cycle-wise. A t-cycle of length M with 3 ∤ M may be an x-cycle by itself (x = t^r, 3r ≡ 1 mod M); otherwise
x groups three t-cycles of equal length into one x-cycle.

**L2 (P3 is self-dual).** x a x^2 = c ⟺ y c y^2 = a with y = x^{-1}. Hence "a involution, c arbitrary" and
"c involution, a arbitrary" are the same problem.

**L3 (gadget B1, intended family).** Block e: c-cycle Γ_e = (g_0 g_1 … g_L) of length L+1 and a-chord
(g_0, g_Δe), label d_e := L − Δ_e ∈ Z_L. Let n = |E| blocks, N = n(L+1). Consider solutions x such that every
block contains a fixed point of t (equivalently: b's transposition in block e is (c(v_e), v_e) for a point v_e).
Then v_e ∈ {g_0, g_Δe} is an a-endpoint (orientation ε_e), the remaining L points of the block form one t-cycle
B_e whose b-endpoint is u_e = c(v_e), and with z the permutation of E induced by x on {B_e}:
 - z has only 3-cycles, plus fixed points allowed only when 3 ∤ L on blocks with effective label ≡ −3^{-1} (mod L);
 - the effective label is d_e (ε_e = 0, v_e = g_Δe) or L−1−d_e (ε_e = 1, v_e = g_0);
 - x is uniquely determined by (ε, z) (x(α_e) = u_{z(e)}, x(v_e) = v_{z(e)}, and x t = t x);
 - x^3 = t holds iff every 3-cycle {e,f,g} of z satisfies  d_e + d_f + d_g ≡ −1 (mod L)  (effective labels).
Proof: x(u_e) = x(t^{−d_e} α_e) = t^{−d_e} u_{z(e)}, hence x^3(u_e) = t^{−(d_e+d_{ze}+d_{z²e})} u_e = t(u_e).
The a-chord (α_e, α'_e) is mapped to (u_{z(e)}, v_{z'(e)}) which must be a b-transposition, forcing the
same block permutation on fixed points and on L-cycles. ∎
Consequence: restricted to this family, P3 on gadget-B1 instances is exactly "signed 3-partition mod L"
(choose per block d_e or −1−d_e, group into triples of sum −1 mod L). With integer ranges
d_e ∈ [D, 2D), L > 6D, 3 | L and a residue trick (d = 9a + s, s ∈ {0,1,4} for the three classes of a
Numerical-3DM instance, L = 9B+6) the flips and class violations are all excluded by counting arguments, so
*if the family were the whole solution set*, P3-inv would be NP-complete (reduction from Numerical 3DM,
N = 3q(9B+7)). This is the target statement; it is NOT established, see §3.

## 2. Verified numerically (complete enumeration)

Grids: L = 2 (N = 9, enumeration of S_9), L = 3 (N = 12) and L = 4 (N = 15) (propagation solver with
max_solutions = 10^7, `exhausted = False` checked on every instance), all label vectors in Z_L^3 (8 + 27 + 64 = 99 instances).
 - Every member of the family of L3 satisfies x a x^2 = c (asserted at construction, all 99 instances); the
   family is contained in the enumerated solution set on all 99 instances (`family_pass` = true everywhere).
 - Every solution outside the family has a block WITHOUT a t-fixed point (0 counter-examples to L3 in 99 instances).
 - Unintended solutions: L = 3: 13 instances × 8 solutions (+1 instance × 33) all with every block split by b
   into two t-cycles of lengths (2,2); L = 4: 56 instances with splits (2,3) in every block; L = 2: cross-block b.
   Gadget passes (solution set = family) on 14/27 instances at L = 3 and exactly on the 8/64 instances with
   all labels in {0, L−1} at L = 4 (chords joining c-adjacent points).
 - Labels (1,1,1) at L = 3 is a NO instance of the 3-partition core (family empty) with 33 solutions of the
   (2,2)-split shape: the reduction is wrong on NO instances as well as inflating YES instances.
 - One probe at N = 30 (L = 4, six blocks) did not finish within 3·10^6 nodes (24 solutions found, family
   size 360) — provisional only, not used for any claim.

## 3. Why gadget B1 fails, and what is conjectured

The chord lets the solver split Γ_e at any position: b = (u, c^j u) produces t-cycles of lengths (j, L+1−j)
with 2 ≤ j ≤ L−1, one a-endpoint in each part, and three blocks with the same j are then grouped exactly as in L3
with two sum conditions (mod j and mod L+1−j) but three sliding parameters, which is generically satisfiable.
Making the intended split (1, L) the only possible one would require either blocks of distinct sizes (which
destroys the equal-length grouping that carries the choice) or number-theoretic conditions on Δ_e for every
j simultaneously (infeasible in general). A variant with the fixed point outside the c-cycle (a = (α_e, φ_e),
φ_e a fixed point of c) forces the b-position (labels restricted to {0, L−1}) and thus carries no numerical
information (computed by hand in this session, not tested by code).

Conjectures / open: (i) the full solution set of gadget-B1 instances is the union over j of "(j, L+1−j)-split
families", each a signed 3-partition mod j and mod L+1−j (consistent with all 99 instances; not proved);
(ii) a rigid version of the gadget needs an extra mechanism to pin the b-transposition to c-adjacent points
while keeping a free label — e.g. a second colour class of blocks of size L'+1 with L' ∉ {L} whose only role is
to absorb fixed points — not tried. No NP-hardness and no algorithm for P3 follow from this track.
