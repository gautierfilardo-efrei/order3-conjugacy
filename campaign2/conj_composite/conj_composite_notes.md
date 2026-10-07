# Track "Conj_n composé" — Conj_n for composite n

Problem. For an integer n ≥ 1, Conj_n: given a, c ∈ Sym(N), is there x ∈ Sym(N) with x^n = 1 and x a x^{-1} = c ?
Conventions as in the paper: image tuples, (pq)(j) = p(q(j)).
Known from the paper (v6): Conj_1 trivial (a = c), Conj_2 ∈ P (orbit pairing of ⟨a, c⟩), Conj_p NP-complete for every prime p ≥ 3
(block lemma + affine criterion for x^p = 1 + NpDM, Thm. thm:conjp).

Scripts / data: `conj_composite_track.py` (stages V1, V1b, V2, V3, V4), `conj_composite_results.json`.

## 0. Summary of the track

| Item | Status | Content |
|---|---|---|
| Theorem C | **PROVED** | For every integer n ≥ 3, Conj_n is NP-complete; hardness already for a = k disjoint L-cycles (L prime or L = 2^r), c ∈ Cent(a), certificates of cycle type n^{N/n} (order exactly n). |
| Corollary | **PROVED** | Conj_n ∈ P iff n ∈ {1, 2} (assuming P ≠ NP). The guess "polynomial iff n is a power of 2" is **refuted**: Conj_4 is NP-complete. Same classification for the exact-order variant Conj_n^{=}. |
| Lemma A | **PROVED** | Affine criterion for x^n = 1 on the blockwise family, any L ≥ 3, any n: solvable iff ∃σ ∈ Sym(k), σ^n = 1, with (∏_{i∈C} m_i)^{n/|C|} ≡ 1 (mod L) on every cycle C of σ. (No translation bookkeeping: t ≡ 0 always suffices.) |
| Lemma B | **PROVED** | The paper's exponent encoding (p → n, same T_0 = M(n+1)^{n−1}) also excludes all σ-cycles of proper-divisor length ℓ | n, ℓ < n — the new case that composite n creates. |
| Lemma D | **PROVED** | Disjoint unions with incompatible cycle lengths: Conj_n(a_1⊔a_2, c_1⊔c_2) ⇔ Conj_n(a_1,c_1) ∧ Conj_n(a_2,c_2). |
| Proposition E | **PROVED** (negative) | The "adjoin a gadget of order m" reduction Conj_p → Conj_{pm} by disjoint union cannot work: Conj_{pm}(a,c) ≠ Conj_p(a,c) in general (a = (1 2 3 4 5), c = a² has conjugators of order 4 only; a 7-cycle with c = a³ has conjugators of order 6 only). Not needed, since Theorem C treats all n directly. |
| V1 | **VERIFIED** | Exhaustive tables of conjugator orders for all pairs (a, c) in Sym(N), N ≤ 9 (all ordered pairs in the same class). |
| V1b | **VERIFIED** | Lemma D on 150 random disjoint unions (brute force over Sym(N1+N2), N1+N2 ≤ 9): 0 failures. |
| V2 | **VERIFIED** | Lemma A against brute force over Sym(N) for N ≤ 9 (all unit vectors, 10 shapes (L,k), 42 instances) and against the full list of L^k k! conjugators (9 shapes, up to 2 016 840 conjugators per instance, N up to 35), n = 1..12: 0 mismatches. |
| V3 | **VERIFIED** | Lemma B by enumeration of all subsets S with |S| | n for n ∈ {3,4,5,6,8,9}, 124 random instances, 931 844 subsets, two values of T each: 0 mismatches. |
| V4 | **VERIFIED** | End-to-end reduction N3DM → Conj_n for n ∈ {3,4,6,9} on 80 random instances (q ≤ 3), answer compared with brute-force N3DM: 0 mismatches; for n = 3, 4, 6 the certificate x was built as an explicit permutation (N up to 13 513 116) and x a x^{-1} = c, x^n = 1, x^d ≠ 1 (d | n, d < n) checked numerically. |

Deviation from the task text: the "natural conjecture" (polynomial iff power of 2) is false; the right statement is the Corollary. The question "does Conj_p reduce to Conj_{pm} by a disjoint gadget" is answered negatively (Prop. E) but becomes moot.

## 1. Setting (blockwise family, as in the paper)

Fix L ≥ 3, k ≥ 1. a consists of k disjoint L-cycles on blocks B_1, …, B_k, each identified with Z_L so that a is z ↦ z + 1.
For units m = (m_1, …, m_k) ∈ (Z_L^*)^k, c = c(m) acts on B_i by z ↦ z + m_i; c ∈ Cent(a), same cycle type as a.
Block lemma (paper, Lemma lem:block, valid for every L ≥ 3): the conjugators of a to c(m) are exactly the maps x sending B_i onto B_{σ(i)} by
z ↦ m_{σ(i)} z + t_i, for σ ∈ Sym(k) and t ∈ Z_L^k.

## 2. Lemma A — affine criterion for x^n = 1 (PROVED)

```latex
\begin{lem}[affine criterion for $x^{n}=1$]\label{lem:affinen}
Let $n\ge1$, $L\ge3$ and $a,c(m)$ be as above. $\mathrm{Conj}_n(a,c(m))$ has a solution iff there is $\sigma\in\Sym(k)$ with
$\sigma^{n}=1$ such that for every cycle $C$ of $\sigma$, writing $R_C=\prod_{i\in C}m_i$,
\[ R_C^{\,n/|C|}\equiv1\pmod L. \]
When the condition holds, $x$ defined by $x(z)=m_{\sigma(i)}z$ on $B_i$ (all translations zero) is a solution; moreover
$x^{d}\neq1$ for every $d<n$ that is not a multiple of the order of $\sigma$.
\end{lem}
\begin{proof}
Let $x$ be a conjugator with block permutation $\sigma$ and translations $t$ (block lemma). Then $x^{j}$ maps $B_i$ onto
$B_{\sigma^{j}(i)}$, so $x^{n}=1$ forces $\sigma^{n}=1$; thus every cycle $C=(i_1\,i_2\cdots i_\ell)$ of $\sigma$ has length
$\ell\mid n$. On $B_{i_1}$, $x^{\ell}$ is the composite of the $\ell$ affine maps $z\mapsto m_{i_2}z+t_{i_1}$,
$z\mapsto m_{i_3}z+t_{i_2}$, \dots, $z\mapsto m_{i_1}z+t_{i_\ell}$, hence an affine map $z\mapsto R_Cz+\tau_C$ with
$R_C=\prod_{i\in C}m_i$; and $x^{n}|_{B_{i_1}}=(x^{\ell}|_{B_{i_1}})^{n/\ell}$ is $z\mapsto R_C^{n/\ell}z+(1+R_C+\dots+R_C^{n/\ell-1})\tau_C$.
An affine map $z\mapsto\alpha z+\beta$ of $\ZZ_L$ ($L\ge3$) is the identity iff $\beta=0$ (take $z=0$) and $\alpha=1$ (take $z=1$).
So $x^{n}=1$ implies $R_C^{n/\ell}\equiv1$ for every cycle: the condition is necessary.
Conversely, if $\sigma^{n}=1$ and $R_C^{n/|C|}\equiv1$ for all cycles, take $t\equiv0$. By the block lemma $x$ conjugates $a$ to $c(m)$;
on each block of a cycle $C$ of length $\ell$, $x^{\ell}$ is the linear map $z\mapsto R_Cz$ (the product is taken in the cyclic order
starting at that block, but all cyclic rotations of a product of elements of the abelian group $\ZZ_L^{*}$ coincide), so
$x^{n}|_{B_i}=(z\mapsto R_C^{n/\ell}z)=\mathrm{id}$ on every block, i.e.\ $x^{n}=1$. Finally $x^{d}$ moves $B_i$ to
$B_{\sigma^{d}(i)}\ne B_i$ whenever $\sigma^{d}\ne1$, so $x^{d}\ne1$ for such $d$.
\end{proof}
```

Remarks. (i) For n = p prime this is the paper's Lemma lem:affinep (cycles of length p: product 1; fixed points: m_i^p = 1).
(ii) The proof needs no case analysis on translations, unlike the paper's version: taking t ≡ 0 suffices for sufficiency, and only the slope
is used for necessity. This simplification can be back-ported to the paper's Lemma lem:affine / lem:affinep.
(iii) For composite n the new feature is the σ-cycles of proper-divisor length ℓ | n, ℓ < n, which are allowed as soon as
R_C^{n/ℓ} = 1 (e.g. for n = 4, L prime: 2-cycles of σ with m_i m_j ≡ −1, on which x² acts as z ↦ −z + τ, an involution). Lemma B
rules these out in the hard instances.

## 3. Lemma B — the exponent encoding also excludes proper-divisor cycles (PROVED)

Numerical n-dimensional matching NnDM (n ≥ 3 fixed): classes C_1, …, C_n of q elements, sizes s(u) ∈ {1, …, B'}, target B' (unary);
decide whether ⋃C_j splits into q n-tuples, one element per class, each of size-sum B'. NnDM is strongly NP-complete for every fixed n ≥ 3
(N3DM = [Garey–Johnson SP16], padded with n − 3 classes of q elements of size 1 and target B + n − 3 — the paper's padding).

Encoding (the paper's, with p replaced by n): M = nB' + 1, w_j = (n+1)^{j−1} (1 ≤ j ≤ n−1), W = Σ_{j<n} w_j = ((n+1)^{n−1} − 1)/n,
T_0 = M (n+1)^{n−1}; for an integer T ≥ T_0 define the integers
e_u = s(u) + M w_j (u ∈ C_j, j < n), e_u = s(u) − MW − B' (u ∈ C_n).

```latex
\begin{lem}[encoding lemma, composite $n$]\label{lem:encn}
Let $n\ge3$, $T\ge T_0$, and let $S\subseteq\bigcup_jC_j$ with $|S|=\ell$, $\ell\mid n$. Then
\[ \frac{n}{\ell}\sum_{u\in S}e_u\equiv0\pmod T \iff \ell=n,\ |S\cap C_j|=1\ \text{for all } j,\ \text{and}\ \sum_{u\in S}s(u)=B'. \]
\end{lem}
\begin{proof}
Put $\beta_j=|S\cap C_j|$, $\Sigma=\sum_{u\in S}s(u)\in[\ell,\ell B']$ and $D=\sum_{j<n}w_j(\beta_j-\beta_n)\in\ZZ$. Summing the
definitions gives the integer identity
\[ \sum_{u\in S}e_u=\Sigma-\beta_nB'+MD. \tag{1} \]
\emph{Size bound.} Every $e_u$ satisfies $|e_u|\le E:=MW+B'$ (for $j<n$: $0<e_u\le B'+Mw_{n-1}\le E$; for $j=n$:
$e_u\in[1-MW-B',\,-MW]$). Now $nE=M\bigl((n+1)^{n-1}-1\bigr)+nB'=T_0-M+nB'=T_0-1<T$, hence
\[ \Bigl|\sum_{u\in S}e_u\Bigr|\le\ell E<\frac{\ell}{n}\,T\le\frac{T}{\gcd(n/\ell,\,T)}. \]
Since $\frac n\ell\,y\equiv0\pmod T$ iff $y\equiv0\pmod{T/\gcd(n/\ell,T)}$, the left-hand side of the lemma holds iff
$\sum_{u\in S}e_u=0$ as an integer.
\emph{(a) $D=0$ iff all $\beta_j$ are equal.} Put $d_j=\beta_j-\beta_n\in[-\ell,\ell]\subseteq[-n,n]$. If some $d_j\ne0$ and $j_0$ is the
least such index, $D\equiv(n+1)^{j_0-1}d_{j_0}\pmod{(n+1)^{j_0}}$, so $D=0$ would force $(n+1)\mid d_{j_0}$, impossible for
$0<|d_{j_0}|\le n$. If all $\beta_j$ equal $\beta$ then $n\beta=\ell\le n$, so $\beta=1$ and $\ell=n$ ($S\ne\emptyset$).
\emph{(b)} If $D=0$ then $\beta_n=1$, $\ell=n$, and (1) reads $\Sigma-B'$, which vanishes iff $\Sigma=B'$.
\emph{(c)} If $D\ne0$ then $|MD|\ge M>M-1=nB'\ge|\Sigma-\beta_nB'|$ (as $\Sigma\le\ell B'\le nB'$ and $\beta_nB'\le nB'$), so (1) is
nonzero. Together: the integer sum vanishes iff $\ell=n$, one element per class, and $\Sigma=B'$.
\end{proof}
```

Remark. The paper's proof of Thm. thm:conjp uses "as Σβ_j = p is prime, equal β_j means β_j = 1"; primality is not needed (n classes,
Σβ_j = |S| ≤ n). The only genuinely new point for composite n is the size bound with the factor n/ℓ, which the paper's T_0 already covers.

## 4. Theorem C — Conj_n is NP-complete for every n ≥ 3 (PROVED)

```latex
\begin{thm}\label{thm:conjn}
For every integer $n\ge3$, $\mathrm{Conj}_n$ is NP-complete. Hardness holds already when $a$ consists of $k$ disjoint cycles of one
length $L$, with $L$ prime (or $L$ a power of two), when $c\in\Cent(a)$, and when the certificate is required to have cycle type
$n^{N/n}$, i.e.\ order exactly $n$. Consequently, unless $\mathrm P=\mathrm{NP}$, $\mathrm{Conj}_n$ is polynomial iff $n\in\{1,2\}$,
and the same holds for the variant asking for a conjugator of order exactly $n$.
\end{thm}
\begin{proof}
Membership in NP is clear. Fix $n\ge3$ and let an $\mathrm{N}n\mathrm{DM}$ instance $(C_1,\dots,C_n;s;B')$ be given. With $M,w_j,W,T_0$
and $e_u$ as in Lemma~\ref{lem:encn}, let $L$ be the least prime with $T:=L-1\ge T_0$ (by Bertrand's postulate $L\le2T_0+2$, found by trial
division in time polynomial in $T_0$) and $g$ a primitive root mod $L$ (found by exhaustive search, $L$ being polynomial in the input);
alternatively $L=2^{r}$ with $T=2^{r-2}\ge T_0$ and $g=5$, of order $T$ in $\ZZ_L^{*}$. Put $m_u=g^{e_u}$, $k=nq$, $N=kL$, and let
$(a,c)=(a,c(m))$ be the blockwise instance; $N=\Theta(qB'(n+1)^{n-1})$ is polynomial in $q+B'$ for fixed $n$. As $g$ has order $T$,
for every set $S$ of blocks with $|S|=\ell$, $R_S:=\prod_{u\in S}m_u=g^{\sum_Se_u}$ satisfies
$R_S^{n/\ell}=1$ iff $\frac n\ell\sum_{u\in S}e_u\equiv0\pmod T$.

By Lemma~\ref{lem:affinen}, $\mathrm{Conj}_n(a,c)$ has a solution iff the blocks can be partitioned into the cycles of some
$\sigma$ with $\sigma^{n}=1$ such that each cycle $S$ (of length $\ell\mid n$) satisfies $R_S^{n/\ell}=1$; by Lemma~\ref{lem:encn} a set
$S$ with $|S|\mid n$ satisfies this iff it is a matching $n$-tuple (one element per class, size-sum $B'$). Hence
$\mathrm{Conj}_n(a,c)$ is solvable iff the blocks can be partitioned into matching $n$-tuples, i.e.\ iff the $\mathrm{N}n\mathrm{DM}$
instance is a yes-instance. In that case every solution $x$ has block permutation $\sigma$ consisting of $n$-cycles only, so every
orbit of $x$ has size exactly $n$ ($x^{d}$ moves every block for $0<d<n$): the certificates have cycle type $n^{N/n}$, and
requiring order exactly $n$ changes nothing. Since $\mathrm{N}n\mathrm{DM}$ is strongly NP-complete for fixed $n\ge3$, so is
$\mathrm{Conj}_n$ (NP-hard, in NP). With $\mathrm{Conj}_1$ trivial and $\mathrm{Conj}_2\in\mathrm P$ (orbit pairing), the classification
follows; for the exact-order variant with $n=2$, if $a\ne c$ every solution of $x^{2}=1$ is an involution, and if $a=c$ the question is
whether $\Cent(a)$ contains an involution, which holds iff $a$ has a cycle of even length or two cycles of the same length (read off the
cycle type).
\end{proof}
```

Why this does not contradict Conj_2 ∈ P: for n = 2 the same construction would encode numerical 2-dimensional matching (pairs of
prescribed sum), which is polynomial; N n DM is NP-complete only for n ≥ 3.

Consequence for the question "Conj_n for n part of the input": NP-complete (contains Conj_3); the reduction above is polynomial only for
fixed n since N = Θ(q B' (n+1)^{n−1}), which is irrelevant for the classification.

## 5. Disjoint unions and the failure of "order-lowering" gadgets

```latex
\begin{lem}[disjoint unions]\label{lem:union}
Let $U=U_1\sqcup U_2$, $a=a_1\sqcup a_2$, $c=c_1\sqcup c_2$ with $a_i,c_i\in\Sym(U_i)$, and write $\Lambda(\cdot)$ for the set of
cycle lengths (fixed points counting as $1$-cycles). If $\Lambda(a_1)\cap\Lambda(c_2)=\emptyset=\Lambda(a_2)\cap\Lambda(c_1)$, then
every conjugator of $a$ to $c$ is of the form $x_1\sqcup x_2$ with $x_ia_ix_i^{-1}=c_i$, and for every $n$,
\[ \mathrm{Conj}_n(a,c)\iff\mathrm{Conj}_n(a_1,c_1)\wedge\mathrm{Conj}_n(a_2,c_2). \]
\end{lem}
\begin{proof}
A conjugator maps each cycle of $a$ onto a cycle of $c$ of the same length. A cycle of $a$ inside $U_1$ has length in $\Lambda(a_1)$,
and no cycle of $c$ inside $U_2$ has such a length, so its image lies in $U_1$; symmetrically $x(U_2)\subseteq U_2$. Hence
$x=x_1\sqcup x_2$, $x^{n}=x_1^{n}\sqcup x_2^{n}$, and $x^{n}=1$ iff $x_1^{n}=x_2^{n}=1$; conversely any pair of such $x_i$ glues.
\end{proof}
```

```latex
\begin{prop}[no order-lowering by disjoint gadgets]\label{prop:nogadget}
Let $p$ be prime and $m\ge2$. There is no pair $(g,g')$ of permutations such that, for all $(a,c)$ with cycle lengths disjoint from those
of $(g,g')$, $\mathrm{Conj}_{pm}(a\sqcup g,c\sqcup g')\iff\mathrm{Conj}_p(a,c)$. More precisely, for every $n\ge1$ there are pairs
$(a,c)$ with $\mathrm{Conj}_n(a,c)$ true and $\mathrm{Conj}_d(a,c)$ false for every proper divisor $d$ of $n$: take a prime
$L\equiv1\pmod n$, $a$ an $L$-cycle and $c=a^{u}$ with $u$ of multiplicative order $n$ modulo $L$.
\end{prop}
\begin{proof}
By Lemma~\ref{lem:union} the left-hand side equals $\mathrm{Conj}_{pm}(a,c)\wedge\mathrm{Conj}_{pm}(g,g')$, so the equivalence would give
$\mathrm{Conj}_{pm}(a,c)\Rightarrow\mathrm{Conj}_p(a,c)$ for all such $(a,c)$. For the second statement, the conjugators of an
$L$-cycle $a$ ($z\mapsto z+1$ on $\ZZ_L$) to $a^{u}$ are the affine maps $z\mapsto uz+t$ (block lemma with $k=1$); for $u\ne1$ and $L$ prime
such a map has the fixed point $z_0=t/(1-u)$ and is conjugate in $\mathrm{AGL}_1(L)$ to $z\mapsto uz$, hence has order exactly
$\mathrm{ord}_L(u)=n$. So $\mathrm{Conj}_{n'}(a,a^{u})$ holds iff $n\mid n'$. With $n=pm$ this contradicts the implication
(e.g.\ $p=2,m=2$: $a=(1\,2\,3\,4\,5)$, $c=a^{2}$, all conjugators of order $4$; $p=3,m=2$: $a$ a $7$-cycle, $c=a^{3}$, all of order $6$).
Primes $L\equiv1\pmod n$ exist by Dirichlet's theorem.
\end{proof}
```

Consequence. The lattice {n : Conj_n(a,c)} is the set of multiples of the orders of conjugators, and it is not determined by its prime
members; a direct hardness proof for each n (Theorem C) is therefore the right route, and it needs no reduction between different Conj's.
(A useful positive by-product of Lemma D: Conj_n hardness instances can be padded with any disjoint part on which an order-n conjugator
exists, e.g. fixed points, without changing the answer.)

## 6. Verification record (script `conj_composite_track.py`, results `conj_composite_results.json`)

V1 — exhaustive profiles, N ≤ 9 (stage `V1`). For each cycle type λ ⊢ N and each of the |class| conjugates c of a representative a, the
orders of all conjugators x (x ranges over Sym(N), c = x a x^{-1}, grouped by c) are tabulated; D(a,c) = {n ≤ 12 : Conj_n(a,c)} is the set
of multiples of those orders. Counts are over all ordered pairs (a, c) in the same class. See the JSON for the per-N tables
(`conj_n_yes_fraction`, `yes_set_distribution`) and the "strict witnesses" (Conj_n true, Conj_d false for all proper d | n) with an example
per (N, n). Fraction of same-class ordered pairs (a, c) with Conj_n(a, c) true (from the JSON):

| N | pairs | n=1 | n=2 | n=3 | n=4 | n=6 | n=8 | n=9 | n=12 | distinct sets D(a,c) |
|---|---|---|---|---|---|---|---|---|---|---|
| 5 | 2 602 | 0.046 | 0.889 | 0.611 | 1.000 | 0.982 | 1.000 | 0.611 | 1.000 | — |
| 6 | 71 412 | 0.010 | 0.568 | 0.519 | 1.000 | 0.976 | 1.000 | 0.519 | 1.000 | — |
| 7 | 2 675 724 | 0.002 | 0.359 | 0.456 | 0.875 | 0.984 | 0.875 | 0.456 | 1.000 | — |
| 8 | 134 269 158 | 0.0003 | 0.192 | 0.273 | 0.815 | 0.974 | 0.943 | 0.273 | 0.999 | — |
| 9 | 8 747 088 662 | 4e-5 | 0.094 | 0.196 | 0.699 | 0.957 | 0.915 | 0.703 | 0.997 | 95 |

Strict witnesses (Conj_n true, Conj_d false for every proper divisor d of n), smallest N and one example (conjugator orders in brackets):
n = 4: N = 5, a = (1 2 3 4 5), c = a² = (1 3 5 2 4), orders [4] (288 pairs at N = 5; 5.29·10⁹ at N = 9);
n = 6: N = 6, a = (1 2 3 4 5 6), c = (1 4 3 6 2 5), orders [4, 5, 6] (11 880 pairs at N = 6);
n = 8: N = 8, a = (1 … 8), c = (1 7 2 4 5 6 3 8), orders [3, 6, 8, 12];
n = 9: N = 9, a = (1 … 9), c = (1 5 6 4 7 3 8 2 9), orders [4, 6, 9, 15] (4.43·10⁹ pairs at N = 9) — Conj_9 true, Conj_3 false;
n = 12: N = 8, a = (1 2 3 4)(5 6), c = (1 7 2 8)(3 4), orders [7, 8, 12, 15].
Empirical observation (VERIFIED, N ≤ 9 only): for N ≤ 6 every same-class pair has a conjugator of order dividing 4 (fraction 1.000);
the first pairs without one appear at N = 7, e.g. the 7-cycle a with c = a³ = (1 4 7 3 6 2 5), whose seven conjugators z ↦ 3z + t all
have order 6 (stage V1b, `seven_cycle_cube_conjugator_orders`). The strict-witness example recorded in the JSON for (N, n) = (7, 6) is a
different pair, c = (1 5 4 3 6 2 7), whose conjugators have orders {6, 10} (797 760 such pairs at N = 7). No general statement is claimed.

V1b — Lemma D (stage `V1b`): 150 random unions (a_1 ⊔ a_2, c_1 ⊔ c_2) with disjoint cycle-length sets, N_1 + N_2 ≤ 9, brute force over
Sym(N_1+N_2): 0 failures. Also: all conjugators of (1 2 3 4 5) → a² have order 4; of a 7-cycle → a³ have order 6 (brute force).

V2 — Lemma A (stage `V2`): (i) brute force over Sym(N) for (L,k) ∈ {(3,1),(4,1),(5,1),(6,1),(7,1),(8,1),(9,1),(3,2),(4,2),(3,3)}, every
unit vector m, n = 1..12: 42 instances, 0 mismatches between D(a,c) and the σ-criterion; (ii) the full list of L^k k! conjugators
(block lemma) with actual permutation powers, (L,k) ∈ {(5,3),(7,3),(11,3),(13,3),(5,4),(7,4),(11,4),(5,5),(7,5)}, 12 or 6 random m per
shape (30 % with a slope equal to 1), 90 instances, up to 2 016 840 conjugators each: 0 mismatches.

V3 — Lemma B (stage `V3`): for n ∈ {3,4,5,6,8,9} and random padded N3DM instances (q ≤ 3, planted and non-planted), all subsets S with
|S| | n, for T = T_0 and T = 2T_0 + 1: 931 844 subsets over 124 instances, 0 mismatches with "ℓ = n, one per class, Σ = B'".

V4 — end-to-end (stage `V4`): random padded N3DM instances, L = least prime with L − 1 ≥ T_0, g a primitive root, m_u = g^{e_u}; the
block partition problem of Lemma A solved by exact cover over all admissible sets (|S| | n and R_S^{n/|S|} = 1); answer compared with
brute-force N3DM: n = 3 (24 inst., N ≤ 4 491), n = 4 (34 inst., N ≤ 55 644), n = 6 (16 inst., N ≤ 13 513 116), n = 9 (6 inst., exponent
level only, N ≈ 2·10^{11}): 0 mismatches, half of the instances planted (yes). For n = 3, 4, 6 the certificate x (t ≡ 0) was built as a
permutation on N points and x a x^{-1} = c, x^n = 1, x^d ≠ 1 (d | n, d < n) were checked.

Limits of the verification. V4 for n = 6 and n = 9 could not be checked against enumeration of all conjugators (L^k k! is astronomically
large); it checks the composition of Lemma A (verified independently in V2) with the exact-cover solver, plus the explicit certificate.
Nothing was verified for n ≥ 13 beyond the proofs.

## 7. Open points

1. Parameterised / restricted versions: Conj_n is polynomial when a has a bounded number of cycles (|Cent(a)| = ∏ L^{k_L} k_L! polynomial,
   enumerate the coset); the hard instances need k = nq → ∞ blocks of equal length. Open: Conj_n when all cycles of a have pairwise
   distinct lengths (Cent(a) abelian, the coset is x_0 · ∏ Z_{L_i}); the slope/translation structure then becomes a system of congruences and
   we did not determine its complexity.
2. Counting conjugators of order dividing n (#Conj_n): #P-hardness is plausible for n ≥ 3 via the same reduction (parsimonious up to the
   n^{q}·(n−1)!^{q} cyclic orderings and the translation families), not written out.
3. Back-port: Lemma A's proof (t ≡ 0) shortens the paper's Lemma lem:affine / lem:affinep; Thm. thm:conjp can be replaced by Thm. C with
   the same figure (fig_reduction) — only the caption's "p prime" changes.
