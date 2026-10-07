# Campagne « prêt pour JSC » — bilan (7 octobre 2026)

Cinq pistes parallèles (sous-agents), toutes terminées ; chaque théorème nouveau re-vérifié dans la session principale par
`verify_round2.py` (code écrit à partir des énoncés) — résultats dans `verify_round2.json`.

## Résultat : p3_article.tex v7 (38 pages, 0 erreur, 0 avertissement, 42 références toutes citées, 5 figures)

| Piste | Résultat | Statut |
|---|---|---|
| Lemmes | Lemme du centralisateur pour **tout** a (critère : chaque classe « liée » (L, ℓ, R) a un effectif divisible par 3 ; formule close du nombre de racines ; O(N)), incluant 3 ∣ L via résultant/Hensel ; lemme du support borné en O(N^{s+1}) avec un nouveau lemme « racines cubiques à images prescrites » | **PROUVÉ** ; vérifié 49 906 paires (N ≤ 8) + 254 136 instances support ; re-vérifié ✔ 6699/6699 paires, 44/44 classes |
| Conj_n composé | **Conj_n NP-complet pour tout n ≥ 3** (critère affine pour xⁿ = 1 valable pour tout n et tout L ≥ 3 ; lemme d'encodage avec le facteur n/ℓ ; réduction depuis NnDM) ; classification complète : polynomial ssi n ∈ {1, 2} ; conjecture « puissance de 2 » réfutée (Conj₄ dur) ; la route « adjoindre un gadget d'ordre m » est impossible (témoins : (01234)→a² n'a que des conjugateurs d'ordre 4) | **PROUVÉ** ; re-vérifié ✔ critère 392/392 instances (L ∈ {4..9}, n ∈ {2,3,4,6,8,9,12}), témoins ordres [4] et [6] |
| Quatre couches | La construction à 4 couches **n'est pas** une réduction Conj₃ → P₃ : à k = 1, L = 5, c = a⁻¹ (Conj₃ NON), P₃(a′, c′) a exactement 10 solutions, aucune ne respectant les couches (type (15, 3, 1, 1)) ; idem pour les 160 vecteurs d'exposants admissibles ; se propage à tout k. Les zéros à m = 2, 3 étaient accidentels (premier moment sur les 20-cycles = 12337775/11639628 ≈ 1,06) | **OBSTRUCTION PROUVÉE** (Prop. prop:fourlayer) ; témoin re-vérifié ✔ (a′, c′ reconstruits depuis la définition, x résout, Conj₃ NON) ; k = 2 mixte non exhaustif (sans objet) |
| Audit | 3 bloquants / 15 à corriger / 12 style : cor:lambda clause 3 contredisait le problème ouvert 2 ; prop:moment affirmait « polynomial » à tort ; lem:support esquisse fautive ; thm:transposition cas ℓ = N − ℓ ; lem:affinep parenthèse fausse pour L composé ; intervalle [3B+3,5B] → [2B+3,5B] ; recensement rigide « 0,8 % » faux (7,1 %, 6,2 %, 3,3 %, 1,1 % pour N = 5…8) ; degrés de dimensionnement non traçables | **Tous appliqués** dans v7 (patches `audit_patches.tex`) ; recensement rigide re-vérifié ✔ N = 5, 6 : 0.071, 0.062 (pondéré par taille de classe) |
| Nouveauté | Conj_p / Conj_n : NOUVEAU (plus proches : Mattes–Ushakov–Weiß LATIN 2024, équations sphériques ; Lohrey–Rosowski–Zetzsche J. Algebra 2025) ; Conj₂ : NOUVEAU comme algorithme (cas c = a⁻¹ = « fortement réel », classique) ; PSL(2,7) : le phénomène est le triplet de Gassmann de degré 7 (Perlis 1977, Sunada 1985, Kammeyer–Kionke 2024), son usage comme obstruction est nouveau ; premier moment : conséquence routinière de comptes connus (Chowla–Herstein–Moore 1951, Moser 1955, Pavlov, Chernoff 1994 ; Goupil 1990, Bédard–Goupil 1992), non énoncé dans la littérature | 12 références ajoutées (DOI résolus via Crossref) et phrases de positionnement insérées ; OpenAlex non interrogé (approbation non revenue) |

## Ce que v7 change par rapport à v6
- Titre : « …the hardness of conjugacy by an element of prescribed order… » (le théorème couvre tout n).
- §3 : quatre lemmes complets (racines cubiques ; images prescrites ; support borné ; centralisateur) remplacent deux esquisses ; prop:moment fusionné dans thm:firstmoment (sans l'affirmation « polynomial ») ; remarque Gassmann après thm:psl (vérifié en session : G₁ et G₂ conjugués dans S₇ par des permutations impaires seulement, φ non induit par S₇, ψ = automorphisme extérieur).
- §4.3 réécrit : « Conjugacy by an element of prescribed order », lemme affine pour xⁿ = 1, lemme d'encodage, Théorème conjn (tout n ≥ 3), corollaire de classification complet, remarque sur le gadget impossible.
- §4.4 : Prop. prop:fourlayer + remarque premier moment après rem:padding.
- §5 : exposant « entre un cinquième et deux cinquièmes », degrés de dimensionnement rattachés aux ajustements de la table.
- §6 : problème 5 remplacé (n dans l'entrée ; L fixé ; types λ ≠ nᴺ/ⁿ).
- Abstract, liste des résultats, 13 patches de style.

## Ce qui reste à faire par un humain
1. Relecture externe des preuves de §4.3 (Conj_n, Conj₂) et des lemmes de §3 — écrits aujourd'hui, vérifiés numériquement, jamais lus par un collègue.
2. Décision sur le titre (`paper/submission_kit/title_options.txt`).
3. Release v0.4.0 sur GitHub (le `.tex` dit « release v0.4.0 »), puis Overleaf avec `P3_JSC_overleaf_v7.zip`.
4. Highlights et lettre de couverture : relire et signer (`submission_kit/`).

## Limites déclarées par les sous-agents (reportées telles quelles)
Lemmes : support vérifié exhaustivement à N ≤ 6 (s = 2, 3, 4) et N = 7 (s = 2), aléatoire N = 7–9. Conj_n : cas positifs à p = 7 non énumérés
(preuve uniforme) ; réduction exponentielle en n. Quatre couches : k = 2 mixte non exhaustif (1500 s). Audit : §5 relue pour la cohérence
seulement ; chiffres de rigidité N = 7, 8 issus de `audit_checks.py` (non re-vérifiés en session, N = 5, 6 l'ont été). Nouveauté : Crossref et
arXiv seulement, pas OpenAlex ; les verdicts « NOUVEAU » valent pour ce qui a été cherché (journal des requêtes dans `novelty_report.md`).
