# PhyloNet

Author: Guillaume Guénard
Institution: Université de Montréal

Développement d'une métrique de distance phylogénétique par induction/déduction
au moyen de réseaux de neurones convolutés.

## Idée sous-jacente

Simuler l'évolution de séquences d'ADN de différentes manières le long de
réseaux phylogénétiques de complexité variées (lignée, arbres, réseaux
faiblement ou fortement convergents) en _gardant la trace des degrés de parenté_
_entre les séquences_.

Les séquences évolueront suivant un processus de Markov aléatoire à partir d'une séquence initiale tirée aléatoirement. Cinq valeurs seront tirées avec
différents niveaux de prévallence: "-", "A", "C", "T", and "G"




pour la
substitution des nucléotides couplé à un processus de délétion/insertion de
nucléotides. Le processus de Markov sera défini par une matrice d'intensité de
mutation (dimension: $5 \times 5$). Chaque ligne de la matrice de mutation
correspond à l'état d'un nucléotide au temps $0$ et chaque colonne à l'état
du même nucléotide après qu'un temps $t$ se soit écoulé. Les six valeurs
d'intensité de mutation se trouvant au-dessus et en dessous de la diagonale sont
des paramètres de simulation; les quatre éléments de la diagonale n'en sont pas
cependant, étant estimés comme la négative de la somme de chaque ligne de la matrice. Cette définition assure que la somme de chaque ligne de la matrice est
bien $0$. En multipliant la matrice d'intensité de mutation par le produit du
pas de temps ($t$) et du taux de mutation moyen du nucléotide (un scalaire) et
en calculant l'exponentiel du produit résultant, on la matrice de probabilité
de mutation. Comme c'était le cas pour la matrice d'intensité de mutation, les
lignes de cette matrice correspondent à l'état des nucléotides au temps $0$,
ses colonnes à l'état des nucléotides au temps $t$ et ses valeurs aux
probabilités conditionnelles d'observer un état donné d'un nucléotide au temps
$t$ compte-tenu de son état au tmps $0$. 




les
valeurs situées sur chaque colonnes d'une ligne donnée correspondent aux
probabilités d'observer un 




Rendu ici...




La ligne de cette matrice correspondant au nucléotide actuel 


Les séquences proviendront toutes d'une séquence initiale générée
alléatoirement. Le taux de mutation sera attribué indépendemment pour chaque
nucléotide sera tiré d'une distribution de gamma. Ce taux suivra le nucléotide
tout au long de son évolution, de sa génération (à initialisation de la séquence
ou lors d'un événement d'insertion) jusqu'à sa délétion ou la fin du processus.
Le taux de mutation sera conservé lors des mutations (transition et
transversions). Lors d'événements de convergence, des sites de recombinaisons
seront générés et les fragments seront échangés. Le niveau de parenté de
l'hybride ainsi formé avec les autres séquences sera calculé comme la moyenne de
celui des parents pondéré par la proportion de matériel hérité de chacun d'eux.
Plusieurs miliers de jeux de données comptant chacun des milions de séquences
pourront ainsi être générés. Dans chacun de ces jeux de données, la distance
entre chaque séquence sera connue, de même que l'historique des divergences et
convergences entre celles-ci. Différents jeux de données pourront être générés
avec différentes matrices de mutation, différents 




Des efforts seront déployées afin de raffiner le processus de génération des
séquences en terme de temps d'exécution et d'utilisation des ressources de
calcul.





En plus des codes IUPAC "standards" qui seront générés (ACGT-),
des codes d'indétermination (RYSWKMBDHVN) pourront optionellement être
substitués aux nucléotides. 





