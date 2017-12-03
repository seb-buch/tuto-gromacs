.. _tuto_lyso:

Tutoriel "Lysozyme"
===================

Le but de ce tutoriel est d'étudier l'influence des ponts disulfures sur la stabilité en température du lysozyme.
Le lysozyme est une protéine présente notamment dans le blanc d'oeuf et qui possède des propriétés antimicrobienne puisque
cette enzyme (EC 3.2.1.17) hydrolyse le peptidoglycane ce qui fragilise la paroi des bactéries à Gram positif:

.. figure:: images/Mecanism_of_action_for_Lysozyme.svg
    :width: 400px
    :align: center

    Liaison hydrolysée par le lysozyme.

Le lysozyme contient aussi 4 ponts disulfure qui lui confère une structure relativement stable puisque sa température de dénaturation élevée (85°C dans le glycérol [Knubovets1999]_).

.. figure:: images/lysozyme_S-S.png
    :width: 400px
    :align: center

    Structure du lysozyme avec ses 4 ponts S-S.

Création de la topologie
------------------------

La structure du lysozyme est bien connue et est disponible sur le site de la `PDB <http://www.rcsb.org/pdb/home/home.do>`_. Le code PDB pour le lysozyme de blanc d'oeuf est: `1AKI <http://www.rcsb.org/pdb/explore/explore.do?structureId=1aki>`_.
La première étape est donc de télécharger les coordonnées atomiques sous le format `.pdb`::

    > wget https://files.rcsb.org/download/1AKI.pdb

.. note::

    Ouvrir le fichier `1AKI.pdb` pour voir toutes les informations qu'il contient.

Le fichier `.pdb` ne contient que les coordonnées des atomes mais ne contient pas la topologie de la molécule qui, elle, est liée à un champ de forces.
Pour générer cette topologie, on peut utiliser l'utilitaire :ref:`pdb2gmx` qui permet de lire un fichier `.pdb` et de générer automatiquement la topologie associée.

.. important::

    Comme le but est d'étudier l'effet des ponts S-S, nous allons générer *2* topologies:

    1. une topologie contenant les ponts S-S

    2. une topologie ne contenant *pas* les ponts S-S

La commande pour générer la topologie du lysozyme est::

    > gmx pdb2gmx -f 1AKI.pdb -ss

L'option `-ss` rend interactif le choix des ponts S-S.

Premièrement :ref:`pdb2gmx` demande de choisir le champ de force::

     1: AMBER03 protein, nucleic AMBER94 (Duan et al., J. Comp. Chem. 24, 1999-2012, 2003)
     2: AMBER94 force field (Cornell et al., JACS 117, 5179-5197, 1995)
     3: AMBER96 protein, nucleic AMBER94 (Kollman et al., Acc. Chem. Res. 29, 461-469, 1996)
     4: AMBER99 protein, nucleic AMBER94 (Wang et al., J. Comp. Chem. 21, 1049-1074, 2000)
     5: AMBER99SB protein, nucleic AMBER94 (Hornak et al., Proteins 65, 712-725, 2006)
     6: AMBER99SB-ILDN protein, nucleic AMBER94 (Lindorff-Larsen et al., Proteins 78, 1950-58, 2010)
     7: AMBERGS force field (Garcia & Sanbonmatsu, PNAS 99, 2782-2787, 2002)
     8: CHARMM27 all-atom force field (CHARM22 plus CMAP for proteins)
     9: GROMOS96 43a1 force field
    10: GROMOS96 43a2 force field (improved alkane dihedrals)
    11: GROMOS96 45a3 force field (Schuler JCC 2001 22 1205)
    12: GROMOS96 53a5 force field (JCC 2004 vol 25 pag 1656)
    13: GROMOS96 53a6 force field (JCC 2004 vol 25 pag 1656)
    14: GROMOS96 54a7 force field (Eur. Biophys. J. (2011), 40,, 843-856, DOI: 10.1007/s00249-011-0700-9)
    15: OPLS-AA/L all-atom force field (2001 aminoacid dihedrals)

Ici, nous allons choisir `GROMOS96 54a7` (choix *14*).

Ensuite, le modèle d'eau doit être choisi::

     1: SPC    simple point charge, recommended
     2: SPC/E  extended simple point charge
     3: None

Il faut choisir `SPC` (choix *1*)

Enfin, le choix des ponts S-S intervient::

    Link CYS-6 SG-48 and CYS-127 SG-981 (y/n) ?

Pour chacun des ponts (4), il faut répondre `y` pour confirmer la présence du pont ou `n` pour générer une topologie sans le pont (= Cystéines libres).

A partir de ces réponses, :ref:`pdb2gmx` génère 3 fichiers:

    1. `conf.gro` qui contient les coordonnées des atomes (similaire au fichier `.pdb` de départ)
    2. `posre.itp` qui contient les données nécessaires pour contraindre les atomes (voir utilité plus loin)
    3. `topol.top` qui contient la topologie du lysozyme et du système.

.. note::
    Ouvrir ces fichiers pour voir ce qu'ils contiennent

.. tip::
    Pour générer les 2 topologies (avec et sans pont S-S), il faut donc exécuter `pdb2gmx` deux fois de suite.

.. important::
    Pour éviter de mélanger les deux topologies, il est fortement conseillé de travailler dans 2 dossiers différents!


Solvatation du lysozyme
-----------------------

Premièrement il est nécessaire de créer la boite de simulation autour de notre protéine grâce à :ref:`editconf`::

    > gmx editconf -f conf.gro -o lyso.gro -d 0.7 -bt dodecahedron

.. note::
    Ouvrir `lyso.gro` avec VMD pour observer la boîte de simulation et sa réplication dans les trois dimensions de l'espace.

Maintenant que la boîte de simulation est définie, on peut rajouter le solvant à l'aide de l'utilitaire :ref:`solvate`::

    > gmx solvate -p topol.top -cp lyso.gro -cs -o lyso_W.gro

.. note::
    Ouvrir `lyso_W.gro` avec VMD pour observer comment l'eau est ajoutée au système.

.. important::
    Il est important de vérifier que :ref:`solvate` a bien rajouté la ligne correspondant au solvant dans `topol.top`


Ajout des ions
--------------

Le lysozyme n'est pas neutre et est chargé positivement, la boîte MD étant répliquée dans les trois dimensions, le système
est virtuellement infini et il est important d'avoir une boîte neutre.
Il faut donc rajouter des ions négatifs (e.g. des Cl-) pour compenser les charges positives du lysozyme.
Cela ce fait en deux étapes:

Premièrement, on crée un fichier de topologie complet (`ions.tpr`) qui contient toutes les informations (coordonnées + topologie) nécessaires.
Cela est fait par l'intermédiaire de l'utilitaire :ref:`grompp` et nécessite, en plus des coordonnées et de la topologie,
un fichier de paramètre de simulation (:download:`ions.mdp <./files/lyso/ions.mdp>`)::

    > gmx grompp -f ions.mdp -c lyso_W.gro -p topol.top -o ions.tpr

.. note::
    Le fichier de paramètre de simulation (`ions.mdp`) n'a pas d'utilité mais il est néanmoins requis par :ref:`grompp`.

    Il peut d'ailleurs être intéressant de l'ouvrir (dans un éditeur de texte) pour voir quelle type de simulation il contient.

Le fichier `ions.tpr` contenant à la fois les coordonnées des atomes et la topologie du système (i.e. les molécules),
on peut maintenant l'utiliser pour remplacer une molécule d'eau par un ion Cl- afin de rendre la boîte MD neutre.
Pour cela on utilise l'utilitaire :ref:`genion`::

    > gmx genion -s ions.tpr -o lyso_WI.gro -p topol.top -neutral

Même si les molécules sont définies dans le fichier `ions.tpr`, :ref:`genion` ne sélectionner que des groupes d'atomes::

    Select a continuous group of solvent molecules
    Group     0 (         System) has 20808 elements
    Group     1 (        Protein) has  1323 elements
    Group     2 (      Protein-H) has  1001 elements
    Group     3 (        C-alpha) has   129 elements
    Group     4 (       Backbone) has   387 elements
    Group     5 (      MainChain) has   517 elements
    Group     6 (   MainChain+Cb) has   634 elements
    Group     7 (    MainChain+H) has   646 elements
    Group     8 (      SideChain) has   677 elements
    Group     9 (    SideChain-H) has   484 elements
    Group    10 (    Prot-Masses) has  1323 elements
    Group    11 (    non-Protein) has 19485 elements
    Group    12 (          Water) has 19485 elements
    Group    13 (            SOL) has 19485 elements
    Group    14 (      non-Water) has  1323 elements

Les molécules d'eau sont contenues dans plusieurs groupes (`System`, `Water` et `SOL`) mais apparaissent comme `SOL` dans le fichier de topologie `topol.top`.

Il faut donc choisir le groupe `13` dans :ref:`genion`.

.. note::
    Ouvrir le fichier `lyso_WI.gro` dans VMD pour visualiser où les ions ont été ajoutés.

.. important::
    Vérifier que le fichier `topol.top` a bien été mis à jour par :ref:`genion`.


.. _minim:

Minimisation énergétique
------------------------

Les molécules ayant été mises dans la boîte aléatoirement, l'énergie du système (calculée à partir du champ de force) est très élevé
et il est nécessaire de la minimiser avant de faire une simulation de dynamique moléculaire.
Pour cela, on va générer une topologie complète (`em.tpr`) qui va contenir les paramètres de minimisation (:download:`minim.mdp <./files/lyso/minim.mdp>`) à l'aide de :ref:`grompp`::

    > gmx grompp -f minim.mdp -c lyso_WI.gro -p topol.top -o em.tpr

.. important::
    Il faut *toujours* bien lire la sortie de :ref:`grompp` afin de vérifier que tout se passe bien (i.e. ni "WARNING", ni "ERROR")

Une fois que la topologie complète est prête, on peut lancer la minimisation à l'aide de :ref:`mdrun`::

    > gmx mdrun -v -deffnm em

.. important::
    Ne pas recopier (et exécuter) cette commande sans réflexion! (:ref:`Pourquoi ?<warning_mdrun>`)

.. note::
    Utiliser :ref:`energy` pour observer la diminution de l'énergie au cours de la minimisation.


.. _equilibration_nvt:

Équilibration NVT
-----------------

Maintenant que la configuration du système est minimisée, on peut équilibrer le système, c'est à dire introduire
les contraintes de température et de pression qui seront appliquées pendant la simulation de production (de données).

On introduit toujours la contrainte de température en premier sous la forme d'une simulation souvent assez courte (inférieure à 1 ns)
dans laquelle:

    1. Le nombre `N` de particles ne varie pas (imposé par les `Conditions Périodiques aux Limites <https://fr.wikipedia.org/wiki/Condition_p%C3%A9riodique_aux_limites>`_)
    2. Le volume `V` de la boîte est fixe (car la pression n'est pas fixée - cf. l'`équation des gaz parfaits <https://fr.wikipedia.org/wiki/Loi_des_gaz_parfaits#.C3.89quation_des_gaz_parfaits>`_)
    3. La température `T` du système est fixé à l'aide d'un thermostat (e.g. `Thermostat de Berendsen <https://en.wikipedia.org/wiki/Berendsen_thermostat>`_)

C'est pourquoi on parle souvent d'équilibration `NVT`.

Comme pour la minimisation énergétique, :ref:`grompp` permet de générer la topologie complète::

    > gmx grompp -f nvt.mdp -c em.gro -p topol.top -o nvt.tpr

Le fichier de paramètre utilisé ici (:download:`nvt.mdp <./files/lyso/nvt.mdp>`) est évidemment adapté pour l'équilibration `NVT` avec en particulier:

    1. la contrainte des positions
    2. l'utilisation d'un thermostat à une température fixée (= couplage en température)
    3. l'absence de couplage en pression (= pression libre)

.. _posres:

Contrainte des positions
++++++++++++++++++++++++

Il est possible de contraindre les atomes afin qu'ils reste fixes. Dans le cas d'une équilibration,
il est souvent intéressant (voire nécessaire) de garder les grosses molécules (e.g. molécules) fixes afin de laisser les petites molécules (e.g. eau et ions), très mobiles, s'équilibrer.

Dans GROMACS, l'introduction de ces contraintes se fait généralement via l'utilisation de fichiers du type `posres.itp`
tels que celui-ci généré par :ref:`pdb2gmx` et dont le contenu ressemble à cela::

    [ position_restraints ]
    ; atom  type  fx    fy    fz
         1     1  1000  1000  1000
         5     1  1000  1000  1000
         6     1  1000  1000  1000
         7     1  1000  1000  1000
         8     1  1000  1000  1000
         9     1  1000  1000  1000
        10     1  1000  1000  1000
    ...

Comme le contenu le laisse supposer, ce fichier donne l'intensité des forces (`fx`, `fy` et `fz`) appliquées aux atomes (identifiés par leur index) pour les contraindre à leur positions.

Ce fichier n'est lu (et considéré pour la simulation) que sous certaines conditions précisées dans la topologie de la molécules (contenue dans `topol.top`)::

    #ifdef POSRES
    #include "posre.itp"
    #endif

Ici, les contraintes de positions ne sont incluses que si `POSRES` est défini. Or le fichier de paramètre `nvt.mdp` contient justement la ligne::

    define          = -DPOSRES

qui active (`-D` suivi de `POSRES`) la contrainte de positions.

.. _temperature_coupling:

Le couplage en température
++++++++++++++++++++++++++

Différents algorithmes sont disponibles afin de coupler le mouvement des atomes (agitation thermique) à la température.
Les paramètres associés doivent évidemment être décrits dans le fichier `.mdp` associé à la simulation.
En l'occurence, le fichier `nvt.mdp` décrit plus haut contient les lignes suivantes::

    tcoupl		= V-rescale
    tc-grps		= Protein Non-Protein
    tau_t		= 0.1	  0.1
    ref_t		= 363 	  363

Cela signifie que le thermostat utilisé est `V-rescale` [Bussi2007]_,
que la température est fixée à 363 K (~90 °C) de manière indépendante pour la protéine et pour le reste du système.


Absence de couplage en pression
+++++++++++++++++++++++++++++++

Le fichier `nvt.mdp` contient le réglage suivant pour le couplage en pression::

    pcoupl		= no

Ce qui confirme bien l'absence de couplage en pression.


Simulation NVT
++++++++++++++

Une fois la topologie complète prête, l'équilibration `NVT` est effectué à l'aide de :ref:`mdrun`::

    > gmx mdrun -v -deffnm nvt

.. important::
    Ne pas recopier (et exécuter) cette commande sans réflexion! (:ref:`Pourquoi ?<warning_mdrun>`)

.. note::
    Vérifier à l'aide de :ref:`energy` que la température est bien stable. Si ce n'est pas le cas, il faut continuer l'équilibration.

.. _equilibration_npt:

Équilibration NPT
-----------------

Une fois que la température est stable, on peut introduire le couplage en pression (toujours en contraignant la position des atomes lourds du lysozyme).

Le fichier de paramètre (:download:`npt.mdp <./files/lyso/npt.mdp>`) contient donc naturellement les réglages du couplage en pression::

    pcoupl		     = Parrinello-Rahman
    pcoupltype	     = isotropic
    tau_p		     = 2.0
    ref_p		     = 1.0
    compressibility  = 4.5e-5

L'algorithme de couplage (`Parinello-Rahman` [Nose1983]_) impose ici une pression isotrope d'un bar avec une compressibilité égale à celle de l'eau.

De ce fait, lors de cette équilibration:

    1. Le nombre `N` de particles reste toujours fixé (cf. `Conditions Périodiques aux Limites <https://fr.wikipedia.org/wiki/Condition_p%C3%A9riodique_aux_limites>`_)
    2. La température `T` du système est fixée (voir équilibration précédente)
    3. La pression `P` est maintenant fixée à l'aide d'un barostat

On parle donc d'équlibration `NPT`.

On peut également noter le réglage suivant dans `npt.mdp`::

    continuation	= yes

Ce paramètre permet d'indiquer que la simulation est la suite d'une simulation précédente: on va "continuer" l'équilibration en activant le couplage en pression. Cela se reflète dans l'appel à :ref:`grompp`::

    > gmx grompp -f npt.mdp -c nvt.gro -p topol.top -t nvt.cpt -o npt.tpr

Le fichier `nvt.cpt` est un fichier :ref:`checkpoint <simul_files>` qui contient toutes les informations nécessaire pour poursuivre la simulation (notamment les vitesses des atomes).
On le fourni à :ref:`grompp` pour permettre de faire l'équilibration `NPT` en se basant sur l'équilibration `NVT`.

On peut alors lancer l'équilibration `NPT` avec :ref:`mdrun`::

    > gmx mdrun -v -deffnm npt

.. important::
    Ne pas recopier (et exécuter) cette commande sans réflexion! (:ref:`Pourquoi ?<warning_mdrun>`)

.. note::
    Utiliser :ref:`energy` pour vérifier que la pression est bien stable. Pour cela, le plus pratique est de regarder la taille de la boîte.

Production & analyse
--------------------

Une fois que le système est équilibré, on peut lancer la "vraie" simulation, sans contraintes, qui va permettre d'acquérir les données que l'on analysera par la suite.

.. tip::
    Comme cette simulation sert à produire des données, on parle souvent de simulation de production.

Simulation de production
++++++++++++++++++++++++

La seule *vraie* différence entre le fichier de paramètre d'équilibration `NPT` (`npt.mdp`) et celui de production (:download:`md.mdp <./files/lyso/md.mdp>`) se situe
dans la disparation de la ligne activant les contraintes de positions. En effet, cette ligne est absente du fichier `md.mdp` puisqu'on veut, justement, que la simulation se fasse sans contraintes.

La génération de la topologie complète se fait comme précédemment::

    > gmx grompp -f md.mdp -c npt.gro -p topol.top -t npt.cpt -o md.tpr


Idem pour lancer la simulation::

    > gmx mdrun -v -deffnm md

.. important::
    Ne pas recopier (et exécuter) cette commande sans réflexion! (:ref:`Pourquoi ?<warning_mdrun>`)

Extension de la simulation
++++++++++++++++++++++++++

La durée simulée est indiquée implicitement dans le fichier `md.mdp` par les lignes::

    nsteps		= 500000
    dt		    = 0.002

`dt` correspond au pas d'intégration (ici 2 fs) et `nsteps` au nombre de pas d'intégration à simuler.
La durée simulée est donc de 500000 x 2 fs soit 1 ns.

Cette durée est malheureusement trop courte pour pouvoir observer une différence de comportement entre le lysozyme natif et le lysozyme sans ponts S-S.
Il faut donc poursuivre la simulation sur une durée plus longue (50 à 100 ns).
L'utilitaire :ref:`convert-tpr` est justement fait pour cela. Pour changer la durée de simulation jusqu'à une durée de 25 ns (soit 25000 ps), il suffit d'exécuter la commande suivante::

    > gmx convert-tpr -s md.tpr -o md.tpr -until 25000


Puis de lancer la suite de la simulation à l'aide de la commande suivante:

    > gmx mdrun -v -deffnm md -cpi md.cpt

.. important::
    Ne pas recopier (et exécuter) cette commande sans réflexion! (:ref:`Pourquoi ?<warning_mdrun>`)

.. important::
    L'ajout de `-cpi md.cpt` permet d'utiliser le fichier :ref:`checkpoint <simul_files>` comme point de départ plutôt
    que de devoir relancer la simulation depuis le début.

Analyse
+++++++

Pour des raisons d'efficacité, GROMACS place toujours tous les atomes dans une boîte de simulation "rectangulaire" pour effectuer les calculs.
Cela croit souvent des artéfacts quand, en particulier, une molécule est à cheval sur deux boîtes (la "vraie" et une réplique).

Il faut donc traiter la trajectoire issue de la simulation de production avant de pouvoir l'analyser.
Pour cela on utilise l'utilitaire :ref:`trjconv`::

    > gmx trjconv -f md.xtc -o md_mol.xtc -pbc mol -ur compact -s md.tpr

* l'option `-pbc mol` permet, fournissant une topologie (`-s md.tpr`), de faire en sorte que les molécules restent entières même si certains de leur atomes "sortent" de la boîte.
* l'option `-ur compact` permet d'utiliser la représentation "compacte" de la boîte (l'eau prendra ici une forme de `dodécaèdre <https://fr.wikipedia.org/wiki/Dod%C3%A9ca%C3%A8dre_r%C3%A9gulier>`_)

On peut maintenant analyser la trajectoire!

Pour rappel, le but est de voir l'effet des ponts S-S sur la stabilité structurale du lysozyme. Un bon indicateur de la stabilité d'une protéine est la
déviation standard des positions atomiques ou `RMSD <https://en.wikipedia.org/wiki/Root-mean-square_deviation_of_atomic_positions>`_.

L'utilitaire :ref:`rms` permet justement de suivre le RMSD au cours de la simulation. Pour cela, il suffit d'exécuter la commande::

    > gmx rms -f md_mol.xtc -s md.tpr -o rmsd.xvg

En premier, il faut spécifier le groupe d'atome de référence::

    Select group for least squares fit
    Group     0 (         System) has 20792 elements
    Group     1 (        Protein) has  1323 elements
    Group     2 (      Protein-H) has  1001 elements
    Group     3 (        C-alpha) has   129 elements
    Group     4 (       Backbone) has   387 elements
    Group     5 (      MainChain) has   517 elements
    Group     6 (   MainChain+Cb) has   634 elements
    Group     7 (    MainChain+H) has   646 elements
    Group     8 (      SideChain) has   677 elements
    Group     9 (    SideChain-H) has   484 elements
    Group    10 (    Prot-Masses) has  1323 elements
    Group    11 (    non-Protein) has 19469 elements
    Group    12 (          Water) has 19461 elements
    Group    13 (            SOL) has 19461 elements
    Group    14 (      non-Water) has  1331 elements
    Group    15 (            Ion) has     8 elements
    Group    16 (             CL) has     8 elements
    Group    17 ( Water_and_ions) has 19469 elements

Ici, on utiliser le squelette protéique (groupe `4`) comme référence.

Ensuite, il faut sélectionner le groupe d'atome dont on veut calculer le RMSD. On choisit également le groupe `4`.

Une fois le calcul effectué, :ref:`rms` génère un graphique (au format texte) du RMSD en fonction du temps de simulation: c'est le fichier `rmsd.xvg`.

Plus le RMSD est faible, plus la protéine est stable (car ses atomes restent proches de leur position initiale)


Les fichiers `.xvg` étant des fichiers texte, il peut être utile de les convertir en image.
Pour cela, un script Python (:download:`xvg2png.py <files/xvg2png.py>` est mis à disposition pour convertir un (ou plusieurs!) fichier `.xvg` en image (format `PNG <https://fr.wikipedia.org/wiki/Portable_Network_Graphics>`_).
Il suffit simplement d'exécuter le script de la façon suivante::

    > python xvg2png.py rmsd.xvg

Cela créera un fichier `rsmd.png` contenant le graphique au format `PNG`.


.. note::
    Tracer ce graphique pour les deux simulations (avec et sans pont S-S) et comparer les pour vérifier si les ponts S-S ont une influence sur la stabilité.

.. rubric:: Références

.. [Knubovets1999] Knubovets T, Osterhout JJ, Connolly PJ, Klibanov AM. Proc Natl Acad Sci. 1999. 96(4):1262–7.
.. [Bussi2007] Bussi G, Donadio D, Parrinello M. J Chem Phys. 2007. 126(1):14101.
.. [Nose1983] Nosé S, Klein ML. Mol Phys. 1983. 50(5):1055–76.
