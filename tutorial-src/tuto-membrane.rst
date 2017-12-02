Tutoriel "Membrane"
===================

.. note::
    Ce tutoriel est inspiré du tutoriel disponible `ici <http://md.chem.rug.nl/index.php/tutorials-general-introduction-gmx5/bilayers-gmx5#Bilayer-self-assembly>`_.

Le but de ce tutoriel est d'observer comment les lipides membranaires s'auto-assemblent pour former une bicouche qui constitue la matrice des membranes biologiques.

En utilisant deux lipides différents (un saturé et un insaturé), le but est aussi d'étudier l'influence de l'insaturation des chaînes sur les propriétés de la membrane.

Création de la topologie
------------------------

Afin d'accélérer les simulations, le champ de force utilisé ici sera de type "gros grains". Il s'agit du champ de force `Martini <http://md.chem.rug.nl/>`_.

Dans ce tutoriel, il n'y a pas besoin d'utiliser :ref:`pdb2gmx` (qui ne marche qu'avec les protéines) car la topologie des lipides est directement disponible.


.. important::

    Dans la suite, le système "DPPC" est pris comme base pour les exemples.

    La transposition du DPPC au DOPC est laissé comme exercice.


Construction du système initial
-------------------------------

Contrairement au :ref:`tutoriel lysozyme <tuto_lyso>`, ici il n'existe pas de fichier :ref:`pdb <coord_files>` à partir duquel le système peut être généré.
Il est donc nécessaire de générer le system *ex nihilo*.
Pour cela, on va utiliser l'utilitaire :ref:`insert-molecules`::

    > gmx insert-molecules -ci DPPC-em.gro -box 15 15 7.5 -nmol 512 -radius 0.21 -try 500 -o system_noW.gro

Cet utilitaire permet, ici, d'insérer aléatoirement 512 (option `-nmol`) molécules de DPPC (option `-ci`) dans une nouvelle boîte MD de taille 15 nm x 15 nm x 7.5nm (option `-box`).
Afin d'éviter les collisions entre les molécules insérées, un rayon de 0.21 nm (option `-radius`) est pris en compte et 500 essais d'insertion (par molécule insérée - option `-try) est effectué``
Quand les 512 molécules sont insérées (ou que les 500x512 essais sont passés), le système est sauvegardé dans `system_noW.gro` (option `-o`).

Construction du fichier de topologie du système
-----------------------------------------------

Malheureusement, :ref:`insert-molecules` ne génère que les coordonnées du système (fichier :ref:`.gro <coord_files>`).
Il faut donc créé le fichier de topologie. Pour simplifier les choses, il est possible de télécharger ce fichier (
:download:`system.top <./files/membrane/system.top>`) dont voici le contenu:

.. code-block:: ahk
    :linenos:

    #include "martini_v2.1.itp"
    #include "martini_v2.0_DPPC_01.itp"

    [ system ]
    DPPC BILAYER SELF-ASSEMBLY

    [ molecules ]
    DPPC 512
    ; W 3072

Ce fichier est très simple à comprendre:

* Le premier `#include` (ligne 1) indique l'emplacement du champ de force :download:`martini_v2.1.itp <./files/membrane/martini_v2.1.itp>`.
* Ensuite (ligne 2), c'est la topologie du lipide :download:`martini_v2.0_DPPC_01.itp <./files/membrane/martini_v2.0_DPPC_01.itp>` qui est incluse.

.. important::
    Pour que `system.top` puisse être interprété par :ref:`grompp`, il faut que *tous* les fichiers inclus soit effectivement à l'emplacement indiqué.
    (Ici, il s'agit du même dossier que `system.top`.

.. note::
    Ouvrir `martini_v2.1.itp` et `martini_v2.0_DPPC_01.itp` pour observer les différences entre un fichier de définition de champ de forces et de topologie d'une molécule.

* la section `[ system ]` contient la description du système qui a simplement une valeur indicative.
* la section `[ molecules ]` contient la liste des molécules présentes dans le système.

.. important::
    Les molécules décrites dans la section `[ molecules ]` doivent:
        1. être obligatoirement définies (directement dans le champ de forces ou dans un fichier de topologie `.itp`).
        2. correspondre à celles présentes dans le fichier `.gro` associé.


Solvation des lipides
---------------------

Comme dans tutoriel :ref:`précédent <tuto_lyso>`, il est nécessaire d'ajouter l'eau puisque seuls les lipides sont présents pour l'instant::

    > gmx solvate -cp system_noW.gro -cs water.gro -o system_W.gro -maxsol 3072 -radius 0.21

Ici, on ajoute 3072 "billes" d'eau; 1 "bille" Martini correspondant à 4 molécules d'eau, cela fait 3072x4=12288 molécules d'eau soit 24 molécules d'eau pour un lipide.

.. important::
    Il faut éditer `system.top` pour y ajouter les 3072 particules d'eau (nom martini: `W`)


Minimisation énergétique
------------------------

Le système étant généré aléatoirement, il est nécessaire de minimiser l'énergie avant toute simulation de dynamique moléculaire::

    > gmx grompp -f minimization.mdp -c system_W.gro -p system.top -o em.tpr

    > gmx mdrun -v -deffnm em

.. important::
    Ne pas recopier (et exécuter) cette commande sans réflexion! (:ref:`Pourquoi ?<warning_mdrun>`)


Simulation de l'auto-assemblage
-------------------------------

La description mésoscopique (i.e. "gros grains") du système permet, cans ce cas, de se passer de la phase d'équilibration `NVT` et `NPT`.
On peut lancer directement la simulation de l'auto-assemblage des lipides en bicouche en utilisant des paramètres (:download:`martini_md.mdp <./files/membrane/martini_md.mdp>`) spécifiques à Martini::

    > gmx grompp -f martini_md.mdp -c em.gro -p system.top -o md.tpr

    > gmx mdrun -v -deffnm md

.. important::
    Ne pas recopier (et exécuter) cette commande sans réflexion! (:ref:`Pourquoi ?<warning_mdrun>`)

.. note::
    Ouvrir le fichier `martini_md.mdp` pour voir quels paramètres sont différents par rapport au fichier utilisé pour simuler le lysozyme (:download:`md.mdp <./files/lyso/md.mdp>`).


Analyse de la bicouche formée
-----------------------------

Avant d'analyser la trajectoire, il est nécessaire de rendre entières les molécules "cassées" par la réplication de la boîte::

    > gmx trjconv -s md.xtc -o md_mol.xtc -pbc mol -s md.tpr

Aire par lipide
+++++++++++++++

Un paramètre important quand on simule une membrane est l'aire occupée par un lipide.
Quand la membrane est plane, il suffit simplement de diviser l'aire du plan XY de la boîte de MD par le nombre de lipides par feuillet.
Le fichier `md.edr` contient justement les valeurs des axes de la boîte que l'on peut extraire à l'aide :ref:`energy`::

    > gmx energy -f md.edr -o box-x.xvg

Sachant que le couplage en pression est anisotrope (cf `martini_md.mdp`), l'axe X et l'axe Y de la boîte sont strictement identiques
et il suffit d'extraire seulement la coordonnée X de la boîte (`Box-X`) pour pouvoir calculer l'aire du plan XY::

    Select the terms you want from the following list by
    selecting either (part of) the name or the number or a combination.
    End your selection with an empty line or a zero.
    -------------------------------------------------------------------
      1  Bond             2  G96Angle         3  LJ-(SR)          4  Coulomb-(SR)
      5  Potential        6  Kinetic-En.      7  Total-Energy     8  Temperature
      9  Pressure        10  Box-X           11  Box-Y           12  Box-Z
     13  Volume          14  Density         15  pV              16  Enthalpy
     17  Vir-XX          18  Vir-XY          19  Vir-XZ          20  Vir-YX
     21  Vir-YY          22  Vir-YZ          23  Vir-ZX          24  Vir-ZY
     25  Vir-ZZ          26  Pres-XX         27  Pres-XY         28  Pres-XZ
     29  Pres-YX         30  Pres-YY         31  Pres-YZ         32  Pres-ZX
     33  Pres-ZY                             34  Pres-ZZ
     35  #Surf*SurfTen                       36  Coul-SR:DOPC-DOPC
     37  LJ-SR:DOPC-DOPC                     38  Coul-SR:DOPC-W
     39  LJ-SR:DOPC-W    40  Coul-SR:W-W     41  LJ-SR:W-W       42  T-DOPC
     43  T-W             44  Lamb-DOPC       45  Lamb-W

Il faut donc sélectionner seulement le terme `10`.

Le fichier ainsi créé (`box-x.xvg`) contient la longueur de l'axe X de la boîte au cours de la trajectoire.

.. note::
    Inspecter le fichier `box-x.xvg` (il s'agit d'un fichier texte)

Il faut donc retraiter les valeurs contenues dans ce fichier pour pouvoir calculer l'aire par lipide au cours de la simulation.
Pour cela, on va utiliser un script Python (:download:`box2apl.py <files/membrane/box2apl.py>`) qui va lire le fichier `box-x.xvg`
et calculer l'aire par lipide::

    > python box2apl.py box-x.xvg 256

.. important::
    Il faut fournir le nombre de lipides par feuillet (ici 256) au script car cette information n'est pas disponible dans `box-x.xvg`.

Le script crée un fichier `PNG` (`box-x_APL.png`) contenant le graphique correspondant.

.. note::
    Comparer l'aire par lipide entre DPPC et DOPC.


Epaisseur de la bicouche
++++++++++++++++++++++++

Pour déterminer l'épaisseur, on va utiliser l'utilitaire :ref:`density` qui permet de tracer la densité des atomes suivant un axe (ici Z).
Classiquement, on détermine l'épaisseur membranaire à partir de la distance en les atomes de phosphore de chacun des feuillets.

La première chose à faire est donc de créer un fichier `.ndx` (cf :ref:`topol_files`) contenant un groupe d'atome correspondant aux atomes de phosophore.
:ref:`make_ndx` permet de faire cela::

    > gmx make_ndx -f md.gro

Par défaut, GROMACS reconnait et crée un certain nombre de groupes::

    0 System              :  9216 atoms
    1 Other               :  9216 atoms
    2 DOPC                :  6144 atoms
    3 W                   :  3072 atoms

Il faut créer un nouveau groupe pour les phosphate. Pour cela on va sélectionner les atomes (`a`) qui s'appellent `PO4` (nom de la bille phosphate dans Martini)::

    > a PO4

    Found 512 atoms with name PO4

      4 PO4                 :   512 atoms

Une fois :ref:`make_ndx` quitté (commande `q`), un fichier `index.ndx` est créé et il contient la sélection d'atome.

On peut alors se servir de :ref:`density` pour calculer le profil de densité::

    > gmx density -f md.xtc -s md.tpr -b 50000 -n index.ndx -d z -o p-density.xvg

.. important::
    Il faut évidemment choisir le groupe `PO4` (numéro `4`) pour avoir la densité des phosphates...

On peut évaluer l'épaisseur membranaire à partir de la distance entre les deux pics de densité correspondant aux deux feuillets.

.. note::
    Comparer les valeurs d'épaisseur membranaire entre le DPPC et le DOPC.


Diffusion latérale
++++++++++++++++++

La diffusion latérale des lipides peut se calculer facilement avec GROMACS à partir du moment où il est possible d'enlever
les "sauts" des molécules::

    > gmx trjconv -f md.xtc -s md.tpr -pbc nojump -o nojump.xtc

La nouvelle trajectoire `nojump.xtc` ne contenant plus aucun saut de molécules, on peut calculer la diffusion des lipides à l'aide de
:ref:`msd`::

    > gmx msd -f nojump.xtc -s md.tpr -rmcomm -lateral z -b 50000

.. note::
    Comparer les diffusion latérale des 2 lipides.


Paramètre d'ordre
+++++++++++++++++

Enfin, le paramètre d'ordre peut être calculé à l'aide du script Python :download:`do-order-gmx5.py <./files/membrane/do-order-gmx5.py>`::

    > python do-order-gmx5.py md.xtc md.tpr 15000 30000 20 0 0 1 512 DPPC


