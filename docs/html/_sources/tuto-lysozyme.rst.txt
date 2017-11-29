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
    :height: 400px
    :align: center

    Structure du lysozyme avec ses 4 ponts S-S.

Création de la topologie
------------------------

La structure du lysozyme est bien connue et est disponible sur le site de la `PDB <http://www.rcsb.org/pdb/home/home.do>`_. Le code PDB pour le lysozyme de blanc d'oeuf est: `1AKI <http://www.rcsb.org/pdb/explore/explore.do?structureId=1aki>`_.
La première étape est donc de télécharger les coordonnées atomiques sous le format `.pdb`::

    wget https://files.rcsb.org/download/1AKI.pdb

.. note::

    Ouvrir le fichier `1AKI.pdb` pour voir toutes les informations qu'il contient.

Le fichier `.pdb` ne contient que les coordonnées des atomes mais ne contient pas la topologie de la molécule qui, elle, est liée à un champ de forces.
Pour générer cette topologie, on peut utiliser l'utilitaire `pdb2gmx` qui permet de lire un fichier `.pdb` et de générer automatiquement la topologie associée.

.. important::
    Comme le but est d'étudier l'effet des ponts S-S, nous allons générer *2* topologies:
        1. une topologie contenant les ponts S-S
        2. une topologie ne contenant *pas* les ponts S-S

La commande pour générer la topologie du lysozyme est::

    gmx pdb2gmx -f 1AKI.pdb -ss

L'option `-ss` rend interactif le choix des ponts S-S.

Premièrement `pdb2gmx` demande de choisir le champ de force::

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

A partir de ces réponses, `pdb2gmx` génère 3 fichiers:

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

Premièrement il est nécessaire de créer la boite de simulation autour de notre protéine::

    gmx editconf -f conf.gro -o lyso.gro -d 0.7 -bt dodecahedron

.. note::
    Ouvrir `lyso.gro` avec VMD pour observer la boîte de simulation.

gmx solvate -p system_noSS.top -cp lyso_noSS.gro -cs -o lyso_noSS_W.gro


Ajout des ions
--------------

gmx grompp -f ../ions.mdp -c lyso_noSS_W.gro -p system_noSS.top -o ions.tpr

gmx genion -s ions.tpr -o lyso_noSS_WI.gro -p system_noSS.top -neutral


Minimisation énergétique
------------------------

gmx grompp -f ../minim.mdp -c lyso_noSS_WI.gro -p system_noSS.top -o em.tpr

gmx mdrun -v -deffnm em


Équilibration NVT
-----------------

gmx grompp -f ../nvt.mdp -c em.gro -p system_noSS.top -o nvt.tpr

gmx mdrun -v -deffnm nvt


Équilibration NPT
-----------------

gmx grompp -f ../npt.mdp -c nvt.gro -p system_noSS.top -t nvt.cpt -o npt.tpr

gmx mdrun -v -deffnm npt


Production & analyse
--------------------


Simulation de production
++++++++++++++++++++++++

gmx grompp -f ../md.mdp -c npt.gro -p system_noSS.top -t npt.cpt -o md.tpr

gmx mdrun -v -deffnm md


Analyse
+++++++

gmx trjconv -f md.xtc -fit rot+trans -o md_fit.xtc -s md.tpr

gmx trjconv -f md_fit.xtc -o md_fit.xtc -pbc mol -ur compact -s md.tpr

gmx rms -f md_fit.xtc -s md.tpr -o rmsd.xvg

.. rubric:: Références

.. [Knubovets1999] Knubovets T, Osterhout JJ, Connolly PJ, Klibanov AM. Proc Natl Acad Sci. 1999. 96(4):1262–7.
