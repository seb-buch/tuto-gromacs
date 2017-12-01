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


Construction du système initial (DPPC)
--------------------------------------

Pour cela, on utilise l'utilitaire :ref:`insert-molecules`::

    > gmx insert-molecules -ci DPPC-em.gro -box 15 15 7.5 -nmol 512 -radius 0.21 -try 500 -o system_noW.gro

Construction du fichier de topologie du système (DPPC)
------------------------------------------------------

Création/édtion du fichier `system.top`.

Solvation des lipides (DPPC)
----------------------------

gmx solvate -cp system_noW.gro -cs water.gro -o system_W.gro -maxsol 3072 -radius 0.21

Minimisation énergétique (DPPC)
-------------------------------

gmx grompp -f minimization.mdp -c system_W.gro -p system.top -o em.tpr

gmx mdrun -v -deffnm em

Simulation de l'auto-assemblage (DPPC)
--------------------------------------

gmx grompp -f martini_md.mdp -c em.gro -p system.top -o md.tpr

gmx mdrun -v -deffnm md


Analyse de la bicouche formée
-----------------------------

gmx trjconv -s md.xtc -o md_mol.xtc -pbc mol -s md.tpr

Aire par lipide
+++++++++++++++

gmx energy -f md.edr -o box-xy.xvg

TODO: script python xvg/512 -> png

Epaisseur de la bicouche
++++++++++++++++++++++++

gmx make_ndx -f md.gro

    > a P*

    > q

$ gmx density -f md.xtc -s md.tpr -b 15000 -n index.ndx -o p-density.xvg -xvg no

    > 4



Diffusion latérale
++++++++++++++++++

gmx trjconv -f md.xtc -s md.tpr -pbc nojump -o nojump.xtc

gmx msd -f nojump.xtc -s md.tpr -rmcomm -lateral z -b 15000

Paramètre d'ordre
+++++++++++++++++

python do-order-gmx5.py md.xtc md.tpr 15000 30000 20 0 0 1 512 DPPC


