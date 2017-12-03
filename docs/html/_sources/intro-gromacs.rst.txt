.. _presentation-gromacs:

Introduction à GROMACS
======================

Toutes les simulations de dynamique moléculaire présentées ici sont effectuées à l'aide de la suite logicielle `GROMACS`_
qui permet de préparer les simulations (`gmx grompp`), de les exécuter (`gmx mdrun`) et également de les analyser (e.g. `gmx energy`).

Utiliser GROMACS
----------------

.. important::
    GROMACS n'a pas d'interface graphique. Il est donc nécessaire de passer par l'invite de commande.

Toutes les commandes propres à GROMACS doivent être exécutées par l'intermédiaire de l'exécutable `gmx`.
Il est par exemple possible d'afficher l'aide de GROMACS grâce à::

    > gmx help

ou en utilisation l'option `-h` de la commande dont on veut l'aide.

Les fichiers GROMACS
--------------------

.. _coord_files:

.. _topol_files:

.. _simul_files:

.. _select_file:

.. _misc_files:

L'ensemble des fichiers supportés par GROMACS est disponible à l'adresse suivante:

    `http://manual.gromacs.org/online/files.html <http://manual.gromacs.org/online/files.html>`_



Outils GROMACS pour lancer une simulation
-----------------------------------------

.. important::
    Cette liste n'est pas exhaustive et ne reprend que les outils utilisés dans ce tutoriel.

.. _insert-molecules:

gmx insert-molecules
++++++++++++++++++++

Voir `http://manual.gromacs.org/programs/gmx-insert-molecules.html <http://manual.gromacs.org/programs/gmx-insert-molecules.html>`_


.. _pdb2gmx:

gmx pdb2gmx
+++++++++++

Voir `http://manual.gromacs.org/programs/gmx-pdb2gmx.html <http://manual.gromacs.org/programs/gmx-pdb2gmx.html>`_


.. _editconf:

gmx editconf
++++++++++++

Voir `http://manual.gromacs.org/programs/gmx-editconf.html <http://manual.gromacs.org/programs/gmx-editconf.html>`_


.. _solvate:

gmx solvate
+++++++++++

Voir `http://manual.gromacs.org/programs/gmx-solvate.html <http://manual.gromacs.org/programs/gmx-solvate.html>`_


.. _genion:

gmx genion
++++++++++

Voir `http://manual.gromacs.org/programs/gmx-genion.html <http://manual.gromacs.org/programs/gmx-genion.html>`_


.. _grompp:

gmx grompp
++++++++++

Voir `http://manual.gromacs.org/programs/gmx-grompp.html <http://manual.gromacs.org/programs/gmx-grompp.html>`_


.. _mdrun:

gmx mdrun
+++++++++

Voir `http://manual.gromacs.org/programs/gmx-mdrun.html <http://manual.gromacs.org/programs/gmx-mdrun.html>`_


.. _warning_mdrun:

.. important::
    La commande :ref:`mdrun` **ne doit jamais** être exécutée seule sur un supercalculateur!

    Par défaut, :ref:`mdrun` utilise *toutes* les ressources disponibles or, sur un supercalculateur, l'attribution des ressources se fait par l'intermédiaire d'un gestionnaire.
    Il faut donc **impérativement** passer par un gestionnaire de tâches (via un script de lancement).


.. _convert-tpr:

gmx convert-tpr
+++++++++++++++

Voir `http://manual.gromacs.org/programs/gmx-convert-tpr.html <http://manual.gromacs.org/programs/gmx-convert-tpr.html>`_


Outils GROMACS pour analyser une simulation
-------------------------------------------

.. important::
    GROMACS offre beaucoup d'outils pour l'analyse, la liste suivante ne contient que ceux utilisés dans ce tutoriel.

.. _energy:

gmx energy
++++++++++

Voir `http://manual.gromacs.org/programs/gmx-energy.html <http://manual.gromacs.org/programs/gmx-energy.html>`_


.. _trjconv:

gmx trjconv
+++++++++++

Voir `http://manual.gromacs.org/programs/gmx-trjconv.html <http://manual.gromacs.org/programs/gmx-trjconv.html>`_


.. _rms:

gmx rms
+++++++

Voir `http://manual.gromacs.org/programs/gmx-rms.html <http://manual.gromacs.org/programs/gmx-rms.html>`_


.. _density:

gmx density
+++++++++++

Voir `http://manual.gromacs.org/programs/gmx-density.html <http://manual.gromacs.org/programs/gmx-density.html>`_


.. _make_ndx:

gmx make_ndx
++++++++++++

Voir `http://manual.gromacs.org/programs/gmx-make_ndx.html <http://manual.gromacs.org/programs/gmx-make_ndx.html>`_


.. _msd:

gmx msd
+++++++

Voir `http://manual.gromacs.org/programs/gmx-msd.html <http://manual.gromacs.org/programs/gmx-msd.html>`_


.. _GROMACS: http://www.gromacs.org