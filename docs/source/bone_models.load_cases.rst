Load Cases
================================

This module contains the load cases used in the bone models. The load cases are used to implement e.g. disease states
or different loading conditions for a certain time interval.

The :ref:`lemaire_load_cases` contains the load cases based on the work of Lemaire et al. (2004) to simulate respective
disease states with increase/decreased concentrations.

The :ref:`pivonka_load_cases` contains the load cases based on the work of Pivonka et al. (2008) to simulate respective
disease states with increase/decreased concentrations.

The :ref:`modiz_load_cases` contains the load cases based on the work of Modiz et al. (2025) to simulate respective
disease states using altered pulse characteristics. It includes both :ref:`lemaire_load_cases` and :ref:`martonova_load_cases`.

The :ref:`martonova_load_cases` contains the load cases based on the work of Martonova et al. (2023) to simulate respective
disease states using altered pulse characteristics. Drug injections are also included.

The :ref:`scheiner_load_cases` contains the load cases based on the work of Scheiner et al. (2021) to simulate disuse
or overuse scenarios based on altered stress tensors. Former load cases (Lemaire et al., Pivonka et al.) can be included.


.. _lemaire_load_cases:
Lemaire Load Cases
----------------------------------------------------

This module contains the load cases for the bone cell population model by Lemaire et al. (2004).
The load cases are used to define disease conditions phenomenologically using external injections of respective
concentrations to either increase or decrease these concentrations.

.. automodule:: bone_models.load_cases.lemaire_load_cases
   :members:
   :undoc-members:
   :show-inheritance:

.. _martonova_load_cases:
Martonova Load Cases
------------------------------------------------------

This module contains the load cases for the two-state receptor model by Martonova et al. (2023).
The load cases are used to define disease conditions based on changed pulsatile characteristics with or without drug
administration.

.. automodule:: bone_models.load_cases.martonova_load_cases
   :members:
   :undoc-members:
   :show-inheritance:

.. _modiz_load_cases:
Modiz Load Cases
--------------------------------------------------

This module contains the load cases for the bone cell population model by Modiz et al. (2025).
The load cases are used to define disease conditions based on changed pulsatile characteristics for the cellular
activity constants. For the bone cell population model by Lemaire et al., (2004) used as a reference the load cases
define disease states phenomenologically. In conclusion, the load cases contain both Lemaire and Martonova load cases.


.. automodule:: bone_models.load_cases.modiz_load_cases
   :members:
   :undoc-members:
   :show-inheritance:

.. _pivonka_load_cases:
Pivonka Load Cases
----------------------------------------------------

This module contains the load cases for the bone cell population model by Pivonka et al. (2008).
The load cases are used to define disease conditions phenomenologically using external injections of respective
concentrations to either increase or decrease these concentrations.

.. automodule:: bone_models.load_cases.pivonka_load_cases
   :members:
   :undoc-members:
   :show-inheritance:


.. _scheiner_load_cases:
Scheiner Load Cases
----------------------------------------------------

This module contains the load cases for the bone cell population model by Scheiner et al. (2013).
The load cases are used to define disuse or overuse conditions during a certain time interval.
The load cases from Lemaire et al. or Pivonka et al., can also be included.

.. automodule:: bone_models.load_cases.scheiner_load_cases
   :members:
   :undoc-members:
   :show-inheritance:

