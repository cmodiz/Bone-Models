Load Cases
================================

This module contains the load cases used in the bone models. The load cases are used to implement e.g. disease states
or different loading conditions for a certain time interval.

.. _lemaire_load_cases:
Lemaire Load Cases
----------------------------------------------------

This module contains the load cases for the bone cell population model by Lemaire et al. (2004).
The load cases are used to define disease conditions phenomenologically using external injections of respective
concentrations to either increase or decrease these concentrations.

.. automodule:: bone_models.bone_cell_population_models.load_cases.lemaire_load_cases
   :members:
   :undoc-members:
   :show-inheritance:

.. _martonova_load_cases:
Martonova Load Cases
------------------------------------------------------

This module contains the load cases for the two-state receptor model by Martonova et al. (2023).
The load cases are used to define disease conditions based on changed pulsatile characteristics with or without drug
administration.

.. automodule:: bone_models.bone_cell_population_models.load_cases.martonova_load_cases
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


.. automodule:: bone_models.bone_cell_population_models.load_cases.modiz_load_cases
   :members:
   :undoc-members:
   :show-inheritance:

.. _pivonka_load_cases:
Pivonka Load Cases
----------------------------------------------------

This module contains the load cases for the bone cell population model by Pivonka et al. (2008).
The load cases are used to define disease conditions phenomenologically using external injections of respective
concentrations to either increase or decrease these concentrations.

.. automodule:: bone_models.bone_cell_population_models.load_cases.pivonka_load_cases
   :members:
   :undoc-members:
   :show-inheritance:


.. _scheiner_load_cases:
Scheiner Load Cases
----------------------------------------------------

This module contains the load cases for the bone cell population model by Scheiner et al. (2013).
The load cases are used to define disuse or overuse conditions during a certain time interval.
The load cases from Lemaire et al. or Pivonka et al., can also be included.

.. automodule:: bone_models.bone_cell_population_models.load_cases.scheiner_load_cases
   :members:
   :undoc-members:
   :show-inheritance:

.. _martinez_reina_load_cases:
Martinez-Reina Load Cases
----------------------------------------------------

This module contains the load cases for the bone cell population model by Martinez-Reina et al. (2019).
The load cases are used to define disuse or overuse conditions, PMO and denosumab treatment during certain time intervals.
The load cases from Lemaire et al. or Pivonka et al., can also be included.

.. automodule:: bone_models.bone_cell_population_models.load_cases.martinez_reina_load_cases
   :members:
   :undoc-members:
   :show-inheritance:

