Parameters
===============================

This module defines the parameters used in the bone cell population models. Each model has a specific set of parameters
governing its behavior, interaction dynamics, and biological processes.
For the bone cell population models, the parameters are named consistently across models to facilitate comparison and understanding.
In each parameter file, the parameter names are matched with the corresponding model nomenclature and units.

The :ref:`base_parameters` define the most general set of parameters for bone cell population models.

The :ref:`lemaire_parameters` correspond to the parameters for the model by Lemaire et al. (2004).

The :ref:`pivonka_parameters` correspond to the parameters for the model by Pivonka et al. (2008).

The :ref:`scheiner_parameters` correspond to the parameters for the model by Scheiner et al. (2013).

The :ref:`modiz_parameters`  correspond to the parameters for the model by Modiz et al. (2025).

The :ref:`martonova_parameters` correspond to the parameters for the model by Martonova et al. (2023).

.. _base_parameters:
Base Parameters
-----------------------------------------------
This module contains the base parameters for the bone cell population models.
The base model is not intended to be used as a standalone model but as a base class for other models.
The base parameters are thus all set to None; a value is assigned to them in the derived models.

.. automodule:: bone_models.parameters.base_parameters
   :members:
   :undoc-members:
   :show-inheritance:

.. _lemaire_parameters:
Lemaire Parameters
--------------------------------------------------
This module contains the parameters for the model by Lemaire et al. (2004).
In each parameter class, the parameter names are matched with the corresponding model nomenclature and units.

.. automodule:: bone_models.parameters.lemaire_parameters
   :members:
   :undoc-members:
   :show-inheritance:

.. _martonova_parameters:
Martonova Parameters
----------------------------------------------------
This module contains the parameters for the model by Martonova et al. (2023).
In each parameter class, the parameter names are matched with the corresponding model nomenclature and units.

.. automodule:: bone_models.parameters.martonova_parameters
   :members:
   :undoc-members:
   :show-inheritance:

.. _modiz_parameters:
Modiz Parameters
------------------------------------------------
This module contains the parameters for the model by Modiz et al. (2025).
In each parameter class, the parameter names are matched with the corresponding model nomenclature and units.

.. automodule:: bone_models.parameters.modiz_parameters
   :members:
   :undoc-members:
   :show-inheritance:

.. _pivonka_parameters:
Pivonka Parameters
--------------------------------------------------
This module contains the parameters for the model by Pivonka et al. (2008).
In each parameter class, the parameter names are matched with the corresponding model nomenclature and units.

.. automodule:: bone_models.parameters.pivonka_parameters
   :members:
   :undoc-members:
   :show-inheritance:

.. _scheiner_parameters:
Scheiner Parameters
---------------------------------------------------
This module contains the parameters for the model by Scheiner et al. (2013).
In each parameter class, the parameter names are matched with the corresponding model nomenclature and units.

.. automodule:: bone_models.parameters.scheiner_parameters
   :members:
   :undoc-members:
   :show-inheritance:
