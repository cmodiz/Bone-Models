Models
===========================

This module contains computational models relevant to bone remodelling.
The :ref:`base_model` is the base class for all models and contains the most general representation of a bone cell population models.

The :ref:`lemaire_model` contains the bone cell population model based on the work of Lemaire et al. (2004) and is a submodule of the Base Model.

The :ref:`pivonka_model` contains the bone cell population model based on the work of Pivonka et al. (2008) and is a
submodule of the Lemaire Model. It adds more details to the intracellular communication pathways.

The :ref:`scheiner_model` contains the bone cell population model based on the work of Scheiner et al. (2013) and is a submodule of the Pivonka Model.
It adds the effect of mechanical loading on the bone cell population model using a homogenisation framework to calculate
microscopic strain energy density from macroscopic stress tensors.

The :ref:`modiz_model` contains the bone cell population model based on the work of Modiz et al. (2025) and is a submodule of the Lemaire Model.
It changes the activation of osteoblasts by parathyroid hormone (PTH) to include details about pulse characteristics and two states of the receptor.
The extension is based on the :ref:`martonova_model`.

The :ref:`martonova_model` calculates cellular activity constants using pulsatile PTH based on the work of Martonova et al. (2023).


.. _base_model:
Base Model
--------------------------------------

.. automodule:: bone_models.models.base_model
   :members:
   :undoc-members:
   :show-inheritance:

.. _lemaire_model:
Lemaire Model
-----------------------------------------

.. automodule:: bone_models.models.lemaire_model
   :members:
   :undoc-members:
   :show-inheritance:

.. _martonova_model:
Martonova Model
-------------------------------------------

.. automodule:: bone_models.models.martonova_model
   :members:
   :undoc-members:
   :show-inheritance:

.. _modiz_model:
Modiz Model
---------------------------------------

.. automodule:: bone_models.models.modiz_model
   :members:
   :undoc-members:
   :show-inheritance:

.. _pivonka_model:
Pivonka Model
-----------------------------------------

.. automodule:: bone_models.models.pivonka_model
   :members:
   :undoc-members:
   :show-inheritance:

.. _scheiner_model:
Scheiner Model
------------------------------------------

.. automodule:: bone_models.models.scheiner_model
   :members:
   :undoc-members:
   :show-inheritance:
