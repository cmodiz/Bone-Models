Models
===========================

This module contains computational models relevant to bone remodelling analysis at macro- and organ-level scales.
The :ref:`lerebours_model` implements the spatial component of the bone remodelling model based on the work of Lerebours et al. (2016).
It initialises a cross-section and simulates the spatial distribution of bone volume fractions over time under various remodelling conditions.
The model calls the respective bone cell population model for each RVE.

.. _lerebours_model:

Lerebours Model
--------------------------------------

.. automodule:: bone_models.bone_spatial_models.models.lerebours_model
   :members:
   :undoc-members:
   :show-inheritance:

