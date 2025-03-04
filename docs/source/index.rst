.. bone_models documentation master file, created by
   sphinx-quickstart on Sun Feb  9 15:21:51 2025.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to Bone Models documentation!
=========================================

This documentation provides an in-depth guide to the **Bone Cell Population Models** used for simulating bone remodeling
and disease progression. The models implemented here are based on established theoretical frameworks from the literature
and aim to describe the complex interactions between osteoblasts, osteoclasts, and biochemical signaling pathways under
various physiological and pathological conditions.

This documentation is continuously updated as new models and features are integrated. If you have any questions or suggestions, feel free to contribute!

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   bone_models

What's included?
--------------

This repository includes multiple computational models, each extending previous work to incorporate additional biological mechanisms:

- **Base Model** (:ref:`base_model`): A general framework for bone cell population dynamics.
- **Lemaire Model** (:ref:`lemaire_model`): Implements the model by Lemaire et al. (2004), describing fundamental osteoblast-osteoclast interactions.
- **Pivonka Model** (:ref:`pivonka_model`): Extends the Lemaire Model with intracellular signaling based on Pivonka et al. (2008).
- **Scheiner Model** (:ref:`scheiner_model`): Further expands the Pivonka Model to include mechanical loading effects (Scheiner et al., 2013).
- **Martinez-Reina Model** (:ref:`martinez_reina_model`): Expands the Scheiner Model to account for trabecular bone disuse, overuse, PMO and denosumab treatment scenarios (Martinez-Reina et al., 2019).
- **Martonova Model** (:ref:`martonova_model`): Computes cellular activity constants based on pulsatile PTH (Martonova et al., 2023).
- **Modiz Model** (:ref:`modiz_model`): Introduces receptor dynamics and pulsatile PTH effects (Modiz et al., 2025), building upon the Lemaire Model and the Martonova Model.

Each model is accompanied by a set of **load cases** that define specific disease states or disuse/ overuse scenarios by altering the
model parameters. The model **parameters** are named in a consistent manner across all models, allowing for easy comparison and calibration.
In the corresponding parameter files, the parameters are matched with the nomenclature used in the original publications.
On a separate page, **examples** are provided to demonstrate how to use the models and load cases to simulate different scenarios and analyze the results.

Getting Started
---------------

To start using the models, install the required dependencies and follow the setup instructions.
You can also find example scripts to help you implement simulations and analyze results.

To install the package, simply use pip:

.. code-block:: bash

    pip install bone_models

Once installed, you can import the available models, load cases, and run simulations.
Each bone cell population model is implemented as a class that automatically initializes the required parameters.
If needed, the parameter files can be accessed and modified.
The models can be imported and instantiated with a specific load case as shown below for the Lemaire and Pivonka models:

.. code-block:: python

    # Import models
    from bone_models.models import Lemaire_Model, Pivonka_Model

    # Optional: if needed, parameters can be imported via the following
    from bone_models.parameters.lemaire_parameters import Lemaire_Parameters

    # Import load cases (for some models there is more than one load case available)
    from bone_models.load_cases.lemaire_load_cases import Lemaire_Load_Case_3
    from bone_models.load_cases.pivonka_load_cases import Pivonka_Load_Case_2

    # Initialize a model instance with a load case
    lemaire_model = Lemaire_Model(Lemaire_Load_Case_3)
    pivonka_model = Pivonka_Model(Pivonka_Load_Case_2)


All bone cell population models are solved using the `solve_bone_cell_population_model` function. This function requires a time span (`tspan`), which specifies the simulation duration.

.. code-block:: python

    # Define time span for simulation
    tspan = [0, 1000]  # Start at time 0, end at time 1000 days

    # Solve the model
    lemaire_solution = lemaire_model.solve_bone_cell_population_model(tspan)
    pivonka_solution = pivonka_model.solve_bone_cell_population_model(tspan)

After solving, the model contains the simulation results, including the evolution of cell populations over time.
More detailed examples, including plots and analysis of the simulation results, can be found on the :ref:`bone_models.examples` page.



