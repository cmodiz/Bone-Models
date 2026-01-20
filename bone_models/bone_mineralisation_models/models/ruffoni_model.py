import numpy as np
import fipy as fp
from bone_models.bone_mineralisation_models.parameters.ruffoni_parameters import Ruffoni_Parameters
import matplotlib.pyplot as plt
import logging
import scipy as sc
from scipy.stats import skewnorm
from scipy.interpolate import PchipInterpolator
from scipy.ndimage import gaussian_filter1d

log = logging.getLogger("fipy")
console = logging.StreamHandler()
console.setLevel(logging.INFO)
log.addHandler(console)


class Ruffoni_Model:
    r"""
    Implements the Ruffoni et al. (2007) model for Bone Mineralization Density Distribution (BMDD)
    using FiPy's Finite Volume Method for robust PDE solving.

    .. note::
       **Source Publication**:
       Ruffoni D., Fratzl P., Roschger P., Klaushofer K., Weinkamer R. (2007).
       *The bone mineralization density distribution as a fingerprint of the mineralization process.*
       **Bone**, 40, 1308–1319.
       :doi:`10.1016/j.bone.2007.01.012`

    The model solves the following partial differential equation:

    .. math::

       \frac{\partial \rho(C, t)}{\partial t} +
       \frac{\partial \left( \rho(C, t) \, v_m(C) \right)}{\partial C}
       = - r(C, t) \, \rho(C, t)

    with boundary condition:

    .. math::

       \rho(0, t) \, v_m(0) = f(t)

    Where:

    - :math:`\rho(C, t)` — bone mineral density distribution (BMDD)
    - :math:`v_m(C)` — mineralization velocity
    - :math:`r(C, t)` — resorption rate at time :math:`t`
    - :math:`f(t)` — formation rate of new bone at time :math:`t`
    """

    def __init__(self, simulation_time, number_of_grid_points, start='BMDD'):
        """
        Initialize the Ruffoni BMDD model:

        There are two options to initialize the steady state, determined by the ``start`` parameter:

        - **'BMDD'** — Initialize from a reference BMDD distribution.
          The mineralization velocity and mineralization law are calculated from this steady-state BMDD using analytical formulas.
        - **'mineralization law'** — Initialize from a predefined mineralization law.
          The BMDD and mineralization velocity are calculated from this mineralization law using analytical formulas.

        :param simulation_time: Total simulation time in years.
        :type simulation_time: float
        :param number_of_grid_points: Number of grid points for the 1D calcium mesh.
        :type number_of_grid_points: int
        :param start: Initialization method, either ``'BMDD'`` or ``'mineralization law'``.
        :type start: str
        :raises ValueError: If the ``start`` parameter is not ``'BMDD'`` or ``'mineralization law'``.

        During initialization, the model sets up a 1D mesh for calcium content,
        initializes the BMDD cell variable, and defines the mineralization law and velocity
        according to the specified starting condition.
        The mesh is defined over the range of calcium content from ``minimum_content`` to
        ``maximum_content`` as specified in the :class:`Ruffoni_Parameters` class.
        """

        if not hasattr(self, "parameters"):
            self.parameters = Ruffoni_Parameters()
        self.simulation_time = simulation_time
        self.nx = number_of_grid_points
        self.mesh = fp.Grid1D(dx=(
                                         self.parameters.calcium.maximum_content - self.parameters.calcium.minimum_content) / number_of_grid_points,
                              nx=number_of_grid_points)
        self.mesh = self.mesh + self.parameters.calcium.minimum_content
        self.BMDD = fp.CellVariable(name="BMDD", mesh=self.mesh, hasOld=True)
        self.mineralization_law = None  # Will be initialized later
        self.inverse_mineralization_law = None  # Will be initialized later
        self.mineralization_velocity = None  # Will be initialized later
        self.start = start

    def initialize_model(self):
        """ Initialize the model based on the specified starting condition. If start is 'BMDD',
        it initializes a reference BMDD - in this case a Gaussian distribution - calculates the mineralization velocity,
        inverse mineralization law, and finally the mineralization law.
        If start is 'mineralization law', it initializes the mineralization law first, calculates the mineralization
        velocity, and finally the steady-state BMDD.

        :raises ValueError: If the start parameter is not 'BMDD' or 'mineralization law'."""
        if self.start == 'BMDD':
            self.initialize_steady_state_BMDD(self.start)
            self.initialize_mineralization_velocity(self.start)
            self.initialize_inverse_mineralization_law(self.start)
            self.initialize_mineralization_law(self.start)
        elif self.start == 'mineralization law':
            self.initialize_mineralization_law(self.start)
            self.initialize_mineralization_velocity(self.start)
            self.initialize_steady_state_BMDD(self.start)
        else:
            raise ValueError("Invalid start parameter. Use 'BMDD' or 'mineralization law'.")

    def solve_for_BMDD(self, save_interval=0.1):
        """
        Solve the time evolution of the Bone Mineralization Density Distribution (BMDD)
        using FiPy's Finite Volume Method (FVM).

        The model is solved iteratively in time, updating velocity, resorption, and formation
        at each step. A transient term models time evolution, a convection term models advection
        by the mineralization velocity, and an implicit source term models resorption.

        *Boundary conditions:*
        - Left boundary (low calcium): formation flux, given by the formed bone volume divided by
          the mineralization velocity at zero calcium.
        - Right boundary (high calcium): absorbing boundary (no-flux condition).

        The PDE is solved using FiPy’s ``sweep`` method with adaptive time-stepping controlled by
        a residual constraint. Results are stored at regular intervals.

        :param save_interval: Time interval (in years) between stored simulation states. Defaults to ``0.1``.
        :type save_interval: float, optional

        :return: (BMDD_evolution, BV_evolution, time_points) as BMDD, bone volume, and time points for each saved time step.
        :rtype: tuple of numpy.ndarray

        :raises ValueError: If the adaptive time step becomes too small during the simulation.

        Notes
        -----
        - The velocity field is defined as a **face variable** because it represents BMDD flux
          across cell boundaries.
        - The resorption coefficient is a **cell variable** since it is a local property.
        - The adaptive solver halves the time step until the desired residual is achieved,
          or terminates if the minimum time step is reached.
        """

        self.initialize_model()

        BMDD_evolution = [self.BMDD.value.copy()]
        BV_evolution = [self.calculate_bone_volume()]
        time_points = [0.0]

        velocity_field = fp.FaceVariable(mesh=self.mesh, rank=1)
        resorption_coeff = fp.CellVariable(mesh=self.mesh)
        equation = (
                fp.TransientTerm(var=self.BMDD)
                + fp.ConvectionTerm(coeff=velocity_field, var=self.BMDD)
                + fp.ImplicitSourceTerm(coeff=resorption_coeff, var=self.BMDD)
        )
        self.BMDD.faceGrad.constrain([0], self.mesh.facesRight)
        # Start time loop
        time = 0.0
        last_save_time = 0.0
        initial_time_step = 1.0e-2
        desired_residual = 1.0e-1  # Desired residual for convergence
        minimum_time_step = 1.0e-4  # Minimum time step to avoid too small dt
        vel_cells = [self.mineralization_velocity(c)
                     for c in self.mesh.cellCenters[0].value]
        velocity_field.setValue([[fp.CellVariable(mesh=self.mesh, value=vel_cells).faceValue]])
        while time < self.simulation_time:
            residual = 1.0  # Initial residual to enter the loop
            # Update values for this time step
            self.BMDD.updateOld()
            resorption_coeff.setValue(self.calculate_resorption_rate(time))
            BMDD_at_zero_calcium = (
                    self.calculate_formation_rate(time) /
                    self.mineralization_velocity(self.parameters.calcium.minimum_content)
            )
            self.BMDD.constrain(BMDD_at_zero_calcium, self.mesh.facesLeft)
            # Use sweep method to solve the PDE with constraint on residual
            dt = initial_time_step
            while residual > desired_residual:
                residual = equation.sweep(var=self.BMDD, dt=dt)
                dt = dt / 2
                if dt < minimum_time_step:
                    log.warning("Time step too small, stopping iteration.")
                    break
            time += dt
            # Save results
            if time - last_save_time >= save_interval:
                BMDD_evolution.append(self.BMDD.value.copy())
                BV_evolution.append(self.calculate_bone_volume())
                time_points.append(time)
                last_save_time = time
        return np.array(BMDD_evolution), np.array(BV_evolution), np.array(time_points)

    def check_mass_conservation(self, formation_rate, resorption_rate, bone_volume):
        """ This function checks the mass conservation in PDE solving - should be close to zero for equilibrium.

        :param formation_rate: Formation rate at current time step.
        :type formation_rate: float
        :param resorption_rate: Resorption rate at current time step.
        :type resorption_rate: float
        :param bone_volume: Current total bone volume.
        :type bone_volume: float
        :return: Change in bone volume (should be close to zero for mass conservation).
        :rtype: float"""

        total_formation = formation_rate
        print(f"Formed Bone Volume: {total_formation:.8f} [microm^3/time step]")
        total_resorption = resorption_rate * bone_volume
        print(f"Resorbed Bone Volume: {total_resorption:.8f} [microm^3/time step]")
        change_in_bone_volume = total_formation - total_resorption
        return change_in_bone_volume

    def initialize_steady_state_BMDD(self, start):
        """
        Initialize BMDD depending on the start parameter.

        If start is 'BMDD', it initializes from a reference BMDD distribution.

        If start is 'mineralization law', it initializes from a mineralization law.

        The BMDD is set to the initial value, and the formation rates are calculated based on the initial BMDD (formation depends on bone volume).

        :param start: Initialization method, either 'BMDD' or 'mineralization law'.
        :type start: str
        :raises ValueError: If start parameter is not 'BMDD' or 'mineralization law'.
        :return: None
        """
        if start == 'BMDD':
            initial_BMDD = self.initialize_BMDD_from_reference()
        elif start == 'mineralization law':
            initial_BMDD = self.initialize_BMDD_from_mineralization_law()
        else:
            raise ValueError("Invalid start parameter. Use 'BMDD' or 'mineralization law'.")
        self.BMDD.setValue(initial_BMDD)
        # For equilibrium: formation_rate = resorption_rate * bone_volume
        calcium_values = self.mesh.cellCenters[0].value
        bone_volume = np.trapz(initial_BMDD, calcium_values)
        self.initialize_formed_bone_volume(bone_volume, start)
        pass

    def initialize_formed_bone_volume(self, bone_volume, start):
        """ Initialize the formed bone volume based on the initial BMDD and initial bone volume. This is used to set the initial condition for the
        BMDD.

        :param bone_volume: The total bone volume calculated from the initial BMDD.
        :type bone_volume: float
        :param start: Initialization method, either 'BMDD' or 'mineralization law
        :type start: str
        :raises ValueError: If start parameter is not 'BMDD' or 'mineralization law'.
        :return: None
        """
        if start == 'BMDD':
            self.parameters.rate.initial_formation = self.parameters.rate.initial_resorption * bone_volume
            self.parameters.rate.final_formation = self.parameters.rate.final_resorption * bone_volume
        elif start == 'mineralization law':
            self.parameters.rate.initial_formation = self.parameters.rate.initial_formation * self.parameters.bone_volume.initial_value
            self.parameters.rate.final_formation = self.parameters.rate.final_formation * self.parameters.bone_volume.initial_value
        else:
            raise ValueError("Invalid start parameter. Use 'BMDD' or 'mineralization law'.")
        pass

    def initialize_BMDD_from_mineralization_law(self, plot=False):
        r"""
        Initialize the Bone Mineralization Density Distribution (BMDD) from a prescribed
        mineralization law. This method is called when the model is initialized with
        ``start='mineralization law'``.

        The initialization is based on the steady-state analytical formulation described in
        the source paper. The BMDD is computed by integrating the ratio of
        the resorption rate to the mineralization velocity across the calcium content range, taking the exponential
        of the integral and multiplying by a factor depending on formed volume and velocity.
        The resulting BMDD is then smoothed using a Gaussian filter and interpolated to match
        the mesh cell centers used in the numerical model.

        If ``plot=True``, the function additionally displays:
        - The exponential of the integrated inverse mineralization law (``exp(-∫ r/v)``),
        - The initial smoothed BMDD distribution over calcium content.

        :param plot: If ``True``, plots the exponential of the inverse mineralization integral
                     and the smoothed initial BMDD distribution. Defaults to ``False``.
        :type plot: bool
        :return: The initialized BMDD distribution mapped to the model mesh.
        :rtype: numpy.ndarray

        **Details:**
        - The BMDD is initialized as:
          :math:`\rho(C) = \frac{f(0)}{v_m(C)} \exp\left(-\int_{C_{min}}^{C} \frac{r(C')}{v_m(C')} \, dC'\right)`
          where:

          - :math:`f(0)` is the formed bone volume at steady state,
          - :math:`r(C)` is the resorption rate,
          - :math:`v_m(C)` is the mineralization velocity.

        - The integral is evaluated numerically using adaptive quadrature.
        - The smoothed BMDD is interpolated with a monotone cubic spline
          (``PchipInterpolator``) and clipped to enforce non-negativity.
        """
        calcium_values = np.linspace(self.parameters.calcium.minimum_content + 1e-6,
                                     self.parameters.calcium.maximum_content - 1e-6, self.nx * 20)
        formed_bone_volume = (self.calculate_formation_rate(t=0) * self.parameters.bone_volume.initial_value)
        resorption_rate = self.calculate_resorption_rate(t=0)

        vel = self.mineralization_velocity(calcium_values)
        vel = np.maximum(vel, 1e-9)

        # Vectorized integration using numpy and list comprehension
        def integrand(c_prime):
            return resorption_rate / np.maximum(self.mineralization_velocity(c_prime), 1e-9)

        # Use numpy vectorization for integration
        from functools import partial
        quad_func = np.vectorize(lambda c:
                                 sc.integrate.quad(partial(integrand), self.parameters.calcium.minimum_content, c,
                                                   epsabs=1e-10, epsrel=1e-8, limit=100)[0])
        integrand_values = quad_func(calcium_values)

        initial_BMDD = (formed_bone_volume / vel) * np.exp(-integrand_values)

        initial_BMDD_smoothed = gaussian_filter1d(initial_BMDD, sigma=5)
        bmdd_interp = PchipInterpolator(calcium_values, initial_BMDD_smoothed)
        mesh_calcium = self.mesh.cellCenters[0].value
        initial_BMDD_mesh = bmdd_interp(mesh_calcium)
        initial_BMDD_mesh = np.clip(initial_BMDD_mesh, 0, None)

        if plot:
            plt.figure(figsize=(10, 6))
            plt.plot(calcium_values, np.exp(-integrand_values), label="exp(-integral)")
            plt.xlabel("Calcium Content [%wt]")
            plt.ylabel("exp(-integral r/v)")
            plt.title("Exponential of the Integral of the Inverse Mineralization Law")
            plt.legend()
            plt.grid(True, alpha=0.3)

            plt.figure(figsize=(10, 6))
            plt.plot(mesh_calcium, initial_BMDD_mesh, label="Initial BMDD (smoothed)")
            plt.xlabel("Calcium Content [%wt]")
            plt.ylabel("BMDD [mm^3/wt%]")
            plt.title("Initial Bone Mineralization Density Distribution (Smoothed)")
            plt.legend()
            plt.grid(True, alpha=0.3)
            plt.show()

        return np.array(initial_BMDD_mesh)

    def initialize_BMDD_from_reference(self, plot=False):
        r"""
        Initialize the Bone Mineralization Density Distribution (BMDD) from a reference
        skewed Gaussian distribution. This method is called when the model is initialized
        with ``start='BMDD'``.

        The reference BMDD parameters are taken from Roschger et al. (2003), providing a
        physiological steady-state distribution representative of healthy trabecular bone.
        The distribution is constructed using a skew-normal probability density function (PDF)
        defined by its peak calcium content and spread, and is scaled to match the reference
        peak amplitude (×25 for unit consistency).

        If ``plot=True``, the method plots the initial BMDD distribution over calcium content.

        :param plot: If ``True``, plots the initial BMDD distribution. Defaults to ``False``.
        :type plot: bool
        :return: The initialized BMDD distribution defined over the mesh calcium content range.
        :rtype: numpy.ndarray

        **Details:**
        - The reference BMDD is computed as:
          :math:`\rho(C) = A \cdot \text{SkewNorm}(C; a, \mu, \sigma)`
          where:

          - :math:`a` is the skewness parameter (negative for left-skew, positive for right-skew),
          - :math:`\mu` is the peak calcium content (from Roschger reference data),
          - :math:`\sigma` is the mean standard deviation of lower and upper percentiles.

        - The resulting distribution is normalized such that its integral equals 1 and
          scaled by a constant factor (25) to match the reference BMDD peak.
        - The function optionally visualizes the initialized BMDD.
        """
        calcium_values = self.mesh.cellCenters[0].value
        shape = -3  # positive value for right-skew
        peak = self.parameters.reference_bmdd.peak
        scale = (self.parameters.reference_bmdd.standard_deviation_lower_percentile +
                 self.parameters.reference_bmdd.standard_deviation_upper_percentile) / 2
        reference_BMDD = skewnorm.pdf(calcium_values, a=shape, loc=peak, scale=scale)
        reference_BMDD = reference_BMDD / np.trapz(reference_BMDD, calcium_values)  # Normalize to ensure integral is 1
        reference_BMDD *= 25  # Scale to desired peak value

        if plot:
            plt.figure(figsize=(10, 6))
            plt.plot(calcium_values, reference_BMDD, label="Initial BMDD")
            plt.xlabel("Calcium Content [%wt]")
            plt.ylabel("BMDD [micro m^3/wt%]")
            plt.title("Initial Bone Mineralization Density Distribution")
            plt.legend()
            plt.grid(True, alpha=0.3)
            plt.show()
        return reference_BMDD

    def calculate_formation_rate(self, t):
        """ Return the formation rate at time t.

        :param t: Time in years.
        :type t: float
        :return: Formation rate at time t.
        :rtype: float
        """
        return self.parameters.rate.initial_formation

    def calculate_resorption_rate(self, t):
        """ Return resorption rate at time t.

        :param t: Time in years.
        :type t: float
        :return: Resorption rate at time t.
        :rtype: float """
        return self.parameters.rate.initial_resorption

    def calculate_bone_volume(self):
        """ Calculate current total bone volume by integrating BMDD.
        The bone volume is calculated as the sum of the product of BMDD and cell volumes across all cells.
        The BMDD is averaged over the mesh cell volumes to ensure correct scaling.

        :return: Total bone volume in micro m^3.
        :rtype: float"""
        return float(self.BMDD.cellVolumeAverage * self.mesh.cellVolumes.sum())

    def initialize_mineralization_velocity_from_BMDD(self):
        """Initialize mineralization velocity from BMDD steady state when start is 'BMDD'. The velocity is calculated
        by integrating the product of BMDD and resorption rate over the calcium content range and dividing the result by
        the corresponding BMDD value.

        The integration is performed using the PchipInterpolator for smooth interpolation of BMDD and resorption rates
        and numpy quad for numerical integration.

        :return: Mineralization velocity values as a numpy array and corresponding calcium values as a numpy array.
        :rtype: tuple of np.ndarray"""
        calcium_values = self.mesh.cellCenters[0].value
        initial_BMDD = self.BMDD.value.copy()
        bmdd_interp = PchipInterpolator(calcium_values, initial_BMDD)
        resorption_rates = np.array([self.calculate_resorption_rate(c) for c in calcium_values])
        resorption_interp = PchipInterpolator(calcium_values, resorption_rates)

        def integrand(c_prime):
            return bmdd_interp(c_prime) * resorption_interp(c_prime)

        quad_func = np.vectorize(lambda c, i:
                                 sc.integrate.quad(integrand, c, self.parameters.calcium.maximum_content,
                                                   epsabs=1e-12, epsrel=1e-9, limit=100)[0] / max(initial_BMDD[i],
                                                                                                  1e-30))
        mineralization_velocity_values = quad_func(calcium_values, np.arange(len(calcium_values)))
        return calcium_values, mineralization_velocity_values

    def initialize_mineralization_velocity_from_mineralization_law(self, start, plot=False):
        """Initialize mineralization velocity from the mineralization law when start is 'mineralization law'.
        The mineralization velocity is calculated as 1 / derivative of the inverse mineralization law.
        This function calls the `initialize_inverse_mineralization_law` method to ensure the inverse law is set up correctly.
        The calcium values are generated over a range defined by the minimum and maximum calcium content.
        The PChipInterpolator.derivative() method is used to calculate the derivative.
        Extrapolation handled by clamping the values to avoid issues with division by zero or negative derivatives.

        :param start: Initialization method, either 'BMDD' or 'mineralization law'.
        :type start: str
        :param plot: If True, plots the derivative of the inverse mineralization law.
        :type plot: bool
        :return: Mineralization velocity values as a numpy array and corresponding calcium values as a numpy array.
        :rtype: tuple of np.ndarray
        """
        calcium_values = np.linspace(
            self.parameters.calcium.minimum_content + 1e-6,
            self.parameters.calcium.maximum_content - 1e-6,
            self.nx * 5
        )
        self.initialize_inverse_mineralization_law(start)
        derivative_values = self.inverse_mineralization_law.derivative()(calcium_values)
        mineralization_velocity_values = 1 / derivative_values
        mineralization_velocity_values = np.clip(mineralization_velocity_values, 0, None)

        if plot:
            plt.figure(figsize=(10, 6))
            plt.plot(calcium_values, derivative_values,
                     label="Derivative of the inverse mineralization law", color='orange')
            plt.xlabel("Calcium Content [%wt]")
            plt.ylabel("Derivative of the inverse mineralization law")
            plt.title("Derivative of the inverse mineralization law")
            plt.legend()
            plt.grid(True, alpha=0.3)
            plt.show(block=False)
        return calcium_values, mineralization_velocity_values

    def initialize_mineralization_velocity(self, start, plot=False):
        """
        Initialize mineralization velocity based on the initial BMDD or mineralization law by calling the respective functions.
        This method sets up the mineralization velocity law as a PchipInterpolator based on the calcium content and
        corresponding mineralization velocity values.

        :param start: Initialization method, either 'BMDD' or 'mineralization law'.
        :type start: str
        :param plot: If True, plots the mineralization velocity law.
        :type plot: bool
        :raises ValueError: If start parameter is not 'BMDD' or 'mineralization law'.
        :return: None
        """
        if start == 'BMDD':
            calcium_values, mineralization_velocity_values = self.initialize_mineralization_velocity_from_BMDD()
        elif start == 'mineralization law':
            calcium_values, mineralization_velocity_values = self.initialize_mineralization_velocity_from_mineralization_law(
                start, plot=False)
        else:
            raise ValueError("Invalid start parameter. Use 'BMDD' or 'mineralization law'.")
        self.mineralization_velocity = PchipInterpolator(calcium_values, mineralization_velocity_values)

        if plot:
            plt.figure(figsize=(10, 6))
            plt.plot(calcium_values, self.mineralization_velocity(calcium_values),
                     label="Mineralization Velocity Law", color='orange')
            plt.xlabel("Calcium Content [%wt]")
            plt.ylabel("Mineralization Velocity")
            plt.title("Mineralization Velocity Law Over Calcium Content")
            plt.legend()
            plt.grid(True, alpha=0.3)
            plt.show(block=False)
        pass

    def initialize_inverse_mineralization_law_from_BMDD(self):
        """ Initialize the inverse mineralization law from BMDD steady state by integrating 1 / mineralization velocity
        from 0 to the calcium content value using numerical integration.

        :return: Tuple of numpy arrays containing calcium values and corresponding inverse mineralization law values.
        :rtype: tuple of np.ndarray
        :raises ValueError: If the initial bone volume is less than or equal to zero."""
        # Initialize inverse mineralization law from mineralization velocity
        calcium_values = self.mesh.cellCenters[0].value
        initial_bone_volume = self.calculate_bone_volume()

        if initial_bone_volume <= 0:
            raise ValueError("Initial bone volume must be greater than zero.")

        def integrand(c_prime):
            return 1 / self.mineralization_velocity(c_prime)

        inverse_values = []
        for c in calcium_values:
            [integral_result, _] = sc.integrate.quad(integrand, 0, c)
            inverse_values.append(integral_result)
        return calcium_values, np.array(inverse_values)

    def initialize_inverse_mineralization_law_from_mineralization_law(self):
        """ Initialize the inverse mineralization law from the mineralization law by generating a range of time values
        and calculating the corresponding calcium values using the mineralization law. The method ensures that the
        calcium values are unique and monotonic, and clamps them to the valid range defined by the minimum and maximum
        calcium content.

        :return: Tuple of numpy arrays containing calcium values and corresponding inverse mineralization law values.
        :rtype: tuple of np.ndarray"""
        time_values = np.linspace(0, 10000, 100000)
        calcium_values = self.mineralization_law(time_values)
        calcium_values, indices = np.unique(calcium_values, return_index=True)
        inverse_values = time_values[indices]
        # Clamp to valid range
        mask = (calcium_values >= self.parameters.calcium.minimum_content) & \
               (calcium_values <= self.parameters.calcium.maximum_content)
        calcium_values = calcium_values[mask]
        inverse_values = inverse_values[mask]
        return calcium_values, inverse_values

    def initialize_inverse_mineralization_law(self, start, plot=False):
        """Initialize the inverse mineralization law based on either BMDD or mineralization law.
        This method sets up the inverse mineralization law as a PchipInterpolator based on the calcium content and
        corresponding inverse mineralization law values.

        :param start: Initialization method, either 'BMDD' or 'mineralization law'.
        :type start: str
        :param plot: If True, plots the inverse mineralization law.
        :type plot: bool
        :raises ValueError: If start parameter is not 'BMDD' or 'mineralization law'.
        :return: None"""
        if start == 'BMDD':
            calcium_values, inverse_values = self.initialize_inverse_mineralization_law_from_BMDD()
        elif start == 'mineralization law':
            calcium_values, inverse_values = self.initialize_inverse_mineralization_law_from_mineralization_law()
        else:
            raise ValueError("Invalid start parameter. Use 'BMDD' or 'mineralization law'.")
        self.inverse_mineralization_law = PchipInterpolator(calcium_values, inverse_values, extrapolate=True)

        if plot:
            plt.figure(figsize=(10, 6))
            plt.plot(calcium_values, self.inverse_mineralization_law(calcium_values),
                     label="Inverse Mineralization Law",
                     color='orange')
            plt.xlabel("Calcium Content [%wt]")
            plt.ylabel("Time since bone formation [years]")
            plt.title("Inverse Mineralization Law Over Calcium Content")
            plt.legend()
            plt.grid(True, alpha=0.3)
            plt.show(block=False)

    def initialize_mineralization_law(self, start, plot=False):
        """ Initialize the mineralization law based on either BMDD or given mineralization law. For the 'BMDD' start,
        it inverts the values of the inverse mineralization law. For 'mineralization law' start, a double hyperbolic or double exponential
        mineralization law is prescribed. Duplicates are removed, and the values are ensured to be monotonic.
        The values are interpolated using PchipInterpolator for smoothness.

        :param start: Initialization method, either 'BMDD' or 'mineralization law'.
        :type start: str
        :param plot: If True, plots the mineralization law. Defaults to False.
        :type plot: bool
        :raises ValueError: If start parameter is not 'BMDD' or 'mineralization law'.
        :return: None
        """
        if start == 'BMDD':
            # Initialize mineralization law from BMDD steady state
            calcium_values = self.mesh.cellCenters[0].value
            time_values = self.inverse_mineralization_law(calcium_values)
        elif start == 'mineralization law':
            time_values = np.linspace(0, 10000, 1000000)
            law_choice = 'double hyperbolic'  # or 'double exponential'
            if law_choice == 'double hyperbolic':
                calcium_values = self.initialize_double_hyperbolic_mineralization_law(time_values)
            elif law_choice == 'double exponential':
                calcium_values = self.initialize_double_exponential_mineralization_law(time_values)
            # Remove duplicates and ensure monotonicity
            calcium_values, indices = np.unique(calcium_values, return_index=True)
            time_values = time_values[indices]
        else:
            raise ValueError("Invalid start parameter. Use 'BMDD' or 'mineralization law'.")

        self.mineralization_law = PchipInterpolator(time_values, calcium_values)

        if plot:
            plt.figure(figsize=(10, 6))
            plt.plot(time_values, self.mineralization_law(time_values), label="Mineralization Law", color='teal')
            plt.xlabel("Time since bone formation [years]")
            plt.ylabel("Calcium content [%wt]")
            plt.title("Mineralization Law Over Time")
            plt.legend()
            plt.grid(True, alpha=0.3)
            plt.show()

    def initialize_double_hyperbolic_mineralization_law(self, time_values):
        r"""
        Initialize a double-hyperbolic mineralization law, which maps tissue age
        (time since formation) to calcium content. This function evaluates the
        mineralization law at the provided ``time_values``.

        The mineralization law is defined as:

        .. math::
            C(t) =
            C_1 \frac{t / r_1}{1 + (t / r_1)} +
            C_2 \frac{t / r_2}{1 + (t / r_2)}

        where:

        - :math:`C_1, C_2` are the asymptotic calcium contents of the fast and slow
          mineralization phases,
        - :math:`r_1, r_2` are characteristic time constants of the two phases,
        - :math:`t` is the time/ tissue age.

        :param time_values: Time array at which the mineralization
            law is evaluated.
        :type time_values: np.ndarray
        :return: Calcium content values corresponding to the provided time points.
        :rtype: np.ndarray
        """
        calcium_values = (self.parameters.mineralization_law.primary_mineral_content *
                          (time_values / self.parameters.mineralization_law.primary_apposition_rate) /
                          (1 + (time_values / self.parameters.mineralization_law.primary_apposition_rate))) + (
                                 self.parameters.mineralization_law.maximum_mineral_content *
                                 (time_values / self.parameters.mineralization_law.secondary_apposition_rate) /
                                 (1 + (
                                         time_values / self.parameters.mineralization_law.secondary_apposition_rate)))
        return calcium_values

    def initialize_double_exponential_mineralization_law(self, time_values):
        r"""
        Initialize a double-exponential mineralization law, which maps tissue age
        (time since formation) to calcium content. This function evaluates the
        mineralization law at the provided ``time_values``.

        The mineralization law is defined as:

        .. math::
            C(t) =
            C_1 (1 - \exp(\frac{t}{r_1}))+
            (C_2 - C_1) (1 - \exp(\frac{t}{r_2}))

        where:

        - :math:`C_1, C_2` are the asymptotic calcium contents of the fast and slow
          mineralization phases,
        - :math:`r_1, r_2` are characteristic time constants of the two phases,
        - :math:`t` is the time/ tissue age.

        :param time_values: Time array at which the mineralization
            law is evaluated.
        :type time_values: np.ndarray
        :return: Calcium content values corresponding to the provided time points.
        :rtype: np.ndarray
        """
        calcium_values = ((self.parameters.mineralization_law.primary_mineral_content *
                   (1 - np.exp(-time_values /
                               self.parameters.mineralization_law.primary_apposition_rate)))
                  + ((self.parameters.mineralization_law.maximum_mineral_content -
                      self.parameters.mineralization_law.primary_mineral_content) *
                     (1 - np.exp(
                         - time_values /
                         self.parameters.mineralization_law.secondary_apposition_rate))))
        return calcium_values

    def plot_results(self, BMDD_evolution, BV_evolution, time_points):
        """
        Plot normalized BMDD and bone volume evolution at saved / specified time points.
        It also plots the mineralization law, mineralization velocity, and the formed and resorbed bone volume over time.

        :param BMDD_evolution: Evolution of BMDD over time as a numpy array.
        :type BMDD_evolution: np.ndarray
        :param BV_evolution: Evolution of bone volume over time as a numpy array.
        :type BV_evolution: np.ndarray
        :param time_points: Time points at which the BMDD and bone volume were saved.
        :type time_points: np.ndarray
        :return: None
        """
        calcium_values = self.mesh.cellCenters[0].value

        # --- Plot BMDD evolution ---
        plt.figure(figsize=(10, 6))
        for i, t in enumerate(time_points):
            plt.plot(calcium_values, BMDD_evolution[i] / BV_evolution[i],
                     label=f't = {t:.2f} years')
        plt.xlabel('Calcium Content [wt%]')
        plt.ylabel('BMDD [1/wt%]')
        plt.title('Bone Mineralization Density Distribution Evolution')
        plt.legend(loc='upper left', fontsize='small', ncol=2)
        plt.grid(True, alpha=0.3)

        # --- Plot BV evolution ---
        plt.figure(figsize=(10, 6))
        plt.plot(time_points, BV_evolution)
        plt.ylabel('BV [micro m^3]')
        plt.xlabel('Time')
        plt.title('Bone Volume Evolution')
        plt.grid(True, alpha=0.3)

        # --- Plot Mineralization Law ---
        plt.figure(figsize=(10, 6))
        mineralization_law = [self.mineralization_law(time) for time in time_points]
        plt.plot(time_points, mineralization_law, color='teal')
        plt.xlabel("Time [years]")
        plt.ylabel("Calcium Content [%wt]")
        plt.title("Mineralization Law Over Time")
        plt.grid(True, alpha=0.3)

        # --- Plot Mineralization Velocity ---
        plt.figure(figsize=(10, 6))
        mineralization_velocity = [self.mineralization_velocity(c) for c in calcium_values]
        plt.plot(calcium_values, mineralization_velocity, color='blue')
        plt.xlabel("Calcium Content [%wt]")
        plt.ylabel("Velocity [dc/dt, %wt/year]")
        plt.title("Mineralization Velocity as a Function of Calcium Content")
        plt.grid(True, alpha=0.3)

        plt.show()
