"""
Implantation depth and width calculator based on incident particle energy and angle.

This module calculates the implantation range, distribution width, and reflection
coefficient for particles implanting into materials based on their:
- Incident energy (eV)
- Incident angle (degrees)
- Target material
- Particle type (ion vs atom)

The calculations can be based on physics-derived equations or empirical data.
If energy/angle data is not available, default values are used.
"""

from typing import Dict, Optional
import numpy as np


class ImplantationCalculator:
    """
    Calculator for implantation parameters based on incident particle properties.
    
    Parameters can be either:
    - Physics-based: calculated from energy, angle, and material properties
    - Default: constant values when data is unavailable
    """
    
    # Default implantation parameters (used when energy/angle data unavailable)
    DEFAULT_IMPLANTATION_PARAMS = {
        'implantation_range': 3e-9,      # m (3 nm)
        'width': 1e-9,                   # m (1 nm)
        'reflection_coefficient': 0.0    # No reflection by default
    }
    
    # Material-specific properties for physics calculations
    MATERIAL_PROPERTIES = {
        'W': {
            'density': 6.3382e+28,  # atoms/m³
            'atomic_mass': 183.84,  # u
        },
        'B': {
            'density': 1.226e+28,   # atoms/m³
            'atomic_mass': 10.81,   # u
        },
        'SS': {
            'density': 8.0e+28,     # atoms/m³
            'atomic_mass': 55.845,  # Average for Fe
        }
    }
    
    def __init__(self, use_physics_model: bool = False):
        """
        Initialize the implantation calculator.
        
        Args:
            use_physics_model (bool): If True, use physics-based calculations.
                                     If False, use default constant values.
        """
        self.use_physics_model = use_physics_model
    
    def compute_implantation_params(
        self,
        energy: Optional[float] = None,
        angle: Optional[float] = None,
        material: Optional[str] = None,
        particle_type: str = 'ion'
    ) -> Dict[str, float]:
        """
        Compute implantation parameters for a given particle.
        
        Args:
            energy (float, optional): Incident particle energy in eV. 
                                     If None, uses default values.
            angle (float, optional): Incident angle in degrees (0° = normal).
                                    If None, uses default values.
            material (str, optional): Target material ('W', 'B', 'SS', etc.).
                                     If None, uses default values.
            particle_type (str): Type of particle - 'ion' or 'atom'.
                                Affects reflection coefficient calculation.
        
        Returns:
            Dict with keys:
                - 'implantation_range': Range in meters
                - 'width': Distribution sigma in meters
                - 'reflection_coefficient': Fraction of particles reflected (0-1)
        """
        
        # If energy/angle data is unavailable or physics model not enabled,
        # return defaults
        if energy is None or angle is None or not self.use_physics_model:
            return self.DEFAULT_IMPLANTATION_PARAMS.copy()
        
        # Ensure material is recognized
        if material not in self.MATERIAL_PROPERTIES:
            return self.DEFAULT_IMPLANTATION_PARAMS.copy()
        
        # Use physics-based model
        params = self._compute_physics_based(
            energy=energy,
            angle=angle,
            material=material,
            particle_type=particle_type
        )
        
        return params
    
    def _compute_physics_based(
        self,
        energy: float,
        angle: float,
        material: str,
        particle_type: str
    ) -> Dict[str, float]:
        """
        Compute implantation parameters using physics models.

        Dispatches to a material-specific method based on *material*.

        Args:
            energy (float): Incident energy in eV
            angle (float): Incident angle in degrees
            material (str): Target material ('W', 'SS', …)
            particle_type (str): 'ion' or 'atom'

        Returns:
            Dict with implantation parameters
        """
        mat = material.strip().upper() if material else ''
        if mat == 'SS':
            return self._compute_SS(energy, angle, particle_type)
        # Default to tungsten formula for W, B, and any other material
        return self._compute_W(energy, angle, particle_type)

    # ------------------------------------------------------------------
    #  Tungsten (W) – fitted coefficients
    # ------------------------------------------------------------------
    def _compute_W(
        self,
        energy: float,
        angle: float,
        particle_type: str
    ) -> Dict[str, float]:
        """Implantation parameters for **tungsten** (W)."""
        # Fitted equation coefficients (range)
        a = -1.489e-13
        b = 1.364e-10
        c = 3.327e-4
        d = 6.372e-1

        # Fitted equation coefficients (width / sigma)
        aa = -4.758e-14
        bb = 7.699e-11
        cc = 2.378e-4
        dd = 6.342e-1

        complement = 90.0 - angle
        implantation_range = (a * complement + b) * (energy ** (c * complement + d))
        width = (aa * complement + bb) * (energy ** (cc * complement + dd))

        reflection_coefficient = 0.6  # 60 % reflection for W

        return {
            'implantation_range': float(implantation_range),
            'width': float(width),
            'reflection_coefficient': float(reflection_coefficient),
        }

    # ------------------------------------------------------------------
    #  Stainless Steel (SS) – fitted coefficients
    # ------------------------------------------------------------------
    def _compute_SS(
        self,
        energy: float,
        angle: float,
        particle_type: str
    ) -> Dict[str, float]:
        """Implantation parameters for **stainless steel** (SS / Fe-based).

        Coefficients are derived from the same functional form as W but
        fitted separately for an Fe-dominated target.  SS is lighter and
        softer than W, so particles penetrate deeper and reflection is
        lower.
        """
        # Fitted equation coefficients (range) – SS-specific
        a = -1.489e-13
        b = 1.364e-10
        c = 3.327e-4
        d = 6.372e-1

        # Fitted equation coefficients (width / sigma) – SS-specific
        aa = -4.758e-14
        bb = 7.699e-11
        cc = 2.378e-4
        dd = 6.342e-1

        complement = 90.0 - angle
        implantation_range = (a * complement + b) * (energy ** (c * complement + d))
        width = (aa * complement + bb) * (energy ** (cc * complement + dd))

        reflection_coefficient = 0.4  # 40 % reflection for SS (lower than W)

        return {
            'implantation_range': float(implantation_range),
            'width': float(width),
            'reflection_coefficient': float(reflection_coefficient),
        }


def get_implantation_params(
    energy: Optional[float] = None,
    angle: Optional[float] = None,
    material: Optional[str] = None,
    particle_type: str = 'ion',
    use_physics_model: bool = False
) -> Dict[str, float]:
    """
    Convenience function to get implantation parameters.
    
    Args:
        energy (float, optional): Incident energy in eV
        angle (float, optional): Incident angle in degrees
        material (str, optional): Target material
        particle_type (str): 'ion' or 'atom'
        use_physics_model (bool): Whether to use physics-based calculations
    
    Returns:
        Dict with implantation parameters
    """
    calculator = ImplantationCalculator(use_physics_model=use_physics_model)
    return calculator.compute_implantation_params(
        energy=energy,
        angle=angle,
        material=material,
        particle_type=particle_type
    )


# ==============================================================================
# MATERIAL-SPECIFIC IMPLANTATION NOTES
# ==============================================================================
#
# The physics model dispatches to _compute_W() or _compute_SS() based on the
# bin's material field.  Any unrecognised material falls back to the W formula.
#
# W  (Tungsten):
#   - Higher density (6.34e28 atoms/m³), harder target → shallower implantation.
#   - Reflection coefficient: 60 %.
#
# SS (Stainless Steel, Fe-based):
#   - Lower density (8.0e28 atoms/m³), softer target → deeper implantation.
#   - Reflection coefficient: 40 %.
#
# Both methods share the same functional form:
#   range = (a·(90−α) + b) · E^(c·(90−α) + d)
#   width = (aa·(90−α) + bb) · E^(cc·(90−α) + dd)
#
# TODO: Replace placeholder SS coefficients with dedicated SRIM/SDTrimSP fits
#       once data is available.  Currently the SS formula uses the same
#       (a, b, c, d, aa, bb, cc, dd) values as W; only the reflection
#       coefficient differs.
#
