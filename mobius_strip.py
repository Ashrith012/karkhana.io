"""
Mobius Strip Visualization and Analysis

This script models a Mobius strip using parametric equations, calculates 
its surface area and edge length, and renders a 3D visualization.

Author: Cascade AI Assistant
Date: May 27, 2025
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D
from scipy import integrate

class MobiusStrip:
    """
    Represents a Mobius strip and provides methods to compute
    its geometric properties and visualize it.

    Parameters:
    ----------
    R : float
        Radius from center to midline of the strip.
    w : float
        Width of the strip.
    n : int
        Resolution of the mesh grid.
    """

    def __init__(self, R, w, n):
        self.R = R
        self.w = w
        self.n = n

        # Generate a mesh grid for parameters u and v
        self.u_values = np.linspace(0, 2 * np.pi, n)
        self.v_values = np.linspace(-w / 2, w / 2, n)
        self.u_mesh, self.v_mesh = np.meshgrid(self.u_values, self.v_values)

        # Compute 3D coordinates of the Mobius surface
        self._compute_surface_points()

    def _compute_surface_points(self):
        """Compute 3D coordinates of the Mobius strip."""
        u, v = self.u_mesh, self.v_mesh
        self.x = (self.R + v * np.cos(u / 2)) * np.cos(u)
        self.y = (self.R + v * np.cos(u / 2)) * np.sin(u)
        self.z = v * np.sin(u / 2)

    def parametric_point(self, u, v):
        """Return a 3D point (x, y, z) on the surface for given (u, v)."""
        x = (self.R + v * np.cos(u / 2)) * np.cos(u)
        y = (self.R + v * np.cos(u / 2)) * np.sin(u)
        z = v * np.sin(u / 2)
        return np.array([x, y, z])

    def _partial_derivatives(self, u, v):
        """
        Compute the partial derivatives of the parametric surface
        with respect to u and v at a point (u, v).
        """
        # ∂r/∂u
        dx_du = -0.5 * v * np.sin(u / 2) * np.cos(u) - (self.R + v * np.cos(u / 2)) * np.sin(u)
        dy_du = -0.5 * v * np.sin(u / 2) * np.sin(u) + (self.R + v * np.cos(u / 2)) * np.cos(u)
        dz_du = 0.5 * v * np.cos(u / 2)
        du = np.array([dx_du, dy_du, dz_du])

        # ∂r/∂v
        dx_dv = np.cos(u / 2) * np.cos(u)
        dy_dv = np.cos(u / 2) * np.sin(u)
        dz_dv = np.sin(u / 2)
        dv = np.array([dx_dv, dy_dv, dz_dv])

        return du, dv

    def _surface_area_element(self, u, v):
        """
        Compute the magnitude of the cross product of partial derivatives,
        which gives the local surface area element.
        """
        du, dv = self._partial_derivatives(u, v)
        return np.linalg.norm(np.cross(du, dv))

    def compute_surface_area(self):
        """
        Approximate the surface area using numerical double integration.

        Returns:
        --------
        float
            Estimated surface area of the Mobius strip.
        """
        def integrand(v, u):
            return self._surface_area_element(u, v)

        area, _ = integrate.dblquad(
            integrand,
            0, 2 * np.pi,  # u limits
            lambda u: -self.w / 2,
            lambda u: self.w / 2  # v limits
        )
        return area

    def compute_edge_length(self):
        """
        Estimate the total edge length of the Mobius strip.

        Since the Mobius strip has one boundary that loops twice around,
        we evaluate the curve at v = w/2 from u = 0 to 4π.
        """
        u_vals = np.linspace(0, 4 * np.pi, 1000)
        edge_points = np.array([self.parametric_point(u, self.w / 2) for u in u_vals])
        distances = np.linalg.norm(edge_points[1:] - edge_points[:-1], axis=1)
        return np.sum(distances)

    def visualize(self, show=True, save_path=None):
        """
        Plot the Mobius strip in 3D.

        Parameters:
        ----------
        show : bool
            Whether to display the plot window.
        save_path : str or None
            Path to save the plot as an image file.
        """
        fig = plt.figure(figsize=(10, 8))
        ax = fig.add_subplot(111, projection='3d')

        # Plot surface
        surf = ax.plot_surface(
            self.x, self.y, self.z,
            cmap=cm.viridis, edgecolor='none', alpha=0.85
        )

        # Plot edge line
        u_edge = np.linspace(0, 4 * np.pi, 1000)
        edge = np.array([self.parametric_point(u, self.w / 2) for u in u_edge])
        ax.plot(edge[:, 0], edge[:, 1], edge[:, 2], color='red', linewidth=2)

        # Axis labels
        ax.set_xlabel('X')
        ax.set_ylabel('Y')
        ax.set_zlabel('Z')
        ax.set_title(f'Mobius Strip (R={self.R}, w={self.w})')
        ax.set_box_aspect([1, 1, 1])

        fig.colorbar(surf, ax=ax, shrink=0.5, aspect=5)

        if save_path:
            plt.savefig(save_path, dpi=300, bbox_inches='tight')

        if show:
            plt.show()

        return fig

def main():
    """Run a demo of the MobiusStrip class."""
    R = 3
    w = 1
    resolution = 100

    strip = MobiusStrip(R, w, resolution)

    area = strip.compute_surface_area()
    edge_len = strip.compute_edge_length()

    print("Mobius Strip Properties:")
    print(f"  Radius       : {R}")
    print(f"  Width        : {w}")
    print(f"  Mesh Points  : {resolution}")
    print(f"  Surface Area : {area:.4f} square units")
    print(f"  Edge Length  : {edge_len:.4f} units")

    # Display and save the plot
    strip.visualize(save_path="mobius_strip_visualization.png")

if __name__ == "__main__":
    main()
