# Möbius Strip Modeling

This project implements a Python script that models a Möbius strip using parametric equations and computes key geometric properties.

## Overview

A Möbius strip is a surface with only one side and one edge. This implementation:

1. Creates a 3D mesh/grid of points on the Möbius strip surface
2. Computes the surface area using numerical integration
3. Calculates the edge length
4. Provides visualization of the Möbius strip

## Requirements

- Python 3.x
- NumPy
- Matplotlib
- SciPy

## Installation

```bash
pip install numpy matplotlib scipy
```

## Usage

Run the script directly:

```bash
python mobius_strip.py
```

Or import the `MobiusStrip` class in your own code:

```python
from mobius_strip import MobiusStrip

# Create a Möbius strip with R=3, w=1, n=100
mobius = MobiusStrip(R=3, w=1, n=100)

# Compute geometric properties
surface_area = mobius.compute_surface_area()
edge_length = mobius.compute_edge_length()

# Visualize
mobius.visualize()
```

## Parameters

- `R`: Radius (distance from the center to the strip)
- `w`: Width of the strip
- `n`: Resolution (number of points in the mesh)

## Mathematical Background

The Möbius strip is parameterized by the equations:
- x(u,v) = (R + v·cos(u/2))·cos(u)
- y(u,v) = (R + v·cos(u/2))·sin(u)
- z(u,v) = v·sin(u/2)

Where:
- u ∈ [0,2π]
- v ∈ [-w/2, w/2]
