# effective_property_codex

Python template for driving Ansys MAPDL / Material Designer to extract
homogenized elastic constants and CTE for a simple fiber-in-matrix
representative volume element.

## Usage
1. Install the `ansys-mapdl-core` package in your Python environment and
   ensure MAPDL is licensed and reachable.
2. Run the script:
   ```bash
   python effective_property.py
   ```
3. The script prints effective elastic moduli, shear moduli, Poisson's
   ratios, and coefficients of thermal expansion for a 1×1×1 mm cube
   containing a centrally located cylindrical fiber sized for a 10%
   volume fraction.

Edit `effective_property.py` to adjust constituent material properties,
mesh size, or boundary conditions to suit your specific study.
