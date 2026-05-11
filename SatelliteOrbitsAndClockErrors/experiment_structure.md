**broadcast_orbits.py**
**_Ephemeris Propagation Parses RINEX navigation files._**
Uses Keplerian propagation (solving Kepler's equation E - e \sin E = M) for GPS/Galileo/BeiDou, and high-precision numerical integration (scipy.integrate.solve_ivp) for GLONASS ECEF differential equations.

**precise_orbits.py**
**_SP3 Loading & Interpolation Parses official IGS combination solutions (SP3c format)._**
Interpolates the sparse 5-minute or 15-minute data onto a dense time grid using piecewise cubic spline interpolation, and renders 3D trajectories.

**compare_broadcast_precise.py**
**_Comparative Validation Aligns both generated trajectories onto an identical time vector._**
Computes the 3D Euclidean error vector and decomposes residuals into local Radial, Along-track, and Cross-track (RAC) components to quantify broadcast accuracy degradation.

**skyplot_precise.py**
**_User Visibility Mapping Transforms ECEF coordinates to local Topocentric (East-North-Up) and Azimuth-Elevation frames relative to a reference station (e.g., Padova)._**
Applies an elevation mask to render polar skyplots and visible satellite counts over time.
