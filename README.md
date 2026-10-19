# MM309-SurfaceEnergies

Exploration of FCC(Ni, Al, Fe) ternary surface energies using LAMMPS

## Assignment

### Objective

Investigate how alloy chemistry, crystallographic orientation, and temperature influence the surface energy of a chemically disordered ternary FCC alloy using atomistic simulations. *The objective is not simply to obtain numerical values, but to establish that the calculated values are physically meaningful and reproducible.*

### Computational Requirements

The minimum dataset is

$$
3 \text{ Alloys} × 3 \text{ Surfaces} × 3 \text{ Temperatures} × 3 \text{ Independent Configurations} = 81 \text{ Surface calculations per group}.
$$

Required temperatures: 
$$
300 \text{ K}, 500 \text{ K}, 700 \text{ K}
$$

The three configurations at each temperature must represent chemically disordered alloy configurations. Describe how these configurations were generated and explain whether the sampling is sufficient. In addition to the 81 primary calculations, perform the calculations needed for bulk reference states and model convergence tests.

### Interatomic Potential

Select an appropriate interatomic potential for the assigned alloy. Briefly document:
- potential formalism and source
- elements and interactions represented
- properties or systems used for validation
- relevance to surface energy calculations
- important limitations

The potential choice must be supported by published literature and/or other reliable technical documentation.

### Model and Convergence

- Determine suitable dimensions for the surface slab and the seperatation between periodic surface images.
- Perform a systematic convergence study for at least one representative surface. Use the results to justify the model used for the remaining calculations.
- The bulk reference energy required for the surface-energy calculation must also be clearly defined and calculated consistently.

### Analysis, Required Results, and Submission

Calculate $\gamma_{111}$, $\gamma_{110}$, and $\gamma_{100}$ for all three assigned alloys at 300K, 500K, and 700K. Report the configuration-to-configuration variation and clearly define the uncertainty reported.

The analysis should address:
- The effect of chemical disorder and alloy chemistry
- The temperature dependence of surface energy
- The sensitivity to model size and configuration sampling

The submitted report must include:
- Potential selection and justification
- Computational methodology and configuration generation
- Bulk reference-state calculation
- Convergence study
- Surface energies and statistical analysis
- Comparison of the three alloy systems
- Literature comparison
- Limitations and uncertainty
- Conclusions and references

