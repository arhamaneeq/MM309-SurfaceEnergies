> [!IMPORTANT]
> 
> Throughout this project, the composition will always be represented as
> 
> $$
> (c_\text{Fe}, c_\text{Ni}, c_\text{Al})
> $$
> 
> I.e., this ordering will be enforced in every generated `.data` file used as an input for a LAMMPS script.
> 
> ```data
> 3 atom types
> 
> Masses
> 
> 1 55.845     # Fe
> 2 58.6934    # Ni
> 3 26.9815385 # Al
> ```
> 
> And the `Atoms # atomic` section *must* use type IDs according to this mapping. Thus,
> 
> ```lammps
> pair_style eam/alloy
> pair_coeff * * potential.eam.alloy Fe Ni Al
> ```
> 
> The element names must appear in precisely *this* order.

