# RichardsDAE

This repository contains the code associated with the short communication:
"Why modified Picard works: mass-conservative and higher order
time integration for Richardson-Richards equation", submitted for
consideration in Advances in Water Resources.

# Reproducing Results

The `cases/` directory contains Julia scripts that reproduce all figures and tables
from the paper. All package dependencies and versions are specified 
in `Project.toml` and `Manifest.toml`, ensuring reproducible results.

To reproduce all results, run:

```console
julia --project cases/haverkamp.jl
julia --project cases/newmexico.jl

julia --project cases/plothead.jl
```

Output figures and data will be saved to `cases/output/`.

## License

This project is licensed under the [MIT License](LICENSE), see the LICENSE 
file for details. See `data/README` for details on the data.

## Citation

If you use this code in your research, please cite the software:

Bootsma, H. P. (2026). Code and analysis for "Why modified Picard works:
mass-conservative and higher order time integration for Richardson-Richards
equation", (Version 1.0) [Software]. Zenodo. 

The associated manuscript has been submitted for consideration at Advances in
Water Resources.

