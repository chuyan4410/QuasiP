# Rotating-frame force coefficients for PDIM

The PDIM-only rotating-frame interface assumes a fixed unit axis `axisR`
and the parameter ordering

```text
mu = [omega; dotOmega] = mu0 + delta.
```

It does not replace or reconstruct the complete `Frot` used by the time
integrator. Aerodynamic-force parameter dependence and gravity are outside
the scope of this interface.

## Homogeneous coefficient evaluator

```matlab
Fpart = myAssembly.genSysFrotPdimPart( ...
    z, axisR, mu0, stateOrder, parameterMultiIndex);
```

`parameterMultiIndex = [ellOmega, ellDotOmega]` selects the coefficient of

```text
deltaOmega^ellOmega * deltaDotOmega^ellDotOmega.
```

The nonzero combinations are:

| State order | Parameter multi-index |
| ---: | :--- |
| 1 | `[0,0]`, `[1,0]`, `[0,1]`, `[2,0]` |
| 3 | `[0,0]`, `[1,0]`, `[0,1]`, `[2,0]` |

The `[1,1]` and `[0,2]` coefficients are identically zero in the current
structural model.

For example, the state-cubic coefficient multiplying `deltaOmega` is

```matlab
F310 = myAssembly.genSysFrotPdimPart(z,axisR,mu0,3,[1,0]);
```

## Multilinear state input

For a cell `inputx` containing state directions, use

```matlab
val = get_frot_pdim_multilinear_form( ...
    myAssembly,axisR,mu0,parameterMultiIndex,inputx);
```

The wrapper calls `get_multilinear_form` for state order 3. The
state-linear terms are evaluated directly.

The rotating force uses the direct time derivative of the interpolated
quaternion rotation matrix, without normalizing the Gauss-point quaternion.
Only its state-linear and state-cubic terms are registered. Structural
elastic fifth-order terms and higher-order ROM expansions are unchanged.
Eigenvalues, eigenvectors and ROMs generated using the old rotating force
must be recomputed before reuse with the updated full-order model.
