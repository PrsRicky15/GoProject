# Classical mechanics and molecular dynamics

[Documentation index](README.md)

## One-coordinate trajectories

The classical page solves:

```text
dx/dt = p/m
dp/dt = -dV/dx
E = p²/(2m) + V(x)
```

It accepts the harmonic, Morse, symmetric double-well and free potentials. The
adapter converts momentum to velocity, divides force by mass, and invokes
velocity-Verlet, leapfrog or Yoshida. Units follow the atomic-unit model convention.
There is no damping, forcing, thermostat or constraint in this mode.

Every step is executed. Output uses stride `ceil(steps/400)` plus the last sample,
producing at most 401 stored times. Charts show x(t), phase space `(x,p)`, and energy.
Metrics report initial/final energy, maximum absolute energy drift over **all**
steps, and final position. Excessive position/energy or non-finite values fail the
job. Smaller error than this hard failure threshold is not automatically accurate.

For harmonic motion with center c, `omega=sqrt(k/m)` and initial velocity p0/m:

```text
x(t) = c + (x0-c)*cos(omega*t) + (p0/(m*omega))*sin(omega*t)
```

This provides a direct reference for checking mass handling, phase error and order.

## Morse bond vibration

The bond mode fixes the potential to Morse:

```text
V(x) = De * [1-exp(-alpha*(x-center))]²
```

Mass is interpreted as reduced mass `mu=m1*m2/(m1+m2)`. Near equilibrium, the
effective spring constant is `2*De*alpha²`, so the small-amplitude angular frequency
is `sqrt(2*De*alpha²/mu)`. Larger excursions are anharmonic. Classical energies
at/above De can dissociate toward the positive-coordinate side.

This is a single bond coordinate, not a multi-atom molecule. No atom identities,
rotations, constraints or chemical units conversion are supplied. The coordinate
is allowed to be negative; choose a center and range consistent with your model.

## Two-dimensional Lennard–Jones fluid

The new `classical.LJSystem` is separate from the old `ModelSystem` skeleton.
It contains arrays `X`, `Y`, `VX`, `VY`, internal forces, box length, cutoff and
current potential energy. Public operations are `NewLJSystem`, `Forces`, `Kinetic`,
`Temperature` and `Step`.

All particles have unit mass in reduced units `m=sigma=epsilon=kB=1`. Coordinates
are two-dimensional and periodic in a square. There are no bonds or atom types.

### Potential, force and cutoff

For pair separation r:

```text
U(r) = 4*(r^-12-r^-6)
F(r) = -dU/dr = 24*(2*r^-13-r^-7)
rc = min(2.5, box/2)
```

Below rc, the code uses a shifted force and matching shifted potential:

```text
Usf(r) = U(r)-U(rc)+(r-rc)*F(rc)
Fsf(r) = F(r)-F(rc)
```

At/above rc the interaction is zero. Both expressions approach zero at rc. Each
pair's vector force is `Fsf(r)*(dx,dy)/r`; equal and opposite forces are added to
the two particles. The reported energy is the shifted/truncated model's energy,
not the infinite-range LJ energy. There is no long-range tail correction.

Minimum-image differences are `dx -= box*round(dx/box)` and similarly for dy.
Each unordered pair is evaluated once. There is no neighbor list; at most 64
particles keeps the all-pairs O(P²) calculation small.

### Initialization

1. Require 4–64 particles, box 5–30, and initial temperature 0–2.
2. Set `side=ceil(sqrt(P))` and `spacing=box/side`; reject spacing below 1.1.
3. Place particles at half-cell offsets on the lattice, filling rows.
4. Draw Gaussian velocities from `math/rand/v2` PCG with seeds `(seed,seed+1)`.
5. Remove mean velocity in each direction, making total momentum approximately zero.
6. Scale velocities to `K=(P-1)*temperature` and calculate initial forces.

In two dimensions with center-of-mass momentum removed, the temperature convention
is `T=2*K/(2*P-2)=K/(P-1)`. The seed reproduces the initialization under the same
implementation, but long chaotic trajectories need not be bitwise identical across
different platforms or future numerical implementations.

### Integration and boundaries

The fluid's velocity-Verlet step is:

```text
v <- v + dt*force/2
x <- x + dt*v
x <- x - box*floor(x/box)     # apply in both directions
recompute all pair forces
v <- v + dt*newForce/2
```

`Step` requires `0<dt<=0.01`; the web minimum is 0.00001. If a pair comes closer
than 0.6, or its distance becomes non-finite, force calculation returns an error.
This is a fail-fast guard, not a repulsive-wall replacement. A failed step can
have partially advanced the public state; the web job discards it. Library callers
must not continue blindly from a failed step.

Public position/velocity arrays should not be resized independently. If a library
caller intentionally edits coordinates, it must recompute forces before stepping.
The application owns each system for one run and never shares it concurrently.

### Results and interpretation

The model evolves without a thermostat: initial temperature is an initial condition,
not a target maintained throughout the run. Potential/kinetic exchange changes
instantaneous temperature while total energy should remain approximately constant.

Output includes total energy and temperature time series, initial/final total
energy, maximum absolute energy drift, and final total momentum magnitude.
Frames store only time and positions, at stride `ceil(steps/150)` including the
final step. The frontend plays those frames at a display interval of 80 ms; this
wall-clock playback speed is unrelated to the physical dt. Periodic wrapping can
make a particle appear to jump from one box edge to the other.

## Model limitations

No thermostat/barostat, pressure estimator, radial distribution function, diffusion
analysis, trajectory export, electrostatics, chemical topology, multiple species,
3D particles, neighbor lists, constraints or equilibration phase is implemented.
An initial lattice and short transient are not automatically an equilibrated fluid.

For a meaningful run, reduce dt and compare energy drift, allow the model to evolve
past initialization transients, and distinguish numerical conservation from physical
validity. Force-gradient and conservation tests are described in [testing](10-testing.md).
