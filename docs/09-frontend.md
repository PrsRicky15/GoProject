# Frontend structure and behavior

[Documentation index](README.md) · [Frontend source](../webpage/web-app/src/App.tsx)

## Entry points and routes

`index.html` loads `src/main.tsx`, which mounts React in `StrictMode`. `App.tsx`
tracks `window.location.hash` and subscribes to `hashchange`:

| Hash | Component |
| --- | --- |
| `#tools` | Lazy-loaded `Workbench` |
| `#calculations` | Lazy-loaded `Calculations` |
| Other/empty | `HomePage` |

There is no React Router dependency. `#features` and `#solutions` refer to landing
page sections. A Suspense fallback appears while a workspace chunk loads.
`HomePage` is imported eagerly; the two workspaces and plotting bundles are deferred.

## Component map

| File | Responsibility |
| --- | --- |
| `pages/HomePage.tsx` | Landing layout, animated feature selection, workspace links |
| `components/SolutionsSection.tsx` | Landing-page use-case content |
| `pages/Workbench.tsx` | Plot drafts, submission, status, visualization and calculator composition |
| `components/Calculator.tsx` | Input, keypad, local answer/history and parser errors |
| `pages/Calculations.tsx` | Model form, numerical job lifecycle and results |
| `components/SimulationChart.tsx` | Convert line/density chart data to Plotly figures |
| `components/ParticlePlayback.tsx` | SVG particles, frame slider and play/pause |
| `components/PlotCanvas.tsx` | Shared Plotly lifecycle, resizing and rendering |

The existing landing page contains placeholder controls such as Sign In and Watch
Demo. Its cloud/GPU/collaboration wording is not a capability contract; those
services are absent. The working paths are the two linked workspaces.

## Library modules

| Module | Responsibility |
| --- | --- |
| `lib/types.ts` | Plot mode, axis, draft, request, result and render-status types |
| `lib/sampling.ts` | Mode metadata, defaults, form normalization; does not sample plots |
| `lib/expression.ts` | Restricted calculator grammar and evaluation |
| `lib/api.ts` | Plot fetch and runtime response decoder |
| `lib/figure.ts` | Plotly trace/layout construction and grid reshaping |
| `lib/simulation.ts` | Calculation types/defaults, API call and response decoder |
| `plotly.d.ts` | Declarations for the partial Plotly distribution modules |

`PlotDraft` permits strings for edited numeric fields; `validateRequest` converts
them into a numeric `PlotRequest`. JSON is treated as unknown until checked.
Plot results become `Float64Array` values; simulation charts currently retain
ordinary numeric arrays. Compile-time interfaces do not validate network input.

## Calculation form state

`Calculations` keeps separate draft, submitted request, result, busy and error
state. Selecting another kind resets that draft to its defaults. It does not
silently reinterpret a previous result: the result carries the settings from its
run and receives a “Previous run · settings changed” badge when they differ.

The form shows only relevant controls. Quantum modes expose basis/grid settings;
evolving modes expose dt/steps; driven quantum adds amplitude/frequency; MD exposes
particle count, box, initial temperature and seed. Bond vibration fixes Morse,
and MD fixes velocity-Verlet. HTML numeric constraints provide early feedback;
Go remains authoritative.

Run snapshots draft values, starts an `AbortController`, and fetches the API.
Cancellation aborts the request; unmount/changed request cleans it up. Aborted
responses are ignored, and the timeout is cleared when the request settles.
The client deadline is 20 seconds. Go checks cancellation between numerical steps.

Reset changes the current model's controls; it does not erase a previous result.
Parameters for a completed run remain inspectable in a disclosure below the plots.
Errors are displayed through an alert region; calculation status/results use live
regions. Results are held only in memory.

## Plot rendering lifecycle

`PlotCanvas` accepts either a typed plot result or a supplied figure/label for a
simulation chart. It creates a dedicated DOM node and serializes asynchronous
Plotly updates through a promise queue. A stale/unmounted render is ignored.
`Plotly.react` updates an existing plot rather than replacing the whole component.

A `ResizeObserver` schedules resize work with `requestAnimationFrame`. Cleanup
disconnects observers, cancels a pending frame, removes the node and purges Plotly
resources after queued work. This matters when moving between routes or modes.

The Cartesian bundle serves line/heatmap/contour charts. The GL3D bundle serves
surfaces and isosurfaces. Spatial rendering checks for WebGL and produces a useful
error if unavailable. The Plotly toolbar exports PNG at scale 2, with filenames
based on the plot mode or `calculation`.

`figure.ts` keeps the isosurface typing bridge in one place because the installed
Plotly type declarations omit some trace details. This cast is a library typing
adaptation, not permission to skip validation of returned arrays.

## Particle playback

The player receives frames and box size. Each frame is drawn as an SVG box with
circles, with a y-axis display flip. Play advances the frame index every 80 ms and
wraps after the last frame. Pause stops the interval. Slider changes pause playback.
Unmount clears the interval. The frame's simulation time is displayed independently
of playback wall-clock time.

All integration remains in Go. A frame does not contain velocity or force, and
the frontend does not interpolate a new physical solution between stored frames.

## Styling and build configuration

`index.css` supplies shared/global styling; `workbench.css` styles the scientific
workspace; `calculations.css` supplies model controls, diagnostics, charts and
responsive layouts. Calculation charts become one column below 1150 px, and the
controls/results layout becomes one column below 700 px.

TypeScript targets ES2022, uses strict mode, bundler module resolution, isolated
modules and noEmit. Vite performs the actual bundling. `skipLibCheck` skips type
checking dependency declarations; it does not disable application strictness.
ESLint applies TypeScript, React Hooks and refresh rules. Tailwind/PostCSS provide
the existing landing-page styling pipeline.

The build command runs TypeScript checks first. Plotly still produces large vendor
chunks despite lazy loading; the warning is documented rather than hidden. No
service worker, offline cache, persistence layer, browser compute worker or frontend
simulation engine is configured.

See [development](12-development.md) for coordinated changes to types, validation,
Go fields, controls and tests.
