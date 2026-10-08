# HTTP API reference

[Documentation index](README.md) · [Handler](../internal/webserver/server.go)
· [Simulation validation](../internal/simulation/simulation.go)

## Endpoints and transport

Default origin: `http://127.0.0.1:8080`. Vite proxies the same `/api` paths.

| Method / path | Response |
| --- | --- |
| `GET /api/health` | `{"engine":"go","status":"ok"}` |
| `POST /api/plot` | Sampled expression grid |
| `POST /api/simulations` | Diagnostics, charts, optional particle frames |
| `GET /calculator/` | 307 redirect to `/#tools` |
| `GET` / `HEAD` static paths | Files below configured `-dist` |
| Other `/api/` path or unsupported API method | JSON 404 |

POST requests require `Content-Type: application/json` (parameters such as charset
are accepted by media-type parsing), exactly one JSON value and no unknown
top-level request fields. The body limit is 8,192 bytes. Go struct fields default
to zero/empty when absent: the server **does not apply frontend defaults**.
Send the full UI request or all fields required by the selected validator.

`http.CrossOriginProtection` protects browser requests. Matching-origin access or
Vite's proxy is the intended browser setup; there is no permissive CORS layer.
This is not authentication: non-browser clients can still call an accessible
service. The default bind is loopback. General headers include `nosniff` and
`Referrer-Policy: same-origin`; calculation/error responses use `Cache-Control: no-store`.

## Plot request

| Field | Type | Meaning / rules |
| --- | --- | --- |
| `mode` | string | `line`, `heatmap`, `surface`, `iso` |
| `expression` | string | Restricted arithmetic; at most 1,024 bytes / 256 tokens |
| `resolution` | integer | Per-axis sample count; limits below |
| `bounds` | object | Active axis keys (`x`, `y`, `z`), each `[min,max]` |
| `level` | number | Finite; relevant to isosurface |
| `style` | string | Must be `heatmap` or `contour`, even for other modes |

Required axes are x for line, x/y for heatmap and surface, x/y/z for iso. Bounds
must be finite and strictly increasing with a representable sample spacing.
Limits are line 20–10,000, heatmap 10–200, surface 10–150, iso 10–60.

Example request:

```json
{
  "mode": "line",
  "expression": "exp(-x^2/2)",
  "resolution": 101,
  "bounds": {"x": [-5, 5]},
  "level": 0,
  "style": "heatmap"
}
```

The response fields are:

| Field | Meaning |
| --- | --- |
| `settings` | Submitted request echoed by Go |
| `coords` | Axis vectors, each of length resolution |
| `values` | Flat numeric/null array, x-fastest |
| `count` | `resolution^dimensions` |
| `invalid` | Number of null samples |
| `min`, `max` | Finite sampled extrema |

All-invalid plots fail. Isosurfaces additionally reject any null and levels not
strictly inside sampled extrema. See [plotting](07-plotter-and-calculator.md) for indexing.

## Simulation request fields

The frontend sends a common struct for every kind. Some fields are irrelevant to
a particular kind but still exist in the request. Use the exact field names below;
Go's JSON decoder can also accept case-insensitive field-name matches. Duplicate
keys are not explicitly rejected, so clients should send each field only once.

| Field | Type | Meaning / accepted range |
| --- | --- | --- |
| `kind` | string | `spectrum`, `quantum`, `classical`, `molecular`, `md` |
| `basis` | string | `dvr` or `fourier` for spectrum/quantum |
| `method` | string | Allowed method for kind, listed below |
| `potential` | string | `harmonic`, `morse`, `double-well`, `free`; bond kind requires `morse` |
| `mass` | number | 0.1–100 for every non-MD kind |
| `strength` | number | 0.01–100 for every non-MD kind; k, De or double-well strength |
| `alpha` | number | 0.05–2 for every non-MD kind; used by Morse |
| `center` | number | -10–10 for every non-MD kind; used by harmonic/Morse |
| `halfWidth` | number | 2–20 for spectrum/quantum; spatial domain [-L,L) |
| `points` | integer | 32–192 for spectrum/quantum |
| `states` | integer | 1–8 for spectrum/quantum; only spectrum displays eigenstates |
| `dt` | number | 0.00001–0.05 for non-MD evolution; MD maximum 0.01 |
| `steps` | integer | 1–2,000 for evolving kinds |
| `x0` | number | -10–10 for non-MD evolution |
| `p0` | number | -20–20 for non-MD evolution |
| `sigma` | number | 0.1–4 for quantum; also `abs(x0)+3*sigma<halfWidth` |
| `driveAmplitude` | number | -5–5 for quantum |
| `driveFrequency` | number | 0–10 for quantum |
| `particles` | integer | 4–64 for MD |
| `box` | number | 5–30 for MD; initial lattice spacing must be at least 1.1 |
| `temperature` | number | 0–2 for MD, initial temperature only |
| `seed` | integer | 0–1,000,000 for MD |

Every checked numeric range excludes NaN and Infinity. JSON itself cannot encode
them. `strength` and `alpha` are validated for all non-MD modes even when their
potential does not use them; `states` is validated for quantum dynamics even though
that job does not diagonalize displayed eigenstates. This explains why omitting an
apparently unused field can still produce 400. The MD branch validates its own
fields and does not use the non-MD potential/basis parameters.

| Kind | Allowed methods | Additional behavior |
| --- | --- | --- |
| `spectrum` | Method ignored | Dense symmetric H(0) eigensystem |
| `quantum` | `split`, `rk4`, `ralston3`, `nystrom5` | `points²*steps<=40,000,000`; explicit dt admission bound |
| `classical` | `verlet`, `leapfrog`, `yoshida` | Force/mass adapter |
| `molecular` | Same classical methods | Morse potential, reduced mass |
| `md` | `verlet` | 2D periodic LJ model |

Sampled static quantum potentials must stay within ±1,000,000. Excessive norm drift,
divergent classical motion and particle overlaps can fail during integration even
after parameter validation succeeds.

## Request examples

### Energy levels

This is a sufficient spectrum body; omitted fields have no effect on this job:

```json
{
  "kind": "spectrum",
  "basis": "dvr",
  "potential": "harmonic",
  "mass": 1,
  "strength": 1,
  "alpha": 0.5,
  "center": 0,
  "halfWidth": 8,
  "points": 64,
  "states": 4
}
```

### Driven quantum packet

```json
{
  "kind": "quantum",
  "basis": "fourier",
  "method": "split",
  "potential": "harmonic",
  "mass": 1,
  "strength": 1,
  "alpha": 0.5,
  "center": 0,
  "halfWidth": 8,
  "points": 64,
  "states": 4,
  "dt": 0.005,
  "steps": 200,
  "x0": 0,
  "p0": 0,
  "sigma": 0.7,
  "driveAmplitude": 0.2,
  "driveFrequency": 1
}
```

### Particle fluid

```json
{
  "kind": "md",
  "method": "verlet",
  "particles": 16,
  "box": 6,
  "temperature": 0.2,
  "seed": 42,
  "dt": 0.002,
  "steps": 400
}
```

Save a body as `request.json`, then call the local server from a Linux/WSL shell:

```sh
curl -sS http://127.0.0.1:8080/api/simulations \
  -H 'Content-Type: application/json' \
  --data-binary @request.json
```

For the plot example use `/api/plot` instead. These endpoints return data, not
saved jobs; requests do not persist settings or results on the server.

## Simulation response

```text
{
  request: submitted Request (Go also emits zero-valued unused fields),
  title: string,
  units: string,
  metrics: [{label: string, value: finite number}],
  charts: [Chart],
  notes: [string],
  frames?: [{time: number, x: number[], y: number[]}]
}
```

A line-chart object has `title`, `xLabel`, `yLabel`, and
`series:[{name,x:[],y:[]}]`. X/Y lengths must match for each series.
A density chart instead has axis vectors `x`, `y` and rectangular matrix `z`.
For quantum density, `x` is position, `y` is time, and `z[timeIndex][positionIndex]`
is density. Particle `frames` contain wrapped coordinates, not velocities.

Classical charts contain at most 401 samples. Quantum/MD temporal sampling produces
at most 151 samples. Frontend validation allows small bounded headroom: 402 series
points and 152 time rows/frames. Spatial quantum axes are limited to 192. At most
8 series per chart, 6 charts, 16 metrics and 12 notes are accepted by the current
decoder. The client associates the validated response with its own submitted
settings snapshot rather than trusting the echoed request blindly.

## Capacity, timeout and errors

| Status | Typical reason |
| --- | --- |
| 200 | Successful result |
| 400 | Malformed JSON, unknown field, extra JSON value, invalid model or numerical failure |
| 403 | Cross-origin protection rejected a browser request |
| 404 | Unknown API path or unsupported API method |
| 408 | Calculation context deadline exceeded |
| 413 | Oversized request detected during initial JSON decoding |
| 415 | Missing/unsupported content type |
| 429 | No free plot/simulation slot |
| 503 | Static frontend missing; returned by the static handler |

Most handler-generated errors have JSON shape `{"error":"message"}`. Protection
middleware and static-file errors may return other formats. The frontend handles
non-JSON responses as an unavailable service. A client-aborted job returns early
rather than guaranteeing a JSON error response.

Plot budget: two simultaneous jobs, five-second context deadline. Simulation
budget: one job, ten-second context deadline. The calculations page has a
twenty-second client timeout. Server timeouts are 5 seconds for request headers,
10 for reading, 20 for writing and 60 idle; headers are limited to 16 KiB.

No authentication, user quota, durable queue, streaming, pagination, API version
prefix, cache of completed jobs or public-deployment hardening is implemented.
