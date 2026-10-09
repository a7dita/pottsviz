# Pottsviz

Interactive memory-retrieval and latching simulations for homogeneous and frontal-posterior Potts associative networks.

## Local development

Use Linux, Node.js 22, pnpm 11.25.0, and a C++ compiler (`g++`).

```bash
pnpm install --frozen-lockfile
pnpm dev
```

The dev and build commands compile all three C++ simulators automatically. Open the local URL printed by Vite, then choose a model, set its parameters, and press **Start Simulation**. **Stop** cancels the active request. Each start creates a new independent process; plots are cleared before rerunning.

```bash
pnpm check
pnpm test
pnpm build
```

Tests cover both live APIs, simultaneous independent homogeneous runs, parameter validation, paired snapshots, and three numerical regression fixtures from the original Ryom implementation.

## Models

The homogeneous model retains the original equations, random seed 1990, 500 units, 75 incoming connections per unit, 100 memories, and a cue to memory 0. Its default 5,000-sweep output was checked against the original executable: all 1,001 overlap rows match exactly at S=7, w=1.2, tau2=200.

The fronto-posterior model uses the original Ryom dynamics already present in this repository, also checked against `simcode_pilot` at commit `19d8b8cf93b0aa0f66fe5e67abb37332a5d5daeb`. It does not incorporate the later modified adaptation equations from the current topology-study code. Each region has 256 units, 50 presynaptic units selected from each source region, and 50 random memories. Defaults are S=7, w=1.1 and tau2=200 in both regions, with lambda=0.5. Other fixed parameters are a=0.25, U=0.1, T1=20, T3A=10, T3B=100000, beta=11, and gammaA=0.5.

Memory mu in each region is paired reciprocally with memory mu in the other region. Patterns are independent between regions (pattern seeds 1 and 2); matching colours identify the association, not identical unit states. Dynamic seed is 1. Patterns are generated in memory for the selected state counts using the original generator, eliminating dependence on a small set of pre-generated files. The browser allows equal adaptation times or slower frontal adaptation, including the undifferentiated reference condition.

Lambda scales intra-area coupling by (1+lambda) and inter-area coupling by (1-lambda), with the original normalization retained. Lambda=1 isolates the regions; lambda=0 gives equal coefficients. The separate many-to-many model adds topology selection. The code skips unmatched memory pairs and zero source-unit blocks while retaining the order of nonzero sums. Regression fixtures at N=256, C=50, p=50 match the Ryom equations exactly for S_p/S_f=7/7 with lambda=0.9, 3/7 with lambda=0.4, and 7/7 with lambda=1.

Choose 1,000, 2,500, or 5,000 sweeps. The inherited termination rules can finish a run early when retrieval disappears. The fronto-posterior model is substantially heavier than the homogeneous model; live plots begin after network initialization and update during computation.


### Many-to-many topology model

`/models/potts_topology` uses 49 memories in each region for all 16 association
conditions: `original`, `random`, `m0`–`m6`, `s1`–`s6`, and `ms7`.
The original reference is a 49 × 49 identity matrix. The four radio options are
original, random, modular and shared-target. Selecting a structural family
reveals its integer m or s slider; the matrix and fresh scores update immediately.
Both sliders select the same `ms7` endpoint at seven. Parameter descriptions
are available on hover and keyboard focus. The shared-target construction
groups posterior memories and controls how many frontal partners they share.
Groups have no assumed semantic meaning.

Rows are frontal target memories and columns are posterior source memories.
Every CSV and SVG retains generator order, with no shuffling or spectral sort.
Both repositories use the same generator and seed. The generator can be run
with `python3 scripts/generate_topologies.py` (NumPy and SciPy required).
`analysis/topology_design.ipynb` reconstructs the design and displays all 16
matrices. `src/lib/topologies/manifest.json` and `metrics.csv` hold freshly
measured R4 redundancy, approximate maximum Barber Q* and two-step mixing gap.
All scores, including the degree-one original reference, use the actual graph
degree. Complete-component modularity values for original and ms7 are exact;
other Q* values use a deterministic greedy search.

Rows are normalized to unit association mass, with reciprocal connections
using the transpose. Degree-seven edges have weight 1/7. The neural dynamics,
activity adaptation and both inhibition updates retain the heterogeneous model's
Ryom equations.

Default values are N=256, C=50, p=49, a=0.25, U=0.1, beta=11, T1=20, S_p=S_f=7, w_p=0.6, w_f=1.1, T2_p=T2_f=200, and lambda_p=lambda_f=0.5. The right controls expose shared U (0–0.6), beta (1–21), a (0.1–0.4), T1 (5–35), and presynaptic connections per source region c_m (13–64). The default c_m=50 corresponds exactly to C/N=50/256; the density sent to the backend is derived from this integer control. Regional S, w, T2 and lambda controls, cue index 0–48, and a Random checkbox are also available. Selecting Random greys the cue label and disables its slider without changing its saved value. All parameter symbols have concise descriptions on hover or keyboard focus; fixed grid columns keep sliders aligned for parameter symbols, and group headings sit inside their cards. T3A, T3B and gammaA stay fixed at the heterogeneous model's values and are listed in the model details. The plot supports 49 traces automatically.

Every Start creates independent random positive pattern and runtime seeds on the server. Posterior and frontal memories use consecutive pattern seeds; dilution and update tables are rebuilt. Seeds appear in the `started` event and the expandable run details. The association topology stays frozen between runs. Only allowlisted topology IDs are accepted; paths and seeds are not supplied by the browser. Each request owns its native process and streaming response, using the same cancellation and timeout architecture as the existing models. Saved CSVs are explicitly bundled beside the native program in Vercel functions.

Tests check exact upstream files, matrix dimensions and degrees, redundancy values, all 12 native simulations, eight corrected dense numerical reference cases, parity with heterogeneous adaptation/inhibition, and independent cue routing, concurrent freshly seeded API runs, global parameter bounds, paired 49-memory output, all parameter tooltips, and every SVG cell of every topology against its CSV coordinate.

## Vercel deployment

The Vercel project uses Node.js 22 and the SvelteKit Vercel adapter. `vercel.json` pins the pnpm install command. `pnpm build` compiles optimized native programs, builds the app, and copies the programs and the compiler's shared C++ libraries into every generated function. Programs and libraries are built together in Vercel's Linux build image; do not upload binaries compiled against another system's libc.

`POST /api/simulate` accepts numerical model parameters and returns newline-delimited JSON events: `started`, `sample`, `done`, or `error`. Each fronto-posterior sample contains both regions at the same sweep. The frontend buffers only complete JSON records and throttles SVG redraws. It shows preparation, progress, completion, and errors.

The simulation and response remain within one invocation. Generated functions enable request cancellation so a disconnected browser can terminate computation. There is no Python runtime dependency, shared output directory, cross-invocation filesystem polling, arbitrary shell-command endpoint, or user-supplied file path. The legacy `/api` and `/api2` routes return HTTP 410. The retained historical Python wrappers are not used by the application.

The simulation function allows up to 300 seconds; its application deadlines are 55 seconds for homogeneous runs and 270 seconds for paired runs. Use a Vercel configuration that permits a 300-second function (Fluid Compute on Hobby, or an appropriate existing plan). Choose a shorter run if the deadline is reached. Each run is independent and ephemeral; results are not stored.

## Acknowledgements

Kwang Il Ryom's original simulation code, insights, and collaboration made this project possible.
