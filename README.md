# Pottsviz

Interactive memory-retrieval and latching simulations for homogeneous and frontal-posterior Potts associative networks.

## Local development

Use Linux, Node.js 22, pnpm 11.25.0, and a C++ compiler (`g++`).

```bash
pnpm install --frozen-lockfile
pnpm dev
```

The dev and build commands compile both C++ simulators automatically. Open the local URL printed by Vite, then choose a model, set its parameters, and press **Start Simulation**. **Stop** cancels the active request. Each start creates a new independent process; plots are cleared before rerunning.

```bash
pnpm check
pnpm test
pnpm build
```

Tests cover both live APIs, simultaneous independent homogeneous runs, parameter validation, paired snapshots, and three numerical regression fixtures from the original Ryom implementation.

## Models

The homogeneous model retains the original equations, random seed 1990, 500 units, 75 incoming connections per unit, 100 memories, and a cue to memory 0. Its default 5,000-sweep output was checked against the original executable: all 1,001 overlap rows match exactly at S=7, w=1.2, tau2=200.

The fronto-posterior model uses the original Ryom dynamics already present in this repository, also checked against `simcode_pilot` at commit `19d8b8cf93b0aa0f66fe5e67abb37332a5d5daeb`. It does not incorporate the later modified adaptation equations from the current topology-study code. Each region has 500 units, 75 incoming connections per unit, and 100 random memories. Defaults are S=7 and w=1.1 in both regions, posterior tau2=100, frontal tau2=400, and lambda=0.9. Other fixed parameters are a=0.25, U=0.1, T1=20, T3A=10, T3B=100000, beta=11, and gammaA=0.5.

Memory mu in each region is paired reciprocally with memory mu in the other region. Patterns are independent between regions (pattern seeds 1 and 2); matching colours identify the association, not identical unit states. Dynamic seed is 1. Patterns are generated in memory for the selected state counts using the original generator, eliminating dependence on a small set of pre-generated files. The browser requires frontal adaptation to remain slower than posterior adaptation.

Lambda scales intra-area coupling by (1+lambda) and inter-area coupling by (1-lambda), with the original normalization retained. Lambda=1 isolates the regions; lambda=0 gives equal coefficients. Many-to-many topology selection is deferred. The code skips unmatched memory pairs and zero source-unit blocks while retaining the order of nonzero sums. Regression fixtures at N=50, p=10 match the original single-threaded implementation exactly for S_p/S_f=7/7 with lambda=0.9, 3/7 with lambda=0.4, and 7/7 with lambda=1.

Choose 1,000, 2,500, or 5,000 sweeps. The inherited termination rules can finish a run early when retrieval disappears. The fronto-posterior model is substantially heavier than the homogeneous model; live plots begin after network initialization and update during computation.

## Vercel deployment

The Vercel project uses Node.js 22 and the SvelteKit Vercel adapter. `vercel.json` pins the pnpm install command. `pnpm build` compiles optimized native programs, builds the app, and copies the programs and the compiler's shared C++ libraries into every generated function. Programs and libraries are built together in Vercel's Linux build image; do not upload binaries compiled against another system's libc.

`POST /api/simulate` accepts numerical model parameters and returns newline-delimited JSON events: `started`, `sample`, `done`, or `error`. Each fronto-posterior sample contains both regions at the same sweep. The frontend buffers only complete JSON records and throttles SVG redraws. It shows preparation, progress, completion, and errors.

The simulation and response remain within one invocation. Generated functions enable request cancellation so a disconnected browser can terminate computation. There is no Python runtime dependency, shared output directory, cross-invocation filesystem polling, arbitrary shell-command endpoint, or user-supplied file path. The legacy `/api` and `/api2` routes return HTTP 410. The retained historical Python wrappers are not used by the application.

The simulation function allows up to 300 seconds; its application deadlines are 55 seconds for homogeneous runs and 270 seconds for paired runs. Use a Vercel configuration that permits a 300-second function (Fluid Compute on Hobby, or an appropriate existing plan). Choose a shorter run if the deadline is reached. Each run is independent and ephemeral; results are not stored.

## Acknowledgements

Kwang Il Ryom's original simulation code, insights, and collaboration made this project possible.
