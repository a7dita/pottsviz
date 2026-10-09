<script lang="ts">
  import { onDestroy } from 'svelte';
  import SliderParam from './ParameterSlider.svelte';
  import CueControl from './CueControl.svelte';
  import TopologyView from './TopologyView.svelte';
  import OverlapPlot from '$lib/OverlapPlot.svelte';
  import { topologies } from '$lib/topologies';
  import { runSimulation, type Sample, type RunInfo } from '$lib/simulation-client';
  import katex from 'katex';
  import 'katex/dist/katex.min.css';
  const math = (value: string) => katex.renderToString(value);
  const S = math('S'), W = math('w'), Tau2 = math('\\tau_2'), L = math('\\lambda');
  let topologyId = 'T00';
  $: selected = topologies.find(item => item.id === topologyId) ?? topologies[0];
  let frontalS = 7, posteriorS = 7, frontalW = 1.1, posteriorW = 1.1;
  let frontalTau = 200, posteriorTau = 200, frontalLambda = 0.5, posteriorLambda = 0.5;
  let U = 0.1, beta = 11, a = 0.25, tau1 = 20, connections = 50;
  $: density = connections / 256;
  let frontalCue = 0, posteriorCue = 0, frontalRandom = false, posteriorRandom = false;
  let steps = 2500;
  let samples: Sample[] = [], isRunning = false, status = 'Ready', error = '';
  let runInfo: RunInfo | undefined;
  let abort: AbortController | undefined;
  let refresh: ReturnType<typeof setTimeout> | undefined;
  $: validTimes = frontalTau >= posteriorTau;
  onDestroy(() => { abort?.abort(); clearTimeout(refresh); });
  async function handleClick() {
    if (isRunning) { abort?.abort(); return; }
    const run = new AbortController(); abort = run;
    isRunning = true; samples = []; runInfo = undefined; error = '';
    status = `Preparing ${topologyId} with fresh memories…`;
    const collected: Sample[] = [];
    try {
      await runSimulation({ model: 'topology', topology: topologyId, steps,
        posterior: { S: posteriorS, w: posteriorW, tau2: posteriorTau, lambda: posteriorLambda, cue: posteriorCue, randomCue: posteriorRandom },
        frontal: { S: frontalS, w: frontalW, tau2: frontalTau, lambda: frontalLambda, cue: frontalCue, randomCue: frontalRandom },
        global: { U, beta, a, tau1, density }
      }, run.signal, sample => {
        collected.push(sample);
        if (!refresh) refresh = setTimeout(() => {
          samples = [...collected]; status = `Running · ${samples[samples.length-1].time} / ${steps} sweeps`;
          refresh = undefined;
        }, 100);
      }, info => { runInfo = info; });
      status = 'Simulation complete';
    } catch (cause) {
      if (run.signal.aborted) status = 'Stopped';
      else { error = cause instanceof Error ? cause.message : 'Simulation failed.'; status = 'Error'; }
    } finally {
      clearTimeout(refresh); refresh = undefined; samples = [...collected];
      isRunning = false; abort = undefined;
    }
  }
</script>

<svelte:head><title>Memory association topologies · Pottsviz</title></svelte:head>
<main class="page">
  <a class="back" href="/models">← All models</a>
  <header><h1>Many-to-many Potts Network</h1><p>Explore memory associations between frontal and posterior networks, with equal adaptation times by default.</p></header>
  <div class="toolbar">
    <label title="Number of full update sweeps over both regions.">Sweeps <select disabled={isRunning} bind:value={steps} aria-label="Simulation length"><option value={1000}>1,000</option><option value={2500}>2,500</option><option value={5000}>5,000</option></select></label>
    <button class="start" disabled={!isRunning && !validTimes} on:click={handleClick}>{isRunning ? 'Stop' : 'Start Simulation'}</button>
    <p role="status" aria-live="polite">{status}</p>
  </div>
  {#if !validTimes}<p class="error">Choose a frontal adaptation time at least as large as the posterior adaptation time.</p>{/if}
  {#if error}<p role="alert" class="error">{error}</p>{/if}
  <div class="workspace">
    <aside class="card topology">
      <h2>Association topology</h2>
      <label for="topology">Select a topology</label>
      <select id="topology" disabled={isRunning} bind:value={topologyId} on:change={() => { samples = []; runInfo = undefined; status = 'Ready'; error = ''; }}>
        {#each topologies as item}<option value={item.id}>{item.id} · {item.name.replaceAll('_', ' ')}</option>{/each}
      </select>
      <TopologyView topology={selected} />
      <p class="note">Purple squares are reciprocal associations. Each row has total weight 1: the many-to-many graphs distribute it across seven neighbours.</p>
    </aside>
    <section class="plots">
      <div class="card"><h2>Frontal</h2><OverlapPlot {samples} region="frontal" /></div>
      <div class="card"><h2>Posterior</h2><OverlapPlot {samples} region="posterior" /></div>
      <p class="note">256 units and 49 fresh random memories per region. Colours identify memory indices; read the matrix to see their associations.</p>
      <details class="card">
        <summary>Model and run details</summary>
        <p class="note">Adaptation tracks state activity σ. Fast and slow activity-dependent inhibition use T₃A = 10, T₃B = 100000 and γA = 0.5, as in the heterogeneous model.</p>
        <p class="note">λ scales intra-area connections by (1 + λ) and incoming inter-area connections by (1 − λ). At λ = 1 the region receives no inter-area input. Each region receives its own transient cue: a selected stored memory or a fresh, unstored random pattern. Cue strength is 5 and decays with a 40-sweep time scale; it is removed after 360 sweeps.</p>
        {#if runInfo}<p class="note">Topology {runInfo.topology} · pattern seed {runInfo.patternSeed} (posterior) / {(runInfo.patternSeed ?? 0)+1} (frontal) · runtime seed {runInfo.runtimeSeed}</p>{/if}
        <p class="note">Every Start regenerates both memory sets, unit connection dilution and update tables, then uses a fresh runtime seed. The selected association graph and its scores remain fixed.</p>
      </details>
    </section>
    <aside class="parameters">
      <section class="card parameter-group">
        <fieldset disabled={isRunning}>
        <legend>Entire network</legend>
        <p class="note">Shared by both regions.</p>
        <div class="parameter"><SliderParam ariaLabel="Global threshold U" labelName={math('U')} minValue={0} maxValue={0.6} bind:value={U} stepSize={0.01} tooltipText="Quiet-state threshold. A larger U favours inactivity."/><output>U = {U.toFixed(2)}</output></div>
        <div class="parameter"><SliderParam ariaLabel="Inverse temperature beta" labelName={math('\\beta')} minValue={1} maxValue={21} bind:value={beta} stepSize={1} tooltipText="Inverse temperature. A larger beta sharpens competition between states."/><output>β = {beta}</output></div>
        <div class="parameter"><SliderParam ariaLabel="Memory sparsity" labelName={math('a')} minValue={0.1} maxValue={0.4} bind:value={a} stepSize={0.01} tooltipText="Fraction of active units in each newly generated memory."/><output>a = {a.toFixed(2)}</output></div>
        <div class="parameter"><SliderParam ariaLabel="Integration time" labelName={math('\\tau_1')} minValue={5} maxValue={35} bind:value={tau1} stepSize={1} tooltipText="Time scale for the integration of incoming fields."/><output>τ₁ = {tau1}</output></div>
        <div class="parameter"><SliderParam ariaLabel="Incoming connections per source region" labelName={math('c_m')} minValue={13} maxValue={64} bind:value={connections} stepSize={1} tooltipText="Presynaptic units selected from each source region. The reference value is 50 out of 256; this is separate from memory association degree."/><output>cₘ = {connections} · C/N = {density.toFixed(5)}</output></div>
      </fieldset>
      </section>
      <section class="card parameter-group">
        <fieldset disabled={isRunning}>
        <legend>Frontal</legend>
        <div class="parameter"><SliderParam ariaLabel="Frontal active states" tooltipText="Number of active states per unit, in addition to the quiet state." labelName={S} minValue={3} maxValue={11} bind:value={frontalS} stepSize={1}/><output>S = {frontalS}</output></div>
        <div class="parameter"><SliderParam ariaLabel="Frontal self reinforcement" tooltipText="Self-reinforcement of active states. Larger w supports persistence and retrieval." labelName={W} minValue={0.6} maxValue={2} bind:value={frontalW} stepSize={0.1}/><output>w = {frontalW.toFixed(1)}</output></div>
        <div class="parameter"><SliderParam ariaLabel="Frontal adaptation time" tooltipText="Time scale of state-specific fatigue. Larger τ₂ slows adaptation and memory switching." labelName={Tau2} minValue={100} maxValue={1600} bind:value={frontalTau} stepSize={100}/><output>τ₂ = {frontalTau}</output></div>
        <div class="parameter"><SliderParam ariaLabel="Frontal coupling lambda" tooltipText="Scales this region’s internal connections by (1 + λ) and incoming cross-region connections by (1 − λ). At λ = 1, cross-region input is off." labelName={L} minValue={0} maxValue={1} bind:value={frontalLambda} stepSize={0.05}/><output>λ = {frontalLambda.toFixed(2)}</output></div>
        <CueControl region="Frontal" bind:cue={frontalCue} bind:random={frontalRandom} />
      </fieldset>
      </section>
      <section class="card parameter-group">
        <fieldset disabled={isRunning}>
        <legend>Posterior</legend>
        <div class="parameter"><SliderParam ariaLabel="Posterior active states" tooltipText="Number of active states per unit, in addition to the quiet state." labelName={S} minValue={3} maxValue={11} bind:value={posteriorS} stepSize={1}/><output>S = {posteriorS}</output></div>
        <div class="parameter"><SliderParam ariaLabel="Posterior self reinforcement" tooltipText="Self-reinforcement of active states. Larger w supports persistence and retrieval." labelName={W} minValue={0.6} maxValue={2} bind:value={posteriorW} stepSize={0.1}/><output>w = {posteriorW.toFixed(1)}</output></div>
        <div class="parameter"><SliderParam ariaLabel="Posterior adaptation time" tooltipText="Time scale of state-specific fatigue. Larger τ₂ slows adaptation and memory switching." labelName={Tau2} minValue={100} maxValue={800} bind:value={posteriorTau} stepSize={100}/><output>τ₂ = {posteriorTau}</output></div>
        <div class="parameter"><SliderParam ariaLabel="Posterior coupling lambda" tooltipText="Scales this region’s internal connections by (1 + λ) and incoming cross-region connections by (1 − λ). At λ = 1, cross-region input is off." labelName={L} minValue={0} maxValue={1} bind:value={posteriorLambda} stepSize={0.05}/><output>λ = {posteriorLambda.toFixed(2)}</output></div>
        <CueControl region="Posterior" bind:cue={posteriorCue} bind:random={posteriorRandom} />
      </fieldset>
      </section>
    </aside>
  </div>
</main>

<style>
  .page { width: min(1540px, 100%); min-width: 0; box-sizing: border-box; padding: 1.5rem; color: #374151; }
  .back { font-size: .85rem; color: #6d28d9; }
  header { text-align: center; margin: 1rem 0 1.5rem; }
  h1 { color: #6d28d9; font-size: 1.9rem; margin-bottom: .4rem; }
  h2, legend { font-size: 1.05rem; font-weight: 600; margin-bottom: .7rem; }
  .toolbar { display: flex; flex-wrap: wrap; align-items: center; justify-content: center; gap: 1rem; margin-bottom: 1.5rem; }
  .toolbar label { display: flex; align-items: center; gap: .5rem; }
  select { border: 1px solid #d1d5db; border-radius: .4rem; padding: .45rem; background: white; max-width: 100%; }
  .topology select { width: 100%; margin: .5rem 0; font-size: .85rem; }
  .topology label { font-size: .8rem; }
  .start { background: #6d28d9; color: white; padding: .6rem 1.2rem; border-radius: .4rem; }
  button:disabled { opacity: .5; }
  .workspace { display: grid; grid-template-columns: minmax(240px, 290px) minmax(360px, 1fr) minmax(275px, 310px); gap: 1.2rem; align-items: start; }
  .card { border: 1px solid #e5e7eb; padding: 1rem; border-radius: .7rem; background: white; }
  .plots, .parameters { display: grid; gap: 1rem; }
  .parameters .card { background: #faf8ff; min-width: 0; }
  fieldset { border: 0; padding: 0; margin: 0; min-width: 0; width: 100%; }
  legend { display: block; float: none; width: 100%; padding: 0; margin: 0 0 .7rem; line-height: 1.4; }
  .parameter { margin-top: .9rem; }
  output { display: block; font-size: .75rem; color: #6d28d9; margin-top: .3rem; text-align: right; }
  .note { font-size: .8rem; line-height: 1.5; color: #6b7280; margin-top: .6rem; }
  .error { color: #b91c1c; text-align: center; margin: .75rem; }
  summary { cursor: pointer; font-size: .85rem; }
  @media (max-width: 1100px) { .workspace { grid-template-columns: 260px minmax(0, 1fr); } .parameters { grid-column: 1 / -1; grid-template-columns: repeat(3, minmax(0, 1fr)); } }
  @media (max-width: 700px) { .page { padding: .75rem; } .workspace, .parameters { grid-template-columns: minmax(0, 1fr); } .parameters { grid-column: auto; } .topology { max-width: 380px; width: 100%; justify-self: center; } }
</style>
