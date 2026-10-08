<script lang="ts">
  import { onDestroy } from 'svelte';
  import SliderParam from './SliderParam.svelte';
  import OverlapPlot from '$lib/OverlapPlot.svelte';
  import { runSimulation, type Sample } from '$lib/simulation-client';
  import katex from 'katex';
  import 'katex/dist/katex.min.css';
  const S = katex.renderToString('S'), W = katex.renderToString('w');
  const Tau2 = katex.renderToString('\\tau_2'), L = katex.renderToString('\\lambda');
  let frontalS = 7, posteriorS = 7, frontalW = 1.1, posteriorW = 1.1;
  let frontalTau = 400, posteriorTau = 100, lambda = 0.9;
  let steps = 5000;
  let samples: Sample[] = [], isRunning = false, status = 'Ready', error = '';
  let abort: AbortController | undefined;
  let refresh: ReturnType<typeof setTimeout> | undefined;
  $: validTimes = frontalTau > posteriorTau;
  onDestroy(() => { abort?.abort(); clearTimeout(refresh); });
  async function handleClick() {
    if (isRunning) { abort?.abort(); return; }
    const run = new AbortController(); abort = run;
    isRunning = true; samples = []; error = ''; status = 'Preparing the two networks…';
    const collected: Sample[] = [];
    try {
      await runSimulation({ model: 'hybrid', topology: 'one-to-one', lambda, steps,
        posterior: { S: posteriorS, w: posteriorW, tau2: posteriorTau },
        frontal: { S: frontalS, w: frontalW, tau2: frontalTau }
      }, run.signal, sample => {
        collected.push(sample);
        if (!refresh) refresh = setTimeout(() => {
          samples = [...collected]; status = `Running · ${samples[samples.length-1].time} / ${steps} sweeps`;
          refresh = undefined;
        }, 100);
      });
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

<svelte:head><title>Fronto-posterior Potts · Pottsviz</title></svelte:head>
<div class="space-y-5 p-4 text-gray-700" style="width:min(1240px,100vw)">
  <h1 class="text-3xl text-purple text-center">Fronto-posterior Potts Network</h1>
  <p class="text-center">One-to-one memory pairs · slower frontal adaptation · reciprocal coupling</p>
  <div class="flex flex-wrap items-center justify-center gap-6 bg-sky-500/[.06] rounded p-4">
    <fieldset disabled={isRunning} style="width:320px">
      <SliderParam ariaLabel="Inter-area coupling lambda" labelName={L} minValue={0} maxValue={1} bind:value={lambda} stepSize={0.1}/>
      <p class="text-center mt-2">{@html L} = {lambda.toFixed(1)}</p>
    </fieldset>
    <p class="text-sm" style="max-width:440px">λ = 1 isolates the networks; λ = 0 gives equal intra-area and inter-area coupling coefficients.</p>
    <label class="flex gap-2 items-center">Sweeps
      <select aria-label="Simulation length" disabled={isRunning} bind:value={steps} class="border rounded p-1">
        <option value={1000}>1,000</option><option value={2500}>2,500</option><option value={5000}>5,000</option>
      </select>
    </label>
    <button class="px-4 py-2 rounded bg-purple-100 disabled:opacity-50" disabled={!isRunning && !validTimes} on:click={handleClick}>{isRunning ? 'Stop' : 'Start Simulation'}</button>
  </div>
  <p role="status" aria-live="polite" class="text-center">{status}</p>
  {#if !validTimes}<p class="text-red-700 text-center">Choose a frontal adaptation time larger than the posterior adaptation time.</p>{/if}
  {#if error}<p role="alert" class="text-red-700 text-center">{error}</p>{/if}
  <div class="grid grid-cols-1 md:grid-cols-2 gap-6">
    <section class="space-y-4">
      <h2 class="text-2xl text-center">Frontal · slow</h2>
      <div class="border-2 border-gray-200 rounded p-2"><OverlapPlot {samples} region="frontal"/></div>
      <fieldset disabled={isRunning} class="space-y-4 bg-sky-500/[.06] rounded p-4">
        <SliderParam ariaLabel="Frontal active states" labelName={S} minValue={3} maxValue={11} bind:value={frontalS} stepSize={1}/>
        <SliderParam ariaLabel="Frontal self reinforcement" labelName={W} minValue={0.6} maxValue={2} bind:value={frontalW} stepSize={0.1}/>
        <SliderParam ariaLabel="Frontal adaptation time" labelName={Tau2} minValue={100} maxValue={1600} bind:value={frontalTau} stepSize={100}/>
        <p>S = {frontalS} · w = {frontalW.toFixed(1)} · τ₂ = {frontalTau}</p>
      </fieldset>
    </section>
    <section class="space-y-4">
      <h2 class="text-2xl text-center">Posterior · fast</h2>
      <div class="border-2 border-gray-200 rounded p-2"><OverlapPlot {samples} region="posterior"/></div>
      <fieldset disabled={isRunning} class="space-y-4 bg-sky-500/[.06] rounded p-4">
        <SliderParam ariaLabel="Posterior active states" labelName={S} minValue={3} maxValue={11} bind:value={posteriorS} stepSize={1}/>
        <SliderParam ariaLabel="Posterior self reinforcement" labelName={W} minValue={0.6} maxValue={2} bind:value={posteriorW} stepSize={0.1}/>
        <SliderParam ariaLabel="Posterior adaptation time" labelName={Tau2} minValue={100} maxValue={800} bind:value={posteriorTau} stepSize={100}/>
        <p>S = {posteriorS} · w = {posteriorW.toFixed(1)} · τ₂ = {posteriorTau}</p>
      </fieldset>
    </section>
  </div>
  <p class="text-sm">500 units and 100 random memories in each region · matching colours identify paired memories · both plots share the same time axis.</p>
</div>
