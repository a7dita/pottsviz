<script lang="ts">
  import { onDestroy } from 'svelte';
  import SliderParam from './SliderParam.svelte';
  import OverlapPlot from '$lib/OverlapPlot.svelte';
  import { runSimulation, type Sample } from '$lib/simulation-client';
  import katex from 'katex';
  import 'katex/dist/katex.min.css';
  const S = katex.renderToString('S');
  const W = katex.renderToString('w');
  const Tau2 = katex.renderToString('\\tau_2');
  let valueS = 7, valueW = 1.2, valueTau2 = 200;
  let steps = 5000;
  let samples: Sample[] = [];
  let isRunning = false;
  let status = 'Ready', error = '';
  let abort: AbortController | undefined;
  let refresh: ReturnType<typeof setTimeout> | undefined;
  onDestroy(() => { abort?.abort(); clearTimeout(refresh); });
  async function handleClick() {
    if (isRunning) { abort?.abort(); return; }
    const run = new AbortController();
    abort = run;
    isRunning = true;
    error = '';
    samples = [];
    status = 'Preparing the network…';
    const collected: Sample[] = [];
    try {
      await runSimulation({ model: 'homo', S: valueS, w: valueW, tau2: valueTau2, steps }, run.signal, sample => {
        collected.push(sample);
        if (!refresh) refresh = setTimeout(() => {
          samples = [...collected];
          status = `Running · ${samples[samples.length-1].time} / ${steps} sweeps`;
          refresh = undefined;
        }, 100);
      });
      status = 'Simulation complete';
    } catch (cause) {
      if (run.signal.aborted) status = 'Stopped';
      else { error = cause instanceof Error ? cause.message : 'Simulation failed.'; status = 'Error'; }
    } finally {
      clearTimeout(refresh); refresh = undefined;
      samples = [...collected];
      isRunning = false;
      abort = undefined;
    }
  }
</script>

<svelte:head><title>Homogeneous Potts · Pottsviz</title></svelte:head>
<div class="space-y-6 p-4 text-gray-700">
  <h1 class="text-3xl text-purple text-center">Potts Associative Network (Homogeneous)</h1>
  <div class="flex flex-wrap gap-8">
    <div class="space-y-4">
      <div class="space-y-2 flex flex-col items-center bg-sky-500/[.06] rounded p-4">
        <p>Your choices:</p><p>{@html S} = {valueS}</p><p>{@html W} = {valueW}</p><p>{@html Tau2} = {valueTau2}</p>
      </div>
      <fieldset disabled={isRunning} class="space-y-4">
        <SliderParam ariaLabel="Number of active states" labelName={S} minValue={3} maxValue={11} bind:value={valueS} stepSize={1}/>
        <SliderParam ariaLabel="Self reinforcement" labelName={W} minValue={0} maxValue={2} bind:value={valueW} stepSize={0.2}/>
        <SliderParam ariaLabel="Adaptation time" labelName={Tau2} minValue={100} maxValue={800} bind:value={valueTau2} stepSize={100}/>
      </fieldset>
      <label class="flex gap-2 items-center">Sweeps
      <select aria-label="Simulation length" disabled={isRunning} bind:value={steps} class="border rounded p-1">
        <option value={1000}>1,000</option><option value={2500}>2,500</option><option value={5000}>5,000</option>
      </select>
    </label>
    <button class="px-4 py-2 rounded bg-purple-100" on:click={handleClick}>{isRunning ? 'Stop' : 'Start Simulation'}</button>
      <p role="status" aria-live="polite">{status}</p>
      {#if error}<p role="alert" class="text-red-700">{error}</p>{/if}
    </div>
    <div class="border-2 border-gray-200 rounded p-4 w-[700px] max-w-full"><OverlapPlot {samples}/></div>
  </div>
  <p class="text-sm">500 units · 100 memories · each colour follows one memory. The plot updates as the simulation runs.</p>
</div>
