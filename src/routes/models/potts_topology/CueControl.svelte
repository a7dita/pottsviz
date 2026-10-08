<script lang="ts">
  import ParameterSlider from './ParameterSlider.svelte';
  import katex from 'katex';
  export let region: 'Frontal' | 'Posterior';
  export let cue = 0;
  export let random = false;
  const symbol = katex.renderToString(String.raw`\mu_{\mathrm{cue}}`);
</script>

<div class="cue-control">
  <div class="heading">
    <span>Cue or random</span>
    <label title="A fresh pattern with the same sparsity and active states, not stored among the 49 memories. Recreated on every Start.">
      <input type="checkbox" bind:checked={random} aria-label={`${region} random cue`} /> Random
    </label>
  </div>
  <ParameterSlider ariaLabel={`${region} cue memory index`} labelName={symbol} minValue={0} maxValue={48} bind:value={cue} stepSize={1} disabled={random}
    tooltipText="Stored memory index receiving the transient cue in this region. Random uses a fresh, unstored pattern instead." />
  <output class:muted={random}>Cue = {cue}</output>
</div>

<style>
  .cue-control { margin-top: 1.1rem; padding-top: .8rem; border-top: 1px solid #e5e7eb; }
  .heading { display: flex; justify-content: space-between; align-items: center; gap: .6rem; margin-bottom: .7rem; font-size: .8rem; }
  label { display: flex; gap: .4rem; align-items: center; cursor: pointer; }
  input { accent-color: #7c3aed; }
  output { display: block; margin-top: .3rem; font-size: .75rem; color: #6d28d9; text-align: right; }
  output.muted { color: #b7bac4; }
</style>
