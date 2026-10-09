<script lang="ts">
  import { createEventDispatcher } from 'svelte';
  import SliderParam from './ParameterSlider.svelte';
  import { topologyForFamily, type TopologyFamily } from '$lib/topologies';
  export let topologyId = 'original';
  export let disabled = false;
  let family: TopologyFamily = 'original';
  let m = 0, s = 1;
  const dispatch = createEventDispatcher<{ change: string }>();
  const families: { value: TopologyFamily; label: string }[] = [
    { value: 'original', label: 'original' }, { value: 'random', label: 'random' },
    { value: 'modular', label: 'modular' }, { value: 'shared-target', label: 'shared-target' }
  ];
  function selectTopology() {
    topologyId = topologyForFamily(family, family === 'modular' ? m : s).id;
    dispatch('change', topologyId);
  }
</script>

<fieldset {disabled}>
  <legend>Topology family</legend>
  <div class="families">
    {#each families as option}
      <label><input type="radio" name="topology-family" value={option.value} bind:group={family} on:change={(event) => { family = event.currentTarget.value as TopologyFamily; selectTopology(); }} />{option.label}</label>
    {/each}
  </div>
  {#if family === 'modular'}
    <div class="family-slider" on:input={(event) => { m = Number((event.target as HTMLInputElement).value); selectTopology(); }}>
      <SliderParam ariaLabel="Modular parameter m" labelName="m" minValue={0} maxValue={7} bind:value={m} stepSize={1} {disabled}
        tooltipText="m is the number of each memory’s seven partners inside its corresponding structural group. The remaining 7 − m partners are outside. Groups need no semantic meaning. At m = 7, ms7 has seven separate association blocks." />
      <output>m = {m} · {m === 7 ? 'ms7' : `m${m}`}</output>
    </div>
  {:else if family === 'shared-target'}
    <div class="family-slider" on:input={(event) => { s = Number((event.target as HTMLInputElement).value); selectTopology(); }}>
      <SliderParam ariaLabel="Shared-target parameter s" labelName="s" minValue={1} maxValue={7} bind:value={s} stepSize={1} {disabled}
        tooltipText="s is the number of frontal partners shared by every pair of posterior memories in a structural group of seven. Every memory keeps seven partners. s = 1 represents the equivalent s = 0 construction. At s = 7, ms7 is the same endpoint as m = 7." />
      <output>s = {s} · {s === 7 ? 'ms7' : `s${s}`}</output>
    </div>
  {/if}
</fieldset>

<style>
  fieldset { border: 0; padding: 0; margin: .5rem 0 1rem; min-width: 0; }
  legend { font-size: .8rem; color: #6b7280; margin-bottom: .5rem; }
  .families { display: grid; gap: .5rem; }
  label { display: flex; align-items: center; gap: .55rem; font-size: .9rem; cursor: pointer; }
  input { accent-color: #7c3aed; }
  .family-slider { margin-top: 1rem; }
  output { display: block; text-align: right; margin-top: .4rem; font-size: .8rem; color: #6d28d9; }
</style>
