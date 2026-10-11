<script lang="ts">
  import { createEventDispatcher } from 'svelte';
  import SliderParam from './ParameterSlider.svelte';
  import { topologyForFamily, type TopologyFamily } from '$lib/topologies';
  export let topologyId = 'original';
  export let disabled = false;
  let family: TopologyFamily = 'original';
  let m = 0, r = 1;
  const dispatch = createEventDispatcher<{ change: string }>();
  const families: { value: TopologyFamily; label: string }[] = [
    { value: 'original', label: 'original' }, { value: 'random', label: 'random' },
    { value: 'modular', label: 'modular' }, { value: 'ring-local', label: 'rings' }
  ];
  function selectTopology() {
    topologyId = topologyForFamily(family, family === 'modular' ? m : r).id;
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
        tooltipText="m is the number of each memory’s seven partners inside its corresponding structural group. The remaining 7 − m partners are randomly chosen outside, while preserving seven partners on both sides. Groups need no semantic meaning. At m = 7, m7 has seven separate association blocks." />
      <output>m = {m} · m{m}</output>
    </div>
  {:else if family === 'ring-local'}
    <div class="family-slider" on:input={(event) => { r = Number((event.target as HTMLInputElement).value); selectTopology(); }}>
      <SliderParam ariaLabel="Ring-local parameter r" labelName="r" minValue={1} maxValue={7} bind:value={r} stepSize={1} {disabled}
        tooltipText="r is the number of each memory’s seven partners in its local ring neighbourhood: the corresponding memory and three neighbours in each direction, wrapping around 49 memories. The remaining 7 − r partners are randomly chosen outside. Both sides keep seven partners. r counts local connections, not jump distance or exactly shared targets. At r = 7, r7 is the fully local connected ring." />
      <output>r = {r} · r{r}</output>
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
