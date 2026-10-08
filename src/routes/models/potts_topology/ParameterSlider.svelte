<script lang="ts">
  export let labelName: string;
  export let minValue: number;
  export let maxValue: number;
  export let value: number;
  export let stepSize: number;
  export let tooltipText: string;
  export let ariaLabel: string;
  export let disabled = false;
  $: helpId = `parameter-${ariaLabel.toLowerCase().replace(/[^a-z0-9]+/g, '-')}-help`;
</script>

<div class="control" class:muted={disabled}>
  <div class="symbol">
    <button type="button" aria-label={`${ariaLabel}: description`} aria-describedby={helpId} title={tooltipText} {disabled}>{@html labelName}</button>
    <div class="help" id={helpId} role="tooltip">{tooltipText}</div>
  </div>
  <span class="bound">{minValue}</span>
  <input aria-label={ariaLabel} type="range" min={minValue} max={maxValue} bind:value step={stepSize} {disabled} />
  <span class="bound">{maxValue}</span>
</div>

<style>
  .control { display: grid; grid-template-columns: 3.25rem 2.1rem minmax(0, 1fr) 2.1rem; gap: .45rem; align-items: center; }
  .symbol { position: relative; }
  button { display: block; width: 100%; text-align: left; cursor: help; color: #374151; }
  button:focus-visible { outline: 2px solid #7c3aed; outline-offset: 3px; border-radius: .2rem; }
  .bound { font-size: .75rem; text-align: center; font-variant-numeric: tabular-nums; color: #6b7280; }
  input { min-width: 0; width: 100%; margin: 0; accent-color: #7c3aed; cursor: pointer; }
  .help { display: none; position: absolute; right: 0; top: 100%; z-index: 10; width: min(20rem, 85vw); padding: .7rem; background: white; color: #374151; border: 1px solid #ddd6fe; box-shadow: 0 6px 20px #0002; border-radius: .4rem; font-size: .8rem; line-height: 1.5; }
  .symbol:hover .help, .symbol:focus-within .help { display: block; }
  .muted button, .muted .bound { color: #b7bac4; }
  .muted input { accent-color: #d1d5db; opacity: .45; cursor: default; }
</style>
