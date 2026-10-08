<script lang="ts">
  export let labelName: string;
  export let minValue: number;
  export let maxValue: number;
  export let value: number;
  export let stepSize: number;
  export let tooltipText = '';
  export let ariaLabel = 'Model parameter';
</script>
<label class="flex gap-3 items-center">
  <span title={tooltipText}>{@html labelName}</span>
  <span>{minValue}</span>
  <input aria-label={ariaLabel} type="range" min={minValue} max={maxValue} bind:value step={stepSize}
    class="w-full h-2 bg-gray-200 rounded-lg appearance-none cursor-pointer" />
  <span>{maxValue}</span>
</label>
