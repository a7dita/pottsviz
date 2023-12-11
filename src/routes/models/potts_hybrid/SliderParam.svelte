<script lang="ts">
	export let labelName: string;
	export let minValue: number;
	export let maxValue: number;
	export let value: number;
	export let stepSize: number;
	export let tooltipText: string = '';
</script>

<div class="flex flex-row space-x-4 space-y-4 place-items-center">
	<!-- FIXME make the sliders equal -->
	<label for="steps-range" class="block mb-2 text-sm font-medium text-gray-900 dark:text-white" />

	<button data-toggle="tooltip" title={tooltipText}>
		{@html labelName}
	</button>
	<p>{minValue}</p>
	<input
		id="steps-range"
		type="range"
		min={minValue}
		max={maxValue}
		bind:value
		step={stepSize}
		class="w-full h-2 bg-gray-200 rounded-lg appearance-none cursor-pointer dark:bg-gray-700"
	/>
	<p>{maxValue}</p>
</div>
