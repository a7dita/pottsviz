<script lang="ts">
	// import necessary modules
	import { onMount } from 'svelte';
	import SliderParam from './SliderParam.svelte';
	import PlotLatch from './PlotLatch.svelte';
	import katex from 'katex';

	// define latex strings for parameters, using katex
	const S = katex.renderToString('S');
	const W = katex.renderToString('w');
	const Tau2 = katex.renderToString('\\tau_2');

	// set default values for the parameters
	let valueS = 7;
	let valueW = 1.2;
	let valueTau2 = 200;

	// set state and state function for the start button
	let isRunning = false;

	const handleClick = () => {
		isRunning = !isRunning;
	};

	// create frontend for calling python scripts at the server
	let result: string | undefined = undefined;
	let command: string;

	const handleRun = async (command: string) => {
		const r = await fetch('/api', {
			method: 'POST',
			body: JSON.stringify({
				command
			})
		});
		result = await r.json();
	};
</script>

<svelte:head>
	<link
		rel="stylesheet"
		href="https://cdn.jsdelivr.net/npm/katex@0.12.0/dist/katex.min.css"
		integrity="sha384-AfEj0r4/OFrOo5t7NnNe46zW/tFgW6x/bCJG8FqQCEo3+Aro6EYUG4+cU+KJWu/X"
		crossorigin="anonymous"
	/>
</svelte:head>

<div class="space-y-16 p-4 text-gray-700">
	<div class="flex flex-row place-content-center">
		<p class="text-2xl">Demo - Potts Associative Network (Homogenous)</p>
	</div>
	<div class="flex space-x-4">
		<div class="space-y-4">
			<!-- display the parameter choices -->
			<div class="space-y-2 flex flex-col place-items-center bg-sky-500/[.06] rounded p-4">
				<p class="">Your choices:</p>
				<p>{@html S} = {valueS}</p>
				<p>{@html W} = {valueW}</p>
				<p>{@html Tau2} = {valueTau2}</p>
			</div>
			<!-- make the sliders -->
			<div>
				<SliderParam labelName={S} minValue={3} maxValue={11} bind:value={valueS} stepSize={1} />
				<SliderParam
					labelName={W}
					minValue={0.6}
					maxValue={2.0}
					bind:value={valueW}
					stepSize={0.2}
				/>
				<SliderParam
					labelName={Tau2}
					minValue={100}
					maxValue={800}
					bind:value={valueTau2}
					stepSize={100}
				/>
			</div>

			<!-- create the button for starting/stopping the simulation -->
			<div class="space-y-2 flex flex-col place-items-center p-4 text-sm">
				<button
					class="px-4 py-2 rounded"
					class:bg-pink-300={isRunning}
					class:bg-gray-200={!isRunning}
					class:text-gray-700={!isRunning}
					class:text-white={isRunning}
					on:click={handleClick}>{isRunning ? 'Stop' : 'Start Simulation'}</button
				>
			</div>
		</div>
		<!-- display the plot -->
		<div
			class="space-y-2 flex flex-col border-4 border-gray-200 place-items-center place-items-center rounded p-4 h-[550px] w-[700px]"
		>

	<div class="space-y-4">
		<input
			bind:value={command}
			type="text"
			class="border-2 border-gray-300 bg-gray-100 round py-1 px-4"
			placeholder="enter command to run..."
		/>
		<button class="px-4 py-1 rounded bg-gray-200" on:click={async () => await handleRun(command)}>
			Run command
		</button>

		<div class="border-2 border-gray-300 bg-blue-200 rounded h-40">
			{result}
			<PlotLatch />
		</div>
	</div>
</div>
