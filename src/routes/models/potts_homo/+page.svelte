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

	// set the command for the backend c++ code
	let command: string;
	let textfileName: string;
	$: {
		command = `python3 automate.py ${valueS} ${valueW} ${valueTau2}`;
		textfileName = `mall_S${valueS}_w${valueW}0_gA0.5_T${valueTau2}.0_cue0`;
	}

	// set state and state function for the start button
	let isRunning = false;
	let text: any = '';
	$: {
		console.log(text);
	}

	const handleClick = () => {
		// if (!isRunning) {
		// runCppProgram(command);
		// }
		if (isRunning) {
			getText(textfileName);
		}
		isRunning = !isRunning;
	};

	// create frontend for calling python scripts at the server
	const getText = async (textfileName: string) => {
		const response = await fetch('/api', {
			method: 'POST',
			body: JSON.stringify({ textfileName })
		});
		const receivedText = await response.text();
		// console.log(receivedText);
		// console.log(typeof receivedText);
		text = receivedText;
	};

	let result: string | undefined = undefined;

	const runCppProgram = async (command: string) => {
		const response = await fetch('/api2', {
			method: 'POST',
			headers: {
				'Content-Type': 'application/json'
			},
			body: JSON.stringify({ command })
		});
		result = await response.json(); // processing the text response
		// Display or use that text file content
		console.log(result);
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
			class="space-y-2 flex flex-col border-4 border-gray-200 place-items-left rounded p-4 h-[550px] w-[700px]"
		>
			<PlotLatch {text} />
		</div>
	</div>
</div>
