<script lang="ts">
	// I choose to use TypeScript,
	// which is almost like JavaScript -  but one needs to define types explicitely, when defining variables,
	//  which significantly helps in debugging.

	// now I import necessary modules
	import SliderParam from './SliderParam.svelte'; // it creates the slider elements for the parameters.
	import PlotLatch from './PlotLatch.svelte'; // it creates the svg our latching, when text is passed as arguments.
	import PlotLatch2 from './PlotLatch2.svelte'; // it creates the svg our latching, when text is passed as arguments.
	import katex from 'katex'; // this library renders plane-text strings to latex elements.

	// I define latex strings for parameters, using katex
	const S = katex.renderToString('S');
	const S_fnet = katex.renderToString('S^f');
	const S_pnet = katex.renderToString('S^p');
	const W = katex.renderToString('w');
	const W_fnet = katex.renderToString('w^f');
	const W_pnet = katex.renderToString('w^p');
	const Tau2 = katex.renderToString('\\tau_2');
	const Tau2_fnet = katex.renderToString('\\tau_2^f');
	const Tau2_pnet = katex.renderToString('\\tau_2^p');
	const L = katex.renderToString('\\lambda');

	// set default values for the parameters
	let valueS_fnet = 7;
	let valueS_pnet = 7;
	let valueW_fnet = 1.1;
	let valueW_pnet = 1.1;
	let valueTau2_fnet = 400.0;
	let valueTau2_pnet = 100.0;
	let valueLambda = 0.9;

	let Z_flag = 0;
	let Z = 1;

	// set the command for the backend python script (that calls the c++ executable)
	// set the textfile name using the user inputs.
	let command: string;
	let textfileName: string;
	$: {
		command = `python3 automate2.py ${valueW_fnet} ${valueW_pnet} ${valueS_fnet} ${valueS_pnet} ${valueTau2_fnet} ${valueTau2_pnet} ${valueLambda} ${Z_flag} ${Z}`;
		textfileName = `backend/data/Zf${Z_flag}_Z${Z}/lambP${valueLambda.toFixed(
			1
		)}_lambF${valueLambda.toFixed(1)}_${valueS_pnet}S${valueS_fnet}_${valueW_pnet.toFixed(
			2
		)}w${valueW_fnet.toFixed(2)}_${valueTau2_pnet.toFixed(1)}T${valueTau2_fnet.toFixed(
			1
		)}_p100_seed1/mall_cue0_seed1`;
	}
	// initialize the variable that stores the textfile content in frontend.
	let text: string = '';

	// create frontend for getting the textfile from the server.
	const getText = async (textfileName: string) => {
		const response = await fetch('/api', {
			method: 'POST',
			body: JSON.stringify({ textfileName })
		});
		const receivedText = await response.text();
		text = receivedText;
	};

	// I create frontend for calling python scripts at the server.
	let result: string | undefined = undefined;

	const runCppProgram = async (command: string) => {
		const response = await fetch('/api2', {
			method: 'POST',
			headers: {
				'Content-Type': 'application/json'
			},
			body: JSON.stringify({ command })
		});
		result = await response.json();
		// console.log(result);
	};

	// here I initialize a helper variable to track the state of an ongoing simulation.
	let intervalId: any = null;

	// I set the state variable and the state function for the start/stop button.
	// I initialize the state variable as false -> that will show 'start'
	let isRunning = false;
	// I want to trigger the state function every time the start/stop button is pressed
	// so that it can call certain other things.
	const handleClick = () => {
		// I call the script on 'start' button pressing only.
		if (!isRunning) {
			runCppProgram(command);
			// FIXME need to kill the backend cpp process upon pressing the 'stop' button.
			// maybe also delete the generated files of the backend?
			// for the moment, I'm manually killing the process.
		}

		// in both starting/stopping cases, I toggle the isRunning state
		isRunning = !isRunning;

		//when the command is running, I keep calling the getText function every 2 sec
		if (isRunning) {
			intervalId = setInterval(() => {
				getText(textfileName);
			}, 2000);
		}

		if (!isRunning) {
			// I stop execution if it is already running, and 'stop' button pressed
			clearInterval(intervalId);
			intervalId = null;
		}
	};
</script>

<!-- from here the script section ends, and the HTML starts -->
<!-- where we put different elements - divisions that have paragraphs, sliders, buttons, images, etc. -->
<!-- interestingly, in Svelte framework, the JavaScript/TypeScript variables can be called directly from HTML part. -->
<!-- when we change the values of variables programmetically, HTML elements also stay "reactive" to the changes -->

<!-- I don't understand fully the following part - but it adds some stylesheet to the HTML header that helps katex to render -->
<svelte:head>
	<link
		rel="stylesheet"
		href="https://cdn.jsdelivr.net/npm/katex@0.12.0/dist/katex.min.css"
		integrity="sha384-AfEj0r4/OFrOo5t7NnNe46zW/tFgW6x/bCJG8FqQCEo3+Aro6EYUG4+cU+KJWu/X"
		crossorigin="anonymous"
	/>
</svelte:head>

<!-- below I am trying to comment what each division does. -->
<div class="space-y-4 p-4 text-gray-700">
	<!-- model title -->
	<div class="flex flex-row place-content-center p-4">
		<p class="text-2xl">Demo - Potts Associative Network (Hybrid)</p>
	</div>
	<!-- rest of the content -->
	<div class="flex space-x-4">
		<div class="py-20">
			<SliderParam
				labelName={L}
				minValue={0.0}
				maxValue={1.0}
				bind:value={valueLambda}
				stepSize={0.1}
			/>
		</div>

		<!-- display the plot -->
		<div class="flex space-x-2">
			<div>
				<div class="flex place-content-center">
					<p class="text-2xl">fnet</p>
				</div>
				<div
					class="space-y-2 flex flex-col border-4 border-gray-200 place-items-left rounded p-4 h-[550px] w-[700px]"
				>
					<PlotLatch {text} />
					<!-- here we put the PlotLatch element with the argument it needs to be provided to render the svg. -->
				</div>
				<div class="px-40 py-4">
					<SliderParam
						labelName={S}
						minValue={3}
						maxValue={11}
						bind:value={valueS_fnet}
						stepSize={1}
					/>

					<SliderParam
						labelName={W}
						minValue={0.6}
						maxValue={2.0}
						bind:value={valueW_fnet}
						stepSize={0.2}
					/>

					<SliderParam
						labelName={Tau2}
						minValue={100.0}
						maxValue={800.0}
						bind:value={valueTau2_fnet}
						stepSize={100.0}
					/>
				</div>
				<div class="space-y-2 flex flex-col place-items-center bg-sky-500/[.06] rounded p-4">
					<p class="">fnet choices:</p>
					<p>{@html S} = {valueS_fnet}</p>
					<p>{@html W} = {valueW_fnet}</p>
					<p>{@html Tau2} = {valueTau2_fnet}</p>
				</div>
			</div>
			<div>
				<div class="flex place-content-center">
					<p class="text-2xl">pnet</p>
				</div>

				<div
					class="space-y-2 flex flex-col border-4 border-gray-200 place-items-left rounded p-4 h-[550px] w-[700px]"
				>
					<PlotLatch2 {text} />
					<!-- here we put the PlotLatch element with the argument it needs to be provided to render the svg. -->
				</div>
				<div class="px-40 py-4">
					<SliderParam
						labelName={S}
						minValue={3}
						maxValue={11}
						bind:value={valueS_pnet}
						stepSize={1}
					/>

					<SliderParam
						labelName={W}
						minValue={0.6}
						maxValue={2.0}
						bind:value={valueW_pnet}
						stepSize={0.2}
					/>

					<SliderParam
						labelName={Tau2}
						minValue={100.0}
						maxValue={800.0}
						bind:value={valueTau2_pnet}
						stepSize={100.0}
					/>
				</div>
				<div class="space-y-2 flex flex-col place-items-center bg-sky-500/[.06] rounded p-4">
					<p class="">pnet choices:</p>
					<p>{@html S} = {valueS_pnet}</p>
					<p>{@html W} = {valueW_pnet}</p>
					<p>{@html Tau2} = {valueTau2_pnet}</p>
				</div>
			</div>
		</div>
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
<!-- what is written inside the quotation marks following 'class=' is the easier way of writing HTML styling using TailwindCSS or UnoCSS. -->
<!-- in this project, I used UnoCSS. -->
