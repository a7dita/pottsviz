<script lang="ts">
	// I choose to use TypeScript,
	// which is almost like JavaScript -  but one needs to define types explicitely, when defining variables,
	// which significantly helps in debugging.

	// now I import necessary modules
	import SliderParam from './SliderParam.svelte'; // it creates the slider elements for the parameters.
	import PlotLatch from './PlotLatch.svelte'; // it creates the svg our latching, when text is passed as arguments.
	import katex from 'katex'; // this library renders plane-text strings to latex elements.

	// I define latex strings for parameters, using katex
	const S = katex.renderToString('S');
	const W = katex.renderToString('w');
	const Tau2 = katex.renderToString('\\tau_2');

	// I set default values for the parameters
	let valueS = 7;
	let valueW = 1.2;
	let valueTau2 = 200;

	// I set the command for the backend python script (that calls the c++ executable)
	// I set the textfile name using the user inputs.
	let command: string;
	let textfileName: string;
	$: {
		command = `python3 automate.py ${valueS} ${valueW} ${valueTau2}`;
		textfileName = `mall_S${valueS}_w${valueW}0_gA0.5_T${valueTau2}.0_cue0`;
		// TODO we can send multiple textfile names here,
		// to pass them to multiple plot elements in separate divisions,
		// to get 3 plots for our 3 cues.
	}
	// I initialize the variable that stores the textfile content in frontend.
	let text: string = '';

	// I create frontend for getting the textfile from the server.
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

	// here I initialize a helper variable to track the state of an ongoing execution - for calling our getText function again and again.
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
		}

		// in both starting/stopping cases, I toggle the isRunning state
		isRunning = !isRunning;

		//when the command is running, I keep calling the getText function.
		// FIXME I need to stop the intervaled execution upon 'stop' button press
		// it is not working and I need to manually refresh the page upon finishing the simulation.
		if (isRunning) {
			// I stop execution if it was already running
			if (intervalId) {
				clearInterval(intervalId);
				intervalId = null;
			} else {
				// if not running, I start execution
				// I call the function every 2 seconds
				intervalId = setInterval(() => {
					getText(textfileName);
				}, 2000);
			}
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
<div class="space-y-16 p-4 text-gray-700">
	<!-- model title -->
	<div class="flex flex-row place-content-center">
		<p class="text-2xl">Demo - Potts Associative Network (Homogenous)</p>
	</div>
	<!-- rest of the content -->
	<div class="flex space-x-4">
		<div class="space-y-4">
			<!-- display the parameter choices -->
			<div class="space-y-2 flex flex-col place-items-center bg-sky-500/[.06] rounded p-4">
				<p class="">Your choices:</p>
				<p>{@html S} = {valueS}</p>
				<p>{@html W} = {valueW}</p>
				<p>{@html Tau2} = {valueTau2}</p>
			</div>
			<!-- create the sliders -->
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
			<!-- here we put the PlotLatch element with the argument it needs to be provided to render the svg. -->
		</div>
	</div>
</div>
<!-- what is written inside the quotation marks following 'class=' is the easier way of writing HTML styling using TailwindCSS or UnoCSS. -->
<!-- in this project, I used UnoCSS. -->
