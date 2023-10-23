<script lang="ts">
	import { onMount } from 'svelte';
	import * as d3 from 'd3';

	export let text: string = '';
	$: {
		// console.log(text);
	}
	interface DataPoint {
		timeStep: number[];
		values: number[][];
	}

	let data: DataPoint = { timeStep: [], values: [] };

	// Fetch data

	// FIXME THE ONLY FINAL STEP MISSING (AND I COULD NOT FIX SO FAR) IS - ONMOUNT FUNCTION IS NOT REACTIVE.
	// HENCE THE DATA IS NOT GETTING UPDATED WHEN I AM TRYING TO LOAD IT FROM THE VARIABLE 'TEXT'.
	// WE WILL HAVE TO FIND ANOTHER WAY.
	onMount(() => {
		// const response = await fetch('../src/routes/data/mall_cue0_seed1_d');
		// const fileContent = await response.text();
		// const lines = fileContent.split('\n');
		console.log(text);
		let lines = text.split('\t');

		lines.forEach((line, i) => {
			const parts = line.split('\t').map(Number);
			data.timeStep[i] = parts[0];

			parts.slice(1).forEach((value, j) => {
				if (!data.values[j]) {
					data.values[j] = [];
				}
				data.values[j][i] = value;
			});
		});
	});

	// FIXME debug the extra line of entry
	$: {
		if (data.timeStep.length && data.values.length) {
			const svg = d3
				.select('#chart')
				.attr('width', Math.max(...data.timeStep) * 1.3)
				.attr('height', 500);
			// console.log(Math.max(...data.timeStep) * 1.3);
			const x = d3
				.scaleLinear()
				.domain([0, Math.max(...data.timeStep)])
				.range([0, Math.max(...data.timeStep) * 1.3]);

			const y = d3.scaleLinear().domain([0, 1]).range([500, 0]);

			data.values.forEach((lineData, j) => {
				const line = d3
					.line<number>()
					.x((d, i) => x(data.timeStep[i]))
					.y((d) => y(d));

				svg
					.append('path')
					.datum(lineData)
					.attr('d', line)
					.attr('stroke-width', 1.2)
					.attr('stroke', `hsl(${j * (360 / data.values.length)}, 100%, 50%)`) // Each line will have a unique hue
					// .attr('stroke', d3.schemeCategory10[j])
					.attr('fill', 'none');
			});
		}
	}
</script>

<div class="">
	<svg id="chart" />
</div>
