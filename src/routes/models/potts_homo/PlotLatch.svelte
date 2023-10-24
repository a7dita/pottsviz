<script lang="ts">
	import { onMount } from 'svelte';
	import * as d3 from 'd3';

	export let text: string = '';

	interface DataPoint {
		timeStep: number[];
		values: number[][];
	}

	let data: DataPoint = { timeStep: [], values: [] };

	// load data
	$: {
		let lines = text.split('\n');
		lines.shift();
		// console.log(lines);

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
	}

	// FIXME debug the extra line of entry
	$: {
		if (data.timeStep.length && data.values.length) {
			const svg = d3
				.select('#chart')
				.attr('width', Math.max(...data.timeStep) * 0.13)
				.attr('height', 500);
			const x = d3
				.scaleLinear()
				.domain([0, Math.max(...data.timeStep)])
				.range([0, Math.max(...data.timeStep) * 0.13]);

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
