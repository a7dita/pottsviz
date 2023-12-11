<script lang="ts">
	import * as d3 from 'd3'; // I import the plotting library

	// here I export the variable, so that the value can be changed from +page.svelte (main page for this directory)
	export let isOdd = false;
	export let text: string = '';

	// I create an object type
	interface DataPoint {
		timeStep: number[];
		values: number[][];
	}

	// ...and an instance of that object type.
	let data: DataPoint = { timeStep: [], values: [] };

	// I read data from the text variable
	$: {
		let lines = text.split('\n');
		// lines.shift(); // remove the first line, which is empty (in the present version of the c++ code)
		if (isOdd) {
			lines = lines.filter((_, i) => i % 2 !== 0); // selecting odd lines only
		} else {
			lines = lines.filter((_, i) => i % 2 === 0);
		}

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

	$: {
		if (data.timeStep.length && data.values.length) {
			const svg = d3
				.select('#chart2') // 'chart is the ID of the svg element below.'
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
					.attr('stroke', `hsl(${j * (360 / data.values.length)}, 100%, 50%)`) // I try to give each line will have a unique hue
					// .attr('stroke', d3.schemeCategory10[j])
					.attr('fill', 'none');
			});
		}
	}
</script>

<div class="">
	<svg id="chart2" />
</div>
